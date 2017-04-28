# -*- coding: UTF-8 -*-

import numpy as np
import math, sys, time, cmath
import multiprocessing as mp
from spin_states import *
#from superoperator import *
try:
	import qutip as qt
	
except:
	print ' Some modules are missing, try installing them using \'installModules.py\''
	sys.exit()
THRUSTY = 1
NPROCS = 8
UPARROW = u'\u2191'
DOWNARROW = u'\u2193'



def alfa(K,I,M,nu,rho, labda):
	c = labda * (I*(I+1) - M*(M+1))**0.5

	return 2*c

def beta(K,I,M,nu,rho, labda):
	a = - 0.5 * K * rho * nu -(nu -1)*M
	b =  0.5 * K * rho * nu -(nu +1)*(M+1)
	c = labda * (I*(I+1) - M*(M+1))**0.5

	return a-b+((a-b)**2+4*c**2)**0.5

def phi(K,I,M,nu,rho, labda):
	return cmath.phase( 2*M+1 - nu *(K * rho -1) + 2* labda * 1j * ( I*(I+1) - M*(M+1) )**0.5 )

def x(M,nu):
	return ( -1-(2* M+1)*nu)

def y(K,I,M,nu,rho, labda):
	return ( ((2*I+1)**2)*labda**2 - ((2*M+1)**2)*(labda**2-1) + nu*(K*rho-1)*(nu*(K*rho-1) -4*M-2) )**0.5
def Ifunc(x, tau_lead):
	threshold = 14
	if x == 0:
		return tau_lead
	
	elif tau_lead == 0:
		if x > 0:
			return 0
		elif x < 0:
			return -x.real
	
	elif x/tau_lead < -threshold:
		return -x.real
	elif x/tau_lead > threshold:
		return 0
	return (x / (math.exp(x/tau_lead)-1)).real



	

class Computation(object):
	def __init__(self, master):
		self.master = master
	
	
	
	
	def defineOperators(self):
		#define all the operators here
		###OPERATOR LIST, if you add any, add it to this list too
		#Sm(.dag()): electron spin destruction (creation) operator
		#Sx, Sy, Sz: electron spin
		#Im(.dag()): nuclear spin destruction (creation) operator
		#Ix, Iy, Iz: nuclear spin
		#H: hyperfine hamiltonian, defined as in the notes
		#M_stg: staggered magnetization, defined as in the notes
		def electronOperator(QObj):
			return qt.tensor([QObj] + [qt.qeye(2) for _ in xrange(self.parameters.K)])
		self.IdentityElec = qt.qeye(2)
		self.Sm = qt.sigmam()
		self.Sx = 0.5 * qt.sigmax()
		self.Sy = 0.5 * qt.sigmay()
		self.Sz = 0.5 * qt.sigmaz()
		

		def nuclearOperator(QObj, n):
			nucleiList = [qt.qeye(2) for _ in xrange(self.parameters.K)]
			nucleiList[n] = QObj
			return qt.tensor(nucleiList)
		s = 0.5
		self.IdentityNucl = qt.qeye(self.parameters.K + 1)
		self.Ix = qt.jmat((self.parameters.K/2.), 'x')
		self.Iy = qt.jmat((self.parameters.K/2.), 'y')
		self.Iz = qt.jmat((self.parameters.K/2.), 'z')
		self.Isq = (self.Ix**2 + self.Iy**2 + self.Iz**2).tidyup()
	
		self.Ip = qt.jmat((self.parameters.K/2.), '+')
		self.Im = qt.jmat((self.parameters.K/2.), '-')
		print 'Evolution Magnetic field: ', [self.parameters.nu_evol for n in xrange(self.parameters.K)]
		
		self.H_elec_evol =  - self.parameters.K * self.parameters.rho_evol * self.parameters.nu_evol * self.Sz
		self.H_nuclei_evol = -self.parameters.nu_evol * self.Iz
		
		self.H_hf_evol = 2 * qt.tensor(self.Sz, self.Iz) + self.parameters.labda * qt.tensor(self.Sm.dag(),  self.Im) + self.parameters.labda * qt.tensor(self.Sm, self.Im.dag() )
		
		self.H_nogradient_evol = qt.tensor(self.H_elec_evol, self.IdentityNucl)  + qt.tensor(self.IdentityElec, self.H_nuclei_evol) + self.H_hf_evol 
		self.H_evol = self.H_nogradient_evol 
		self.H_elec_init =  - self.parameters.K * self.parameters.rho_init * self.parameters.nu_init * self.Sz
		self.H_nuclei_init = -self.parameters.nu_init * self.Iz
		self.H_hf_init = 2 * qt.tensor(self.Sz, self.Iz) + qt.tensor(self.Sm.dag(),  self.Im) + qt.tensor(self.Sm, self.Im.dag() )
		self.H_init = qt.tensor(self.H_elec_init, self.IdentityNucl)  + qt.tensor(self.IdentityElec, self.H_nuclei_init) + self.H_hf_init 
		self.M_stg = qt.tensor(self.Sz,  self.IdentityNucl) -  qt.tensor(self.IdentityElec, self.Iz) / self.parameters.K
		self.S2 = qt.tensor(self.Sx,  self.IdentityNucl) **2 + qt.tensor(self.Sy,  self.IdentityNucl) **2 + qt.tensor(self.Sz,  self.IdentityNucl) **2		
		self.I2 = qt.tensor(self.IdentityElec, self.Ix) **2 + qt.tensor(self.IdentityElec, self.Iy) **2 + qt.tensor(self.IdentityElec, self.Iz) **2
		self.H_init = qt.Qobj(self.H_init.full())
		self.H_nogradient_evol = qt.Qobj(self.H_nogradient_evol.full())
		self.H_evol = qt.Qobj(self.H_evol.full())


		if self.H_init == self.H_evol:
			self.isHConserved = True
		else:
			self.isHConserved = False
	def getEigenstate(self,ps,I,M,nu,rho, labda,p, ham,ns1,ns2 = False):
		
		if ps == 'u':
			state = qt.tensor(qt.basis(2, 0), ns1['state'])
			te = -0.5*self.parameters.K*rho*nu-I*(nu-1)
			rep = '|%s >|%g %g %d>'%(UPARROW, ns1['M'], ns1['I'], ns1['p'])
		elif ps == 'd':
			state = qt.tensor(qt.basis(2, 1), ns1['state'])	
			te = 0.5*self.parameters.K*rho*nu+I*(nu+1)
			rep = '|%s >|%g %g %d>'%(DOWNARROW, ns1['M'], ns1['I'], ns1['p'])	
		elif ps == '+':
			f = phi(self.parameters.K,I,M,nu,rho, labda)
			state = math.cos(f/2.) * qt.tensor(qt.basis(2, 0),ns1['state']) + math.sin(f/2.) * qt.tensor(qt.basis(2, 1),ns2['state'])
			te = 0.5* x(M,nu)  + 0.5*  y(self.parameters.K,I,M,nu,rho, labda)
			rep = '%+.3f|%s >|%g %g %d> %+.3f|%s >|%g %g %d>'%(math.cos(f/2.), UPARROW, ns1['M'], ns1['I'], ns1['p'], math.sin(f/2.), DOWNARROW, ns2['M'], ns2['I'], ns2['p'])	
		elif  ps == '-':
			f = phi(self.parameters.K,I,M,nu,rho, labda)
			state = -math.sin(f/2.) * qt.tensor(qt.basis(2, 0),ns1['state']) + math.cos(f/2.) * qt.tensor(qt.basis(2, 1),ns2['state'])
			te = 0.5* x(M,nu)  - 0.5* y(self.parameters.K,I,M,nu,rho, labda)
		 	rep = '%+.3f|%s >|%g %g %d> %+.3f|%s >|%g %g %d>'%(-math.sin(f/2.), UPARROW, ns1['M'], ns1['I'], ns1['p'], math.cos(f/2.), DOWNARROW, ns2['M'], ns2['I'], ns2['p'])	
		return {'ps' : ps,'I' : I,'M' : M,'p' : p,'state' : qt.Qobj(state.full()).tidyup().unit(), 'theorEnergy': te,'energy' : ham.matrix_element(state.dag(),state).real, 'representation' : rep} 
	
		
	def calculateEigenstates(self,nu, rho, labda, ham):
		eigenstates = []
		for ns in self.nucleiStates:
			I, M, p = ns['I'], ns['M'], ns['p']
			if round(M,2) == round(I,2):
				eigenstates.append(self.getEigenstate('u',I,M,nu, rho, labda,p, ham,ns))
			
			else:
				for ps in ('-', '+'):
					ns2 = [elem for elem in self.nucleiStates if (elem['I'] == I and elem['p'] == p and elem['M'] == M + 1)][0]
					
					eigenstates.append(self.getEigenstate(ps,I,M,nu, rho, labda,p, ham,ns,ns2))
			if M == -I:
				eigenstates.append(self.getEigenstate('d',I,M,nu, rho, labda,p, ham,ns))
				
		return 	eigenstates

	
	def checkEigenstates(self, states, Ham):
		for i in xrange(len(states)):
			es = states[i]
			state = es['state']
			if round(es['energy'], 12) != round(es['theorEnergy'], 12):
				print 'Wrong check ',es['I'],es['M'],es['ps'],es['p']
				print es['energy'], es['theorEnergy']	
				
			for j in xrange(i+1, len(states)):
			
				es2 = states[j]
				ol = es['state'].overlap(es2['state'])
				if ol > 10**(-10):
					print 'non zero overlap', es, es2,ol
			
	def defineStates(self):	
		self.eigenstatesInit = self.calculateEigenstates(self.parameters.nu_init, self.parameters.rho_init, 1, self.H_init)#[{'ps' : 0,'I' : 0,'M' : 0,'p' : 0,'state' : self.H_init.eigenstates()[1][i], 'theorEnergy': self.H_init.eigenstates()[0][i],'energy' : self.H_init.eigenstates()[0][i], 'representation' : 0} for i in xrange(len(self.H_init.eigenstates()[0]))]#
		self.eigenstatesEvol = self.calculateEigenstates(self.parameters.nu_evol, self.parameters.rho_evol, self.parameters.labda, self.H_nogradient_evol)# the eigenstates must be calculated with no MF gradient
		if not THRUSTY:
			print 'Check Init'
			self.checkEigenstates(self.eigenstatesInit, self.H_init)
			print 'Check Evol'
			self.checkEigenstates(self.eigenstatesEvol, self.H_nogradient_evol)
		 
		self.eigenstatesInit = sorted([s for s in self.eigenstatesInit], key = lambda s : s['energy'])
		groundStates = [s for s in self.eigenstatesInit if round(s['energy'], 15) == round(self.eigenstatesInit[0]['energy'], 15)]
		if not THRUSTY or 1:
			print 'Initial Hamiltonian: ',self.H_init
			print 'Eigenstates: (energy)(theorEnergy) State', len(self.eigenstatesInit)	
			for gs in self.eigenstatesInit:
				print '(%.3f)(%.3f) |%s %g %+g %d> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['I'],gs['M'],gs['p'], gs['representation'])
		print 'Groundstates: (energy)(theorEnergy) State'			
		for gs in groundStates:
			print '(%.3f)(%.3f) |%s %g %+g %d> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['I'],gs['M'],gs['p'], gs['representation'])
		if self.parameters.tau_dot > 0:
			dm = (-self.H_init/self.parameters.tau_dot).expm().matrix_element(self.eigenstatesInit[0]['state'].dag(), self.eigenstatesInit[0]['state']) * qt.ket2dm(self.eigenstatesInit[0]['state'])
			for es in self.eigenstatesInit[1:]:
				dm += (-self.H_init/self.parameters.tau_dot).expm().matrix_element(es['state'].dag(), es['state']) * qt.ket2dm(es['state'])
				

		else:
			dm = qt.ket2dm(groundStates[0]['state'])
			for gs in groundStates[1:]:
				dm += qt.ket2dm(gs['state'])
				
		return dm/dm.tr()
	def getHamiltonian(self, K, nu, rho, labda):
		H_elec =  - K * rho * nu * self.Sz
		H_nuclei = -nu * self.Iz
		
		H_hf = 2 * qt.tensor(self.Sz, self.Iz) + labda * qt.tensor(self.Sm.dag(),  self.Im) + labda * qt.tensor(self.Sm, self.Im.dag() )
		
		H = qt.Qobj((qt.tensor(H_elec, self.IdentityNucl)  + qt.tensor(self.IdentityElec, H_nuclei) + H_hf ).full())

		return H

		
	def initializeDissipators(self,states, meanEnergy, Ham, printTrans = False):
		baseRate = (self.parameters.gamma0**2/ (2*math.pi*(self.parameters.mu-meanEnergy)**2))**0.5
		intermediate = lambda ns : qt.tensor(qt.basis(2,0), ns) +  qt.tensor(qt.basis(2,1), ns)
		if printTrans:
			self.transition_rules = []
			self.state_index = 1
			f = open('transitions.txt' , 'w+')
			f.write('state1 state2 rate Iz_exp_value')
			f.flush()
			g = open('rules.txt', 'w+')
			
			
			g.write('nu %f, Krho %f, gammatilde %f \n'%(self.parameters.nu_evol, self.parameters.K * self.parameters.rho_evol, baseRate))
			g.flush()
			#f.write('Transition |ps M I p> --> |ps M I p> -- baserate * overlap * Ifunction\n')
			#f.flush()
		def getStateIndex(state):
			for s in self.transition_rules:
				if state == s['state']:
					index = s['index']
					break
			else:
				index = int(self.state_index)
				self.state_index += 1
				rule = '%s %+g %g %d = %d \n'%(state['ps'], state['M'],state['I'],state['p'], index)
				self.transition_rules.append({'state' : state, 'index' : index})
				#print index, self.state_index, rule
				g.write(rule)
				g.flush()
					
			return index
	

		def addDissipator(state1, state2, Ham):
			dissipators12 = [] 
			initial = state1['state']
			final = state2['state']
			overlap = math.fabs( sum([initial.overlap(intermediate(ns['state']))  *  final.overlap(intermediate(ns['state'])) for ns in self.nucleiStates]) )**2
			if overlap > 10**(-15):
				ifunction = Ifunc(Ham.matrix_element(final.dag(), final) - Ham.matrix_element(initial.dag(), initial), self.parameters.tau_lead)
				rate = baseRate * (overlap * ifunction)**0.5
				print 'Transition |%s %+g %g %d> --> |%s %+g %g %d> -- rate %.3f'%(state1['ps'], state1['M'],state1['I'],state1['p'],state2['ps'], 
						state2['M'],state2['I'],state2['p'],rate)
				if printTrans:				
					#print state1['state'], qt.tensor(self.IdentityElec, self.H_nuclei_evol+self.H_nuclei_gradient)
					m_value = qt.tensor(self.IdentityElec, self.Iz).matrix_element(state1['state'].dag(), state1['state'])
					f.write('%d %d %f %f\n'%(getStateIndex(state1),getStateIndex(state2),rate, m_value))
					f.flush()
				if ifunction > 10**(-15):
					dissipators12.append(rate * final * initial.dag())
				if state1 != state2:
					ifunction = Ifunc(Ham.matrix_element(initial.dag(), initial) - Ham.matrix_element(final.dag(), final), self.parameters.tau_lead)
					rate = baseRate * (overlap * ifunction)**0.5
					print 'Transition |%s %+g %g %d> --> |%s %+g %g %d> -- rate %.3f'%(state2['ps'], state2['M'],state2['I'],state2['p'],state1['ps'], 
						state1['M'],state1['I'],state1['p'],rate)
					if printTrans:			
						m_value = qt.tensor(self.IdentityElec, self.Iz).matrix_element(state2['state'].dag(), state2['state'])		
						f.write('%d %d %f %f\n'%(getStateIndex(state2), getStateIndex(state1), rate, m_value ))
						f.flush()						
					if ifunction > 10**(-15):					
						dissipators12.append(rate * initial * final.dag())
			
			return dissipators12
		

		def getDisList(chunk, meanEnergy, Ham, out_q):
			dissipators = []
			for i,j in chunk:	
				state1 = states[i]
				state2 = states[j]

			
				if not THRUSTY or (state1['I'] == state2['I'] and state1['p']== state2['p'] and -1 <= state1['M'] - state2['M'] <= 1):#this is just to speed things up. only transitions allowed conserve I and p
					 dissipators += addDissipator(state1, state2, Ham)
							
			out_q.put(dissipators)

		
		out_q = mp.Queue()
		
		procs = []
		DisList = []
		indeciesPairs = []
		#chunkedStates is a list of pairs of indecies, we cutted the original double for loop to divide the job between multiple threads
		for i in xrange(len(states)):
			for j in xrange(i, len(states)):
				indeciesPairs.append((i,j))	
		size = max(1, int(round(1.*len(indeciesPairs)/NPROCS)))
		chunkedStates =[indeciesPairs[s:s + size] for s in xrange(0, len(indeciesPairs), size)]
		
		for cs in chunkedStates:
			p = mp.Process(target = getDisList, args = (cs, meanEnergy, Ham, out_q))
			procs.append(p)
			p.start()
		for p in procs:
			DisList.extend(out_q.get())
		for p in procs:
			p.join()
		if printTrans:
			f.close()	
			g.close()	
		return DisList
	

	def initializeLiouvillians(self, name, Ham, collapseOperators):
		l = qt.liouvillian(Ham, collapseOperators)
		setattr(self, name, l)
	def get_operators_compute(self):
		operators = []
		if self.master.IzBool.get():
			operators.append(qt.tensor(self.IdentityElec, self.Iz))
		
			#op = qt.tensor(self.IdentityElec, self.Iz) * qt.tensor(self.IdentityElec, self.Iz)
			#operators.append(op)
		if self.master.SzBool.get():
			operators.append(qt.tensor(self.Sz, self.IdentityNucl))
			
		if self.master.M_stg_bool.get():
			operators.append(self.M_stg)
		if self.master.HBool.get():
			operators.append(self.H_init)
			operators.append(self.H_evol)
		if self.master.I2Bool.get():
			operators.append(self.I2)
		if self.master.S2Bool.get():
			operators.append(self.S2)
		if self.master.CohBool.get():
			operators.append(qt.tensor(self.IdentityElec, self.Ip))
			operators.append(qt.tensor(self.IdentityElec, self.Ip * self.Im))
		if self.master.IzNormBool.get():
			operators.append(qt.tensor(self.IdentityElec, self.Iz)/(self.parameters.K/2.))
		return operators

	def compute_phase_diagram(self,parameters):
		self.parameters = parameters
		self.defineOperators()
		N = 30
		
		xmin, xmax = min(0,self.parameters.nu_init,self.parameters.nu_evol), max(self.parameters.nu_init,self.parameters.nu_evol)+1
		ymin, ymax = 0, max(self.parameters.rho_init,self.parameters.rho_evol)+1#min(self.parameters.rho_init,self.parameters.rho_evol)
		
		
		x_values = np.linspace(xmin, xmax, num=N, endpoint=True)
		y_values = np.linspace( ymax,ymin, num=N, endpoint=True)
		pdbuffer = []
		T = max(self.parameters.tau_dot, 0.1)
		iz1 = qt.Qobj( qt.tensor( qt.tensor(self.IdentityElec, self.Iz), qt.qeye(2 * (self.parameters.K + 1))).full() )
		iz2 = qt.Qobj( qt.tensor( qt.qeye(2 * (self.parameters.K + 1))).full(), qt.tensor(self.IdentityElec, self.Iz) )
		
		for nu in x_values:
			for rho in y_values:
				h1 = qt.tensor(self.getHamiltonian(self.parameters.K, nu, rho , 0), qt.qeye(2 * (self.parameters.K + 1)) ) 	
				h2 = qt.tensor( qt.qeye(2 * (self.parameters.K + 1)), self.getHamiltonian(self.parameters.K, -0.5 * nu, 2 * rho , 0) ) 
				h_int = - 5 * qt.tensor( qt.tensor( self.Sz, qt.qeye( (self.parameters.K + 1)) ) , qt.tensor( self.Sz, qt.qeye( (self.parameters.K + 1)) ) )
				#print h1, h2, h_int
				h =qt.Qobj( h1.full() ) + qt.Qobj(h2.full() ) + qt.Qobj(h_int.full() )
				ens = zip(h.eigenstates()[0], h.eigenstates()[1])
				pdbuffer.append(sum( [math.exp(-s[0]/T) * qt.expect(iz1 ,s[1]) for s in ens] )/sum( [math.exp(-s[0]/T) for s in ens] ))
				

		



		pd = np.ndarray( shape=(N,N), dtype=float,  buffer=np.array(pdbuffer) )
		print pd
		return (y_values,x_values,  pd.transpose() )


	def compute_dynamics(self, parameters):
		startTime = time.time()
		self.parameters = parameters
		self.nucleiStates = [{'I' : self.parameters.K/2.,'M' : - i + self.parameters.K/2.,'p' : 1,'state' : qt.basis(self.parameters.K + 1, i)} for i in xrange(self.parameters.K + 1)]#decompose(self.parameters.K, self.parameters.K/2.)
		self.nucleiStates = sorted([el for el in self.nucleiStates], key = lambda el: (el['I'], el['M']))
		if not THRUSTY:
			represent(self.nucleiStates)
		self.defineOperators()
		self.density_matrix = self.defineStates()
		statesTime = time.time()
		print 'States and operators got in %f'%(statesTime-startTime)
		self.operators_compute = self.get_operators_compute()
		
		
		self.meanEnergyInit = sum([es['energy'] for es in self.eigenstatesInit])/len(self.eigenstatesInit)
		self.meanEnergyEvol = sum([es['energy'] for es in self.eigenstatesEvol])/len(self.eigenstatesEvol)
		if self.operators_compute:
			if self.parameters.gamma0 > 0:
				print 'Transition |ps M I p> --> |ps M I p> -- rate [baseRate * overlap * Ifunction]'
				if self.parameters.dis_time >= self.parameters.evo_time:
					DisListEvol = self.initializeDissipators(self.eigenstatesEvol, self.meanEnergyEvol, self.H_evol, printTrans = True)
					DisListInit = []
				else:
					DisListInit = self.initializeDissipators(self.eigenstatesInit, self.meanEnergyInit, self.H_init)
					if not self.isHConserved:
						DisListEvol = self.initializeDissipators(self.eigenstatesEvol, self.meanEnergyEvol, self.H_evol, printTrans = True)
					elif self.isHConserved:
						DisListEvol = DisListInit
			else:
				DisListInit = []
				DisListEvol = []
			dis_time = time.time()
			print 'Dissipators (%d) got in %f'%(len(DisListEvol), dis_time-statesTime)
			if self.parameters.gamma0 > 0:

				if self.parameters.dis_time <= self.parameters.startTime:
					LInitDiss = qt.liouvillian(self.H_init, DisListInit)	
					LInitNoDiss = LInitDiss
				else:
					LInitDiss = qt.liouvillian(self.H_init, DisListInit)	
					LInitNoDiss = qt.liouvillian(self.H_init, [])
				if not self.isHConserved:
					if self.parameters.dis_time >= self.parameters.evo_time:
						LEvolNoDiss = qt.liouvillian(self.H_evol, [])
						LEvolDiss = qt.liouvillian(self.H_evol, DisListEvol)
					else:
						LEvolDiss = qt.liouvillian(self.H_evol, DisListEvol)
						LEvolNoDiss = LEvolDiss
				else:
					LEvolNoDiss = LInitNoDiss
					LEvolDiss = LInitDiss
		
			else:
				LInitNoDiss = qt.liouvillian(self.H_init, [])
				LInitDiss = LInitNoDiss
				if not self.isHConserved:
					LEvolNoDiss = qt.liouvillian(self.H_evol, [])
					LEvolDiss = LEvolNoDiss
				else:
					LEvolNoDiss = LInitNoDiss
					LEvolDiss = LInitDiss
			liouTime = time.time()
			print 'Liouvillians got in %f'%(liouTime-dis_time)		
		
		
		
			opts = qt.Odeoptions(nsteps = 10000000000, rhs_reuse=True, num_cpus = 4, rtol = 0.001, atol = 0.00001)
			tlist = np.linspace(0,self.parameters.tInterval,num=self.parameters.tSteps)

		
		
			def L_func(t, args):
				if t <= self.parameters.evo_time:
					if t < self.parameters.dis_time:
						L = LInitNoDiss
					elif t >= self.parameters.dis_time:
						L = LInitDiss
				elif t > self.parameters.evo_time:
					if t < self.parameters.dis_time:
						L = LEvolNoDiss
					elif t >= self.parameters.dis_time:
						L = LEvolDiss
				return L.data 

			try:
				self.expt_list = qt.mesolve(L_func, self.density_matrix, tlist, [], self.operators_compute, progress_bar = True)
				with open('results.txt', 'w+') as h:	
					for r in list(self.expt_list.expect[0]):
						h.write('%f \n'%r)
				
				#self.expt_list = qt.mesolve([self.H_init, [(self.H_evol- self.H_init), H_evol_coeff]] , self.density_matrix, tlist, collapseList, self.operators_compute, options = opts)
				print 'States and operators got in %f'%(statesTime-startTime)
				print 'Dissipators (%d) got in %f'%(len(DisListEvol), dis_time-statesTime)			
				print 'Liouvillians got in %f'%(time.time()-dis_time)				
				print 'Finished in %f'%(time.time()-startTime)	

				
			except Exception as e:
				print e
				self.expt_list = 'Error'

		return self.expt_list
		
				
		




