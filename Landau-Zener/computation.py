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
THRUSTY = 0
HYSTERESIS = 0
NPROCS = 1
UPARROW = u'\u2191'
DOWNARROW = u'\u2193'





def phi(N,I,M,b,gammaI, alfa):
	return cmath.phase( 2*M+1 - b *(N * gammaI -1) + 2* alfa * 1j * ( I*(I+1) - M*(M+1) )**0.5 )

def x(M,b):
	return -1-(2* M+1)*b

def y(N,I,M,b,gammaI, alfa):
	if alfa == 0:
		return 2*M + 1 - b*(N*gammaI-1)

	return ( (2*M + 1 - b*(N*gammaI-1))**2 + 4 * alfa**2 * (I-M) * (I+1+M) )**0.5
def Ifunc(x, tau_lead):
	x = x.real
	threshold = 14
	if x == 0:
		return tau_lead
	
	elif tau_lead == 0:
		if x > 0:
			return 0
		elif x < 0:
			return -x
	
	elif x/tau_lead < -threshold:
		return -x
	elif x/tau_lead > threshold:
		return 0
	return (x / (math.exp(x/tau_lead)-1))



	

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
			return qt.tensor([QObj] + [qt.qeye(2) for _ in xrange(self.parameters.N)])
		self.IdentityElec = qt.qeye(2)
		self.Sm = qt.sigmam()
		self.Sx = 0.5 * qt.sigmax()
		self.Sy = 0.5 * qt.sigmay()
		self.Sz = 0.5 * qt.sigmaz()
		

		def nuclearOperator(QObj, n):
			nucleiList = [qt.qeye(2) for _ in xrange(self.parameters.N)]
			nucleiList[n] = QObj
			return qt.tensor(nucleiList)
		s = 0.5
		self.IdentityNucl = qt.qeye(self.parameters.N + 1)
		self.Ix = qt.jmat((self.parameters.N/2.), 'x')
		self.Iy = qt.jmat((self.parameters.N/2.), 'y')
		self.Iz = qt.jmat((self.parameters.N/2.), 'z')
		self.Isq = (self.Ix**2 + self.Iy**2 + self.Iz**2).tidyup()
		self.Ip = qt.jmat((self.parameters.N/2.), '+')
		self.Im = qt.jmat((self.parameters.N/2.), '-')

		self.S2 = qt.tensor(self.Sx,  self.IdentityNucl) **2 + qt.tensor(self.Sy,  self.IdentityNucl) **2 + qt.tensor(self.Sz,  self.IdentityNucl) **2		
		self.I2 = qt.tensor(self.IdentityElec, self.Ix) **2 + qt.tensor(self.IdentityElec, self.Iy) **2 + qt.tensor(self.IdentityElec, self.Iz) **2
		
		self.H_elec_init =  - self.parameters.N * self.parameters.gammaI_init * self.parameters.b_init * self.Sz
		self.H_nuclei_init = -self.parameters.b_init * self.Iz
		self.H_hf_init = 2 * qt.tensor(self.Sz, self.Iz) + self.parameters.alfa * qt.tensor(self.Sm.dag(),  self.Im) +  self.parameters.alfa * qt.tensor(self.Sm, self.Im.dag() )
		self.H_init = qt.tensor(self.H_elec_init, self.IdentityNucl)  + qt.tensor(self.IdentityElec, self.H_nuclei_init) + self.H_hf_init 

		self.H_elec_final =  - self.parameters.N * self.parameters.gammaI_final * self.parameters.b_final * self.Sz
		self.H_nuclei_final = -self.parameters.b_final * self.Iz
		self.H_hf_final = 2 * qt.tensor(self.Sz, self.Iz) + self.parameters.alfa * qt.tensor(self.Sm.dag(),  self.Im) + self.parameters.alfa * qt.tensor(self.Sm, self.Im.dag() )
		self.H_final = qt.tensor(self.H_elec_final, self.IdentityNucl)  + qt.tensor(self.IdentityElec, self.H_nuclei_final) + self.H_hf_final

		self.M_stg = qt.tensor(self.Sz,  self.IdentityNucl) -  qt.tensor(self.IdentityElec, self.Iz) / self.parameters.N
		self.H_init = qt.Qobj(self.H_init.full())
		self.H_final = qt.Qobj(self.H_final.full())
		
		self.identity = qt.tensor(self.IdentityElec, self.IdentityNucl)
		if HYSTERESIS:	
			self.b_med1 = self.parameters.b_init/10.
			self.b_med2 = -self.parameters.b_init/10.
			self.H_med1 = self.getHamiltonian(self.parameters.N, self.b_med1, self.parameters.gammaI_init, self.parameters.alfa)
			self.H_med2 = self.getHamiltonian(self.parameters.N, self.b_med2, self.parameters.gammaI_init, self.parameters.alfa)
	def getHamiltonian(self, N, b, gammaI, alfa):
		H_elec =  - N * gammaI * b * self.Sz
		H_nuclei = -b * self.Iz
		
		H_hf = 2 * qt.tensor(self.Sz, self.Iz) + alfa * qt.tensor(self.Sm.dag(),  self.Im) + alfa * qt.tensor(self.Sm, self.Im.dag() )
		
		H = qt.Qobj((qt.tensor(H_elec, self.IdentityNucl)  + qt.tensor(self.IdentityElec, H_nuclei) + H_hf ).full())
		return H	
	def getEigenstate(self,ps,I,M,b,gammaI, alfa,p, ham,ns1,ns2 = False):
		
		if ps == 'u':
			state = qt.tensor(qt.basis(2, 0), ns1['state'])
			te = -0.5*self.parameters.N*gammaI*b-I*(b-1)
			rep = '|%s >|%g %g %d>'%(UPARROW, ns1['M'], ns1['I'], ns1['p'])
		elif ps == 'd':
			state = qt.tensor(qt.basis(2, 1), ns1['state'])	
			te = 0.5*self.parameters.N*gammaI*b+I*(b+1)
			rep = '|%s >|%g %g %d>'%(DOWNARROW, ns1['M'], ns1['I'], ns1['p'])	
		elif ps == '+':
			f = phi(self.parameters.N,I,M,b,gammaI, alfa)
			state = math.cos(f/2.) * qt.tensor(qt.basis(2, 0),ns1['state']) + math.sin(f/2.) * qt.tensor(qt.basis(2, 1),ns2['state'])
			te = 0.5* x(M,b)  + 0.5*  y(self.parameters.N,I,M,b,gammaI, alfa)
			rep = '%+.3f|%s >|%g %g %d> %+.3f|%s >|%g %g %d>'%(math.cos(f/2.), UPARROW, ns1['M'], ns1['I'], ns1['p'], math.sin(f/2.), DOWNARROW, ns2['M'], ns2['I'], ns2['p'])	
		elif  ps == '-':
			f = phi(self.parameters.N,I,M,b,gammaI, alfa)
			state = -math.sin(f/2.) * qt.tensor(qt.basis(2, 0),ns1['state']) + math.cos(f/2.) * qt.tensor(qt.basis(2, 1),ns2['state'])
			te = 0.5* x(M,b)  - 0.5* y(self.parameters.N,I,M,b,gammaI, alfa)
		 	rep = '%+.3f|%s >|%g %g %d> %+.3f|%s >|%g %g %d>'%(-math.sin(f/2.), UPARROW, ns1['M'], ns1['I'], ns1['p'], math.cos(f/2.), DOWNARROW, ns2['M'], ns2['I'], ns2['p'])	
		return {'ps' : ps,'I' : I,'M' : M,'p' : p,'state' : qt.Qobj(state.full()).tidyup().unit(), 'theorEnergy': te,'energy' : ham.matrix_element(state.dag(),state).real, 'representation' : rep} 
	
		
	def calculateEigenstates(self,b, gammaI, alfa, ham):
		eigenstates = []
		for ns in self.nucleiStates:
			I, M, p = ns['I'], ns['M'], ns['p']
			if round(M,2) == round(I,2):
				eigenstates.append(self.getEigenstate('u',I,M,b, gammaI, alfa,p, ham,ns))
			
			else:
				for ps in ('-', '+'):
					ns2 = [elem for elem in self.nucleiStates if (elem['I'] == I and elem['p'] == p and elem['M'] == M + 1)][0]
					
					eigenstates.append(self.getEigenstate(ps,I,M,b, gammaI, alfa,p, ham,ns,ns2))
			if M == -I:
				eigenstates.append(self.getEigenstate('d',I,M,b, gammaI, alfa,p, ham,ns))
				
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
		self.eigenstates_init = self.calculateEigenstates(self.parameters.b_init, self.parameters.gammaI_init, self.parameters.alfa, self.H_init)#[{'ps' : 0,'I' : 0,'M' : 0,'p' : 0,'state' : self.H_init.eigenstates()[1][i], 'theorEnergy': self.H_init.eigenstates()[0][i],'energy' : self.H_init.eigenstates()[0][i], 'representation' : 0} for i in xrange(len(self.H_init.eigenstates()[0]))]#
		self.eigenstates_final = self.calculateEigenstates(self.parameters.b_final, self.parameters.gammaI_final, self.parameters.alfa, self.H_final)# the eigenstates must be calculated with no MF gradient	
		self.eigenstates_init = sorted([s for s in self.eigenstates_init], key = lambda s : s['energy'])
		self.eigenstates_final = sorted([s for s in self.eigenstates_final], key = lambda s : s['energy'])		
		if HYSTERESIS:
			self.eigenstates_med1 = self.calculateEigenstates(self.b_med1, self.parameters.gammaI_init, self.parameters.alfa, self.H_med1)
			self.eigenstates_med2 = self.calculateEigenstates(self.b_med2, self.parameters.gammaI_init, self.parameters.alfa, self.H_med2)
			self.eigenstates_med1 = sorted([s for s in self.eigenstates_med1], key = lambda s : s['energy'])
			self.eigenstates_med2 = sorted([s for s in self.eigenstates_med2], key = lambda s : s['energy'])
		if not THRUSTY:
			print 'Check Init'
			self.checkEigenstates(self.eigenstates_init, self.H_init)
			print 'Check _final'
			self.checkEigenstates(self.eigenstates_final, self.H_final)
		 
		
		self.groundStates = [s for s in self.eigenstates_init if round(s['energy'], 15) == round(self.eigenstates_init[0]['energy'], 15)]
		print 'Initial Hamiltonian: b= %f  g= %f'%(self.parameters.b_init, self.parameters.gammaI_init)
		print 'Eigenstates: (energy)(theorEnergy) State |ps M I>', len(self.eigenstates_init)	
		for gs in self.eigenstates_init:
			print '(%.3f)(%.3f) |%s %+g %g> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['M'],gs['I'], gs['representation'])
		print 'Final Hamiltonian: b= %f  g= %f'%(self.parameters.b_final, self.parameters.gammaI_final)
		print 'Eigenstates: (energy)(theorEnergy) State |ps M I>', len(self.eigenstates_final)	
		for gs in self.eigenstates_final:
			print '(%.3f)(%.3f) |%s %+g %g> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['M'],gs['I'], gs['representation'])
		if HYSTERESIS:
			print 'Med 1 Hamiltonian: b= %f  g= %f'%(self.b_med1, self.parameters.gammaI_init)
			print 'Eigenstates: (energy)(theorEnergy) State |ps M I>', len(self.eigenstates_init)	
			for gs in self.eigenstates_med1:
				print '(%.3f)(%.3f) |%s %+g %g> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['M'],gs['I'], gs['representation'])
			print 'Med 2 Hamiltonian: b= %f  g= %f'%(self.b_med2, self.parameters.gammaI_init)
			print 'Eigenstates: (energy)(theorEnergy) State |ps M I>', len(self.eigenstates_init)	
			for gs in self.eigenstates_med2:
				print '(%.3f)(%.3f) |%s %+g %g> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['M'],gs['I'], gs['representation'])
		print 'Initial ground states: (energy)(theorEnergy) State'			
		for gs in self.groundStates:
			print '(%.3f)(%.3f) |%s %g %+g %d> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['I'],gs['M'],gs['p'], gs['representation'])
		if self.parameters.tau_dot > 0:
			dm = (-self.H_init/self.parameters.tau_dot).expm().matrix_element(self.eigenstates_init[0]['state'].dag(), self.eigenstates_init[0]['state']) * qt.ket2dm(self.eigenstates_init[0]['state'])
			for es in self.eigenstates_init[1:]:
				dm += (-self.H_init/self.parameters.tau_dot).expm().matrix_element(es['state'].dag(), es['state']) * qt.ket2dm(es['state'])
				

		else:
			dm = qt.ket2dm(self.groundStates[0]['state'])
			for gs in self.groundStates[1:]:
				dm += qt.ket2dm(gs['state'])
				
		return dm/dm.tr()
	

		
	def initialize_dissipators(self,states, mean_energy, Ham):	
		def addDissipator(state1, state2, Ham):
			dissipators12 = [] 
			initial = state1['state']
			final = state2['state']
			overlap = math.fabs( sum([initial.overlap(intermediate(ns['state']))  *  final.overlap(intermediate(ns['state'])) for ns in self.nucleiStates]) )**2
			if overlap > 10**(-15):
				ifunction = Ifunc(Ham.matrix_element(final.dag(), final) - Ham.matrix_element(initial.dag(), initial), self.parameters.tau_lead)
				rate = baseRate * (overlap * ifunction)**0.5
				if rate >0:
					print 'Transition |%s %+g %g> --> |%s %+g %g> -- rate %.3f * %.3f * %.3f -- Delta E %.3f'%(state1['ps'], state1['M'],state1['I'],state2['ps'], 
							state2['M'],state2['I'],baseRate ,overlap , ifunction, (Ham.matrix_element(final.dag(), final) - Ham.matrix_element(initial.dag(), initial)).real)
					
					dissipators12.append(rate * final * initial.dag())
				if state1 != state2:
					ifunction = Ifunc(Ham.matrix_element(initial.dag(), initial) - Ham.matrix_element(final.dag(), final), self.parameters.tau_lead)
					rate = baseRate * (overlap * ifunction)**0.5
					if rate >0:
						print 'Transition |%s %+g %g> --> |%s %+g %g> -- rate %.3f * %.3f * %.3f -- Delta E %.3f'%(state2['ps'], state2['M'],state2['I'],state1['ps'], 
						state1['M'],state1['I'],baseRate ,overlap , ifunction, (Ham.matrix_element(initial.dag(), initial) - Ham.matrix_element(final.dag(), final)).real)
						dissipators12.append(rate * initial * final.dag())
			
			return dissipators12
		baseRate = (self.parameters.eta**2/ (2*math.pi*(self.parameters.mu-mean_energy)**2))**0.5
		intermediate = lambda ns : qt.tensor(qt.basis(2,0), ns) +  qt.tensor(qt.basis(2,1), ns)
		dissipators = []
		for i in xrange(len(states)):
			for j in xrange(i, len(states)):
				state1, state2 = (states[i], states[j])
				#this is just to speed things up. only transitions allowed conserve I and p
				if not THRUSTY or (state1['I'] == state2['I'] and state1['p']== state2['p'] and -1 <= state1['M'] - state2['M'] <= 1):
					dissipators += addDissipator(state1, state2, Ham)
		
		return dissipators
	
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
			operators.append(self.H_final)
		if self.master.I2Bool.get():
			operators.append(self.I2)
		if self.master.S2Bool.get():
			operators.append(self.S2)
		if self.master.CohBool.get():
			operators.append(qt.tensor(self.IdentityElec, self.Ip))
			operators.append(qt.tensor(self.IdentityElec, self.Ip * self.Im))
		if self.master.IzNormBool.get():
			operators.append(qt.tensor(self.IdentityElec, self.Iz)/(self.parameters.N/2.))
		return operators

	def compute_phase_diagram(self):
		
		N = 20
		
		xmin, xmax = min(0,self.parameters.b_init,self.parameters.b_final), max(self.parameters.b_init,self.parameters.b_final)+1
		ymin, ymax = 0, max(self.parameters.gammaI_init,self.parameters.gammaI_final)+1#min(self.parameters.gammaI_init,self.parameters.gammaI_final)
		
		
		x_values = np.linspace(xmin, xmax, num=N, endpoint=True)
		y_values = np.linspace( ymax,ymin, num=N, endpoint=True)
		pdbuffer = []
		T = max(self.parameters.tau_dot, 0.1)
		iz = qt.Qobj(qt.tensor(self.IdentityElec, self.Iz).full())
		for b in x_values:
			for gammaI in y_values:
				h = self.getHamiltonian(self.parameters.N, b, gammaI , self.parameters.alfa)	
				c, d = h.eigenstates()
				ens = zip(c, d)
				pdbuffer.append(sum( [math.exp(-s[0]/T) * qt.expect(iz ,s[1]) for s in ens] )/sum( [math.exp(-s[0]/T) for s in ens] ))
				


		pd = np.ndarray( shape=(N,N), dtype=float,  buffer=np.array(pdbuffer) )
		print pd
		return (y_values,x_values,  pd.transpose() )


	def compute_dynamics(self):
		def L_func(t, args):
			if t <= self.parameters.evo_time:
				if t < self.parameters.dis_time:
					return self.L_init.data
				elif t >= self.parameters.dis_time:
					return self.L_init.data + self.L_dissipators_init.data
			elif t > self.parameters.evo_time:
				if t < self.parameters.dis_time:
					return self.L_final.data
				elif t >= self.parameters.dis_time:
					return self.L_final.data + self.L_dissipators_final.data
		def L_cycle(t, args):
			for i, tist in enumerate(split_tlist):
				if np.amin(tist) <= t <= np.amax(tist):
					if i%2:
						return self.L_init.data + self.L_dissipators_init.data
					else:
						return self.L_final.data + self.L_dissipators_final.data
			else:
				
				return self.L_final.data + self.L_dissipators_final.data
		def L_hysteresis(t, args):
			#print t
			for i, tist in enumerate(split_tlist):
				if np.amin(tist) <= t <= np.amax(tist):
					
					if i%8==0 or i%8==7:
						return self.L_init.data + self.L_dissipators_init.data
					elif i%8==1 or i%8==6:
						return self.L_med1.data + self.L_dissipators_med1.data
					elif i%8==2 or i%8==5:
						return self.L_med2.data + self.L_dissipators_med2.data
					elif i%8==3 or i%8==4:
						return self.L_final.data + self.L_dissipators_final.data
					
			#terrible hack because of unespected behaviour of mesolve, the times used are larger than the list given			
			
			else:		
				if t < np.amax(self.tlist)+1:
					return self.L_init.data + self.L_dissipators_init.data	
				else:
					print 'extreme problem  ', t#, split_tlist
					return self.L_final.data + self.L_dissipators_final.data
			
		self.operators_compute = self.get_operators_compute()
		if HYSTERESIS:
			split_tlist = np.array_split(self.tlist, 8)
		
			try:
				self.expt_list = qt.mesolve(L_hysteresis, self.density_matrix, self.tlist, [], self.operators_compute, progress_bar = True) 
				
			except Exception as e:
				print e
				self.expt_list = 'Error'
		else:
			split_tlist = np.array_split(self.tlist, self.parameters.cycles)
		
			try:
				self.expt_list = qt.mesolve(L_cycle, self.density_matrix, self.tlist, [], self.operators_compute, progress_bar = True) 
				with open('results.txt', 'w+') as h:	
					for r in list(self.expt_list.expect[0]):
						h.write('%f \n'%r)	
			except Exception as e:
				print e
				self.expt_list = 'Error'
		
		
		return self.expt_list
	def compute_density_matrix(self):
		def L_func(t, args):
			v = args[0]
			if 0 <= v * t <=1:
				return self.L_init.data * (1 - v * t) + self.L_final.data * v * t
			else:
				return self.L_final.data
		
		
		print 'Calculating density matrix'
		a, b = self.H_init.eigenstates()
		c, d = self.H_final.eigenstates()
		l = b + d
		
		#self.expt_density_matrix_list = qt.mesolve(L_func, self.density_matrix, self.tlist, [], [], progress_bar = True, args=[self.parameters.quench_v]) 
		#pgt = [[qt.expect(rho, l[i]) for rho in self.expt_density_matrix_list.states] for i in xrange(len(l))]
		#return (self.tlist, pgt)

		dtlist = np.linspace(0.00001, 1./self.parameters.quench_v, num=min(80, self.parameters.tSteps))
		pgt = []
		for dt in dtlist:
			self.expt_density_matrix_list = qt.mesolve(L_func, self.density_matrix, dtlist, [], [], args=[1./dt]) 
			rho = self.expt_density_matrix_list.states[-1]
			pgt.append([qt.expect(rho, d[i]) for i in xrange(len(d))])
		return (dtlist, zip(*pgt))



	def compute_energy_levels(self):
		def H(t):		
			if 0 <= self.parameters.quench_v * t <=1:
				return self.H_init * (1 - self.parameters.quench_v * t) + self.H_final * self.parameters.quench_v * t
			else:
				return self.H_final
		print 'Calculating energy levels'
		energy_levels = zip(*[[a for a in H(t).eigenstates()[0]] for t in self.tlist])
		
		return  (self.tlist, energy_levels)
	def prepare_computation(self,parameters):		
		startTime = time.time()
		
		self.parameters = parameters
		p = self.parameters
		self.nucleiStates = [{'I' : p.N/2.,'M' : - i + p.N/2.,'p' : 1,'state' : qt.basis(p.N + 1, i)} for i in xrange(p.N + 1)]#decompose(p.N, p.N/2.)
		self.nucleiStates = sorted([el for el in self.nucleiStates], key = lambda el: (el['I'], el['M']))
		if not THRUSTY:
			represent(self.nucleiStates)
		self.defineOperators()
		self.density_matrix = self.defineStates()
		self.mean_energy_init = sum([es['energy'] for es in self.eigenstates_init])/len(self.eigenstates_init)
		self.mean_energy_final = sum([es['energy'] for es in self.eigenstates_final])/len(self.eigenstates_final)
		statesTime = time.time()
		print 'States and operators got in %f'%(statesTime-startTime)

		if p.eta > 0 and p.dissipators_kind.get() == 'fermi':
			print 'Transition (Fermi golden rule) |ps M I> --> |ps M I> -- rate [baseRate * overlap * Ifunction] Delta E'
			print 'Final dissipators'
			dissipators_list_final = self.initialize_dissipators(self.eigenstates_final, self.mean_energy_final, self.H_final)
			print 'Initial dissipators'
			dissipators_list_init = self.initialize_dissipators(self.eigenstates_init, self.mean_energy_init, self.H_init)
			
		elif p.eta > 0 and p.dissipators_kind.get() == 'central':	
			baserate = (p.eta**2/ (2*math.pi*(p.mu-self.mean_energy_final)**2))**0.5
			
			rate_plus = baserate * Ifunc(-p.b_final * p.N * p.gammaI_final, p.tau_lead)**0.5
			rate_minus = baserate * Ifunc(p.b_final * p.N * p.gammaI_final, p.tau_lead)**0.5
			dissipators_list_final = [rate_plus * qt.tensor(self.Sm.dag(), self.IdentityNucl), rate_minus * qt.tensor(self.Sm, self.IdentityNucl)]
			print 'Transition (central spin) final rate plus = %g, rate minus = %g'%(rate_plus, rate_minus)
			baserate = (p.eta**2/ (2*math.pi*(p.mu-self.mean_energy_init)**2))**0.5
			rate_plus = baserate * Ifunc(-p.b_init * p.N * p.gammaI_init, p.tau_lead)**0.5
			rate_minus = baserate * Ifunc(p.b_init * p.N * p.gammaI_init, p.tau_lead)**0.5
			dissipators_list_init = [rate_plus * qt.tensor(self.Sm.dag(), self.IdentityNucl), rate_minus * qt.tensor(self.Sm, self.IdentityNucl)]
			print 'Transition (central spin) initial rate plus = %g, rate minus = %g'%(rate_plus, rate_minus)
			
				
		else:
			dissipators_list_init = []
			dissipators_list_final = []
			dissipators_list_med1 = []
			dissipators_list_med2 = []
			print 'No dissipation'
		
		dis_time = time.time()
		print 'Dissipators (%d) got in %f'%(len(dissipators_list_final), dis_time-statesTime)

		self.L_init = qt.Qobj(qt.liouvillian(self.H_init, []).full())
		self.L_final = qt.Qobj(qt.liouvillian(self.H_final, []).full())
		self.L_dissipators_init = qt.Qobj(qt.liouvillian(self.identity, dissipators_list_init).full())
		self.L_dissipators_final = qt.Qobj(qt.liouvillian(self.identity, dissipators_list_final).full())
		
		if HYSTERESIS:
			self.mean_energy_med1 = sum([es['energy'] for es in self.eigenstates_med1])/len(self.eigenstates_med1)
			self.mean_energy_med2 = sum([es['energy'] for es in self.eigenstates_med2])/len(self.eigenstates_med2)
			if p.eta > 0 and p.dissipators_kind.get() == 'fermi':
				print 'Med 1 dissipators'
				dissipators_list_med1 = self.initialize_dissipators(self.eigenstates_med1, self.mean_energy_med1, self.H_med1)
				print 'Med 2 dissipators'
				dissipators_list_med2 = self.initialize_dissipators(self.eigenstates_med2, self.mean_energy_med2, self.H_med2)
			
			elif p.eta > 0 and p.dissipators_kind.get() == 'central':	
				baserate = (p.eta**2/ (2*math.pi*(p.mu-self.mean_energy_med1)**2))**0.5
			
				rate_plus = baserate * Ifunc(-self.b_med1 * p.N * p.gammaI_init, p.tau_lead)**0.5
				rate_minus = baserate * Ifunc(self.b_med1 * p.N * p.gammaI_init, p.tau_lead)**0.5
				dissipators_list_med1 = [rate_plus * qt.tensor(self.Sm.dag(), self.IdentityNucl), rate_minus * qt.tensor(self.Sm, self.IdentityNucl)]
				print 'Transition (central spin) Med 1 rate plus = %g, rate minus = %g'%(rate_plus, rate_minus)
				baserate = (p.eta**2/ (2*math.pi*(p.mu-self.mean_energy_med2)**2))**0.5
				rate_plus = baserate * Ifunc(-self.b_med2 * p.N * p.gammaI_init, p.tau_lead)**0.5
				rate_minus = baserate * Ifunc(self.b_med2 * p.N * p.gammaI_init, p.tau_lead)**0.5
				dissipators_list_med2 = [rate_plus * qt.tensor(self.Sm.dag(), self.IdentityNucl), rate_minus * qt.tensor(self.Sm, self.IdentityNucl)]
				print 'Transition (central spin) Med 2 rate plus = %g, rate minus = %g'%(rate_plus, rate_minus)
			else:
				dissipators_list_med1 = []
				dissipators_list_med2 = []
			
			self.L_med1 = qt.Qobj(qt.liouvillian(self.H_med1, []).full())
			self.L_med2 = qt.Qobj(qt.liouvillian(self.H_med2, []).full())
			self.L_dissipators_med1 = qt.Qobj(qt.liouvillian(self.identity, dissipators_list_med1).full())
			self.L_dissipators_med2 = qt.Qobj(qt.liouvillian(self.identity, dissipators_list_med2).full())	

		
		liouTime = time.time()
		print 'Liouvillians got in %f'%(liouTime-dis_time)		
		opts = qt.Odeoptions(nsteps = 10000000000, rhs_reuse=True, num_cpus = 4, rtol = 0.001, atol = 0.00001)
		self.tlist = np.linspace(0,p.tInterval,num=p.tSteps, endpoint=True)
		
	
			#self.expt_list = qt.mesolve([self.H_init, [(self.H_final- self.H_init), H_final_coeff]] , self.density_matrix, tlist, collapseList, self.operators_compute, options = opts)
			#print 'States and operators got in %f'%(statesTime-startTime)
			#print 'Dissipators (%d) got in %f'%(len(dissipators_list_final), dis_time-statesTime)			
			#print 'Liouvillians got in %f'%(time.time()-dis_time)				
			#print 'Finished in %f'%(time.time()-startTime)	
		

