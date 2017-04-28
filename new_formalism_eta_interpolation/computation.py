# -*- coding: UTF-8 -*-

import numpy as np
import math, sys, time, cmath, os
import multiprocessing as mp
from spin_states import *
#from superoperator import *
try:
    import qutip as qt
    
except:
    print ' Some modules are missing, try installing them using \'installModules.py\''
    sys.exit()
THRUSTY = 0
NPROCS = 1
UPARROW = u'\u2191'
DOWNARROW = u'\u2193'



def HST(x):
    if x==0:
        return 0.5 
    else:
        return int(x>0)

def phi(N,I,M,b,gamma, alfa):
    return cmath.phase( 2*M + 1 + b - N*gamma + 2* alfa * 1j * ( I*(I+1) - M*(M+1) )**0.5 )

def x(M,b):
    return -1-(2* M+1)*b

def y(N,I,M,b,gamma, alfa):
    if alfa == 0:
        return 2*M + 1 + b - N*gamma

    return ( (2*M + 1 + b - N*gamma)**2 + 4 * alfa**2 * (I-M) * (I+1+M) )**0.5
    
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
    def __init__(self, master = None):
        self.master = master
    
    def define_operators(self):
        self.identity_central = qt.qeye(2)
        self.Sm = qt.sigmam()
        self.Sp = qt.sigmap()
        self.Sx = 0.5 * qt.sigmax()
        self.Sy = 0.5 * qt.sigmay()
        self.Sz = 0.5 * qt.sigmaz()
        
        self.identity_ancilla = qt.qeye(self.parameters.N + 1)
        self.Ix = qt.jmat((self.parameters.N/2.), 'x')
        self.Iy = qt.jmat((self.parameters.N/2.), 'y')
        self.Iz = qt.jmat((self.parameters.N/2.), 'z')
        self.Isq = (self.Ix**2 + self.Iy**2 + self.Iz**2).tidyup()
        self.Ip = qt.jmat((self.parameters.N/2.), '+')
        self.Im = qt.jmat((self.parameters.N/2.), '-')
        self.m_dm = lambda m : qt.ket2dm(qt.basis(self.parameters.N + 1, m))
        
        self.S2 = qt.tensor(self.Sx,  self.identity_ancilla) **2 + qt.tensor(self.Sy,  self.identity_ancilla) **2 + qt.tensor(self.Sz,  self.identity_ancilla) **2        
        self.I2 = qt.tensor(self.identity_central, self.Ix) **2 + qt.tensor(self.identity_central, self.Iy) **2 + qt.tensor(self.identity_central, self.Iz) **2
        

        self.M_stg = qt.tensor(self.Sz, self.identity_ancilla)- 1/(2*self.parameters.N/2.) * qt.tensor(self.identity_central, self.Iz)
        self.identity = qt.tensor(self.identity_central, self.identity_ancilla)
        
    def get_hamiltonian(self, N, b, gamma, alfa):
        H_elec =  - N * gamma * self.Sz
        H_nuclei = -b * self.Iz
        
        H_hf = 2 * qt.tensor(self.Sz, self.Iz) + alfa * qt.tensor(self.Sp,  self.Im) + alfa * qt.tensor(self.Sm, self.Ip )
        
        H = qt.Qobj((qt.tensor(H_elec, self.identity_ancilla)  + qt.tensor(self.identity_central, H_nuclei) + H_hf ).full())
        return H    


    
    def calculate_eigenstates(self,b, gamma, alfa, ham):
        def get_eigenstate(ps,I,M,b,gamma, alfa,p, ham,ns1,ns2 = False):
            if ps == 'u':
                state = qt.tensor(qt.basis(2, 0), ns1['state'])
                te = -0.5*self.parameters.N*gamma-I*(b-1)
                rep = '|%s >|%g %g %d>'%(UPARROW, ns1['M'], ns1['I'], ns1['p'])
            elif ps == 'd':
                state = qt.tensor(qt.basis(2, 1), ns1['state'])    
                te = 0.5*self.parameters.N*gamma+I*(b+1)
                rep = '|%s >|%g %g %d>'%(DOWNARROW, ns1['M'], ns1['I'], ns1['p'])    
            elif ps == '+':
                f = phi(self.parameters.N,I,M,b,gamma, alfa)
                state = math.cos(f/2.) * qt.tensor(qt.basis(2, 0),ns1['state']) + math.sin(f/2.) * qt.tensor(qt.basis(2, 1),ns2['state'])
                te = 0.5* x(M,b)  + 0.5*  y(self.parameters.N,I,M,b,gamma, alfa)
                rep = '%+.3f|%s >|%g %g %d> %+.3f|%s >|%g %g %d>'%(math.cos(f/2.), UPARROW, ns1['M'], ns1['I'], ns1['p'], math.sin(f/2.), DOWNARROW, ns2['M'], ns2['I'], ns2['p'])    
            elif  ps == '-':
                f = phi(self.parameters.N,I,M,b,gamma, alfa)
                state = -math.sin(f/2.) * qt.tensor(qt.basis(2, 0),ns1['state']) + math.cos(f/2.) * qt.tensor(qt.basis(2, 1),ns2['state'])
                te = 0.5* x(M,b)  - 0.5* y(self.parameters.N,I,M,b,gamma, alfa)
                rep = '%+.3f|%s >|%g %g %d> %+.3f|%s >|%g %g %d>'%(-math.sin(f/2.), UPARROW, ns1['M'], ns1['I'], ns1['p'], math.cos(f/2.), DOWNARROW, ns2['M'], ns2['I'], ns2['p'])    
            return {'ps' : ps,'I' : I,'M' : M,'p' : p,'state' : qt.Qobj(state.full()).tidyup().unit(), 'theorEnergy': te,'energy' : ham.matrix_element(state.dag(),state).real, 'representation' : rep} 
    
        eigenstates = []
        for ns in self.nucleiStates:
            I, M, p = ns['I'], ns['M'], ns['p']
            if round(M,2) == round(I,2):
                eigenstates.append(get_eigenstate('u',I,M,b, gamma, alfa,p, ham,ns))
            
            else:
                for ps in ('-', '+'):
                    ns2 = [elem for elem in self.nucleiStates if (elem['I'] == I and elem['p'] == p and elem['M'] == M + 1)][0]
                    
                    eigenstates.append(get_eigenstate(ps,I,M,b, gamma, alfa,p, ham,ns,ns2))
            if M == -I:
                eigenstates.append(get_eigenstate('d',I,M,b, gamma, alfa,p, ham,ns))
                
        eigenstates = sorted([s for s in eigenstates], key = lambda s : s['energy'])
        return     eigenstates
    

    
    
    def get_dm(self, eigenstates):    
        self.groundStates = [s for s in eigenstates if round(s['energy'], 15) == round(eigenstates[0]['energy'], 15)]
        h_init = self.get_hamiltonian(self.parameters.N, self.parameters.b_init, self.parameters.gamma_init, self.parameters.alfa)
        if self.parameters.tau_dot > 0:
            dm = (-h_init/self.parameters.tau_dot).expm().matrix_element(eigenstates[0]['state'].dag(), eigenstates[0]['state']) * qt.ket2dm(eigenstates[0]['state'])
            for es in eigenstates[1:]:
                dm += (-h_init/self.parameters.tau_dot).expm().matrix_element(es['state'].dag(), es['state']) * qt.ket2dm(es['state'])
                

        else:
            dm = qt.ket2dm(self.groundStates[0]['state'])
            for gs in self.groundStates[1:]:
                dm += qt.ket2dm(gs['state'])
                
        return dm/dm.tr()
    

        
    def get_fermi_dissipators_full(self,states, Ham):    
        def addDissipator(state1, state2, Ham):
            dissipators12 = [] 
            initial = state1['state']
            final = state2['state']
            overlap = math.fabs( sum([initial.overlap(intermediate(ns['state']))  *  final.overlap(intermediate(ns['state'])) for ns in self.nucleiStates]) )**2
            ifunction = Ifunc(Ham.matrix_element(final.dag(), final) - Ham.matrix_element(initial.dag(), initial), self.parameters.tau_lead)
            rate = (baseRate * overlap * ifunction)**0.5
            if rate >0:
                print 'Transition |%s %+g %g> --> |%s %+g %g> -- rate %.3f * %.3f * %.3f -- Delta E %.3f'%(state1['ps'], state1['M'],state1['I'],state2['ps'], 
                        state2['M'],state2['I'],baseRate ,overlap , ifunction, (Ham.matrix_element(final.dag(), final) - Ham.matrix_element(initial.dag(), initial)).real)
                dissipators12.append(rate * final * initial.dag())
            if state1 != state2:
                ifunction = Ifunc(Ham.matrix_element(initial.dag(), initial) - Ham.matrix_element(final.dag(), final), self.parameters.tau_lead)
                rate = (baseRate * overlap * ifunction)**0.5
                if rate >0:
                    print 'Transition |%s %+g %g> --> |%s %+g %g> -- rate %.3f * %.3f * %.3f -- Delta E %.3f'%(state2['ps'], state2['M'],state2['I'],state1['ps'], 
                    state1['M'],state1['I'],baseRate ,overlap , ifunction, (Ham.matrix_element(initial.dag(), initial) - Ham.matrix_element(final.dag(), final)).real)
                    dissipators12.append(rate * initial * final.dag())
            
            return dissipators12
        baseRate = self.parameters.eta
        intermediate = lambda ns : qt.tensor(qt.basis(2,0), ns) +  qt.tensor(qt.basis(2,1), ns)
        dissipators = []
        for i in xrange(len(states)):
            for j in xrange(i, len(states)):
                state1, state2 = (states[i], states[j])
                #this is just to speed things up. only transitions allowed conserve I and p
                if not THRUSTY or (state1['I'] == state2['I'] and state1['p']== state2['p'] and -1 <= state1['M'] - state2['M'] <= 1):
                    dissipators += addDissipator(state1, state2, Ham)
        
        return dissipators
    def get_fermi_dissipators(self,states, Ham, baseRate):    
        def addDissipator(state1, state2, Ham):
            dissipators12 = [] 
            initial = state1['state']
            final = state2['state']
            overlap = math.fabs( initial.overlap(qt.Qobj(qt.tensor(self.Sp,  self.identity_ancilla).full()) * final ) )**2 + math.fabs( initial.overlap(qt.Qobj(qt.tensor(self.Sm,  self.identity_ancilla).full()) * final ) )**2
            ifunction = Ifunc(Ham.matrix_element(final.dag(), final) - Ham.matrix_element(initial.dag(), initial), self.parameters.tau_lead)
            rate = (baseRate * overlap * ifunction)**0.5
            if rate >0:
                print 'Transition |%s %+g %g> --> |%s %+g %g> -- rate %.3f * %.3f * %.3f -- Delta E %.3f'%(state1['ps'], state1['M'],state1['I'],state2['ps'], 
                        state2['M'],state2['I'],baseRate ,overlap , ifunction, (Ham.matrix_element(final.dag(), final) - Ham.matrix_element(initial.dag(), initial)).real)
                dissipators12.append(rate * final * initial.dag())
            if state1 != state2:
                ifunction = Ifunc(Ham.matrix_element(initial.dag(), initial) - Ham.matrix_element(final.dag(), final), self.parameters.tau_lead)
                rate = (baseRate * overlap * ifunction)**0.5
                if rate >0:
                    print 'Transition |%s %+g %g> --> |%s %+g %g> -- rate %.3f * %.3f * %.3f -- Delta E %.3f'%(state2['ps'], state2['M'],state2['I'],state1['ps'], 
                    state1['M'],state1['I'],baseRate ,overlap , ifunction, (Ham.matrix_element(initial.dag(), initial) - Ham.matrix_element(final.dag(), final)).real)
                    dissipators12.append(rate * initial * final.dag())
            
            return dissipators12

        
        dissipators = []
        for i in xrange(len(states)):
            for j in xrange(i, len(states)):
                state1, state2 = (states[i], states[j])
                #this is just to speed things up. only transitions allowed conserve I and p
                if not THRUSTY or (state1['I'] == state2['I'] and state1['p']== state2['p'] and -1 <= state1['M'] - state2['M'] <= 1):
                    dissipators += addDissipator(state1, state2, Ham)
        
        return dissipators


    def get_central_dissipators(self,  N, b, gamma, baseRate):
        rate_plus = (baseRate * Ifunc(- N * gamma, self.parameters.tau_dot) )**0.5
        rate_minus = (baseRate * Ifunc(N * gamma, self.parameters.tau_dot) )**0.5
        print 'r+={0} r-={1}'.format(rate_plus, rate_minus)
        dissipators = [rate_plus * qt.tensor(self.Sp, self.identity_ancilla), rate_minus * qt.tensor(self.Sm, self.identity_ancilla)]
        return dissipators
    
    def get_operators_compute(self):
        operators = []
        if self.master.IzBool.get():
            operators.append(qt.tensor(self.identity_central, self.Iz))
        if self.master.SzBool.get():
            operators.append(qt.tensor(self.Sz, self.identity_ancilla))
            
        if self.master.M_stg_bool.get():
            operators.append(self.M_stg)
        
        if self.master.I2Bool.get():
            operators.append(self.I2)
        if self.master.S2Bool.get():
            operators.append(self.S2)
        if self.master.CohBool.get():
            operators.append(qt.tensor(self.identity_central, self.Ip))
            operators.append(qt.tensor(self.identity_central, self.Ip * self.Im))
        if self.master.IzNormBool.get():
            operators.append(qt.tensor(self.identity_central, self.Iz)/(self.parameters.N/2.))
        if self.master.pBool.get():
            operators.append(0.5 * (1- qt.tensor(self.identity_central, self.Iz)/(self.parameters.N/2.)))
        return operators

    
    def compute_dynamics(self):
        def L_eta(t, args):
            #print t
            return self.L0.data + self.L1.data
            
            
        #self.operators_compute = self.get_operators_compute()
        split_tlist = np.array_split(self.tlist, 2)
        self.mean_time = 0
        spin_flip = qt.tensor(qt.sigmay(), qt.sigmay())
        p = 0.5 * (1- qt.tensor(self.identity_central, self.Iz)/(self.parameters.N/2.))
        ETA_POINTS = 10
        
        #for i in xrange(ETA_POINTS):
            #baserate = (0.0 + self.parameters.eta*i/(ETA_POINTS-1))
        for baserate in [20]:#[x/10. for x in xrange(0, 10, 3)] + [x for x in xrange(10, 50, 10)]:
            rate1 = baserate
            rate2 = baserate
            g = self.parameters.gamma_final
            if self.parameters.dissipators_kind == 'fermi':
                all_dissipators = self.get_fermi_dissipators(self.all_eigenstates[-1], self.all_hamiltonians[-1], baserate**0.5)
            elif self.parameters.dissipators_kind == 'central':
                all_dissipators = []
                for n in xrange(self.parameters.N + 1):
                    m = self.parameters.N/2. - n
                    Delta_e = g - 2 * m / self.parameters.N#Eup,m - Edown, m
                    rp = baserate**0.5 * HST(Delta_e)
                    rm = baserate**0.5 * HST(-Delta_e)
                    #rpu = baserate**0.5 * HST(g - 1)
                    #rpd = baserate**0.5 * HST(g + 1)
                    #rmu = baserate**0.5 * HST(-g+ 1)
                    #rmd = baserate**0.5 * HST(-g- 1)
                    all_dissipators += [rp * qt.tensor(self.Sp, self.m_dm(n)), rm * qt.tensor(self.Sm, self.m_dm(n))]
                                   
                                   
                                   
                                   
                #print 'r+u={} r+d={}'.format(rpu, rpd)
                #print 'r-u={} r-d={}'.format(rmu, rmd)
            self.L1 = qt.Qobj(qt.liouvillian(self.identity, all_dissipators).full())
            operators_compute = [qt.Qobj(qt.tensor(self.Sz, self.identity_ancilla).full()), qt.Qobj(qt.tensor(self.identity_central, self.Iz).full())]     
            
            try:
                opts = qt.Odeoptions(nsteps = 10000000000, rhs_reuse=True, num_cpus = 4, rtol = 0.001, atol = 0.00001, store_states = True)
                self.expt_list = qt.mesolve(L_eta, self.density_matrix, self.tlist, [], operators_compute, progress_bar = True, options=opts)
                #print  self.expt_list.states, len(self.expt_list.states), len(self.expt_list.times)
                #raw_input()
                #self.mean_time = sum([self.expt_list.times[i] * self.expt_list.expect[-1][i] for i in xrange(len(self.expt_list.times))])/sum(self.expt_list.expect[-1])
                #with open('archive/central/mean_times_N=%d'%(int(self.parameters.N)), 'a') as f:
                #    f.write('%f %f %f'%(self.mean_time, rate1, rate2 ) + '\n') 
                
            except Exception as e:
                print e
                self.expt_list = 'Error'
                
            folder = '_'.join(['{}={}'.format(k, v) for k, v in self.parameters.__class__.__dict__.iteritems() if (not k.startswith('_') and k != 'get_parameters')])
        
            print folder
            if not os.path.exists('archive/double_{}'.format(folder)):
                os.makedirs('archive/double_{}'.format(folder))
            for j, r in enumerate(self.expt_list.expect): 
                with open('archive/double_{}/expectation_{}'.format(folder, j, ), 'w') as f:
                    for t, res in zip(self.expt_list.times, r):
                        f.write('{} {}'.format(t, res) + '\n') 
                
            #with open('archive/concurrence/eta={}'.format(rate1), 'a') as f:
            #    for i in xrange(len(self.expt_list.times)):
                    #concurrence = self.expt_list.states[i] * spin_flip * self.expt_list.states[i].conj() * spin_flip
                    
            #        conc1 = qt.concurrence(self.expt_list.states[i].tidyup())#its a little sloppy here as it takes the abs value of eigenvalue, still i have problems for negative dm!
                     
                    
                    
            #        for ee in self.expt_list.states[i].tidyup().eigenenergies():
            #            if isinstance(ee, complex) or ee.real < -0.001:
            #                print self.expt_list.states[i].eigenenergies(), self.expt_list.times[i], baserate
            #                print self.expt_list.states[i]
                            #raw_input()
            #        f.write('{} {} {}'.format(rate1, self.expt_list.times[i], conc1) + '\n')  
                    #for l in concurrence.eigenenergies():
                    #    f.write('{} '.format(l)) 
                    #f.write('\n')  
                

        
        
        
        return (self.expt_list,  self.mean_time)
    
    def prepare_computation(self,parameters):    
        startTime = time.time()
        
        self.parameters = parameters
        p = self.parameters
        self.nucleiStates = [{'I' : p.N/2.,'M' : - i + p.N/2.,'p' : 1,'state' : qt.basis(p.N + 1, i)} for i in xrange(p.N + 1)]
        self.nucleiStates = sorted([el for el in self.nucleiStates], key = lambda el: (el['I'], el['M']))
        if not THRUSTY:
            represent(self.nucleiStates)

        self.define_operators()
        self.all_b = [p.b_init, p.b_final]
        self.all_gamma = [p.gamma_init, p.gamma_final]
        self.all_hamiltonians = [self.get_hamiltonian(p.N, self.all_b[i], self.all_gamma[i], p.alfa) for i in xrange(2)]
        self.all_eigenstates = [self.calculate_eigenstates(self.all_b[i], self.all_gamma[i], p.alfa, self.all_hamiltonians[i]) for i in xrange(2)]
        
        
        

        self.density_matrix = qt.ket2dm(qt.tensor(qt.basis(2, 1), qt.basis(p.N + 1, p.N)))#I impose to start in the fully polarized state. self.get_dm(self.all_eigenstates[0])


        print 'Initial Hamiltonian: b= %f  g= %f N=%d'%(self.parameters.b_init, self.parameters.gamma_init, self.parameters.N)
        print 'Eigenstates: (energy)(theorEnergy) State |ps M I>', len(self.all_eigenstates[0])    
        for gs in self.all_eigenstates[0]:
            print '(%.3f)(%.3f) |%s %+g %g> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['M'],gs['I'], gs['representation'])
        
        print 'Final Hamiltonian: b= %f  g= %f N=%d'%(self.parameters.b_final, self.parameters.gamma_final, self.parameters.N)
        print 'Eigenstates: (energy)(theorEnergy) State |ps M I>', len(self.all_eigenstates[-1])    
        for gs in self.all_eigenstates[-1]:
            print '(%.3f)(%.3f) |%s %+g %g> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['M'],gs['I'], gs['representation'])

        statesTime = time.time()
        print 'States and operators got in %f'%(statesTime-startTime)
        

        self.L0 = qt.Qobj(qt.liouvillian(self.all_hamiltonians[-1], []).full())
        #self.all_L_dissipators = [qt.Qobj(qt.liouvillian(self.identity, all_dissipators[i]).full()) for i in xrange(2)]
        liouTime = time.time()
        print 'Liouvillians got in %f'%(liouTime-statesTime)        
        
        self.tlist = np.linspace(0,p.points * p.wait_time,num=p.tSteps, endpoint=True)# np.linspace(0,p.tInterval,num=p.tSteps, endpoint=True)
    
