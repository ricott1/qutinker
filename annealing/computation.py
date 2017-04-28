# -*- coding: UTF-8 -*-

import numpy as np
np.set_printoptions(threshold=np.nan)
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
    
    def get_ancilla_states(self, N):
        i = int(N)
        while i >= 0:
            state = {'I' : N/2.,'M' : - i + N/2.,'p' : 1,'state' : qt.basis(N + 1, i)} 
            
            yield state
            i -= 1
    def define_operators(self):
        self.identity_central = qt.qeye(2)
        self.Sm = qt.sigmam()
        self.Sp = qt.sigmap()
        self.Sx = 0.5 * qt.sigmax()
        self.Sy = 0.5 * qt.sigmay()
        self.Sz = 0.5 * qt.sigmaz()
        
        self.identity_ancilla = lambda N : qt.qeye(N + 1)
        self.Ix = lambda N : qt.jmat((N/2.), 'x')
        self.Iy = lambda N : qt.jmat((N/2.), 'y')
        self.Iz = lambda N : qt.jmat((N/2.), 'z')
        self.Isq = lambda N : (self.Ix(N)**2 + self.Iy(N)**2 + self.Iz(N)**2).tidyup()
        self.Ip = lambda N : qt.jmat((N/2.), '+')
        self.Im = lambda N : qt.jmat((N/2.), '-')

        self.S2 = lambda N : qt.tensor(self.Sx(N),  self.identity_ancilla) **2 + qt.tensor(self.Sy(N),  self.identity_ancilla) **2 + qt.tensor(self.Sz(N),  self.identity_ancilla) **2        
        self.I2 = lambda N : qt.tensor(self.identity_central, self.Ix(N)) **2 + qt.tensor(self.identity_central, self.Iy(N)) **2 + qt.tensor(self.identity_central, self.Iz(N)) **2
        self.m_dm = lambda m, N : qt.ket2dm(qt.basis(N + 1, m))
        self.s_dm = lambda s : qt.ket2dm(qt.basis(2, s))

        self.M_stg = lambda N : qt.tensor(self.Sz, self.identity_ancilla(N))- 1/(N) * qt.tensor(self.identity_central, self.Iz(N))
        self.identity = lambda N1, N2 : qt.tensor(self.identity_central, self.identity_ancilla(N1),self.identity_central, self.identity_ancilla(N2))
        
    def get_local_hamiltonian(self, N, b, gamma, alfa):
        H_elec =  - N * gamma * self.Sz
        H_nuclei = -b * self.Iz(N)
        
        H_hf = 2 * qt.tensor(self.Sz, self.Iz(N)) + alfa * qt.tensor(self.Sp,  self.Im(N)) + alfa * qt.tensor(self.Sm, self.Ip(N) )
        
        H = qt.Qobj((qt.tensor(H_elec, self.identity_ancilla(N))  + qt.tensor(self.identity_central, H_nuclei) + H_hf ).full())
        return H
        
    def get_transverse_hamiltonian(self, N, Delta):
        H_elec =  Delta * self.Sx
        
        H = qt.Qobj((qt.tensor(H_elec, self.identity_ancilla(N))).full())
        return H
        
    def get_ising_hamiltonian(self, N1, N2, J):
        H = J * qt.Qobj((qt.tensor(self.Sz, self.identity_ancilla(N1),self.Sz, self.identity_ancilla(N2))).full())
        return H    

    def get_total_hamiltonian(self, N1, N2, v, t, Delta, J):
        H =  max(0, 1- v*t) * qt.Qobj(qt.tensor( self.get_transverse_hamiltonian(N1, Delta), qt.tensor(self.identity_central, self.identity_ancilla(N2)) ).full()) + \
             max(0, 1- v*t) * qt.Qobj(qt.tensor( qt.tensor(self.identity_central, self.identity_ancilla(N1)), self.get_transverse_hamiltonian(N2, Delta) ).full()) + \
             min(1, v*t) * qt.Qobj(qt.tensor( self.local_hamiltonian_1, qt.tensor(self.identity_central, self.identity_ancilla(N2)) ).full()) + \
             min(1, v*t) * qt.Qobj(qt.tensor( qt.tensor(self.identity_central, self.identity_ancilla(N1)), self.local_hamiltonian_2 ).full()) + \
             min(1, v*t) * self.get_ising_hamiltonian(N1, N2, J)
        return H 
    
    def calculate_eigenstates(self, N, b, gamma, alfa, ham):
        def get_eigenstate(ps,I,M,N, b,gamma, alfa,p, ham,ns1,ns2 = False):
            if ps == 'u':
                state = qt.tensor(qt.basis(2, 0), ns1['state'])
                te = -0.5*N*gamma-I*(b-1)
                rep = '|%s >|%g %g %d>'%(UPARROW, ns1['M'], ns1['I'], ns1['p'])
            elif ps == 'd':
                state = qt.tensor(qt.basis(2, 1), ns1['state'])    
                te = 0.5*N*gamma+I*(b+1)
                rep = '|%s >|%g %g %d>'%(DOWNARROW, ns1['M'], ns1['I'], ns1['p'])    
            elif ps == '+':
                f = phi(N,I,M,b,gamma, alfa)
                state = math.cos(f/2.) * qt.tensor(qt.basis(2, 0),ns1['state']) + math.sin(f/2.) * qt.tensor(qt.basis(2, 1),ns2['state'])
                te = 0.5* x(M,b)  + 0.5*  y(N,I,M,b,gamma, alfa)
                rep = '%+.3f|%s >|%g %g %d> %+.3f|%s >|%g %g %d>'%(math.cos(f/2.), UPARROW, ns1['M'], ns1['I'], ns1['p'], math.sin(f/2.), DOWNARROW, ns2['M'], ns2['I'], ns2['p'])    
            elif  ps == '-':
                f = phi(N,I,M,b,gamma, alfa)
                state = -math.sin(f/2.) * qt.tensor(qt.basis(2, 0),ns1['state']) + math.cos(f/2.) * qt.tensor(qt.basis(2, 1),ns2['state'])
                te = 0.5* x(M,b)  - 0.5* y(N,I,M,b,gamma, alfa)
                rep = '%+.3f|%s >|%g %g %d> %+.3f|%s >|%g %g %d>'%(-math.sin(f/2.), UPARROW, ns1['M'], ns1['I'], ns1['p'], math.cos(f/2.), DOWNARROW, ns2['M'], ns2['I'], ns2['p'])    
            return {'ps' : ps,'I' : I,'M' : M,'p' : p,'state' : qt.Qobj(state.full()).tidyup().unit(), 'theorEnergy': te,'energy' : ham.matrix_element(state.dag(),state).real, 'representation' : rep} 
    
        eigenstates = []
        for ns in self.get_ancilla_states(N):
            I, M, p = ns['I'], ns['M'], ns['p']
            
            if round(M,2) == round(I,2):
                eigenstates.append(get_eigenstate('u',I,M,N,b, gamma, alfa,p, ham,ns))
            
            else:
                for ps in ('-', '+'):
                    states2 = [elem for elem in self.get_ancilla_states(N) if (elem['I'] == I and elem['p'] == p and elem['M'] == M + 1)]
                    if states2:
                        ns2 = states2[0]
                    else:
                        print '{} {}'.format(ns['I'], ns['M'])
                    eigenstates.append(get_eigenstate(ps,I,M,N,b, gamma, alfa,p, ham,ns,ns2))
            if M == -I:
                eigenstates.append(get_eigenstate('d',I,M,N,b, gamma, alfa,p, ham,ns))
                
        eigenstates = sorted([s for s in eigenstates], key = lambda s : s['energy'])
        return     eigenstates
    

    
    
    def get_dm(self, hamiltonian): 
        states = sorted(zip(hamiltonian.eigenstates()[1], hamiltonian.eigenstates()[0]), key= lambda x: x[1])
        gs_energy = states[0][1]
        
        if self.parameters.tau_dot > 0:
            dm = (-h_init/self.parameters.tau_dot).expm().matrix_element(eigenstates[0]['state'].dag(), eigenstates[0]['state']) * qt.ket2dm(eigenstates[0]['state'])
            for es in eigenstates[1:]:
                dm += (-h_init/self.parameters.tau_dot).expm().matrix_element(es['state'].dag(), es['state']) * qt.ket2dm(es['state'])
                

        else:
            dm = qt.ket2dm(states[0][0])
            for state in states[1:]:
                ham1 = qt.Qobj((qt.tensor(self.identity_central, self.Iz(self.parameters.N1), self.identity_central, self.identity_ancilla(self.parameters.N2))).full())
                iz1 = ham1.matrix_element(state[0].dag(),state[0]).real
                ham2 = qt.Qobj((qt.tensor(self.identity_central, self.identity_ancilla(self.parameters.N1), self.identity_central, self.Iz(self.parameters.N2))).full())
                iz2 = ham2.matrix_element(state[0].dag(),state[0]).real
                print 'Energy={} iz1={} iz2={}'.format(state[1], iz1, iz2)
                if round(state[1], 9) == round(gs_energy, 9):
                    dm += qt.ket2dm(state[0])
                    #raw_input()
                    
        print 'Ground state energy: {}'.format(gs_energy)
        print '''Initial Density matrix:
                {}'''.format(dm)       
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
    
    def get_central_dissipators(self,  N, b, gamma):
        
        baseRate = self.parameters.eta
        
        rate_plus = (baseRate * Ifunc(- N * gamma, self.parameters.tau_dot) )**0.5
        rate_minus = (baseRate * Ifunc(N * gamma, self.parameters.tau_dot) )**0.5
        print 'r+={0} r-={1}'.format(rate_plus, rate_minus)
        dissipators = [rate_plus * qt.Qobj((qt.tensor(self.Sp, self.identity_ancilla, self.identity_central, self.identity_ancilla) + qt.tensor( self.identity_central, self.identity_ancilla, self.Sp, self.identity_ancilla)).full()), rate_minus * qt.Qobj((qt.tensor(self.Sm, self.identity_ancilla, self.identity_central, self.identity_ancilla) + qt.tensor( self.identity_central, self.identity_ancilla, self.Sm, self.identity_ancilla)).full())]
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
        return operators

    
    def compute_dynamics(self):
        
        def L(t, args):
            L0 = qt.Qobj(qt.liouvillian(self.get_total_hamiltonian(p.N1, p.N2, p.annealing_v, t, p.J, p.Delta), []).full())
            
            #i = int(t)
            #print t, 
            return L0.data + L1.data
            
        p = self.parameters
        
        all_dissipators = []
        
        for n in xrange(p.N1 + 1):
            m = p.N1/2. - n             
            for s in (-1,1):
                Delta_e = p.gamma1 - 2. * m / p.N1 - 0.5 * p.J * s
                rp = p.eta1**0.5 * HST(Delta_e)
                rm = p.eta1**0.5 * HST(-Delta_e)
                
                all_dissipators += [rp * qt.tensor(self.Sp, self.m_dm(n, p.N1), self.s_dm((1-s)/2), self.identity_ancilla(p.N2) ), rm * qt.tensor(self.Sm, self.m_dm(n, p.N1), self.s_dm((1-s)/2), self.identity_ancilla(p.N2))]
        for n in xrange(p.N2 + 1):
            m = p.N2/2. - n             
            for s in (-1,1):                
                Delta_e = p.gamma2 - 2 * m / p.N2 - 0.5 * p.J * s
                rp = p.eta2**0.5 * HST(Delta_e)
                rm = p.eta2**0.5 * HST(-Delta_e)
                
                all_dissipators += [rp * qt.tensor(self.s_dm((1-s)/2), self.identity_ancilla(p.N1), self.Sp, self.m_dm(n, p.N2) ), rm * qt.tensor(self.s_dm((1-s)/2), self.identity_ancilla(p.N1), self.Sm, self.m_dm(n, p.N2))]
        print "There are {} dissipators".format(len(all_dissipators))                       
        L1 = qt.Qobj(qt.liouvillian(self.identity(p.N1, p.N2), all_dissipators).full())
        
        #self.operators_compute = self.get_operators_compute()
        operators_compute = [qt.Qobj((qt.tensor(self.Sz, self.identity_ancilla(p.N1),self.identity_central, self.identity_ancilla(p.N2))).full()), 
                            qt.Qobj((qt.tensor(self.identity_central, self.identity_ancilla(p.N1),self.Sz, self.identity_ancilla(p.N2))).full())]     
        try:
            #print self.tlist
            opts = qt.Odeoptions(nsteps = 10000000, rhs_reuse=True, min_step=0.001, max_step=0, num_cpus = 4, rtol = 0.001, atol = 0.00001, store_states = True)
            self.expt_list = qt.mesolve(L, self.density_matrix, self.tlist, [], operators_compute, progress_bar = True, options=opts)
            
        except Exception as e:
            print e
            self.expt_list = 'Error'
            
        folder = '_'.join(['{}={}'.format(k, v) for k, v in p.__class__.__dict__.iteritems() if (not k.startswith('_') and k != 'get_parameters')])
        
        print folder
        if not os.path.exists('archive/annealing_{}'.format(folder)):
            os.makedirs('archive/annealing_{}'.format(folder))
        for i, r in enumerate(self.expt_list.expect): 
            with open('archive/annealing_{}/expectation_{}'.format(folder, i), 'w') as f:
                for t, res in zip(self.expt_list.times, r):
                    f.write('{} {}'.format(t, res) + '\n') 
        #with open('archive/annealing_{}/states'.format(folder), 'w') as f:
        #    for t, st in zip(self.expt_list.times, self.expt_list.states):
        #        f.write('{} {}'.format(t, st.full()) + '\n')   
        
        return self.expt_list
    
    def prepare_computation(self,parameters):    
        startTime = time.time()
        
        self.parameters = parameters
        p = self.parameters
        
        #if not THRUSTY:
        #    represent(self.nucleiStates)

        self.define_operators()

        self.local_hamiltonian_1 = self.get_local_hamiltonian(p.N1, p.b1, p.gamma1, p.alfa1)
        self.local_hamiltonian_2 = self.get_local_hamiltonian(p.N2, p.b2, p.gamma2, p.alfa2)
        
        h0 = self.get_total_hamiltonian(p.N1, p.N2, p.annealing_v, 0, p.J, p.Delta)
        
        self.density_matrix = self.get_dm(h0)
        
        self.tlist = np.linspace(0, p.tInterval, num=p.tSteps, endpoint=True)
       
        self.eigenstates_1 = self.calculate_eigenstates(p.N1,p.b1, p.gamma1, p.alfa1, self.local_hamiltonian_1)
        self.eigenstates_2 = self.calculate_eigenstates(p.N2,p.b2, p.gamma2, p.alfa2, self.local_hamiltonian_2)

        print 'Hamiltonian 1: b= {}  g= {}'.format(p.b1, p.gamma1)
        print 'Eigenstates: (energy)(theorEnergy) State |ps M I>', len(self.eigenstates_1)    
        for gs in self.eigenstates_1:
            print '(%.3f)(%.3f) |%s %+g %g> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['M'],gs['I'], gs['representation'])
        
        print 'Hamiltonian 2: b= {}  g= {}'.format(p.b2, p.gamma2)
        print 'Eigenstates: (energy)(theorEnergy) State |ps M I>', len(self.eigenstates_2)    
        for gs in self.eigenstates_2:
            print '(%.3f)(%.3f) |%s %+g %g> = %s'%(gs['energy'],gs['theorEnergy'],gs['ps'], gs['M'],gs['I'], gs['representation'])

        statesTime = time.time()
        print 'States and operators got in %f'%(statesTime-startTime)
        
        print 'J={} Delta={}'.format(p.J, p.Delta)
   
        #print 'Total Hamiltonian at t=0:{}'.format(h0)
        
        #self.all_L0 = [qt.Qobj(qt.liouvillian(self.get_total_hamiltonian(p.annealing_v, t, p.J, p.Delta), []).full()) for t in self.tlist]
        #self.all_L1 = [qt.Qobj(qt.liouvillian(self.identity, self.all_dissipators[-1]).full())  for t in self.tlist]
        #self.L = [self.all_L0[i].data + self.all_L1[i].data for i in xrange(len(self.tlist))]
        #liouTime = time.time()
        #print 'Liouvillians got in %f'%(liouTime-dis_time)        
        
        
    
