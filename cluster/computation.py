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
    print " Some modules are missing, try installing them using \"installModules.py\""
    sys.exit()
THRUSTY = 0
NPROCS = 1
UPARROW = u"\u2191"
DOWNARROW = u"\u2193"


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

class Cluster(object):
    def __init__(self, p):
        self.parameters = {"N" : p[0], "gamma" : p[1], "b" : p[2], "alfa" : p[3], "beta" : p[4],"eta" : p[5], "tau_dot" : p[6], "tau_lead" : p[7]}
        for p, v in self.parameters.iteritems():
            setattr(self, p, v)
        self.define_operators()
        self.local_hamiltonian = self.get_local_hamiltonian()
        
    def get_ancilla_states(self):
        N = self.N
        i = int(N)
        while i >= 0:
            state = {"I" : N/2.,"M" : - i + N/2.,"p" : 1,"state" : qt.basis(N + 1, i)} 
            
            yield state
            i -= 1
    def define_operators(self):
        N = self.N
        self.identity_central = qt.qeye(2)
        self.Sm = qt.sigmam()
        self.Sp = qt.sigmap()
        self.Sx = 0.5 * qt.sigmax()
        self.Sy = 0.5 * qt.sigmay()
        self.Sz = 0.5 * qt.sigmaz()
        
        self.identity_ancilla = qt.qeye(N + 1)
        self.Ix = qt.jmat((N/2.), "x")
        self.Iy = qt.jmat((N/2.), "y")
        self.Iz = qt.jmat((N/2.), "z")
        self.Isq = (self.Ix**2 + self.Iy**2 + self.Iz**2).tidyup()
        self.Ip = qt.jmat((N/2.), "+")
        self.Im = qt.jmat((N/2.), "-")

        self.S2 = qt.tensor(self.Sx,  self.identity_ancilla) **2 + qt.tensor(self.Sy,  self.identity_ancilla) **2 + qt.tensor(self.Sz,  self.identity_ancilla) **2        
        self.I2 = qt.tensor(self.identity_central, self.Ix) **2 + qt.tensor(self.identity_central, self.Iy) **2 + qt.tensor(self.identity_central, self.Iz) **2
        self.m_dm = lambda m : qt.ket2dm(qt.basis(N + 1, m))
        self.s_dm = lambda s : qt.ket2dm(qt.basis(2, s))
        self.identity = qt.tensor(self.identity_central, self.identity_ancilla)
        self.M_stg = qt.tensor(self.Sz, self.identity_ancilla)- 1/(N) * qt.tensor(self.identity_central, self.Iz)
        
        
    def get_local_hamiltonian(self):
        N, b, gamma, alfa, beta = self.N, self.b, self.gamma, self.alfa, self.beta     
    
        H_elec =  - N * gamma * self.Sz
        H_nuclei = -b * self.Iz
        
        H_hf = beta * (2 * qt.tensor(self.Sz, self.Iz) + alfa * qt.tensor(self.Sp,  self.Im) + alfa * qt.tensor(self.Sm, self.Ip ) )
        
        H = qt.Qobj((qt.tensor(H_elec, self.identity_ancilla)  + qt.tensor(self.identity_central, H_nuclei) + H_hf ).full())
        return H
    
    def get_transverse_hamiltonian(self, Delta):
        H_elec =  Delta * self.Sx
        
        H = qt.Qobj((qt.tensor(H_elec, self.identity_ancilla)).full())
        return H  
    
    
    def calculate_eigenstates(self, N, b, gamma, alfa, ham):
        def get_eigenstate(ps,I,M,N, b,gamma, alfa,p, ham,ns1,ns2 = False):
            if ps == "u":
                state = qt.tensor(qt.basis(2, 0), ns1["state"])
                te = -0.5*N*gamma-I*(b-1)
                rep = "|%s >|%g %g %d>"%(UPARROW, ns1["M"], ns1["I"], ns1["p"])
            elif ps == "d":
                state = qt.tensor(qt.basis(2, 1), ns1["state"])    
                te = 0.5*N*gamma+I*(b+1)
                rep = "|%s >|%g %g %d>"%(DOWNARROW, ns1["M"], ns1["I"], ns1["p"])    
            elif ps == "+":
                f = phi(N,I,M,b,gamma, alfa)
                state = math.cos(f/2.) * qt.tensor(qt.basis(2, 0),ns1["state"]) + math.sin(f/2.) * qt.tensor(qt.basis(2, 1),ns2["state"])
                te = 0.5* x(M,b)  + 0.5*  y(N,I,M,b,gamma, alfa)
                rep = "%+.3f|%s >|%g %g %d> %+.3f|%s >|%g %g %d>"%(math.cos(f/2.), UPARROW, ns1["M"], ns1["I"], ns1["p"], math.sin(f/2.), DOWNARROW, ns2["M"], ns2["I"], ns2["p"])    
            elif  ps == "-":
                f = phi(N,I,M,b,gamma, alfa)
                state = -math.sin(f/2.) * qt.tensor(qt.basis(2, 0),ns1["state"]) + math.cos(f/2.) * qt.tensor(qt.basis(2, 1),ns2["state"])
                te = 0.5* x(M,b)  - 0.5* y(N,I,M,b,gamma, alfa)
                rep = "%+.3f|%s >|%g %g %d> %+.3f|%s >|%g %g %d>"%(-math.sin(f/2.), UPARROW, ns1["M"], ns1["I"], ns1["p"], math.cos(f/2.), DOWNARROW, ns2["M"], ns2["I"], ns2["p"])    
            return {"ps" : ps,"I" : I,"M" : M,"p" : p,"state" : qt.Qobj(state.full()).tidyup().unit(), "theorEnergy": te,"energy" : ham.matrix_element(state.dag(),state).real, "representation" : rep} 
    
        eigenstates = []
        for ns in self.get_ancilla_states(N):
            I, M, p = ns["I"], ns["M"], ns["p"]
            
            if round(M,2) == round(I,2):
                eigenstates.append(get_eigenstate("u",I,M,N,b, gamma, alfa,p, ham,ns))
            
            else:
                for ps in ("-", "+"):
                    states2 = [elem for elem in self.get_ancilla_states(N) if (elem["I"] == I and elem["p"] == p and elem["M"] == M + 1)]
                    if states2:
                        ns2 = states2[0]
                    else:
                        print "{} {}".format(ns["I"], ns["M"])
                    eigenstates.append(get_eigenstate(ps,I,M,N,b, gamma, alfa,p, ham,ns,ns2))
            if M == -I:
                eigenstates.append(get_eigenstate("d",I,M,N,b, gamma, alfa,p, ham,ns))
                
        eigenstates = sorted([s for s in eigenstates], key = lambda s : s["energy"])
        return     eigenstates
      
    

    

class Computation(object):
    def __init__(self, parameters):
        self.parameters = parameters
        
        self.folder = "try4"
        self.cluster_parameters = []
        self.clusters = [Cluster(self.parameters.cluster_parameters[c]) for c in xrange(self.parameters.clusters)]
        self.identity = qt.tensor([c.identity for c in self.clusters])
        self.local_H = self.get_local_hamiltonian()
        self.transverse_H = self.get_transverse_hamiltonian()
        self.ising_H = self.get_ising_hamiltonian()
    
    def get_total_hamiltonian(self, t):
        v = self.parameters.annealing_v
        H =  max(0, 1- v*t) * self.transverse_H + min(1, v*t) * self.local_H + min(1, v*t) * self.ising_H     
        return H     

    def get_local_hamiltonian(self):
        for i, c in enumerate(self.clusters):
            vector = [g.local_hamiltonian if i==j else g.identity for j, g in enumerate(self.clusters) ]
            if i==0:
                H = qt.Qobj(qt.tensor(vector).full())
            else:
                H += qt.Qobj(qt.tensor(vector).full())
            
        return H
    
    def get_operators_compute(self):
        ops = []
        for i, c in enumerate(self.clusters):
            vector = [qt.tensor(g.Sz, g.identity_ancilla) if i==j else g.identity for j, g in enumerate(self.clusters) ]
            ops.append(qt.Qobj(qt.tensor(vector).full()) )
            
        return ops
         
    def get_transverse_hamiltonian(self):
        for i, c in enumerate(self.clusters):
            vector = [g.get_transverse_hamiltonian(self.parameters.Delta)  if i==j else g.identity for j, g in enumerate(self.clusters)]
            if i==0:
                H = qt.Qobj(qt.tensor(vector).full())
            else:
                H += qt.Qobj(qt.tensor(vector).full())
            
        return H
        
    def get_ising_hamiltonian(self):
        
        for i, c in enumerate(self.clusters):
            
            for j, g in enumerate(self.clusters):
                if i<j: 
                    vector = [qt.tensor(d.Sz, d.identity_ancilla) if k in (i, j) else d.identity for k, d in enumerate(self.clusters) ]
                    if i==0 and j==1:
                        H = self.parameters.J[i,j] * qt.Qobj(qt.tensor(vector).full())
                    else:
                        H += self.parameters.J[i,j] * qt.Qobj(qt.tensor(vector).full())
            
        return H
    
    
    def get_dm(self, hamiltonian): 
        states = sorted(zip(hamiltonian.eigenstates()[1], hamiltonian.eigenstates()[0]), key= lambda x: x[1])
        gs_energy = states[0][1]
        tau_dot = 0
        if tau_dot > 0:
            dm = (-h_init/tau_dot).expm().matrix_element(eigenstates[0]["state"].dag(), eigenstates[0]["state"]) * qt.ket2dm(eigenstates[0]["state"])
            for es in eigenstates[1:]:
                dm += (-h_init/tau_dot).expm().matrix_element(es["state"].dag(), es["state"]) * qt.ket2dm(es["state"])
                

        else:
            dm = qt.ket2dm(states[0][0])
            for state in states[1:]:
                if round(state[1], 9) == round(gs_energy, 9):
                    dm += qt.ket2dm(state[0])
                    #raw_input()
                    
        print "Ground state energy: {}".format(gs_energy)
        #print """Initial Density matrix:
        #        {}""".format(dm)       
        return dm/dm.tr()
    
    def get_dissipators(self):
        dissipators = []
        done = []
        
        #here im double counting cause i replace both 0 and 1 with operators on that cluster.
        for i in xrange(2**len(self.clusters)):#cycle over all possible ising states
            b = format(i, "0{}b".format(len(self.clusters)) )
            
            base = []
            for j, s in enumerate(b):
                base += [self.clusters[j].s_dm(int(s)), self.clusters[j].identity_ancilla ]
            print b#, base
            #now we have a base vector containing this particular b state
            for j, c in enumerate(self.clusters):
                q = list(b)
                q[j] = 'd'
                    
                if q not in done:
                    done.append(q)   
                    for n in xrange(c.N + 1):
                        m = c.N/2. - n
                        Delta_e = c.gamma - 2 * c.beta * m / c.N - 0.5 * sum([self.parameters.J[j,l] * np.sign(0.5 - int(b[l])) for l in xrange(len(self.clusters))])
                        vec = list(base)
                        
                        vec[2*j], vec[2*j+1] = c.Sp, c.m_dm(n)
                        rate = (c.eta * HST(Delta_e))**0.5
                        dissipators.append( rate * qt.tensor(vec) )
                        
                        vec[2*j], vec[2*j+1] = c.Sm, c.m_dm(n)
                        rate = (c.eta * HST(-Delta_e))**0.5
                        dissipators.append( rate * qt.Qobj(qt.tensor(vec).full()) )
                                
                
        print len(dissipators)#, dissipators[0]
        return dissipators
    
    
    def prepare_computation(self):    
        startTime = time.time()
        p = self.parameters
        self.local_hamiltonians = []
        for c in self.clusters:
            self.local_hamiltonians.append(c.local_hamiltonian)
        
        h0 = self.get_total_hamiltonian(0)
        
        self.density_matrix = self.get_dm(h0)
        
        self.tlist = np.linspace(0, p.tInterval, num=p.tSteps, endpoint=True)
        print "J={} Delta={}".format(p.J, p.Delta)
   
        #print "Total Hamiltonian at t=0:{}".format(h0)
    
    def compute_dynamics(self):
        
        def L(t, args):
            L0 = qt.Qobj(qt.liouvillian(self.get_total_hamiltonian(t), []).full())
            return L0.data + L1.data
            
        p = self.parameters
        
        dissipators = self.get_dissipators()   
        print "There are {} dissipators".format(len(dissipators))     
        L1 = qt.Qobj(qt.liouvillian(self.identity, dissipators).full())
        
        operators_compute = self.get_operators_compute()
        try:
            #print self.tlist
            opts = qt.Odeoptions(nsteps = 10000000, rhs_reuse=True, min_step=0.001, max_step=0, num_cpus = 4, rtol = 0.001, atol = 0.00001, store_states = True)
            self.expt_list = qt.mesolve(L, self.density_matrix, self.tlist, [], operators_compute, progress_bar = True, options=opts)
            
        except Exception as e:
            print e
            self.expt_list = "Error"
        
        if not os.path.exists("archive/annealing_{}".format(self.folder)):
            os.makedirs("archive/annealing_{}".format(self.folder))
        for i, r in enumerate(self.expt_list.expect): 
            with open("archive/annealing_{}/expectation_{}".format(self.folder, i), "w") as f:
                for t, res in zip(self.expt_list.times, r):
                    f.write("{} {}".format(t, res) + "\n") 
                    
        with open("archive/annealing_{}/parameters".format(self.folder), "w") as f:
            for l in ["{} = {}".format(k, v) for k, v in self.parameters.__class__.__dict__.iteritems() ]:
                f.write(l + "\n")    
        
        return self.expt_list
    
     
     
     
     
            
    def show_states(self):
        self.eigenstates_1 = self.calculate_eigenstates(p.N1,p.b1, p.gamma1, p.alfa1, self.local_hamiltonian_1)
        self.eigenstates_2 = self.calculate_eigenstates(p.N2,p.b2, p.gamma2, p.alfa2, self.local_hamiltonian_2)

        print "Hamiltonian 1: b= {}  g= {}".format(p.b1, p.gamma1)
        print "Eigenstates: (energy)(theorEnergy) State |ps M I>", len(self.eigenstates_1)    
        for gs in self.eigenstates_1:
            print "(%.3f)(%.3f) |%s %+g %g> = %s"%(gs["energy"],gs["theorEnergy"],gs["ps"], gs["M"],gs["I"], gs["representation"])
        
        print "Hamiltonian 2: b= {}  g= {}".format(p.b2, p.gamma2)
        print "Eigenstates: (energy)(theorEnergy) State |ps M I>", len(self.eigenstates_2)    
        for gs in self.eigenstates_2:
            print "(%.3f)(%.3f) |%s %+g %g> = %s"%(gs["energy"],gs["theorEnergy"],gs["ps"], gs["M"],gs["I"], gs["representation"])

        statesTime = time.time()
        print "States and operators got in %f"%(statesTime-startTime)
        
        
    
