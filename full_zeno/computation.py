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
        
    def define_operators(self):
        def nuclearOperator(QObj, n):
            nucleiList = [qt.qeye(2) for _ in xrange(N)]
            nucleiList[n] = QObj
            return qt.tensor(nucleiList)
        N, b, gamma, alfa, beta = self.N, self.b, self.gamma, self.alfa, self.beta
        self.identity_central = qt.qeye(2)
        self.Sm = qt.sigmam()
        self.Sp = qt.sigmap()
        self.Sx = 0.5 * qt.sigmax()
        self.Sy = 0.5 * qt.sigmay()
        self.Sz = 0.5 * qt.sigmaz()
        
        self.identity_ancilla = qt.tensor([qt.qeye(2) for _ in xrange(N)])
        self.Im = sum([nuclearOperator(qt.sigmam(), n) for n in xrange(N)])
        self.Ip = sum([nuclearOperator(qt.sigmap(), n) for n in xrange(N)])
        self.Ix = sum([nuclearOperator(0.5 * qt.sigmax(), n) for n in xrange(N)])
        self.Iy = sum([nuclearOperator(0.5 * qt.sigmay(), n) for n in xrange(N)])
        self.Iz = sum([nuclearOperator(0.5 * qt.sigmaz(), n) for n in xrange(N)])

        self.S2 = qt.tensor(self.Sx,  self.identity_ancilla) **2 + qt.tensor(self.Sy,  self.identity_ancilla) **2 + qt.tensor(self.Sz,  self.identity_ancilla) **2        
        self.I2 = qt.tensor(self.identity_central, self.Ix) **2 + qt.tensor(self.identity_central, self.Iy) **2 + qt.tensor(self.identity_central, self.Iz) **2
        
        self.identity = qt.tensor(self.identity_central, self.identity_ancilla)
        self.M_stg = qt.tensor(self.Sz, self.identity_ancilla)- 1/(N) * qt.tensor(self.identity_central, self.Iz)
        
        self.H_central =  - N * gamma * self.Sz
        
        self.gradient = lambda n : beta  * (n/ (N-1.) -0.5) if N > 1 else 0
        print 'ancilla Magnetic field: ', [-b + self.gradient(n) for n in xrange(N)]
        
        self.H_ancilla = sum([(-b + self.gradient(n)) * nuclearOperator(0.5 * qt.sigmaz(), n) for n in xrange(N)])
        
        self.H_hf = (2 * qt.tensor(self.Sz, self.Iz) + alfa * qt.tensor(self.Sp,  self.Im) + alfa * qt.tensor(self.Sm, self.Ip ) )
        
        self.H = qt.Qobj((qt.tensor(self.H_central, self.identity_ancilla)  + qt.tensor(self.identity_central, self.H_ancilla) + self.H_hf ).full())
        print self.H

class Computation(object):
    def __init__(self, parameters):
        self.parameters = parameters
        self.folder = "dephasing"
        self.cluster = Cluster(self.parameters.cluster_parameters)
        
        self.H = self.cluster.H

    
    def get_operators_compute(self):
        op = 1./(self.cluster.N/2.) * qt.tensor(self.cluster.identity_central, self.cluster.Iz)
        print op
        return op
         

    def get_dm(self): 
        vec = [qt.basis(2, 1) for n in xrange(self.cluster.N+1)]
        state = qt.tensor(vec) 
        dm = qt.ket2dm(state)
        
        print "Initial energy: {}".format(self.H.matrix_element(state.dag(),state).real)
        print "Initial magnetization: {}".format(qt.tensor(self.cluster.identity_central, self.cluster.Iz).matrix_element(state.dag(),state).real)
        print """Initial Density matrix:
                {}""".format(dm)       
        return dm/dm.tr()
    
    def prepare_computation(self):    
        startTime = time.time()
        p = self.parameters       
        self.density_matrix = self.get_dm()        
        self.tlist = np.linspace(0, p.tInterval, num=p.tSteps, endpoint=True)        
    
    def compute_dynamics(self):        
        p = self.parameters
        operators_compute = self.get_operators_compute()
        ETA_POINTS = 200
        results = []
        for i in xrange(ETA_POINTS):
            baserate = (0.0 + self.cluster.eta*i/(ETA_POINTS-1.))
            dissipator = baserate * qt.tensor(self.cluster.Sp, self.cluster.identity_ancilla)  
        
            L = qt.liouvillian(self.H, [dissipator])   
       
        
            try:
                #print self.tlist
                opts = qt.Odeoptions(nsteps = 10000000, rhs_reuse=True, min_step=0.001, max_step=0, num_cpus = 4, rtol = 0.001, atol = 0.00001, store_states = True)
                self.expt_list = qt.mesolve(L, self.density_matrix, self.tlist, [], operators_compute, progress_bar = True, options=opts)
                results.append(list(self.expt_list.expect[0]))
            except Exception as e:
                print e
                self.expt_list = "Error"
        folder = "archive/zeno_beta={}".format(self.cluster.beta)      
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        with open("{}/expectation".format(folder), "w") as f:
            for i in xrange(len(self.expt_list.times)):
                string = "{} {}".format(self.expt_list.times[i], " ".join([str(r[i]) for r in results]))
                f.write("{}".format(string) + "\n") 
                        
        with open("{}/parameters".format(folder), "w") as f:
            for l in ["{} = {}".format(k, v) for k, v in self.parameters.__class__.__dict__.iteritems() ]:
                f.write(l + "\n")    
        
        return self.expt_list
    

    
