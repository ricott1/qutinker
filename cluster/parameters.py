import random
import numpy as np
class Parameters(object):
    clusters = 3
    N = [1, 1, 1]
    tau_dot = [0, 0,0]
    tau_lead = [0, 0,0]
    gamma = [-10, 10,0]
    b = [N[i] * gamma[i] for i in xrange(clusters)]
    
    alfa = [1, 1,1]
    beta = [0, 0,0]
    eta = [4, 4,4]
    
    cluster_parameters = zip(N, gamma, b, alfa, beta, eta, tau_dot, tau_lead)
    
    tInterval = 20
    tSteps = 600
    annealing_v = 0.1
    startTime = 0
    wait_time = 1
    dissipators_kind = 'central'
    
    
    J = np.matrix('0 1 1; 1 0 1; 1 1 0')
    print J
    print "Initial random Ising network: {}".format(J)
    Delta = 10
    
    def get_parameters(self, master):
        self.N1 = max(1, int(master.N1_entry.get()))
        self.N2 = max(1, int(master.N2_entry.get()))
        self.N = self.N1
        self.tau_dot = float(master.tau_dot_entry.get())
        self.tau_lead = max(float(master.tau_lead_entry.get()), 0)
        self.gamma1 = float(master.gamma1_entry.get())
        self.b1 = float(master.b1_entry.get())
        self.gamma2 = float(master.gamma2_entry.get())
        self.b2 = float(master.b2_entry.get())
        self.alfa1 = float(master.alfa1_entry.get())
        self.eta1 = max(float(master.eta1_entry.get()), 0)
        self.alfa2 = float(master.alfa2_entry.get())
        self.eta2 = max(float(master.eta2_entry.get()), 0)
        self.tInterval = float(master.tInterval_entry.get())
        self.tSteps = int(master.tSteps_entry.get())
        self.annealing_v = float(master.annealing_v_entry.get())
        self.startTime = 0
        self.wait_time = float(master.wait_time_entry.get())
        self.dissipators_kind = master.dissipators_kind
        self.J = float(master.J_entry.get())
        self.Delta = float(master.Delta_entry.get())
