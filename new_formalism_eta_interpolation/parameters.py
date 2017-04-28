class Parameters(object):
    N = 1
    
    tau_dot = 0
    tau_lead = 0
    gamma_init = -2
    b_init = N * gamma_init
    gamma_final = -1.55
    b_final = 2
    alfa = 1
    eta = 10
    tInterval = 20
    tSteps = 600
    quench_v = 0.1
    startTime = 0
    wait_time = 10
    dissipators_kind = 'central'
    points = 2
    
    def get_parameters(self, master):
        self.N = max(1, int(master.N_entry.get()))
        
        self.tau_dot = float(master.tau_dot_entry.get())
        self.tau_lead = max(float(master.tau_lead_entry.get()), 0)
        self.gamma_init = float(master.gamma_init_entry.get())
        self.b_init = float(master.b_init_entry.get())
        self.gamma_final = float(master.gamma_final_entry.get())
        self.b_final = float(master.b_final_entry.get())
        self.alfa = float(master.alfa_entry.get())
        self.eta = max(float(master.eta_entry.get()), 0)
        self.tInterval = float(master.tInterval_entry.get())
        self.tSteps = int(master.tSteps_entry.get())
        self.quench_v = float(master.annealing_v_entry.get())
        self.startTime = 0
        self.wait_time = float(master.wait_time_entry.get())
        self.dissipators_kind = master.dissipators_kind
        
