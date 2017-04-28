# -*- coding: UTF-8 -*-
        
if __name__ == '__main__':
    import warnings, sys
    from computation import *
    from parameters import *
    computation = Computation(Parameters())
    computation.prepare_computation()
    dynamics_result = computation.compute_dynamics()
        
        
        


