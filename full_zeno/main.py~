# -*- coding: UTF-8 -*-
        
if __name__ == '__main__':
    import warnings, sys
    warnings.filterwarnings("ignore")
    if len(sys.argv) > 1 and sys.argv[-1] in ('-w', '--window'):
        from window_gui import *
        root = Root()
        root.iconbitmap('@icon.xbm')
        mainWindow = MainWindow(root)    
        root.mainloop()
    else:
        from computation import *
        from parameters import *
        computation = Computation()
        parameters = Parameters()
        computation.prepare_computation(parameters)
        dynamics_result = computation.compute_dynamics()
        
        
        


