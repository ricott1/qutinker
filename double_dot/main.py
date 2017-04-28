# -*- coding: UTF-8 -*-
###
#there are two main modules in this package: the quTinkerGraphic provides a high level user interface. It sends data via the getData function to the quTinkerComputation module where the actual computation take place using qutip. Finally the results are sent back to get displayed.
###





from mainWindow import *			
rescaled = True		
		
if __name__ == '__main__':
	import warnings
	warnings.filterwarnings("ignore")
	root = Root()
	root.iconbitmap('@icon.xbm')
	mainWindow = MainWindow(root)	
		
	
	
	root.mainloop()


