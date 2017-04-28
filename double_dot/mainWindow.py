# -*- coding: UTF-8 -*-
try:
	import matplotlib
	import matplotlib.pyplot as plt
	import matplotlib.figure as mplfig
	from mpl_toolkits.mplot3d import Axes3D

	import matplotlib.backends.backend_tkagg as tkagg
	matplotlib.use('TkAgg')	
	from Tkinter import *
	import numpy as np
except:
	print ' Some modules are missing, try inp.arange(I+1)nstalling them using \'installModules.py\''
	import sys
	sys.exit()
import os
from computation import *
class Root(Tk):
	def __init__(self):
		
		Tk.__init__(self)
		self.wm_title('QuTinker')
		self.geometry("1520x900+0+0")



class Parameters(object):
	def __init__(self, master):
		self.master = master
		self.get_parameters()
	def get_parameters(self):
		self.K = max(1, int(self.master.K_entry.get()))
		self.tau_dot = float(self.master.tau_dot_entry.get())
		self.tau_lead = max(float(self.master.tau_lead_entry.get()), 0)
		self.rho_init = float(self.master.rho_init_entry.get())/self.K
		self.nu_init = float(self.master.nu_init_entry.get())
		self.rho_evol = float(self.master.rho_evol_entry.get())/self.K
		self.nu_evol = float(self.master.nu_evol_entry.get())
		self.labda = float(self.master.labda_entry.get())
		self.gamma0 = max(float(self.master.gamma0_entry.get()), 0)
		self.mu = max(float(self.master.mu_entry.get()), 0)
		self.alfa = float(self.master.alfa_entry.get()) 
		self.tInterval = float(self.master.tInterval_entry.get())
		self.tSteps = int(self.master.tSteps_entry.get())
		self.startTime = 0
		self.dis_time = float(self.master.dis_time_entry.get())
		self.evo_time = float(self.master.evo_time_entry.get())
class MainWindow(Frame):
	'''
	This class contains the main widget as a frame which we fill up with fancy buttons and entry fields. 
	There are many exhaustive documentations of Tkinter online, e.g. http://effbot.org/tkinterbook/tkinter-index.htm
	'''
	def __init__(self, parent):
		Frame.__init__(self, parent)
		self.computation = Computation(self)
		self.dynamics_result = []
		self.oldResults = []
		self.saveText = ''	
		self.build_borders()
		self.build_initial_state_frame()
		self.build_evolution_frame()
		self.build_dissipator_frame()
		self.build_time_schedule_frame()
		self.build_plot_parameters_frame()
		self.build_operators_bools()
		self.buildMatPlotLibFrame()
		self.build_phys_parameters_frame()
		self.build_buttons()
		self.grid()
	#this two custom functions (my_entry and my_check_button) allow an easier creation of entries and checkbutton. They just keep track of the variables associated with them.
	def my_entry(self, text,pos,default=1.0,master='self', comm = False):
		if master == 'self':
			master = self
		Label(master,text= text).grid(row=pos[0],column=pos[1],columnspan=2,sticky=W)
		if comm:
			entry = Entry(master, width=8, command = comm)
		else:
			entry = Entry(master, width=8)

		entry.insert(0,default)
		#entry.bind('<FocusIn>',lambda e: self.update_rescaled_parameters(e))
		entry.grid(row=pos[0],column=pos[1]+2)
		return entry	

	def my_check_button(self, text,pos,master='self', select=True, command=False):
		if master == 'self':
			master = self	
		Label(master,text= text).grid(row=pos[0],column=pos[1], sticky=W)
		booleanVar = BooleanVar()
		if command:
			cb = Checkbutton(master,variable = booleanVar, command = command)
		else:
			cb = Checkbutton(master,variable = booleanVar)
		cb.grid(row=pos[0],column=pos[1]+1)
		if select:
			cb.select()
		return booleanVar

	def build_borders(self):	
		#border frames to have some distance from the border of the widget
		Frame(self, height=20).grid(row=0, column=0, columnspan=10)
		Frame(self, width=35).grid(row=0, column=0, rowspan=7)
	def build_buttons(self):
		#we create the buttons to plot and save the figure
		self.plotButton = Button(self,text='PLOT',width=30, fg='red',height=8,command=self.plot)
		self.plotButton.grid(row = 6, column = 5, columnspan =3,sticky=E+W)		
	
		#Button(self,text='Quit',width=30, fg='red',command=sys.exit).grid(row=9,column=5,columnspan=3)
		
		self.saveButton = Button(self,text='SAVE',width=30, fg='blue',height=8,command=self.save, state=DISABLED)
		self.saveButton.grid(row=6,column=2, columnspan =3,sticky=E+W)
	def buildMatPlotLibFrame(self):
		#here we build the actual matplotlib frame. 
		Frame(self,height=500,width = 1000).grid(row=1,column=2, columnspan=6,rowspan=4, sticky=W+E+N+S)
		#we initialize the figure, setting also the title. There are a couple of tricks to make matplotlib not popup a widget automatically and embed it instead on tkinter
		self.figure = mplfig.Figure(figsize=(14,6))
		#hamiltonian = '$H = K \\rho \\nu S^z - \\nu I^z + 2 S^z I^z + \lambda (S^+ I^-+ S^- I^+)$'
		#self.figure.suptitle(hamiltonian, fontsize=24)
		
		self.canvas = tkagg.FigureCanvasTkAgg(self.figure, master=self)
		
		self.canvas.get_tk_widget().grid(row=1,column=2, columnspan=6,rowspan=4, sticky=W+E+N+S)
		self.canvas._tkcanvas.grid()
	def build_time_schedule_frame(self):
		self.timeScheduleFrame = LabelFrame(self, text='Time Schedule', padx=5, pady=5, relief= RIDGE)
		self.timeScheduleFrame.grid(row = 4, column = 1, sticky=W+E)
		self.evo_time_entry = self.my_entry('t evol: evolution starting time',(0,0),default=0, master= self.timeScheduleFrame)
		self.dis_time_entry = self.my_entry('t diss: dissipation starting time',(1,0),default=0, master= self.timeScheduleFrame)
		self.tInterval_entry = self.my_entry('Ending time',(2,0), default=60, master=self.timeScheduleFrame)
		self.tSteps_entry = self.my_entry('Time steps',(3,0), default=300, master=self.timeScheduleFrame)
	def build_plot_parameters_frame(self):
		##labelFrame for the plot parameters
		self.plotFrame = LabelFrame(self, text='Plot parameters', padx=5, pady=5, relief= RIDGE)
		self.plotFrame.grid(row = 5, column = 2, columnspan =3,sticky=E)	
		self.LegendBool = self.my_check_button(text='Legends',pos=(0,0),master=self.plotFrame, select=True)
		self.LegendManual = Entry(self.plotFrame, width=20)
		self.LegendManual.grid(row=0,column=2, columnspan = 4)
		self.keepButton = Button(self.plotFrame,text='Keep it!',width=10,command=self.keepResults, state=DISABLED)
		self.keepButton.grid(row=1,column=0,columnspan=2)
		self.clearButton = Button(self.plotFrame,text='Clear past',width=10,command=self.clearOldResults, state=DISABLED)
		self.clearButton.grid(row=1,column=2,columnspan=2)
		self.theoLinesButton = Button(self.plotFrame,text='Add Theo Line',width=10,command=self.addTheoLine)
		self.theoLinesButton.grid(row = 2,column = 0)
	def build_operators_bools(self):
		#To add an operator use the my_check_button function, copy paste the following line and change the appropriate parameters. Then you have to define the operator in the main module, here basically you can just add a checkButton which keep tracks of the preferences.
		self.operatorsFrame = LabelFrame(self, text='Operators', padx=5, pady=5, relief= RIDGE)
		self.operatorsFrame.grid(row = 5, column = 5, columnspan =3,sticky=W)
		self.SzBool = self.my_check_button(text='Sz',pos=(0,0),master=self.operatorsFrame, select=True)
		self.IzBool = self.my_check_button(text='Iz',pos=(0,2),master=self.operatorsFrame, select=True)
		self.M_stg_bool = self.my_check_button(text='M_stg',pos=(1,0),master=self.operatorsFrame, select=False)
		self.HBool = self.my_check_button(text='H',pos=(1,2),master=self.operatorsFrame, select=False)
		self.I2Bool = self.my_check_button(text='I2',pos=(2,0),master=self.operatorsFrame, select=False)
		self.S2Bool = self.my_check_button(text='S2',pos=(2,2),master=self.operatorsFrame, select=False)
		self.CohBool = self.my_check_button(text='I+',pos=(3,0),master=self.operatorsFrame, select=False)
		self.IzNormBool = self.my_check_button(text='Iz/I',pos=(3,2),master=self.operatorsFrame, select=True)	
	def build_initial_state_frame(self):
		##labelframe for the initial state. It looks messy, it's because of some triggered command to make it fancy.		
		self.initialFrame = LabelFrame(self, text='Initial Hamiltonian', padx=5, pady=5, relief=RIDGE)
		self.initialFrame.grid(row = 1, column = 1, sticky=W+E+N)
		
		Label(self.initialFrame,text= 'K: # of nuclei').grid(row=0,column=0,columnspan = 2, sticky=W)
		self.K_entry = Spinbox(self.initialFrame, width=3, from_=1, to=40)
		self.K_entry.delete(0, END)
		defaultK = 10
		self.K_entry.insert(0,defaultK)
		self.K_entry.grid(row=0,column=2)

		self.tau_dot_entry = self.my_entry('τd: temperature of the dot',(1,0),default=0, master=self.initialFrame)
		self.rho_init_entry = self.my_entry('Kρ: electron / nucleus Zeeman',(2,0),default=15.0, master= self.initialFrame)
		self.nu_init_entry = self.my_entry('ν: nuclei Zeeman',(3,0),default=0.3, master= self.initialFrame)
		
	def build_evolution_frame(self):
		##physical parameters on the left, with labelframe for the evolution
		#To add an entry use the my_entry function		
		self.evolutionFrame = LabelFrame(self, text='Evolution Hamiltonian', padx=5, pady=5, relief= RIDGE)
		self.evolutionFrame.grid(row = 2, column = 1, sticky=W+E)
		self.rho_evol_entry = self.my_entry('Kρ: electron / nucleus Zeeman',(0,0),default=15.0, master= self.evolutionFrame)
		self.nu_evol_entry = self.my_entry('ν: nuclei Zeeman',(1,0),default=2.0, master= self.evolutionFrame)
		
		self.alfa_entry = self.my_entry('α: magnetic field gradient',(2,0),default=0, master= self.evolutionFrame)
		self.labda_entry = self.my_entry('λ: flip-flop term',(3,0),default=1, master= self.evolutionFrame)
		
	def build_dissipator_frame(self):
		self.dissipatorFrame = LabelFrame(self, text='Dissipating Mechanism', padx=5, pady=5, relief= RIDGE)
		self.dissipatorFrame.grid(row = 3, column = 1, sticky=W+E)
		self.gamma0_entry = self.my_entry('Γ0: base transition rate',(0,0),default=0, master= self.dissipatorFrame)
		self.mu_entry = self.my_entry('μ: chemical potential of the lead',(1,0),default=100, master= self.dissipatorFrame)
		self.tau_lead_entry = self.my_entry('τl: temperature of the lead',(2,0),default=0, master= self.dissipatorFrame)
		

	def build_phys_parameters_frame(self):
		##labelframe for the initial state. It looks messy, it's because of some triggered command to make it fancy.		
		self.physParametersFrame = LabelFrame(self, text='Real parameters', padx=5, pady=5, relief=RIDGE)
		self.physParametersFrame.grid(row = 5, column = 1, sticky=W+E+N)
		self.B_init_entry = self.my_entry('B0: initial magnetic field [T]',(1,0),default=-0.05, master= self.physParametersFrame)
		self.B_evol_entry = self.my_entry('B1: evolution magnetic field [T]',(2,0),default=0.05, master= self.physParametersFrame)
		self.T_dot_entry = self.my_entry('Td: temperature of the dot [K]',(3,0),default=0, master=self.physParametersFrame)
		self.T_lead_entry = self.my_entry('Tl: temperature of the lead [K]',(4,0),default=0, master=self.physParametersFrame)
		self.gE_entry = self.my_entry('g*: electron effective g-factor',(5,0),default=2, master= self.physParametersFrame)
		self.gN_entry = self.my_entry('gN: nuclear g-factor',(6,0),default=2, master= self.physParametersFrame)
		self.Ahf_entry = self.my_entry('Ahf: total hyperfine coupling [μeV]',(7,0),default=2, master= self.physParametersFrame)
		self.defaultMaterial = StringVar()
		rb = Radiobutton(self.physParametersFrame, text='GaAs', variable=self.defaultMaterial, value='GaAs', command = self.update_material)
		rb.grid(row = 0, column = 0)
		
		Radiobutton(self.physParametersFrame, text='InAs', variable=self.defaultMaterial, value='InAs', command = self.update_material).grid(row = 0, column = 1)
		rb.select()
		#rb.invoke()
	
	def update_material(self):
		if self.defaultMaterial.get() == 'GaAs':
			self.gE_entry.delete(0, END)
			self.gE_entry.insert(0,-0.44)
			self.gN_entry.delete(0, END)
			self.gN_entry.insert(0,1.84)
			self.Ahf_entry.delete(0, END)
			self.Ahf_entry.insert(0,84)
		elif self.defaultMaterial.get() == 'InAs':
			self.gE_entry.delete(0, END)
			self.gE_entry.insert(0,-14.9)
			self.gN_entry.delete(0, END)
			self.gN_entry.insert(0,3.48)
			self.Ahf_entry.delete(0, END)
			self.Ahf_entry.insert(0,98)
		self.update_rescaled_parameters(event=None)
	def update_rescaled_parameters(self,event):	
		MUN = 3.15 * 10**(-8)
		MUB = 5.79 * 10**(-5)
		KB = 8.617 * 10**(-5)
		G = 2

		K = int(self.K_entry.get())
		gn = float(self.gN_entry.get())
		ge = float(self.gE_entry.get())
		Binit = float( self.B_init_entry.get())
		Bevol = float( self.B_evol_entry.get())
		Tdot = float(self.T_dot_entry.get())
		Tlead = float(self.T_lead_entry.get())
		Ahf = float(self.Ahf_entry.get())
 


		self.tau_dot_entry.delete(0, END)
		self.tau_lead_entry.delete(0, END)
		self.rho_init_entry.delete(0, END)
		self.nu_init_entry.delete(0, END)
		self.rho_evol_entry.delete(0, END)
		self.nu_evol_entry.delete(0, END)
		
		self.tau_dot_entry.insert(0,KB * Tdot /  (Ahf / (2. * K)))
		self.tau_lead_entry.insert(0,KB * Tlead /  (Ahf / (2. * K)))
		self.rho_init_entry.insert(0,math.fabs(ge) * MUB  / (gn * MUN ))
		self.nu_init_entry.insert(0,gn * MUN * Binit / (Ahf / (2. * K)))
		self.rho_evol_entry.insert(0,math.fabs(ge) * MUB  / (gn * MUN))
		self.nu_evol_entry.insert(0,gn * MUN * Bevol / (Ahf / (2. * K)))


	
	def keepResults(self):
		if (self.dynamics_result,self.legends) not in self.oldResults:
			self.oldResults.append((self.dynamics_result,self.legends,self.saveText))
		self.keepButton['state'] = DISABLED
		self.clearButton['state'] = NORMAL
	def clearOldResults(self):
		self.oldResults = []	
		self.saveText = ''
		self.clearButton['state'] = DISABLED
		self.keepButton['state'] = NORMAL
	def save(self):
		self.saveText = '''

		\\begin{tabular}{ccc}
			&	Initial Hamiltonian	&	Evol Hamiltonian	\\\\ \\hline
		$\\rho$ &%g	&	%g\\\\
		$\\nu$	&%g 	&%g	\\\\
		$\\lambda$&1&%g\\\\
		$\\alpha$&0&%g
		\\end{tabular}\\vline
		\\begin{tabular}{cc}
		$K$ &%g\\\\
		$\\tau_{dot}$ &%g\\\\
		$\\tau_{lead}$	&%g\\\\
		$\\Gamma_0$&%g\\\\
		$\\mu$&%g
		\\end{tabular}\\vline

		\\begin{tabular}{cc}
		$t_{evol}$ &%g\\\\
		$t_{diss}$	&%g\\\\
		$t_{end}$&%g\\\\
		$t_{steps}$&%g
		\\end{tabular}
		'''%(self.parameters.rho_init, self.parameters.rho_evol, self.parameters.nu_init, self.parameters.nu_evol, self.parameters.labda, self.parameters.alfa, self.parameters.K, self.parameters.tau_dot, 
			 self.parameters.tau_lead, self.parameters.gamma0, self.parameters.mu, self.parameters.evo_time, self.parameters.dis_time, self.parameters.tInterval, self.parameters.tSteps)
		def saveAll():
			dirName = saveName.get().strip()
			
			if not dirName:
				dirName = 'QuTinkerPlot'
			print 'Saving in %s'%dirName
			if not os.path.exists('archive/'+dirName):
	    			os.makedirs('archive/'+dirName)
			self.figure.savefig('archive/'+dirName + '/figure.pdf')
			self.figure.savefig('archive/'+dirName + '/figure.png')
			self.figure.savefig('archive/'+dirName + '/figure.jpg')
			
			with open('archive/'+dirName+'/'+dirName+'.tex','w') as parFile:
				parFile.write(incipit)
				for oR in self.oldResults:
					parFile.write( oR[2])
				parFile.write(self.saveText)
				parFile.write(figure)
				parFile.write(excipit)
			with open('archive/'+dirName+'/rawData.txt','w') as dataFile:
				for oR in self.oldResults:
					for r in oR[0].expect:
						for n in r:
							dataFile.write(str(n))
							dataFile.write(' ')
						dataFile.write('\n')
				for r in self.dynamics_result.expect:
					for n in r:
						dataFile.write(str(n))
						dataFile.write(' ')
					dataFile.write('\n')
				
			top.destroy()
		
		incipit = '''
	\\documentclass[10pt,a4paper]{article}
	\\usepackage[utf8]{inputenc}
	\\usepackage{amsmath}
	\\usepackage{amsfonts}
	\\usepackage{amssymb}
	\\usepackage{graphicx}
	\\begin{document}
	\\begin{equation*}
		H = K \\rho \\nu S^z - \\sum\\limits_{l=1}^K\\left(\\nu+\\frac{\\alpha}{2} \\frac{2l-(K+1)}{K-1}\\right) I_l^z+ 2 S^z I^z + \\lambda (S^+ I^-+ S^- I^+)
		\\end{equation*}'''
		figure = '''
	\\begin{figure}[h]
		\\includegraphics[scale=0.7]{figure.pdf}
		\\end{figure}'''
		excipit = '''
	\\end{document}
	'''
		top = Toplevel()
		top.geometry("300x100+0+0")
		top.title('Save the figure')
		Frame(top, height=20).grid(row=0, column=0, columnspan=3)
		Frame(top, width=30).grid(row=0, column=0, rowspan=3)
		saveName = Entry(top, width=24)
		saveName.insert(0,'QuTinkerPlot')
		saveName.grid(row=1, column=1,columnspan=2 ,sticky=W+E)
		Button(top, text='Close',width=12, command= top.destroy).grid(row=2, column=2, sticky=W+E)
		Button(top, text='Ok',width=12, command= saveAll).grid(row=2, column=1, sticky=W+E)
		
		
	def addTheoLine(self):
		self.legends.append('$Coherent \\;\\left<m_I\\right>$')
		self.legends.append('$Incoherent \\;\\left<m_I\\right>$')
		#print self.parameters
		I,m0,gamma = self.parameters.K/2., -self.parameters.K/2., (self.parameters.gamma0**2/(2*math.pi*(self.parameters.mu-self.computation.meanEnergyEvol)**2))*(self.parameters.nu_evol-1)/(self.parameters.nu_evol**2*(self.parameters.K*self.parameters.rho_evol-1)**2)
		print 'gamma = ',gamma
		tlist = np.linspace(0,self.parameters.tInterval,num=self.parameters.tSteps)
		theoCoherent = [I - (2*I+1)*(I-m0)/((I+m0+1)*math.exp(gamma*(2*I+1)*t)+I-m0) if gamma*(2*I+1)*t < 15 else I for t in tlist]
		theoIncoherent = [0.5 * self.parameters.K * (1-2*math.exp(-gamma*t)) if gamma*t < 15 else 0.5 * self.parameters.K for t in tlist]
		#print theoIncoherent
		self.ax.plot(tlist * gamma, theoCoherent,'-')
		self.ax.plot(tlist * gamma, theoIncoherent,'-')
		self.canvas.show()

	def plot(self):
		self.parameters = Parameters(self)
		self.figure.clf()
		self.dynamics_ax = self.figure.add_subplot(121)
		self.pdiagram_ax = self.figure.add_subplot(122)

		def show_dynamics():
			self.dynamics_result = self.computation.compute_dynamics(self.parameters)
		
			if self.dynamics_result == 'Error':
				top = Toplevel()
				top.geometry("400x100+0+0")
				top.title('Error in the OdeSolver')
				Label(top, text = 'Try change the options').pack()
				Button(top, text='Close', command= top.destroy).pack()
				return

		

			self.saveButton['state'] = NORMAL
			self.keepButton['state'] = NORMAL
			self.clearButton['state'] = NORMAL
			self.legends = []
		
			if self.LegendManual.get():
				self.legends.extend([l.strip() for l in self.LegendManual.get().split(',')])
			else: 
				for oR in self.oldResults:
					self.legends += oR[1]
		
				if self.IzBool.get():
					self.legends.append('$I^z$')
					#self.ax.set_ylabel('$\left< I_z(t)\\right>$')
				if self.SzBool.get():
					self.legends.append('$S^z$')
				if self.M_stg_bool.get():
					self.legends.append('$M_{stg}$')
				if self.HBool.get():
					self.legends.append('$H_{in}$')
					self.legends.append('$H_{ev}$')
				if self.I2Bool.get():
					self.legends.append('$I^2$')
				if self.S2Bool.get():
					self.legends.append('$S^2$')
				if self.CohBool.get():
					self.legends.append('$I+$')
					self.legends.append('$I+I-$')
				if self.IzNormBool.get():
					self.legends.append('$I^z/I$')
			
					
		
			self.dynamics_ax.set_aspect('auto')
			gamma = (self.parameters.gamma0**2/(2*math.pi*(self.parameters.mu-self.computation.meanEnergyEvol)**2))*(self.parameters.nu_evol-1)/(self.parameters.nu_evol**2*(self.parameters.K*self.parameters.rho_evol-1)**2)
		


			print 'gamma ', gamma
			for oR in self.oldResults:
				for r in oR[0].expect:
					self.dynamics_ax.plot(oR[0].times * gamma, r)
			if self.dynamics_result:
				for r in self.dynamics_result.expect:
					self.dynamics_ax.plot(self.dynamics_result.times * gamma, r)
		
				
	
			if self.LegendBool.get():
				self.dynamics_ax.legend(self.legends)
			#if not self.computation.isHConserved:
			#	self.dynamics_ax.axvline(linewidth=2, color='r',x = self.parameters.evo_time, ls='--')
			#if self.parameters.gamma0:
			#	self.dynamics_ax.axvline(linewidth=2, color='b',x = self.parameters.dis_time, ls='--')
			#self.dynamics_ax.axhline(linewidth=4)        
			#self.dynamics_ax.axvline(linewidth=4)  
			self.dynamics_ax.xaxis.set_tick_params(labelsize='large')
			self.dynamics_ax.yaxis.set_tick_params(labelsize='large')
			self.dynamics_ax.set_xlabel('$\eta \, t$')
			self.dynamics_ax.set_xlim([0,max(self.dynamics_result.times) * gamma])
			self.dynamics_ax.set_ylim([-1,1])
		
		def show_phase_diagram():
			self.pdiagram_result = self.computation.compute_phase_diagram(self.parameters)
			
		
			y, x, z = self.pdiagram_result
		
			c = self.pdiagram_ax.imshow(z, extent=(np.amin(x), np.amax(x),np.amin(y), np.amax(y) ))
		
			line_xmin, line_ymin = min(self.parameters.nu_init,self.parameters.nu_evol), min(self.parameters.rho_init,self.parameters.rho_evol)
			line_xmax, line_ymax = max(self.parameters.nu_init,self.parameters.nu_evol), max(self.parameters.rho_init,self.parameters.rho_evol)
			self.pdiagram_ax.plot([line_xmin, line_xmax],  [line_ymin, line_ymax ],'k-', lw=2)
			self.pdiagram_ax.set_aspect('auto')
			self.figure.colorbar(c, ticks=[-self.parameters.K/2., 0, self.parameters.K/2.], orientation='vertical')


		show_phase_diagram()
		self.canvas.show()	


















	
								
