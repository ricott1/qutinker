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
		self.geometry("1520x750+0+0")



class Parameters(object):
	def __init__(self, master):
		self.master = master
		self.get_parameters()
	def get_parameters(self):
		self.N = max(1, int(self.master.N_entry.get()))
		self.tau_dot = float(self.master.tau_dot_entry.get())
		self.tau_lead = max(float(self.master.tau_lead_entry.get()), 0)
		self.gammaI_init = float(self.master.gammaI_init_entry.get())
		self.b_init = float(self.master.b_init_entry.get())
		self.gammaI_final = float(self.master.gammaI_final_entry.get())
		self.b_final = float(self.master.b_final_entry.get())
		self.alfa = float(self.master.alfa_entry.get())
		self.eta = max(float(self.master.eta_entry.get()), 0)
		self.mu = max(float(self.master.mu_entry.get()), 0)
		self.tInterval = float(self.master.tInterval_entry.get())
		self.tSteps = int(self.master.tSteps_entry.get())
		self.quench_v = float(self.master.quench_v_entry.get())
		self.startTime = 0
		self.dis_time = float(self.master.dis_time_entry.get())
		self.evo_time = float(self.master.evo_time_entry.get())
		self.dissipators_kind = self.master.dissipators_kind
		self.cycles = max(1, int(self.master.cycles_entry.get()))
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
		self.build_finalution_frame()
		self.build_dissipator_frame()
		self.build_time_schedule_frame()
		self.build_plot_parameters_frame()
		self.build_operators_bools()
		self.buildMatPlotLibFrame()
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
		Frame(self,height=500,width = 1000).grid(row=1,column=2, columnspan=6,rowspan=5, sticky=W+E+N+S)
		#we initialize the figure, setting also the title. There are a couple of tricks to make matplotlib not popup a widget automatically and embed it instead on tkinter
		self.figure = mplfig.Figure(figsize=(14,6))
		#hamiltonian = '$H = K \\rho \\nu S^z - \\nu I^z + 2 S^z I^z + \lambda (S^+ I^-+ S^- I^+)$'
		#self.figure.suptitle(hamiltonian, fontsize=24)
		
		self.canvas = tkagg.FigureCanvasTkAgg(self.figure, master=self)
		
		self.canvas.get_tk_widget().grid(row=1,column=2, columnspan=6,rowspan=5, sticky=W+E+N+S)
		self.canvas._tkcanvas.grid()
	def build_time_schedule_frame(self):
		self.timeScheduleFrame = LabelFrame(self, text='Time Schedule', padx=5, pady=5, relief= RIDGE)
		self.timeScheduleFrame.grid(row = 4, column = 1, sticky=W+E)
		self.evo_time_entry = self.my_entry('t evol: evolution starting time',(0,0),default=0, master= self.timeScheduleFrame)
		self.dis_time_entry = self.my_entry('t diss: dissipation starting time',(1,0),default=0, master= self.timeScheduleFrame)
		self.tInterval_entry = self.my_entry('Ending time',(2,0), default=30, master=self.timeScheduleFrame)
		self.tSteps_entry = self.my_entry('Time steps',(3,0), default=300, master=self.timeScheduleFrame)
		self.quench_v_entry = self.my_entry('Quench V',(4,0), default=0.1, master=self.timeScheduleFrame)

	def build_plot_parameters_frame(self):
		##labelFrame for the plot parameters
		self.plotFrame = LabelFrame(self, text='Plot parameters', padx=5, pady=5, relief= RIDGE)
		self.plotFrame.grid(row = 5, column = 1,sticky=W+E)	
		self.LegendBool = self.my_check_button(text='Legends',pos=(0,0),master=self.plotFrame, select=True)
		self.reverse_dynamics_bool = self.my_check_button(text='Reverse',pos=(4,0),master=self.plotFrame, select=False)
		self.LegendManual = Entry(self.plotFrame, width=20)
		self.LegendManual.grid(row=0,column=2, columnspan = 4)
		self.keepButton = Button(self.plotFrame,text='Keep it!',width=10,command=self.keepResults, state=DISABLED)
		self.keepButton.grid(row=1,column=0)
		self.clearButton = Button(self.plotFrame,text='Clear past',width=10,command=self.clearOldResults, state=DISABLED)
		self.clearButton.grid(row=1,column=1)
		self.theoLinesButton = Button(self.plotFrame,text='Add Theo Line',width=10,command=self.addTheoLine)
		self.theoLinesButton.grid(row = 1,column = 2)
		self.energy_levels_bool = self.my_check_button(text='EL',pos=(2,0),master=self.plotFrame, select=False)
		self.phase_diagram_bool = self.my_check_button(text='PD',pos=(2,2),master=self.plotFrame, select=False)
		self.dynamics_bool = self.my_check_button(text='DY',pos=(3,0),master=self.plotFrame, select=True)
		self.probabilities_bool = self.my_check_button(text='DM',pos=(3,2),master=self.plotFrame, select=False)
	def build_operators_bools(self):
		#To add an operator use the my_check_button function, copy paste the following line and change the appropriate parameters. Then you have to define the operator in the main module, here basically you can just add a checkButton which keep tracks of the preferences.
		self.operatorsFrame = LabelFrame(self, text='Operators', padx=5, pady=5, relief= RIDGE)
		self.operatorsFrame.grid(row = 6, column = 1,sticky=W+E)
		self.SzBool = self.my_check_button(text='Sz',pos=(0,0),master=self.operatorsFrame, select=False)
		self.IzBool = self.my_check_button(text='Iz',pos=(0,2),master=self.operatorsFrame, select=False)
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
		
		Label(self.initialFrame,text= 'N: # of nuclei').grid(row=0,column=0,columnspan = 2, sticky=W)
		self.N_entry = Spinbox(self.initialFrame, width=3, from_=1, to=40)
		self.N_entry.delete(0, END)
		defaultN = 4
		self.N_entry.insert(0,defaultN)
		self.N_entry.grid(row=0,column=2)

		self.tau_dot_entry = self.my_entry('τd: temperature of the dot',(1,0),default=0, master=self.initialFrame)
		self.gammaI_init_entry = self.my_entry('γ/2I: electron / nucleus Zeeman',(2,0),default=2, master= self.initialFrame)
		self.b_init_entry = self.my_entry('b: nuclei Zeeman',(3,0),default=2, master= self.initialFrame)
		
	def build_finalution_frame(self):
		##physical parameters on the left, with labelframe for the evolution
		#To add an entry use the my_entry function		
		self.evolutionFrame = LabelFrame(self, text='Evolution Hamiltonian', padx=5, pady=5, relief= RIDGE)
		self.evolutionFrame.grid(row = 2, column = 1, sticky=W+E)
		self.gammaI_final_entry = self.my_entry('γ/2I: electron / nucleus Zeeman',(0,0),default=2, master= self.evolutionFrame)
		self.b_final_entry = self.my_entry('b: nuclei Zeeman',(1,0),default=0.3, master= self.evolutionFrame)
		
		self.alfa_entry = self.my_entry('α: flip-flop term',(3,0),default=1, master= self.evolutionFrame)
		
	def build_dissipator_frame(self):
		self.dissipatorFrame = LabelFrame(self, text='Dissipating Mechanism', padx=5, pady=5, relief= RIDGE)
		self.dissipatorFrame.grid(row = 3, column = 1, sticky=W+E)
		self.eta_entry = self.my_entry('η: base transition rate',(0,0),default=0, master= self.dissipatorFrame)
		
		self.mu_entry = self.my_entry('μ: chemical potential of the lead',(1,0),default=100, master= self.dissipatorFrame)
		self.tau_lead_entry = self.my_entry('τl: temperature of the lead',(2,0),default=0, master= self.dissipatorFrame)
		self.cycles_entry = self.my_entry('C: number of cycles',(3,0),default=1, master= self.dissipatorFrame)
		
	
		self.dissipators_kind = StringVar()

		a = Radiobutton(self.dissipatorFrame, text='Fermi', variable=self.dissipators_kind, value='fermi')
		a.grid(row=4, column=0)
		a.select()
		b = Radiobutton(self.dissipatorFrame, text='Central', variable=self.dissipators_kind, value='central')
		b.grid(row=4, column=2)
		b.deselect()

	def keepResults(self):
		self.oldResults += [(self.dynamics_result,self.legends,self.saveText)]


		self.keepButton['state'] = DISABLED
		self.clearButton['state'] = NORMAL
	def clearOldResults(self):
		self.oldResults = []	
		self.saveText = ''
		self.clearButton['state'] = DISABLED
		self.keepButton['state'] = NORMAL
	def save(self):
		
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
			
			top.destroy()
		
		
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
		I,m0,gamma = self.parameters.K/2., -self.parameters.K/2., (self.parameters.eta**2/(2*math.pi*(self.parameters.mu-self.computation.mean_energy_final)**2))*(self.parameters.b_final-1)/(self.parameters.b_final**2*(self.parameters.K*self.parameters.gammaI_final-1)**2)
		print 'baserate eta = ',gamma
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
		
		def show_dynamics():
			self.dynamics_result = self.computation.compute_dynamics()
			if self.dynamics_result == 'Error':
				top = Toplevel()
				top.geometry("400x100+0+0")
				top.title('Error in the OdeSolver')
				Label(top, text = 'Try change the options').pack()
				Button(top, text='Close', command= top.destroy).pack()
				return

		

			
		
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
			
		
			for oR in self.oldResults:
				for r in oR[0].expect:
					if self.reverse_dynamics_bool.get():
						split_tlist = np.array_split(self.dynamics_result.times, self.parameters.cycles)
						
						split_result = np.array_split(r, self.parameters.cycles)
						self.dynamics_ax.set_xlim([0,max(split_tlist[0]) * gamma])
						for i, result in enumerate(split_result):
							if i%2:
								self.dynamics_ax.plot(split_tlist[0] * gamma, list(reversed(result)))
							else:
								self.dynamics_ax.plot(split_tlist[0] * gamma, result)
								
					else:
						self.dynamics_ax.plot(self.dynamics_result.times * gamma, r)
						self.dynamics_ax.set_xlim([0,max(self.dynamics_result.times) * gamma])
			if self.dynamics_result:
				for r in self.dynamics_result.expect:
					if self.reverse_dynamics_bool.get():
						split_tlist = np.array_split(self.dynamics_result.times, self.parameters.cycles)
						
						split_result = np.array_split(r, self.parameters.cycles)
						self.dynamics_ax.set_xlim([0,max(split_tlist[0]) * gamma])
						for i, result in enumerate(split_result):
							if i%2:
								self.dynamics_ax.plot(split_tlist[0] * gamma, list(reversed(result)))
							else:
								self.dynamics_ax.plot(split_tlist[0] * gamma, result)
						#for t in xrange(1, 4):
						#	self.dynamics_ax.axvline(linewidth=2, color='b',x = max(split_tlist[0]) * gamma * t/4., ls='--') 	
					else:
						self.dynamics_ax.plot(self.dynamics_result.times * gamma, r)
						self.dynamics_ax.set_xlim([0,max(self.dynamics_result.times) * gamma])
					
				
				
	
			if self.LegendBool.get():
				self.dynamics_ax.legend(self.legends)
			#if not self.computation.isHConserved:
			#	self.dynamics_ax.axvline(linewidth=2, color='r',x = self.parameters.evo_time, ls='--')
			#if self.parameters.eta:
			#	self.dynamics_ax.axvline(linewidth=2, color='b',x = self.parameters.dis_time, ls='--')
			#self.dynamics_ax.axhline(linewidth=4)        
			#self.dynamics_ax.axvline(linewidth=4)  
			self.dynamics_ax.xaxis.set_tick_params(labelsize='large')
			self.dynamics_ax.yaxis.set_tick_params(labelsize='large')
			self.dynamics_ax.set_xlabel('$\eta \, t$')
			
			self.dynamics_ax.set_ylim([-1.2,1.2])
		
		def show_phase_diagram():
			self.pdiagram_result = self.computation.compute_phase_diagram()
			
		
			y, x, z = self.pdiagram_result
		
			c = self.pdiagram_ax.imshow(z, extent=(np.amin(x), np.amax(x),np.amin(y), np.amax(y) ))
		
			line_xmin, line_ymin = min(self.parameters.b_init,self.parameters.b_final), min(self.parameters.gammaI_init,self.parameters.gammaI_final)
			line_xmax, line_ymax = max(self.parameters.b_init,self.parameters.b_final), max(self.parameters.gammaI_init,self.parameters.gammaI_final)
			self.pdiagram_ax.plot([line_xmin, line_xmax],  [line_ymin, line_ymax ],'k-', lw=2)
			self.pdiagram_ax.set_aspect('auto')
			self.figure.colorbar(c, ticks=[-self.parameters.N/2., 0, self.parameters.N/2.], orientation='vertical', ax = self.pdiagram_ax)
		def show_density_matrix():
			times, pg = self.computation.compute_density_matrix()
			
			for p in pg:
				self.density_matrix_ax.plot(times, p)
		def show_energy_levels():
			times, energy_levels = self.computation.compute_energy_levels()
			for el in energy_levels:
				self.levels_ax.plot(times , el)
			#self.dynamics_ax.set_ylim([-10,0])
		self.dynamics_ax = self.figure.add_subplot(221)
		self.dynamics_ax.set_title('Dynamics')
		self.pdiagram_ax = self.figure.add_subplot(224)
		self.pdiagram_ax.set_title('Phase Diagram')
		self.levels_ax = self.figure.add_subplot(223)
		self.levels_ax.set_title('Energy levels')
		self.density_matrix_ax = self.figure.add_subplot(222)
		self.density_matrix_ax.set_title('Probabilities')
		self.computation.prepare_computation(self.parameters)
		print  (self.parameters.eta**2/(2*math.pi*(self.parameters.mu-self.computation.mean_energy_final)**2))
		print  (self.parameters.eta**2/(2*math.pi*(self.parameters.mu-self.computation.mean_energy_init)**2))
		gamma = (self.parameters.eta**2/(2*math.pi*(self.parameters.mu-self.computation.mean_energy_final)**2))#*(self.parameters.b_final-1)/(self.parameters.b_final**2*(self.parameters.N*self.parameters.gammaI_final-1)**2)	

		self.saveButton['state'] = NORMAL
		self.keepButton['state'] = NORMAL
		self.clearButton['state'] = NORMAL
		self.legends = []
		if self.energy_levels_bool.get():
			show_energy_levels()
		if self.phase_diagram_bool.get():
			show_phase_diagram()
		if self.dynamics_bool.get():
			show_dynamics()
		if self.probabilities_bool.get():
			show_density_matrix()
		
		self.canvas.show()	


















	
								
