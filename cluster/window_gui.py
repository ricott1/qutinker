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
from parameters import *
class Root(Tk):
    def __init__(self):
        
        Tk.__init__(self)
        self.wm_title('QuTinker')
        self.geometry("1520x750+0+0")



class MainWindow(Frame):
    '''
    This class contains the main widget as a frame which we fill up with fancy buttons and entry fields. 
    There are many exhaustive documentations of Tkinter online, e.g. http://effbot.org/tkinterbook/tkinter-index.htm
    '''
    def __init__(self, parent):
        Frame.__init__(self, parent)
        self.computation = Computation(self)
        self.dynamics_result = []
        self.old_results = []
        self.build_borders()
        self.build_state_1_frame()
        self.build_state_2_frame()
        self.build_ising_frame()
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
        self.wait_time_entry = self.my_entry('t: waiting time at point',(1,0),default=20, master= self.timeScheduleFrame)
        self.tInterval_entry = self.my_entry('Ending time',(2,0), default=12, master=self.timeScheduleFrame)
        self.tSteps_entry = self.my_entry('Time steps',(3,0), default=20, master=self.timeScheduleFrame)
        self.annealing_v_entry = self.my_entry('Annealing V',(4,0), default=0.1, master=self.timeScheduleFrame)
        self.tau_dot_entry = self.my_entry('τd: temperature of the system',(5,0),default=0, master=self.timeScheduleFrame)
        self.tau_lead_entry = self.my_entry('τl: temperature of the lead',(6,0),default=0, master= self.timeScheduleFrame)
        
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
        self.clearButton = Button(self.plotFrame,text='Clear past',width=10,command=self.clearold_results, state=DISABLED)
        self.clearButton.grid(row=1,column=1)
        
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
        
    def build_state_1_frame(self):    
        self.state1_frame = LabelFrame(self, text='Hamiltonian 1', padx=5, pady=5, relief=RIDGE)
        self.state1_frame.grid(row = 1, column = 1, sticky=W+E+N)
        Label(self.state1_frame,text= 'N1: # of bath spins').grid(row=0,column=0,columnspan = 2, sticky=W)
        self.N1_entry = Spinbox(self.state1_frame, width=3, from_=1, to=40)
        self.N1_entry.delete(0, END)
        defaultN = 4
        self.N1_entry.insert(0,defaultN)
        self.N1_entry.grid(row=0,column=2)
        self.gamma1_entry = self.my_entry('γ1: central spin Zeeman',(2,0),default=2, master= self.state1_frame)
        self.b1_entry = self.my_entry('b1: bath spins Zeeman',(3,0),default=2, master= self.state1_frame)
        self.alfa1_entry = self.my_entry('α1: flip-flop term',(4,0),default=1, master= self.state1_frame)
        
        self.eta1_entry = self.my_entry('η1: base transition rate',(5,0),default=10, master= self.state1_frame)
    
        self.dissipators_kind = StringVar()

        a = Radiobutton(self.state1_frame, text='Fermi', variable=self.dissipators_kind, value='fermi')
        a.grid(row=6, column=0)
        a.deselect()
        b = Radiobutton(self.state1_frame, text='Central', variable=self.dissipators_kind, value='central')
        b.grid(row=6, column=2)
        b.select()
        
    def build_state_2_frame(self):    
        self.state2_frame = LabelFrame(self, text='Hamiltonian 2', padx=5, pady=5, relief=RIDGE)
        self.state2_frame.grid(row = 2, column = 1, sticky=W+E+N)
        Label(self.state2_frame,text= 'N2: # of bath spins').grid(row=0,column=0,columnspan = 2, sticky=W)
        self.N2_entry = Spinbox(self.state2_frame, width=3, from_=1, to=40)
        self.N2_entry.delete(0, END)
        defaultN = 4
        self.N2_entry.insert(0,defaultN)
        self.N2_entry.grid(row=0,column=2)
        self.gamma2_entry = self.my_entry('γ2: central spin Zeeman',(2,0),default=2, master= self.state2_frame)
        self.b2_entry = self.my_entry('b2: bath spins Zeeman',(3,0),default=2, master= self.state2_frame)
        self.alfa2_entry = self.my_entry('α2: flip-flop term',(4,0),default=1, master= self.state2_frame)
        
        self.eta2_entry = self.my_entry('η2: base transition rate',(5,0),default=10, master= self.state2_frame)
    
        self.dissipators_kind = StringVar()

        a = Radiobutton(self.state2_frame, text='Fermi', variable=self.dissipators_kind, value='fermi')
        a.grid(row=6, column=0)
        a.deselect()
        b = Radiobutton(self.state2_frame, text='Central', variable=self.dissipators_kind, value='central')
        b.grid(row=6, column=2)
        b.select()
        
    def build_ising_frame(self):
        self.ising_frame = LabelFrame(self, text='Ising and Transverse Hamiltonian', padx=5, pady=5, relief=RIDGE)
        self.ising_frame.grid(row = 3, column = 1, sticky=W+E+N)
        self.J_entry = self.my_entry('J: Ising coupling',(1,0),default=10, master= self.ising_frame)
        self.Delta_entry = self.my_entry('Delta: transverse field',(26,0),default=10, master= self.ising_frame)        

    def keepResults(self):
        result = {'dynamics' : self.dynamics_result, 'legends' : self.legends}
        self.old_results.append(result)
        self.keepButton['state'] = DISABLED
        self.clearButton['state'] = NORMAL
    def clearold_results(self):
        self.old_results = []
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
            result = {'dynamics' : self.dynamics_result, 'legends' : self.legends}
            i = 0
            for r in self.old_results + [result]:
                
                with open('archive/'+dirName +'/hysteresis_data{}'.format(i), 'w+') as f:
                    lines = []
                    lines = zip(r['hysteresis']['b'],r['hysteresis']['gamma'], r['hysteresis']['expect'])
                    for l in lines:
                        f.write(' '.join([str(x) for x in l]) + '\n')
                eta = 1#self.parameters.eta if  self.parameters.eta > 0 else 1.
                with open('archive/'+dirName +'/dynamics_data{}'.format(i), 'w+') as f:
                    for exp in r['dynamics'].expect:
                        lines = zip(r['dynamics'].times * eta, exp)
                        for l in lines:
                            f.write(' '.join([str(x) for x in l]) + '\n')
                    
                
                i += 1
            
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


    def plot(self):
        self.parameters = Parameters()
        self.parameters.get_parameters(self)
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
                for oR in self.old_results:
                    self.legends += oR['legends']
        
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
            
                    
        
            COLORS = ['r','r','b','b','g','g','y','y']
            STYLES = ['-','--','-','--','-','--','-','--']
            i = 0
            
            for oR in self.old_results:
                for r in oR['dynamics'].expect:
                    if self.reverse_dynamics_bool.get():
                        split_tlist = np.array_split(oR['dynamics'].times, self.parameters.cycles)
                        
                        split_result = np.array_split(r, self.parameters.cycles)
                        self.dynamics_ax.set_xlim([0,max(split_tlist[0]) * eta])
                        for i, result in enumerate(split_result):
                            if i%2:
                                self.dynamics_ax.plot(split_tlist[0] * eta, list(reversed(result)))
                            else:
                                self.dynamics_ax.plot(split_tlist[0] * eta, result)
                                
                    else:
                        self.dynamics_ax.plot(oR['dynamics'].times * eta, r)
                        self.dynamics_ax.set_xlim([0,max(self.dynamics_result.times + oR['dynamics'].times) * eta])
                
                
                i += 1
            if self.dynamics_result:
                                
                for r in self.dynamics_result.expect:
                    self.dynamics_ax.plot(self.dynamics_result.times * eta, r)
                    self.dynamics_ax.set_xlim([0,max(self.dynamics_result.times) * eta])
            
            #for t in xrange(self.parameters.points):
            #    self.dynamics_ax.axvline(linewidth=0.5, color='b',x = max(self.dynamics_result.times) * eta * t/self.parameters.points, ls='--')
                
            #self.dynamics_ax.plot( self.b)
            if self.LegendBool.get():
                self.dynamics_ax.legend(self.legends)
            
            
            
        self.dynamics_ax = self.figure.add_subplot(121)
        self.dynamics_ax.set_title('Dynamics')
        self.dynamics_ax.xaxis.set_tick_params(labelsize='large')
        self.dynamics_ax.yaxis.set_tick_params(labelsize='large')
        self.dynamics_ax.set_xlabel('$A/(2N) \, t$')
        self.dynamics_ax.set_ylim([-1.2,1.2])
        self.dynamics_ax.set_aspect('auto')
        
        self.computation.prepare_computation(self.parameters)
        
        eta = 1#self.parameters.eta if  self.parameters.eta > 0 else 1.    

        self.saveButton['state'] = NORMAL
        self.keepButton['state'] = NORMAL
        self.clearButton['state'] = NORMAL
        self.legends = []
        show_dynamics()
        
        self.canvas.show()    


















    
                                
