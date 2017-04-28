import matplotlib, sys, os
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimization


def t(g, eta, N=1):
    i = N/2.
    return (1. + eta**2 + 4 *eta**4)/(2*eta + 4 *eta**3)
    
try:
    option = sys.argv[1]
except IndexError:
    option = '-t'
    
    
try:
    dis = sys.argv[2]
    print dis
except IndexError:
    dis = '-central'
COLORS = ['b','r','y','g','black','white','magenta','purple', 'cyan','grey','brown','pink','b','r','y','g','black','white','magenta','purple', 'cyan','grey','brown','pink']
t_scale = []
eta_scale = []
short_times = []
long_times = []

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(111)
STYLES = ['-','--','-.',':','-','--','-','--','-','--']
x, y = [], []
with open('archive/central/mean_times_N=1', 'r+') as f:
        for l in f.readlines():
            y.append(float(l.split()[0]))
            x.append(float(l.split()[2])/4)
        ax1.plot(x, y, lw=2)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', which='minor', labelsize=20)
ax1.set_ylabel(r'$\bar{t}$', fontsize=32)
ax1.set_xlabel(r'$\eta_0/\alpha$', fontsize=28)
ax1.set_xlim([0,12])
ax1.set_ylim([0,80])

#


plt.show()
