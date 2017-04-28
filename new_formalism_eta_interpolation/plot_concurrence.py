import matplotlib, sys, os
import matplotlib.pyplot as plt
import numpy as np
import collections
import scipy.optimize as optimization
from matplotlib.widgets import Slider, Button, RadioButtons
import math
def t(g, eta, N=1):
    i = N/2.
    return (1. + eta**2 + 4 *eta**4)/(2*eta + 4 *eta**3)
    
def ent(x):
    if x == 0:
        return 0
    elif x ==1:
        return 0
    else:
        return -x * np.log2(x) - (1 - x) * np.log2(1 - x)
    
COLORS = ['b','eta','concurrence','g','black','white','magenta','purple', 'cyan','grey','brown','pink','b','eta','concurrence','g','black','white','magenta','purple', 'cyan','grey','brown','pink']


fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
STYLES = ['-','--','-.',':','-','--','-','--','-','--']
files = os.listdir('archive/concurrence/')
max_x = 2000
r = []
for fil in files:
    try:
        with open('archive/concurrence/{}'.format(fil), 'r+') as f:
            x, y, e = [], [], []
            for l in f.readlines():
                x.append(float(l.split()[1]))
                try:
                    c = sorted([ complex(i).real**0.5 for i in l.split()[2:]])
                except ValueError:
                    print l.split()
                    raw_input()
                #print conc, l.split()[2:], sorted(conc)
                #print conc[-1] - sum(conc[:-1])
                #raw_input()
                conc = max(0, c[0] - sum(c[1:]))
                y.append(conc)
                enta = ent((1+(1-conc**2)**0.5)/2.)
                if math.isnan(enta):
                    print conc,enta
                e.append(enta)
            
            r.append({'eta' : float(l.split()[0]), 'time' : x, 'concurrence' : y, 'entanglement' : e})
    except IOError:
        continue
r = sorted(r, key= lambda i : i['eta'])
print [r[i]['eta'] for i in xrange(len(r))] 
t = r[0]['time'][:max_x]
idx = 0
l, = ax1.plot(t, r[idx]['concurrence'][:max_x], label='Concurrence', lw=2, ls='-')
m, = ax1.plot(t, r[idx]['entanglement'][:max_x], label='Entanglement', lw=2, ls='-')


axeta = plt.axes([0.2, 0.025, 0.7, 0.03], axisbg='r')
idx_slider = Slider(axeta, '$\eta$', 0, len(r), valfmt='%1.0f',valinit=0)

def update(val):
    idx = int(idx_slider.val)
    l.set_xdata(r[idx]['time'][:max_x])
    l.set_ydata(r[idx]['concurrence'][:max_x])
    m.set_xdata(r[idx]['time'][:max_x])
    m.set_ydata(r[idx]['entanglement'][:max_x])
    print '$\eta={}$'.format(r[idx]['eta'])
    #print r[idx]['concurrence'] 
    plt.draw()
    
    
ax1.set_ylim([0,0.75]) 
idx_slider.on_changed(update)
#ax1.legend(bbox_to_anchor=(0.2, 0.6), loc=3, borderaxespad=0.,prop={'size':20})
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', which='minor', labelsize=20)
ax1.legend()
ax1.legend(bbox_to_anchor=(0.65, 0.4), loc=2, borderaxespad=0.,prop={'size':22})
ax1.set_xlabel('$t$', fontsize=32)


#print [r[i]['entanglement'] for i in xrange(len(r))]
#print r[0]['entanglement']
ax2.scatter([r[i]['eta'] for i in xrange(len(r))], [sum(r[i]['concurrence'])/r[i]['time'][-1] for i in xrange(len(r))], color='r')
ax2.scatter([r[i]['eta'] for i in xrange(len(r))], [sum(r[i]['entanglement'])/r[i]['time'][-1] for i in xrange(len(r))])
ax1.set_xlabel('$\eta$', fontsize=32)

plt.show()
