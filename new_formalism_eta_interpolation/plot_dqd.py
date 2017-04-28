import matplotlib, sys, os
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimization
import math
from matplotlib.widgets import Slider, Button, RadioButtons


def iz(eta, t):
    if eta > 1:
        d = math.sqrt(eta**2-1)
        return 1 + math.exp(-2*eta*t)/(eta**2-1)*(1-eta*d*math.sinh(2*d*t)-eta**2*math.cosh(2*d*t))
    elif eta < 1:
        d = math.sqrt(1-eta**2)
        return 1 + math.exp(-2*eta*t)/(eta**2-1)*(1+eta*d*math.sin(2*d*t)-eta**2*math.cos(2*d*t))
    elif eta == 1:
        return  math.exp(-2 * t) *(-1 +  math.exp(2 * t) - 2 * t - 2 *t**2)

fig = plt.figure(figsize=(10,8)) 
ax1 = fig.add_subplot(111)

t = np.linspace(0, 10, num=1000)
eta = 1.2
l, = ax1.plot(t, [iz(0.1, s) for s in t], label='$\eta=0.1$', lw=2, ls='-')
l, = ax1.plot(t, [iz(1./2**0.5, s) for s in t], label=r'$\eta=\frac{1}{\sqrt{2}}$', lw=2, ls='--')
l, = ax1.plot(t, [iz(4, s) for s in t], label='$\eta=4$', lw=2, ls='-.')

#axeta = plt.axes([0.2, 0.025, 0.7, 0.03], axisbg='r')
#eta_slider = Slider(axeta, '$\eta$', 0.001, 10.0, valinit=1.1)

def update(val):
    eta = eta_slider.val
    l.set_ydata([iz(eta, s) for s in t])
    plt.draw()
    
#eta_slider.on_changed(update)
ax1.set_xlabel('$t$', fontsize=32)
ax1.set_ylabel(r'$\left< I_z\right>/I$', fontsize=32)
ax1.set_ylim([0,1.1])
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', which='minor', labelsize=20)
ax1.legend()
ax1.legend(bbox_to_anchor=(0.65, 0.4), loc=2, borderaxespad=0.,prop={'size':22})
plt.show()
