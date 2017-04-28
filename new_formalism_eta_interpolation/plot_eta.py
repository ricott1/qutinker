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
if dis in ('-c', '-central'):
    print 'ale;'
    folder = 'central'
    y_max = 50
else:
    folder = 'fermi_eta'
    y_max = 100
#/%s_eta
MAX = 2
for i in xrange(1, MAX):
    n = i#2**i
    x, y = [], []
    try:
        with open('archive/{}/mean_times_N={}'.format(folder,  n), 'r+') as f:
            for l in f.readlines():
                y.append(float(l.split()[0]))
                x.append(float(l.split()[2]))
            x_max = max(x)
            if option == '-t':        
                ax1.plot(x, y, label='N=%d'%n, lw=2, ls=STYLES[i])
                #ax1.plot(x, [t(0, eta/4, n) for eta in x], label='analytical N=%d'%n)
            t_scale.append((float(min(y)), n))
            eta_scale.append((float(x[y.index(min(y))]), n))
            short_times.append((y[2], n))
            try:
                long_times.append((y[90], n))
            except:
                continue
    except IOError:
        continue


tt = zip(*t_scale)
ee = zip(*eta_scale)
st = zip(*short_times)
lt = zip(*long_times)


if option == '-t':    
    ax1.set_ylabel(r'$\bar{t}$', fontsize=32)
    ax1.set_xlabel(r'$\eta$', fontsize=32)
    ax1.set_xlim([0, x_max])
    ax1.set_ylim([0,y_max])
    ax1.legend()
    
    #g = 1
if option == '-an':    
    x, y = [], []
    with open('archive/central/mean_times_N=1', 'r+') as f:
        for l in f.readlines():
            y.append(float(l.split()[0]))
            x.append(float(l.split()[2]))
        #ax1.scatter(x, y, label='N=1', color=COLORS[1])
    x = np.arange(0,4.1,0.01)   
    for n in xrange(1,2):
        ax1.plot(x, [t(0, (eta), n) for eta in x], lw=2)#, label='analytical N=%d'%n)
    ax1.set_ylabel(r'$\bar{t}$', fontsize=40)
    ax1.set_xlabel(r'$\Gamma_0/\alpha$', fontsize=40)
    ax1.set_xlim([0,max(x)])
    ax1.set_ylim([0,10])
    ax1.legend()


elif option == '-eta':    
    ax1.scatter(ee[1], ee[0])
    ax1.set_ylabel(r'$\eta_{min}$', fontsize=28)
    ax1.set_xlabel(r'$N$', fontsize=28)
    ax1.set_xlim([0,MAX])
    ax1.set_ylim([0,8])
    def func(x, a, b):
        return (4-2**0.5)/2**0.5*x+1# + b#a + b * x
    #x0 = np.array([0.0, 0.0])
    #(a,b), corr = optimization.curve_fit(func, np.array(ee[1]), np.array(ee[0]), x0)
    #print a,b
    #ax1.annotate(r'$f(N) = +%.3f N +%.3f$'%(2*a,b), xy=(15, 3))#%.3f  %+.3f N$'%(a,b), xy=(15, 3))
    #ax1.plot(tt[1], [func(n, a, b) for n in tt[1]], label='analytical N=%d'%n)
    ax1.plot(tt[1], [4./(2**0.5) for n in tt[1]], label='analytical N=%d'%n)
elif option == '-tau':    
    ax1.scatter(tt[1], tt[0])
    ax1.set_ylabel(r'$\bar{t}_{min}$', fontsize=32)
    ax1.set_xlabel(r'$N$', fontsize=28)
    ax1.set_xlim([0,MAX])
    def func(x, a, b, c, d):
        return (b*x/2.*(c*x/2.+d))**0.5 + a
    x0 = np.array([0.0, 1.0, 2.0, 1.0])
    (a,b,c, d), corr = optimization.curve_fit(func, np.array(tt[1]), np.array(tt[0]), x0)
    
    #ax1.annotate(r'$f(N) = %.3f + \sqrt{%.3f\frac{N}{2}(%.3f\frac{N}{2}%+.3f)}$'%(a,b,c, d), xy=(15, 1))
    #ax1.plot(tt[1], [func(n, a,  b,  c, d) for n in tt[1]], label='analytical N=%d'%n)
    
elif option == '-small':    
    ax1.scatter(st[1], st[0],label=r'$\eta \ll 1$')
    ax1.set_ylabel(r'$t_{I}$', fontsize=28)
    ax1.set_xlabel(r'$N$', fontsize=28)
    ax1.legend()
    ax1.set_xlim([0,MAX])
    def func(x, a, b, c):
        return a + (x/2.*(b*x/2.+c))**0.5
    x0 = np.array([0.0, 0.0, 0.0])
    (a,b,c), corr = optimization.curve_fit(func, np.array(st[1]), np.array(st[0]), x0)
    ax1.annotate(r'$f(N) = %.3f + \sqrt{\frac{N}{2}(%.3f\frac{N}{2}%+.3f)}$'%(a,b,c), xy=(15, 3))
    ax1.plot(st[1], [func(n, a,  b,  c) for n in st[1]], label='analytical N=%d'%n)
elif option == '-large':    
    ax1.scatter(lt[1], lt[0],label=r'$\eta \gg 1$')
    ax1.set_ylabel(r'$t_{I}$', fontsize=28)
    ax1.set_xlabel(r'$N$', fontsize=28)
    ax1.legend()
    ax1.set_xlim([0,MAX])
    def func(x, a, b, c):
        return a + (x/2.*(b*x/2.+c))**-0.5
    x0 = np.array([0.0, 0.0, 0.0])
    (a,b,c), corr = optimization.curve_fit(func, np.array(lt[1]), np.array(lt[0]), x0)
    ax1.annotate(r'$f(N) = %.3f + \frac{1}{\sqrt{\frac{N}{2}(%.3f\frac{N}{2}%+.3f)}}$'%(a,b,c), xy=(15, 3))
    ax1.plot(tt[1], [func(n, a,  b,  c) for n in lt[1]], label='analytical N=%d'%n)

#ax1.legend(bbox_to_anchor=(0.2, 0.6), loc=3, borderaxespad=0.,prop={'size':20})
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', which='minor', labelsize=20)

#ax1.set_ylim([0,20])
#ax1.set_ylim([1,20])

#


plt.show()
