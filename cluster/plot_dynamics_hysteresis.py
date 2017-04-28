import matplotlib, sys, os
import matplotlib.pyplot as plt
import numpy as np 
try:
	folder_name = sys.argv[1]
except IndexError:
	folder_name = raw_input('>Open: ')
try:
	number = int(sys.argv[2])
except IndexError, ValueError:
	number = int(raw_input('>Number: '))

data = [[] for _ in xrange(number)]
to_delete = []
for i in xrange(number):
	try:
		with open('archive/' + folder_name + '/hysteresis_data{}'.format(i), 'r+') as f:
			data[i] = f.readlines()

	except IOError:
		print 'Invalid input file'
		to_delete.append(i)
data = [data[i] for i in xrange(number) if i not in to_delete]


COLORS = ['r','r','b','b','g','g','y','y','o','o']
STYLES = ['-','-','--','--','-','--','-','--','-','--']
LEGENDS = [r'$\left< I_z\right>/I$', 'Hysteresis', 'b', '_nolegend_', '$t_{EQ}=5$', '$t_{EQ}=5$']




fig = plt.figure(figsize=(10, 8)) 
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
b_init, b_final = -2, 4
points = 31
x = np.linspace(0,points * 100,num=24000, endpoint=True)

split_tlist = np.array_split(x, points)
y = []

for i in xrange(0,points,1):
    y += [b_init * (1.- i/(points-1.)) + b_final * i/(points-1.) for s in split_tlist[i]]
 

ax2.plot(x, y, COLORS[4], linewidth=2.0, ls='-', label=LEGENDS[2])

ax1.axvline(linewidth=0.5, color='black',x = 0, ls='-')
ax1.axhline(linewidth=0.5, color='black',y = 0, ls='-') 
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', which='minor', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='minor', labelsize=20)

i = 0
with open('archive/' + folder_name + '/dynamics_data{}'.format(i), 'r+') as f:
	dyn = f.readlines()
x = [float(v.split()[0]) for v in dyn]
y = [float(v.split()[-1]) for v in dyn]
#print x, y
ax1.plot(x, y, COLORS[0], linewidth=2.0, ls='--', label=LEGENDS[0])

labels = np.linspace(0,1,num=11, endpoint=True)
tick_locs = [i*points*10 for i in xrange(11)]
tick_lbls = labels
plt.xticks(tick_locs, tick_lbls)
#ax1.set_xticks(tick_locs, tick_lbls)
for d in data:
	x = [p * 3100/points for p in xrange(1, points+1)]
	y = [float(v.split()[-1]) for v in d]
	print x, y
	ax1.scatter(x, y)# COLORS[i],
	i += 1
	
	
	
	


ax1.set_xlim(0, 3100)
ax1.set_ylim(-1.2, 1.2)
ax2.set_ylim(-2.2, 4.2)

ax1.legend(bbox_to_anchor=(0.7, 0.2), loc=3, borderaxespad=0.,prop={'size':22})
ax2.legend(bbox_to_anchor=(0.7, 0.1), loc=3, borderaxespad=0.,prop={'size':22})
all_b = ['' if i%5!=0 else round(b_init * (1.- i/(points-1.)) + b_final * i/(points-1.),1) for i in xrange(0,points,1)]
#ax2.set_xticks(xrange(0, points, 1))
#ax2.set_xticklabels(all_b)

ax1.set_xlabel('$t/(p\,t_{EQ})$', fontsize=32)
ax2.set_ylabel('$b$', fontsize=32)
ax1.set_ylabel(r'$\left< I_z\right>/I$', fontsize=32)
	
plt.show()
