import matplotlib, sys, os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np

try:
	folder_name = sys.argv[1]
except IndexError:
	folder_name = raw_input('>Open: ')
try:
	number = int(sys.argv[2])
except IndexError, ValueError:
	number = int(raw_input('>Number: '))




COLORS = ['r','r','b','b','g','g','y','y','o','o']
STYLES = ['-','-','--','--','-','--','-','--','-','--']
LEGENDS = ['Golden rule', '_nolegend_', 'Central spin', '_nolegend_', '$t_{EQ}=5$', '$t_{EQ}=5$']
fig = plt.figure(figsize=(10, 20)) 
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', which='minor', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='minor', labelsize=20)





ax1.axvline(linewidth=0.5, color='black',x = 0, ls='-')
ax1.axhline(linewidth=0.5, color='black',y = 0, ls='-') 
ax2.axvline(linewidth=0.5, color='black',x = 0, ls='-')
ax2.axhline(linewidth=0.5, color='black',y = 0, ls='-') 
trans1 = ax1.get_xaxis_transform()
trans2 = ax2.get_xaxis_transform()  
ax1.annotate(r'(a)', xy=(-2.9, 0.9), size=26, xycoords=trans1)
ax2.annotate(r'(b)', xy=(-2.9, 0.9), size=26, xycoords=trans2)





ax1.set_ylim(-1.2, 1.2)
ax1.set_xlabel('$b$', fontsize=32)
ax1.set_ylabel(r'$\left< I_z\right>/I$', fontsize=32)
ax2.set_ylim(-1.2, 1.2)
ax2.set_xlabel('$\gamma$', fontsize=32)
ax2.set_ylabel(r'$\left< I_z\right>/I$', fontsize=32)




folder_name = 'bplots'
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
i = 0
for d in data:
	x = [float(v.split()[0]) for v in d]
	y = [float(v.split()[-1]) for v in d]
	ax1.plot(x, y, COLORS[i], linewidth=2.0, ls=STYLES[i], label=LEGENDS[i])
	
	i += 1

folder_name = 'gammaplots'
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
i = 0
for d in data:
	x = [float(v.split()[1]) for v in d]
	y = [float(v.split()[-1]) for v in d]
	ax2.plot(x, y, COLORS[i], linewidth=2.0, ls=STYLES[i], label=LEGENDS[i])
	
	i += 1






ax1.legend(bbox_to_anchor=(0.6, 0.26), loc=2, borderaxespad=0.,prop={'size':22})
ax2.legend(bbox_to_anchor=(0.6, 0.26), loc=2, borderaxespad=0.,prop={'size':22})
plt.show()
