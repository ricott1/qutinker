import matplotlib, sys, os
import matplotlib.pyplot as plt
try:
	folder_name = sys.argv[1]
except IndexError:
	folder_name = raw_input('>Open: ')
try:
	number = int(sys.argv[2])
except IndexError, ValueError:
	number = raw_input('>Number: ')

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
LEGENDS = ['$t_{EQ}=50$', '_nolegend_', '$t_{EQ}=5$', '_nolegend_', '$t_{EQ}=5$', '$t_{EQ}=5$']




plt.ylim(-1.2, 1.2)
#plt.xlabel('$b$', fontsize=24)
#plt.ylabel(r'$\left< I_z/I\right>$', fontsize=24)



plt.axvline(linewidth=0.5, color='black',x = 0, ls='-')
plt.axhline(linewidth=0.5, color='black',y = 0, ls='-') 



i = 0
for d in data:
	x = [float(v.split()[0]) for v in d]
	y = [float(v.split()[1]) for v in d]
	plt.plot(x, y, COLORS[i], linewidth=2.0, ls=STYLES[i], label=LEGENDS[i])
	i += 1

plt.legend(bbox_to_anchor=(0.75, 0.25), loc=2, borderaxespad=0.)


	
plt.show()
