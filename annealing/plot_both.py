import matplotlib, sys, os
import matplotlib.pyplot as plt
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
for i in xrange(0, number, 2):
	try:
		with open('archive/' + folder_name + '/dynamics_data{}'.format(i), 'r+') as f:
			data[i] = f.readlines()
		with open('archive/' + folder_name + '/hysteresis_data{}'.format(i), 'r+') as f:
			data[i+1] = f.readlines()

	except IOError:
		print 'Invalid input file'
		to_delete.append(i)
data = [data[i] for i in xrange(number) if i not in to_delete]


COLORS = ['b','r','g','y','o','r','b','g','y','o']
LEGENDS = ['golden rule', 'central', '_nolegend_', '_nolegend_', '$t_{EQ}=5$', '$t_{EQ}=5$']
STYLES = ['-','-','o','o','-','--','-','--','-','--']
DASHES = [[8,0],[8,0],[8,0],[8,0],[8,0]]

fig, ax1 = plt.subplots()
ax2 = ax1.twiny()
plt.ylim(-1.2, 1.2)
plt.xlim(0, max([float(v.split()[0]) for v in data[0]+data[1]]))
#plt.xlabel('$\eta\, t$', fontsize=24)
#plt.ylabel(r'$\left< I_z/I\right>$', fontsize=24)
#MAX = 600
#plt.xticks([500 + i for i in xrange(0,MAX,100)], [ i for i in xrange(0,MAX,100)])

Min = 0#len(d)/2-4
x_avg = sum([float(v.split()[0]) for v in data[0][Min:]])/len([float(v.split()[0]) for v in data[0][Min:]])

plt.axvline(linewidth=0.5, color='black',x = x_avg, ls='-')
plt.axhline(linewidth=0.5, color='black',y = 0, ls='-') 

#ax2.set_xticks(new_tick_locations)
#ax2.set_xticklabels(tick_function(new_tick_locations))

i = 0
for d in data:
	x = [float(v.split()[0]) for v in d[Min:]]
	y = [float(v.split()[1]) for v in d[Min:]]
	print data.index(d)
	ax1.plot(x, y, c=COLORS[i], label=LEGENDS[i])
	i += 1

TEQ = 50
for d in data[2:]:
	x = [(d.index(v)+1)*TEQ for v in d]
	y = [float(v.split()[1]) for v in d[Min:]]
	print data.index(d)
	ax2.scatter(x, y, c=COLORS[i], label=LEGENDS[i])
	i += 1
#, dashes=[8, 4, 2, 4, 2, 4]
plt.legend(bbox_to_anchor=(0.65, 0.25), loc=2, borderaxespad=0.)


plt.show()
