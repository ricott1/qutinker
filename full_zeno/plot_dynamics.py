import matplotlib, sys, os
import matplotlib.pyplot as plt

for name in [d for d in os.listdir('.') if os.path.isdir(os.path.join('.', d))]:
	try:
		with open('{}/expectation_0'.format(name), 'r+') as f:
			t1, sz1 = zip(*[l.split() for l in f.readlines()])
	    with open('{}/expectation_1'.format(name), 'r+') as f:
			t2, sz2 = zip(*[l.split() for l in f.readlines()])

	except IOError:
		print 'Invalid input file'
		

plt.ylim(-1.2, 1.2)
plt.xlim(0, max([float(v.split()[0]) for v in data[0]+data[1]]))
plt.plot(t1, sz1, linewidth=2.0, label='$S_1^z$')
plt.plot(t2, sz2, linewidth=2.0, label='$S_2^z$')
plt.legend(bbox_to_anchor=(0.65, 0.25), loc=2, borderaxespad=0.)


plt.show()
