# -*- coding: UTF-8 -*-
import sys, time, math
import qutip as qt
import numpy as np
import copy

def decompose(N, s=0.5):
	global Ix, Iy, Iz, Isq, Ip
	Ix = sum([nuclearOperator(qt.jmat(s, 'x'), n, N, s) for n in xrange(N)])#sum([nuclearOperator(0.5 * qt.sigmax(), n, N) for n in xrange(N)])
	Iy = sum([nuclearOperator(qt.jmat(s, 'y'), n, N, s) for n in xrange(N)])#sum([nuclearOperator(0.5 * qt.sigmay(), n, N) for n in xrange(N)])
	Iz = sum([nuclearOperator(qt.jmat(s, 'z'), n, N, s) for n in xrange(N)])#sum([nuclearOperator(0.5 * qt.sigmaz(), n, N) for n in xrange(N)])
	Isq = (Ix**2 + Iy**2 + Iz**2).tidyup()
	
	Ip = sum([nuclearOperator(qt.jmat(s, '+'), n, N, s) for n in xrange(N)])#sum([nuclearOperator(qt.sigmap(), n, N) for n in xrange(N)])
	Im = sum([nuclearOperator(qt.jmat(s, '-'), n, N, s) for n in xrange(N)])#sum([nuclearOperator(qt.sigmam(), n, N) for n in xrange(N)])
	ALLSTATES = []
	stuff = qt.simdiag([Iz, Isq])
	tot = int((2*s+1))**N
	
	#print stuff
	p = 0
	for i in xrange(tot):
		M = round(stuff[0][0][i] + 0,2)#to be sure that 0 is always positive
		I = round((-1. + (1+ 4* stuff[0][1][i])**0.5)/2. +0 ,2)
		state = stuff[1][i].tidyup().unit()
		if M == +I:	
			p += 1		
			ALLSTATES.append({'I' : I,'M' : M,'p' : p,'state' : state})
			#print I,M,p,state
	
	ALLSTATES = orthogonalize(ALLSTATES)
	
	checkOverlap(ALLSTATES)
	#return ALLSTATES
	#print ALLSTATES
	for elem in ALLSTATES:		
		if elem['M'] - 1 >= -elem['I']:
			Im2 = qt.Qobj(Im.full())
			state2 = qt.Qobj(elem['state'].full())
			elem['state'] = state2
			newState= (Im2 * state2)
			newElem = dict(elem)
			newElem['M'] -= 1
			#print newState
			try:	
				newElem['state'] = newState.unit().tidyup()
			except IndexError, ZeroDivisionError:
				print elem['state'],elem['I'],elem['M'], newState
				pass
				#print elem['state'],elem['I'],elem['M'], newState
				
			ALLSTATES.append(newElem)
		

	ALLSTATES = sorted([el for el in ALLSTATES], key = lambda el: (el['I'], el['M']))
	
	return ALLSTATES
def orthogonalize(allStates):
	def gramSchmidt(states):
		orthostates = []
		for v1 in states:
			for v2 in orthostates:
				v1 = v1 - v2.overlap(v1.dag()) * v2
			orthostates.append(v1.unit().tidyup())
		return orthostates

	orthostates = gramSchmidt([s['state'] for s in allStates])
	for i in xrange(len(allStates)):
		allStates[i]['state'] = orthostates[i]
	#print allStates
	return allStates


def represent(allstates):
	tot = len(allstates)
	print 'There are %d nuclear states'%tot
	N = int(math.log(tot,2))
	for elem in allstates:
		print 'I = %g  M = %g p = %g'%(elem['I'],elem['M'],elem['p'])
		
		#create a list with all possible permutations of 0 and 1 for each spin
		spinsList = [[0 for _ in xrange(N)]]
		for t in xrange(1,tot):
			binT = bin(t)[2:].zfill(N)
			spinsList.append([int(binT[i]) for i in xrange(N)])
		#to each entry in the state correspond a precise spins configuration.
		for i,d in enumerate(elem['state'].full()):
			#print elem['state'].full()
			if d.real != 0:
				print '%+.2f %2s'%(d.real , '|'+' '.join([UPARROW if s==0 else DOWNARROW for s in spinsList[i]])+' >'),
		
		print '\n'


def nuclearOperator(QObj, n, N, s=1./2):
	
	nucleiList = [qt.qeye(2*s + 1) for _ in xrange(N)]
	nucleiList[n] = QObj
	return qt.tensor(nucleiList)	


def checkOverlap(allstates):
	for i in xrange(len(allstates)):
		for j in xrange(i+1,len(allstates)):
			ol = allstates[i]['state'].overlap(allstates[j]['state']).real
			if ol < 0.000000001:
				ol = 0
			if ol:
				print 'Overlap < %g %g %d|%g %g %d > = %g'%(allstates[i]['I'],allstates[i]['M'],allstates[i]['p'],allstates[j]['I'],allstates[j]['M'],allstates[j]['p'],ol)
				
UPARROW = u'\u2191'
DOWNARROW = u'\u2193'
CHECKOVERLAP = 1
ORTHOGONAL = 0
if __name__ == '__main__':
	try:	
		K = int(sys.argv[1])
	except:
		K = 4
	
	iz = qt.jmat((K/2.), 'z')	
	#print  iz
	for i in xrange(K+1):
		state = qt.basis(K + 1, i)
		#print state,iz.matrix_element(state.dag(), state).real
	izc = sum([nuclearOperator(0.5 * qt.sigmaz(), n, K, 0.5) for n in xrange(K)])

	for statec in decompose(K):
		if statec['I'] == K/2.:
			st = statec['state']
			#print st,izc.matrix_element(st.dag(), st).real
		
	print [izc.eigenstates()[0][_],





	
