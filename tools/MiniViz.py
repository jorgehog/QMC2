from scitools.std import *
import os

path = os.path.environ "/NetBeansProjects/nbQMC2/"
filename = "MiniOut.dat"



fil = open(path + filename, 'r')

N = 5

plotData = [[] for i in range(N)]
titles = ["AlphaGrad", "BetaGrad", "alpha", "beta", "E"]

for line in fil:	
	raw = line.split()
	for i in range(len(raw)):
		plotData[i].append(float(raw[i]))

if len(raw) == N:
	for i in range(N):
		figure(i+1)
		plot(plotData[i], title=titles[i])
else:
	print "hei"
	maping = [0,2,4]
	maping2 = [0,2,3];
	for i in range(len(raw)-1):
		figure(i+1)
		plot(plotData[maping2[i]], title=titles[maping[i]])

raw_input()
