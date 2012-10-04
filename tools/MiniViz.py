from scitools.std import plot, figure

from pyLibQMC import paths

filename = "alpha.dat"



fil = open(paths.IDEPath + "/" + filename, 'r')

N = 7

plotData = [[] for i in range(N)]
titles = ["AlphaGrad", "BetaGrad","aStep", "alpha", "beta", "E", "Em"]

for line in fil:	
	raw = line.split()
	for i in range(len(raw)):
         plotData[i].append(float(raw[i]))

if len(raw) == N-1:
    N = 5;
    titles.remove("BetaGrad");
    titles.remove("beta");
    plotData.pop(1);
  
for i in range(N):
    if (titles[i] == "Em"):
        figure(i)
    else:
        figure(i+1)
    
    if (titles[i]=="E"):
        plot(plotData[i], legend=titles[i], hold="on")
    else:
        plot(plotData[i], legend=titles[i])

raw_input()