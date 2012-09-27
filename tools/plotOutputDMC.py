from scitools.std import *
import os

fileName = "DMC_out.dat"
filePath = os.path.expanduser("~") + "/NetBeansProjects/nbQMC2/" + fileName

dynamicAxis = True

cmd = ""

if len(sys.argv) > 1:
	exactE = float(sys.argv[1])
else:
	exactE = 0;

while (cmd != "q"):

	DMCout = open(filePath , "r")

	E = []
	Eavg = []
	N = []
	Navg = []
	ET = []

	for line in DMCout:
		raw = line.split()
		E.append(float(raw[0]))
		Eavg.append(float(raw[1]))
		N.append(float(raw[2]))
		Navg.append(float(raw[3]))
		ET.append(float(raw[4]))

	DMCout.close()

	n= int(len(E)*(9./10))
	k = 10
	tmpEavg = sum(E[-n:])/n
	dE = abs(max(E[-n:]) - tmpEavg) 
	tmpAxis = [n,len(E), tmpEavg - k*dE, tmpEavg + k*dE] 

	if len(E) == len(Eavg) == len(N) == len(Navg) == len(ET):

		figure(1)
		plot(E, 'r')
		hold("on")

		plot(Eavg, 'b')
		title("E")
		xlabel("cycle")
		ylabel("E")
		legend(["E", "Eavg"])
		if (exactE != 0):
			plot(zeros(len(E)) + exactE, 'g')
		
		if dynamicAxis:
			axis(tmpAxis)
		

		figure(2)
		plot(N, 'r')
		hold("on")

		plot(Navg, 'b')
		title("N/N_0")
		xlabel("cycle")
		ylabel("N/N_0")
		legend(["N", "Navg"])

		figure(3)
		plot(ET, 'b')
		title("E_t")
		xlabel("cycle")
		ylabel("Et")
		if (exactE != 0):
			hold("on")
			plot(zeros(len(E)) + exactE, 'g')

		cmd = raw_input("press Enter to plot again or 'q' to end")
	else:
		cmd = "rerun"
