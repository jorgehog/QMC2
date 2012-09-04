from scitools.std import *
import os

fileName = "DMC_out.dat"
filePath = os.path.expanduser("~") + "/NetBeansProjects/nbQMC2/" + fileName

cmd = ""
while (cmd != "q"):

	DMCout = open(filePath , "r")

	E = []
	Eavg = []
	N = []
	Navg = []

	for line in DMCout:
		raw = line.split()
		E.append(float(raw[0]))
		Eavg.append(float(raw[1]))
		N.append(float(raw[2]))
		Navg.append(float(raw[3]))

	DMCout.close()

	if len(E) == len(Eavg) == len(N):

		figure(1)
		plot(E, 'r')
		hold("on")

		plot(Eavg, 'b')
		title("E")
		xlabel("cycle")
		ylabel("E")
		legend(["E", "Eavg"])

		figure(2)
		plot(N, 'r')
		hold("on")

		plot(Navg, 'b')
		title("N/N_0")
		xlabel("cycle")
		ylabel("N/N_0")
		legend(["N", "Navg"])


	cmd = raw_input("press Enter to plot again or 'q' to end")
