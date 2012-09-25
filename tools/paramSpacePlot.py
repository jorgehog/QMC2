# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 16:48:27 2012

@author: jorgehog
"""

from scitools.std import plot, figure
import os

two = False


filename = "alphaVsE.dat"
filenameBF = "alphaVsE_BF.dat"
path = os.path.expanduser("~") + "/NetBeansProjects/nbQMC2/"

inFile_IS = open(path + filename, 'r')
if two:
    inFile_BF = open(path + filenameBF, 'r')


alpha = []
E = []
for line in inFile_IS:
    alphai, Ei = line.split()
    alpha.append(float(alphai))
    E.append(float(Ei))
inFile_IS.close()

if two:
    alpha2 = []
    E2 = []
    for line in inFile_BF:
        alphai, Ei = line.split()
        alpha2.append(float(alphai))
        E2.append(float(Ei))
        inFile_BF.close()
    
def get_grad_1d(E, alpha):
    grad = []

    h = alpha[1] - alpha[0]
    
    for i in range(1, len(alpha)-1):
        deriv = 1/(2*h)*(-E[i-1] + E[i+1])
        grad.append(deriv);
          
    return grad
    
aGrad = get_grad_1d(E, alpha)
if two: 
    aGrad2 = get_grad_1d(E2, alpha2)    
    
    
figure(1)
plot(alpha, E, xlabel="alpha", legend="IS", hold="on")
#figure(2)
#plot(alpha, E, xlabel="alpha",  legend="IS", hold="on")
#figure(3)
#plot(alpha, E, xlabel="alpha", legend="IS", hold="on")
figure(4)
plot(alpha[1:-1], aGrad, legend="IS", hold="on")                

if two:
    figure(1)
    plot(alpha2, E2, xlabel="alpha", ylabel="E", title="2p Qdot noInteraction",\
                                            legend="BF", hardcopy=path+"full.eps")
                                            #figure(2)
                                            #plot(alpha2, E2, xlabel="alpha", ylabel="E", title="2p Qdot noInteraction",\
                                            #    legend="BF", axis=[1,1.5, 1,3], hardcopy=path+"closeToOne.eps")
                                            #figure(3)
                                            #plot(alpha2, E2, xlabel="alpha", ylabel="E", title="2p Qdot noInteraction",\
                                            #                                hardcopy=path+"closeToTwoTwo.eps", legend="BF", axis=[1,1.5, 1,3])
    figure(4)
    plot(alpha2[1:-1], aGrad2, xlabel="alpha", ylabel="dE/dalpha", title="2p Qdot noInteraction",\
                                            hardcopy=path+"grad.eps", legend="BF")

raw_input()