# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 13:05:51 2012

@author: jorgehog
"""

from pyLibQMC import paths

validation = True
fullruns = False

if validation:

    wList = [0.5,1]
    npList = [2,6,12, 20, 30]

    rawFile = """general:
n_p = __NP__
w = __W__
use_coulomb = 0
doVMC = 1
doDMC = 1
doMIN = 1

VMC:
n_c = 1E3

DMC:
n_b=10
therm=100
n_c=100

MIN:
SGDsamples=1000
n_c_SGD=100

"""

    for w in wList:
        for np in npList:
                filename = "test_np%d_w%s.ini" % (np, str(w).replace(".", ""))
                iniFile = rawFile.replace("__NP__", str(np))
                iniFile = iniFile.replace("__W__", str(w))
            
                outFile = open(paths.iniFilePath + "/" + filename, 'w')
                outFile.write(iniFile)
                outFile.close()
                
                
if fullruns:
    
    wList = [1, 0.5, 0.28]
    npList = [2, 6, 12]
    
    rawFile = """general:
n_p = __NP__
w = __W__
doMIN = 1
doVMC = 1
doDMC = 1

DMC:
dist_in = 1
therm=5000

output:
dist_out = 1

MIN:
SGDsamples=10000


"""

    for w in wList:
        for np in npList:
            filename = "fullRun_np%d_w%s.ini" % (np, str(w).replace(".", ""))
            iniFile = rawFile.replace("__NP__", str(np))
            iniFile = iniFile.replace("__W__", str(w))
            
            outFile = open(paths.iniFilePath + "/" + filename, 'w')
            outFile.write(iniFile)
            outFile.close()
