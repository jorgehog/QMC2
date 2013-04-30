# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 13:05:51 2012

@author: jorgehog
"""

from pyLibQMC import paths

validation = True
fullruns = False
test = False
abel = False

if test:
    wList = [0.5,1]
    npList = [2]
    
    rawFile = """general:
n_p = __NP__
systemConstant = __W__
use_coulomb = 0
doDMC = 1

DMC:
n_b=1
n_c=100

"""
    for w in wList:
        for np in npList:
                filename = "SCRIPTTEST_np%d_w%s.ini" % (np, str(w).replace(".", ""))
                iniFile = rawFile.replace("__NP__", str(np))
                iniFile = iniFile.replace("__W__", str(w))
            
                outFile = open(paths.iniFilePath + "/" + filename, 'w')
                outFile.write(iniFile)
                outFile.close()

    
    
    


if validation:

    npList = range(2, 29, 2)

    rawFile = """general:
n_p = __NP__
use_coulomb = 0
system=Atoms
doVMC = 1
doDMC = 1

VMC:
n_c = 1E5

DMC:
n_b=10
therm=10
n_c=10


MIN:
SGDsamples=2000
n_c_SGD=400
therm=10000
"""

  
    for np in npList:
            filename = "test_Atoms_np%d.ini" % (np)
            iniFile = rawFile.replace("__NP__", str(np))
            
            outFile = open(paths.iniFilePath + "/" + filename, 'w')
            outFile.write(iniFile)
            outFile.close()
                
                
if fullruns:
    
    wList = [0.1, 0.28, 0.5, 1]
    npList = [42]
    
    rawFile = """general:
n_p = __NP__
systemConstant = __W__
doVMC = 1
doDMC = 1
do_blocking=1

VMC:
n_c=1E3

DMC:
n_c=3000
n_w=100
"""

    for w in wList:
        for np in npList:
            filename = "fullRun_np%d_w%s.ini" % (np, str(w).replace(".", ""))
            iniFile = rawFile.replace("__NP__", str(np))
            iniFile = iniFile.replace("__W__", str(w))
            
            outFile = open(paths.iniFilePath + "/" + filename, 'w')
            outFile.write(iniFile)
            outFile.close()

if abel:

    raw = """general:
doVMC=1
doDMC=1
do_blocking=1
systemConstant=__w__
n_p=__N__

VMC:
n_c=1E8
DMC:
n_w=256000
n_c=3000
"""

    Nl = [2,6,12]
    wl = [0.01, 0.001]

    for N in Nl:
        for w in wl:
            filename = "lowWabelRun_N%d_w%s.ini" % (N, str(w).replace(".", ""))
            iniFile = raw.replace("__N__", str(N))
            iniFile = iniFile.replace("__w__", str(w))
            
            
            outFile = open(paths.iniFilePath + "/" + filename, 'w')
            outFile.write(iniFile)
            outFile.close()
