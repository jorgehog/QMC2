# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 13:05:51 2012

@author: jorgehog
"""

from pyLibQMC import paths

validation = False
fullruns = True

if validation:

    wList = [0.5,1]
    npList = [2,6,12]
    useColJast = [0,1]

    rawFile = """general:
n_p = __NP__
w = __W__
use_coulomb = __USECJ__
use_jastrow = __USECJ__
doVMC = 1

VMC:
n_c = 1E4

"""

    for w in wList:
        for np in npList:
            for cj in useColJast:
                filename = "test_np%d_w%s_cj%d.ini" % (np, str(w).replace(".", ""), cj)
                iniFile = rawFile.replace("__NP__", str(np))
                iniFile = iniFile.replace("__W__", str(w))
                iniFile = iniFile.replace("__USECJ__", str(cj))
            
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
