# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 16:29:12 2012

@author: jorgehog
"""

import sys, os, re

sys.path.append(os.getcwd() + "/tools")

from pyLibQMC import paths, variables



def dumpStrList(aList):
    i = 0
    for element in aList:
        print "[%d] %s" % (i, element)
        i+=1

def selectFile(Files):
    while True:
        dumpStrList(Files)
    
        action = raw_input("type 'display/run 0' to view content of / run first file.") 
#        action = "run 1"

        if re.findall("run \d+|display \d+", action):
            cmd, n = action.split()
            n = int(n)
            
            if n >= len(Files):
                print "Invalid file"
            else:
                if cmd == "display":
                    os.system("cat %s" % paths.iniFilePath + "/" + Files[n])
                else:
                    return Files[n]
                    
        else:
            print "Invalid option"

def setTag(arglist, line):
   
   if line.startswith('general'):

       if "-g" not in arglist:
           arglist.append("-g")
        
        
   elif line.startswith('VMC'):
       if "-v" not in arglist:
           arglist.append("-v")
    
        
   elif line.startswith('DMC'):
       if "-d" not in arglist:
           arglist.append("-d")
        
   elif line.startswith('MIN'):
       if "-m" not in arglist:
           arglist.append("-m")

   elif line.startswith('output'):
        if "-o" not in arglist:
           arglist.append("-o")
    
   elif line.startswith('variational'):
        if "-vp" not in arglist:
           arglist.append("-vp")


def parseFiles():
    iniDirContent = os.listdir(paths.iniFilePath)
    iniFiles = []
    
    for content in iniDirContent:
        if (re.findall("(.+).ini$", content)):
            iniFiles.append(content)

    if len(iniFiles) == 0:
        print "No iniFiles found in %s" % paths.iniFilePath
    else:
        print "Found iniFile(s):\n-----------------------------"
    
    iniFile = open(paths.iniFilePath + "/" + selectFile(iniFiles), 'r')
    
    
    arglist = []
    for line in iniFile:

        setTag(arglist, line)
        raw = re.findall(".+\s*=\s*.+", line)
            
        if raw:
            arglist.append(raw[0].replace(" ", ""))
        
                
            
            
    return convertToCMLargs(arglist)

def varParameterMap(n_p, dim, w, system):
    
    if dim == 2 and system == "QDots":
        alpha = 0
        beta = 0
        
        if n_p == 2:
            if w == 1.0:
                alpha = 0.987;
                beta = 0.398;
            elif w == 0.5:
                alpha = 0.9826;
                beta = 0.3098;
            elif w == 0.28:
                alpha = 0.9806;
                beta = 0.2483;
            elif w == 0.01:
                alpha = 0.824;
                beta = 0.081;
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 

        elif n_p == 6:
            if w == 1.0:
                alpha = 0.92;
                beta = 0.565;
            elif w == 0.5:
                alpha = 0.9004;
                beta = 0.4099;
            elif w == 0.28:
                alpha = 0.88;
                beta = 0.33;
            elif w == 0.01:
                alpha = 0.64;
                beta = 0.09;
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 

        elif n_p == 12:
            if w == 1.0:
                alpha = 0.87;
                beta = 0.68;
            elif w == 0.5:
                alpha = 0.8453;
                beta = 0.4813;
            elif w == 0.28:
                alpha = 0.8662;
                beta = 0.3346;
            elif w == 0.01:
                alpha = 0.37;
                beta = 0.55;
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 

        elif n_p == 20:
            if w == 1:
                alpha = 0.8351;
                beta = 0.7451;
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 

        elif n_p == 30:
            if w == 1:
                alpha = 0.78;
                beta = 0.85;
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 


        else:
            print "No saved parameters for n_p=", n_p 
            
        return alpha, beta

    else:
        print "Unknown type ", system, "with dim=", dim 


def convertToCMLargs(arglist):
    nInputs = 41;
    print arglist
    cmlArgs = ["def"]*nInputs
    
    cmlMAPo = {"blocking_out":  0,
               "dist_out":      1,
               "outputSuffix":  2,
               "outputPath":    3}    
    
    cmlMAPg = {"n_p":           4,
               "dim":           5,
               "w":             6,
               "random_seed":   7,
               "h":             8,
               "D":             9,
               "parallell":     10,
               "doMIN":         11,
               "doVMC":         12,
               "doDMC":         13,
               "use_coulomb":   14,
               "use_jastrow":   15,
               "sampling":      16,
               "kinetics":      17,
               "system":        18}
    
    cmlMAPv = {"n_c":           19,
               "dt":            20}
    
    cmlMAPd = {"dt":            21,
               "E_T":           22,
               "n_b":           23,
               "n_w":           24,
               "n_c":           25,
               "therm":         26,
               "dist_in":       27}
               
    cmlMAPm = {"max_step":      28,
               "f_max":         29,
               "f_min":         30,
               "omega":         31,
               "A":             32,
               "a":             33,
               "SGDsamples":    34,
               "n_w":           35,
               "therm":         36,
               "n_c":           37,
               "n_c_SGD":       38,
               "alpha":         39,
               "beta":          40}
               
    
    key = ""
    for arg in arglist:
        if arg.startswith("-"):
            key = arg[1:];
        else:
            name, val = arg.split("=");
            index = eval("cmlMAP" + key + "[name]")
            cmlArgs[index] = val

    cmlArgs = " ".join(cmlArgs)
    return [cmlArgs]
            
def initRuns(CMLargs):
    for CMLarg in CMLargs:
        os.system(paths.programPath + "/" + variables.QMC2programName + " " + CMLarg)
        

def main():
    if len(sys.argv) == 1:
        CMLargs = parseFiles()
    else:
        CMLargs = [convertToCMLargs(sys.argv[1:])]
    print CMLargs

    initRuns(CMLargs)
    
main()