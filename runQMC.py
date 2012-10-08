# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 16:29:12 2012

@author: jorgehog
"""

import sys, os, re

sys.path.append(os.getcwd() + "/tools")

from pyLibQMC import paths, variables, add_date



def dumpStrList(aList):
    i = 0
    for element in aList:
        print "[%d] %s" % (i, element)
        i+=1

def selectFiles(Files):
    while True:
        dumpStrList(Files)
    
#        action = raw_input("type 'display/run 0' to view content of / run first file.") 
        action = "run 0 1"

        if action.split()[0] == "run" or action.split()[0] == "display":

            if action.split()[0] == "run":
                fileIDs = [int(ID) for ID in action.split()[1:]]

                legal = True
                for ID in fileIDs:
                    if fileIDs.count(ID) > 1:
                        print "Several instances of file#", ID
                        legal = False
                
                
                for ID in fileIDs:
                    if ID >= len(Files) or ID < 0:
                        print "Invalid file", ID
                        legal = False
                    
                if legal:
                    return [Files[n] for n in fileIDs]
        
        
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


def initializeDir(path, filename, date=True):
    
    dirName = filename.split(".")[0]
    
    if date:
        dirName = add_date(dirName) 
    
    PATH = path + "/" + dirName

    os.mkdir(PATH)

    return PATH    
    

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
    
    fileNames = selectFiles(iniFiles)

    superDir = initializeDir(paths.scratchPath, "QMCrun")
    
    dirs = []
    parsedFiles = []
    for fileName in fileNames:
        #Initialize a new directory for the run
        dirPath = initializeDir(superDir, fileName, date=False)
        dirs.append(dirPath)
        filePath = paths.iniFilePath + "/" + fileName        
        
        #Copy the iniFile to the runDir
        os.system("cp %s %s" % (filePath, dirPath + "/" + fileName))
        iniFile = open(filePath, 'r')
    
        #Read the iniFile
        arglist = []        
        outPathSet = False
        for line in iniFile:

            setTag(arglist, line)
            raw = re.findall(".+\s*=\s*.+", line)
        
            if not outPathSet:
                if arglist[-1] == "-o":
                    arglist.append("outputPath=%s/" % dirPath)
                    outPathSet = True
        
            if raw:
                arglist.append(raw[0].replace(" ", ""))
    
        if not outPathSet:
            arglist.append("-o")
            arglist.append("outputPath=%s/" % dirPath)
            
        parsedFiles.append(convertToCMLargs(arglist))
    
    return parsedFiles, dirs

def valConvert(val):
    
    if re.findall("\d+\.?\d*[E/e][+\-]?\d+", val):
        suff, expo = re.findall("(\d+\.?\d*)[E/e]([+\-]?\d+)", val)[0]
        
        val = str(float(suff)*10**int(expo));
        
    return val
        
        

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


def consistencyCheck(cmlArgs):
    #no minimization initialized -> get param set
    if (cmlArgs[11]=="0" or cmlArgs[11]=="def"):

        n_p = 2;
        dim = 2;
        w = 1;
        system = "QDots"

        if (cmlArgs[4] != "def"):
            n_p = int(cmlArgs[4]);
        if (cmlArgs[5] != "def"):
            dim = int(cmlArgs[5]);
        if (cmlArgs[6] != "def"):
            w = float(cmlArgs[6])
        if (cmlArgs[18] != "def"):
            system = cmlArgs[18]
            
        
        alpha, beta = varParameterMap(n_p, dim, w, system);

        if cmlArgs[41] == "def":
            cmlArgs[41] = str(alpha);
        if cmlArgs[42] == "def":
            cmlArgs[42] = str(beta);

            
    return cmlArgs
            
        


def convertToCMLargs(arglist):
    nInputs = 43;
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
               
    cmlMAPvp = {"alpha":        41,
                "beta":         42}
               
    
    key = ""
    for arg in arglist:
        if arg.startswith("-"):
            key = arg[1:];
        else:
            name, val = arg.split("=");
            index = eval("cmlMAP" + key + "[name]")
            val = valConvert(val)
            cmlArgs[index] = val

    cmlArgs = consistencyCheck(cmlArgs)
    cmlArgs = " ".join(cmlArgs)
    
    return cmlArgs
    
            
def initRuns(CMLargs, stdoutFileFlag, dirs):

    i = 0
    for CMLarg in CMLargs:
        stdout = (" > %s/stdout.txt" % dirs[i])*stdoutFileFlag
        os.system(paths.programPath + "/" + variables.QMC2programName + " " + CMLarg + stdout)
        i+=1

def main():
    stdoutFileFlag = False
    if "stdoutToFile" in sys.argv:
        stdoutFileFlag = True
        sys.argv.remove("stdoutToFile")

    
    if len(sys.argv) == 1:
        CMLargs, dirs = parseFiles()
    else:
        CMLargs = [convertToCMLargs(sys.argv[1:])]
        dirs = [paths.scratchPath]
    

    initRuns(CMLargs, stdoutFileFlag, dirs)
    
main()