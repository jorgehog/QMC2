# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 16:29:12 2012

@author: jorgehog
"""

import sys, os, re, shutil, threading

from pyLibQMC import paths, misc, add_date

from PySide.QtCore import *
from PySide.QtGui import *

openGUI = True
stdoutToFile = False
mpiFlag = True
n_cores = 4

class guiThread(threading.Thread):
    def __init__(self, masterDir):
        super(guiThread, self).__init__()
        self.masterDir = masterDir
        
    def run(self):
        
        os.system("python %s %s > %s" % (os.path.join(paths.toolsPath, 'qmcGUI.py'),\
                                    self.masterDir, \
                                    os.path.join(self.masterDir, "GUI_out.txt")))

            
cmlMAPo = {"dist_out"      : 0,
           "outputPath"    : 1,
           "dmc_out"       : 2,
           "ASGD_out"      : 3}

cmlMAPg = {"n_p"           : 4,
           "dim"           : 5,
           "systemConstant" : 6,
           "random_seed"   : 7,
           "h"             : 8,
           "doMIN"         : 9,
           "doVMC"         : 10,
           "doDMC"         : 11,
           "use_coulomb"   : 12,
           "use_jastrow"   : 13,
           "sampling"      : 14,
           "system"        : 15,
           "do_blocking"   : 16}

cmlMAPv = {"n_c"           : 17,
           "dt"            : 18}

cmlMAPd = {"dt"            : 19,
           "E_T"           : 20,
           "n_b"           : 21,
           "n_w"           : 22,
           "n_c"           : 23,
           "therm"         : 24,
           "dist_in"       : 25,
           "dist_in_path"  : 26}

cmlMAPm = {"max_step"      : 27,
           "f_max"         : 28,
           "f_min"         : 29,
           "omega"         : 30,
           "A"             : 31,
           "a"             : 32,
           "SGDsamples"    : 33,
           "n_w"           : 34,
           "therm"         : 35,
           "n_c"           : 36,
           "n_c_SGD"       : 37,
           "alpha"         : 38,
           "beta"          : 39}

cmlMAPvp = {"alpha"         : 40,
            "beta"          : 41}


def dumpStrList(aList):
    i = 0
    for element in aList:
        print "[%d] %s" % (i, element)
        i+=1

def selectFilesRaw():
    
    iniDirContent = os.listdir(paths.iniFilePath)
    iniFiles = []
    
    for content in iniDirContent:
        if (re.findall("(.+).ini$", content)):
            iniFiles.append(content)

    if len(iniFiles) == 0:
        print "No iniFiles found in %s" % paths.iniFilePath
    else:
        print "Found iniFile(s):\n-----------------------------"    
    
    while True:
        dumpStrList(iniFiles)
    
        action = raw_input("type 'display/run 0' to view content of / run first file:") 
#        action = "run 0 1"

        if action.split()[0] == "run" or action.split()[0] == "display":

            if action.split()[1] == "all":
                fileIDs = range(len(iniFiles))
            else:
                try:
                    fileIDs = [int(ID) for ID in action.split()[1:]]
                except:
                     print "Invalid file", ID
                     legal = False
                     continue
            

            legal = len(fileIDs) > 0
            for ID in fileIDs:
                if fileIDs.count(ID) > 1:
                    print "Several instances of file#", ID
                    legal = False
                    break
                
                
            for ID in fileIDs:
                if ID >= len(iniFiles) or ID < 0:
                    print "Invalid file", ID
                    legal = False
    

            if legal:
                
                
                if action.split()[0] == "run":
                    return [iniFiles[ID] for ID in fileIDs]
                    
                else:
                    for ID in fileIDs:
                        print "----------------------\n" + iniFiles[ID] + ":\n"
                        os.system("cat %s" % paths.iniFilePath + "/" + iniFiles[ID])
                        print "----------------------\n"
         
        else:
            print "Invalid option"


def selectFilesGUI():

    dialog = QFileDialog()
    dialog.setDirectory(paths.iniFilePath)
    dialog.setNameFilter("All ini files (*.ini)")
    dialog.setFileMode(QFileDialog.ExistingFiles)
        
    if dialog.exec_():          
        return [os.path.basename(file_) for file_ in dialog.selectedFiles()]
    
    
    
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
    
    PATH = os.path.join(path, dirName)

    os.mkdir(PATH)

    return PATH    
    

def parseFiles():
    
    
    if openGUI:
        fileNames = None
        app = QApplication(sys.argv)
        while fileNames is None:
            fileNames = selectFilesGUI()
            if fileNames is None:
                print "No files chosen."
        app.exit()
        
    else:
        fileNames = selectFilesRaw()

    #The master directory of runs are set to the scratchPath variable
    superDir = initializeDir(paths.scratchPath, "QMCrun")
    
    dirs = []
    parsedFiles = []
    
    for fileName in fileNames:
        
        #Initialize a new directory for the run
        dirPath = initializeDir(superDir, fileName, date=False)
        dirs.append(dirPath)
        
        #Original ini-file
        filePath = os.path.join(paths.iniFilePath, fileName)        
        
        #Copy the original iniFile to the runDir
        shutil.copy(filePath, os.path.join(dirPath, fileName))
        iniFile = open(filePath, 'r')
    
        #Read the iniFile
        arglist = []        

        for line in iniFile:

            setTag(arglist, line)
            raw = re.findall(".+\s*=\s*.+", line)
      
            if raw:
                arglist.append(raw[0].replace(" ", ""))
    
        #In case the ini file contains no output flags
        try:
            arglist.insert(arglist.index("-o") + 1, "outputPath=%s/" % dirPath)
        except ValueError:
            arglist.append("-o")
            arglist.append("outputPath=%s/" % dirPath)
        
        if "dist_out=0" not in arglist:
            initializeDir(dirPath, "walker_positions", date=False)
        
        parsedFiles.append(convertToCMLargs(arglist))

    return parsedFiles, dirs, superDir

def valConvert(val):
    
    if re.findall("\d+\.?\d*[E/e][+\-]?\d+", val):
        suff, expo = re.findall("(\d+\.?\d*)[E/e]([+\-]?\d+)", val)[0]
        
        val = str(float(suff)*10**int(expo));
        
    return val
        
        

def varParameterMap(n_p, dim, systemConstant, system):
    
    if dim == 2 and system == "QDots":
        
        w = systemConstant
        
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
    if (cmlArgs[cmlMAPg['doMIN']]=="0" or cmlArgs[cmlMAPg['doMIN']]=="def"):

        n_p = 2
        dim = 2
        systemConstant = 1
        system = "QDots"

        if (cmlArgs[cmlMAPg['n_p']] != "def"):
            n_p = int(cmlArgs[cmlMAPg['n_p']])
        if (cmlArgs[cmlMAPg['dim']] != "def"):
            dim = int(cmlArgs[cmlMAPg['dim']])
        if (cmlArgs[cmlMAPg['systemConstant']] != "def"):
            systemConstant = float(cmlArgs[cmlMAPg['systemConstant']])
        if (cmlArgs[cmlMAPg['system']] != "def"):
            system = cmlArgs[cmlMAPg['system']]
            
        
        alpha, beta = varParameterMap(n_p, dim, systemConstant, system)

        if cmlArgs[cmlMAPvp['alpha']] == "def":
            cmlArgs[cmlMAPvp['alpha']] = str(alpha)
        if cmlArgs[cmlMAPvp['beta']] == "def":
            cmlArgs[cmlMAPvp['beta']] = str(beta)
            
     
    #No col -> no jast, alpha=1
    if (cmlArgs[cmlMAPg['use_coulomb']] == "0"):
        cmlArgs[cmlMAPg['use_jastrow']] = "0"
        cmlArgs[cmlMAPvp['alpha']] = "1"
        
    return cmlArgs
            
        


def convertToCMLargs(arglist):
    nInputs = len(cmlMAPo) + len(cmlMAPg) + len(cmlMAPv) \
            + len(cmlMAPd) + len(cmlMAPm) + len(cmlMAPvp) 

    cmlArgs = ["def"]*nInputs

    for arg in arglist:
        
        #Get the correct key
        if arg.startswith("-"):
            key = arg[1:]
        
        #Key is always obtained first. Then we retrieve the correct value
        #from the correct dictionary cmlMAP__key__[name]
        else:
            name, val = arg.split("=");
            index = eval("cmlMAP" + key + "[name]")
            val = valConvert(val)
            cmlArgs[index] = val

    cmlArgs = consistencyCheck(cmlArgs)
    cmlArgs = " ".join(cmlArgs)
    
    return cmlArgs
    
def sendVersion(superDir):
    os.system("git branch -v > " + os.path.join(superDir, "version.txt"))

            
def initRuns(CMLargs, dirs, superDir):

    if openGUI:
        if superDir is None:
            job = guiThread(paths.scratchPath)
        else:            
            job = guiThread(superDir)
        job.start()

    i = 0
    for CMLarg in CMLargs:
        print "Running job ", dirs[i]
        stdout = (" > %s/stdout.txt" % dirs[i])*stdoutToFile
        MPIrun = ("mpiexec -n %d " % n_cores)*mpiFlag        
        
        os.system(MPIrun + os.path.join(paths.programPath, misc.QMC2programName) \
                    + " " + CMLarg + stdout)
    
        i+=1

    print "All queued jobs finished."        
        
    if superDir:
        sendVersion(superDir)
        shutil.copy(os.path.join(paths.toolsPath, "output2tex.py"), superDir)
    if openGUI:
        job.join()
        
    
def getTupleString(pre, suff):
        if pre.startswith('var'):
    	      pre = 'vp'
        else:
            pre = pre[0]
        return '("-%s", "%s")' % (pre, suff)

def getCppMap(raw):
    
    printCppMap = False

    p = "^\s*if \(def\.compare\(argv\[\d+\]\) != 0\) (.+) = "
    vars_ = re.findall(p, raw, re.MULTILINE)
    s = "cppNameMap = {"
    l = len(s)
    i = 0
    for var in vars_:
        s += ''.ljust(l)*(i!=0) + ('"%s"' % var).ljust(32) + " : " \
            + getTupleString(*var.split(".")) + ","*(i!=len(vars_)-1) \
            + "}"*(i==len(vars_)-1) + "\n"
            
        i+=1
    
    exec(s)
    
    if printCppMap:
        for key in cppNameMap.keys():
            print "%s : %s" % (key.ljust(32), cppNameMap[key])
        
    return cppNameMap    
    
def consistentMap():
    
    mainCPPFile = open(os.path.join(paths.CODE, 'src', 'QMCmain.cpp'), 'r')
    p = "^\s*if \(def\.compare\(argv\[(\d+)\]\) != 0\) (.+) = "

    raw = mainCPPFile.read()  
    
    cppOrder = re.findall(p, raw, re.MULTILINE)
    mainCPPFile.close()
    
    cppNameMap = getCppMap(raw)
    
    consistent = True
    for (indexCpp, variable) in cppOrder:
        indexCpp = int(indexCpp)
        key, pyvar = cppNameMap[variable]
        indexPy = eval("cmlMAP" + key[1:] + "[pyvar]") + 1
        
        if indexCpp != indexPy:
            print "---------------------------------------"
            print indexCpp, indexPy
            s = len(str(variable))
            print "MISMATCH: [cpp]    '%s' (index %s) with" % (variable, indexCpp)
            print "          [python] %s   (index %s / %s in code) " % \
            (("'" + key + " " + pyvar + "'").ljust(s) , indexPy, indexPy-1)
            
            consistent = False
    
    return consistent
      

def main():
    
    global stdoutToFile, mpiFlag, openGUI, n_cores
    
    if not consistentMap():
        print "The map is inconsistent with the C++ reader."
        sys.exit(1)
    else:
        print "Map consistent"
    
    #Checking flags
    if "stdoutToFile" in sys.argv:
        stdoutToFile = True
        sys.argv.remove("stdoutToFile")
    if "noMPI" in sys.argv:
        mpiFlag = False
        sys.argv.remove("noMPI")
    if "noGUI" in sys.argv:
        openGUI = False
        sys.argv.remove("noGUI")
    if re.findall("\-n \d+", " ".join(sys.argv)):
        n_cores = int(re.findall("\-n (\d+)", " ".join(sys.argv))[0])
        sys.argv.remove("-n")
        sys.argv.remove(str(n_cores))

    print "MPI nodes: ", n_cores
    
    if len(sys.argv) == 1:
        CMLargs, dirs, superDir = parseFiles()
    else:
        CMLargs = [convertToCMLargs(sys.argv[1:])]
        dirs = [os.path.join(paths.scratchPath, "QMC_SCRATCH")]
        superDir = None
    
    initRuns(CMLargs, dirs, superDir)
  
if __name__ == "__main__":  
    main()