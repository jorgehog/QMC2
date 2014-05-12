# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 16:29:12 2012

@author: jorgehog
"""

import sys, os, re, shutil, subprocess
from inspect import getfile, currentframe
from os.path import join as pjoin, dirname, abspath

try:
	from pyLibQMC import paths, misc, add_date, parseCML
except ImportError:
	print "\n\nYou need to add the QMC paths your PythonPath for some scripts to work."
	_dir = pjoin(dirname(abspath(getfile(currentframe()))), "lib")
	print "\nexport PYTHONPATH=$PYTHONPATH:%s \n" % _dir
	raw_input("Press return to continue.")
	sys.path.append(_dir)
	from pyLibQMC import paths, misc, add_date, parseCML 
except:
     print "Error in pyLibQMC. Run 'python pyLibQMC.py' to debug."
     sys.exit(1)

try:
    from PySide.QtCore import *
    from PySide.QtGui import *
except:
    print "Unable to detect pyside. GUI not supported."


cmlMAPg = {"runpath"       : 0,
           "n_p"           : 1,
           "dim"           : 2,
           "systemConstant" : 3,
           "R"             : 4,
           "random_seed"   : 5,
           "doMIN"         : 6,
           "doVMC"         : 7,
           "doDMC"         : 8,
           "use_coulomb"   : 9,
           "use_jastrow"   : 10,
           "do_blocking"   : 11,
           "sampling"      : 12,
           "system"        : 13,
           "deadlock_x"    : 14}

cmlMAPv = {"n_c"           : 15,
           "dt"            : 16}

cmlMAPd = {"dt"            : 17,
           "n_b"           : 18,
           "n_w"           : 19,
           "n_c"           : 20,
           "therm"         : 21}

cmlMAPm = {"SGDsamples"    : 22,
           "n_w"           : 23,
           "therm"         : 24,
           "n_c_SGD"       : 25,
           "max_step"      : 26,
           "alpha"         : 27,
           "beta"          : 28}

cmlMAPvp = {"alpha"         : 29,
            "beta"          : 30}
            
SYSTEMS = {

    "QDots" : "0",
    "DoubleWell" : "1",
    "Atoms" : "2",
    "Diatom" : "3",
    "QDots3D" : "4"
}

SAMPLING = {
    "IS" : "0",
    "BF" : "1"
}            
            
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
    
    PATH = pjoin(path, dirName)

    os.mkdir(PATH)

    return PATH    

def setOutPath(arglist, dirPath):
    #In case the ini file contains no output flags
    try:
        arglist.insert(arglist.index("-g") + 1, "runpath=%s/" % dirPath)
    except ValueError:
        arglist.append("-g")
        arglist.append("runpath=%s/" % dirPath)

def parseFiles(openGUI, codename): 
    
    if openGUI:
        fileNames = None
        app = QApplication(sys.argv)
        
        fileNames = selectFilesGUI()
            
        app.exit()
        
    else:
        fileNames = selectFilesRaw()

    if fileNames is None:
        print "No files chosen."
        sys.exit(0)

    #The master directory of runs are set to the scratchPath variable
    superDir = initializeDir(paths.scratchPath, codename)
    
    dirs = []
    parsedFiles = []

    for fileName in fileNames:
        
        #Initialize a new directory for the run
        dirPath = initializeDir(superDir, fileName, date=False)
        dirs.append(dirPath)
        
        #Original ini-file
        filePath = pjoin(paths.iniFilePath, fileName)        
        
        #Copy the original iniFile to the runDir
        shutil.copy(filePath, pjoin(dirPath, fileName))
        iniFile = open(filePath, 'r')
    
        #Read the iniFile
        arglist = []        

        for line in iniFile:

            setTag(arglist, line)
            raw = re.findall(".+\s*=\s*.+", line)
      
            if raw:
                arglist.append(raw[0].replace(" ", ""))
    
        
        setOutPath(arglist, dirPath)
        
        parsedFiles.append(convertToCMLargs(arglist))

    return parsedFiles, dirs, superDir

def valConvert(val):
    
    if re.findall("\d+\.?\d*[Ee][+\-]?\d+", val):
        suff, expo = re.findall("(\d+\.?\d*)[Ee]([+\-]?\d+)", val)[0]
        
        val = str(float(suff)*10**int(expo));
        
    return val
        
        

def varParameterMap(n_p, systemConstant, system, R):
    
    if system ==  SYSTEMS["QDots"]:
        
        w = systemConstant
        
        alpha = 0
        beta = 0
        
        if n_p == 2:
            if w == 1.0:
                alpha = 0.98831;
                beta = 0.398664;
            elif w == 0.5:
                alpha = 0.9808;
                beta = 0.30955;
            elif w == 0.28:
                alpha = 0.9716;
                beta = 0.2532;
            elif w == 0.1:
                alpha = 0.9381;
                beta = 0.1808;
            elif w == 0.01:
                alpha = 0.86;
                beta = 0.079;
            elif w == 0.001:
                alpha=0.71;
                beta=0.031;
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 

        elif n_p == 6:
            if w == 1.0:
                alpha = 0.9243;
                beta = 0.5571;
            elif w == 0.5:
                alpha = 0.90004;
                beta = 0.4121;
            elif w == 0.28:
                alpha = 0.8759;
                beta = 0.3257;
            elif w == 0.1:
                alpha = 0.8086;
                beta = 0.2162;
            elif w == 0.01:
                alpha = 0.6;
                beta = 0.096;
            elif w == 0.001:
                alpha = 0.39
                beta = 0.036
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 

        elif n_p == 12:
            if w == 1.0:
                alpha = 0.8756;
                beta = 0.66;
            elif w == 0.5:
                alpha = 0.8432;
                beta = 0.4841;
            elif w == 0.28:
                alpha = 0.8044;
                beta = 0.3734;
            elif w == 0.1:
                alpha = 0.7186;
                beta = 0.2571;
            elif w == 0.01:
                alpha = 0.48; 
                beta = 0.11; 
            elif w == 0.001:
                alpha = 0.26
                beta = 0.05
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 

        elif n_p == 20:
            if w == 1:
                alpha = 0.8361;
                beta = 0.7332;
            elif w == 0.5:
                alpha = 0.8034;
                beta = 0.5356;
            elif w == 0.28:
                alpha = 0.7682;
                beta = 0.4121;
            elif w == 0.1:
                alpha = 0.6756;
                beta = 0.2722;
            elif w == 0.01:
                alpha = 0.39;
                beta = 0.13;
            elif w == 0.001:
                alpha = 0.05
                beta = 0.05
                
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 

        elif n_p == 30:
            if w == 1:
                alpha = 0.8085;
                beta = 0.7944;
            elif w == 0.5:
                alpha = 0.7671;
                beta = 0.5741;
            elif w == 0.28:
                alpha = 0.7229;
                beta = 0.4482;
            elif w == 0.1:
                alpha = 0.6154;
                beta = 0.3001;
            elif w == 0.01:
                alpha = 0.37;
                beta = 0.55;
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 
                
        elif n_p==42:
            if w==1:
                alpha = 0.782778
                beta = 0.84400
            elif w == 0.5:
                alpha = 0.74023;
                beta = 0.610209;
            elif w == 0.28:
                alpha = 0.698021;
                beta = 0.476245;
            elif w == 0.1:
                alpha = 0.59907;
                beta = 0.320888;
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 
            
        elif n_p==56:
            if w==1:
                alpha = 0.76
                beta = 0.886972
            elif w == 0.5:
                alpha = 0.718;
                beta = 0.63;
            elif w == 0.28:
                alpha = 0.6755;
                beta = 0.501;
            elif w == 0.1:
                alpha = 0.5689;
                beta = 0.333102;
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 

        else:
            print "\n\nNo saved parameters for n_p=", n_p , "\n"
            
        return alpha, beta, 0
        
        
    elif system ==  SYSTEMS["Atoms"]:

        Z = systemConstant

        if Z == n_p == 2:
            alpha = 0.922925
            beta = 0.3477
        elif Z == n_p == 4:
            alpha = 0.9915
            beta = 0.95
            
        elif Z == n_p == 10:
            alpha = 1.03
            beta = 0.09

        elif Z == n_p == 12:
            alpha= 1.03
            beta= 0.01024
        
        elif Z == n_p == 18:
            alpha=1.04
            beta=0.15
            
        elif Z == n_p == 36:
            alpha = 1.06
            beta=0.34
        
        else:
            print "\n\nNo saved parameters for n_p=", n_p , "\n"
            return 0, 0
        
        return alpha, beta, 0
        
    elif system ==  SYSTEMS["Diatom"]:
    
        if n_p == 2:
             alpha = 1.35618
             beta = 0.28
             R = 1.4
        elif n_p == 6:
            R = 5.051
            alpha = 0.95
            beta = 0.3568
            
        elif n_p == 8:
            R= 4.63    
            alpha = 1.026
            beta = 0.207

        elif n_p==10:
            R = 3.005
            alpha = 1.032
            beta = 0.2565
            
        elif n_p==12:
            R = 2.3481
            alpha = 1.0367
            beta = 0.2565
            
        elif n_p==14:
            R = 2.068
            alpha = 1.0367
            beta = 0.324            
        
        elif n_p == 16:
            R = 2.282  
            alpha = 1.020
            beta = 0.3638
        elif n_p == 18:
            R = 2.68
            alpha = 1.023
            beta = 0.3638

        

        else:
            print "\n\nNo saved parameters for n_p=", n_p , "\n"
            return 0, 0, 0
        
        return alpha, beta, R
    
    elif system ==  SYSTEMS["DoubleWell"]:
        
        if n_p== 2 and R == 3:
            alpha = 1.36139
            beta = 0.206428
        
        elif n_p == 2 and R == 2:
            alpha = 1.18491
            beta = 0.268434
        
        else:
            print "\n\nNo saved parameters for n_p=", n_p , "R=", R, "\n"
            return 1, 0.4, R            
            
        return alpha, beta, R
        
    elif system == SYSTEMS["QDots3D"]:
        
        w = systemConstant

        alpha = 1
        beta = 0.4        
        
        if n_p == 2:
            if w == 1:
                alpha = 1
                beta = 0.25
            elif w == 0.5:
                alpha=0.995197
                beta=0.202636
            elif w == 0.28:
                alpha=0.988581
                beta=0.163395
            elif w == 0.1:
                alpha=0.963068
                beta = 0.111726
            elif w == 0.01:
                alpha=	0.872512
                beta=0.0474507
            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 
        elif n_p == 8:
            if w == 1:
                alpha = 0.956322  
                beta = 0.377917
            elif w == 0.5:
                alpha=0.944336
                beta=0.275492
            elif w == 0.28:
                alpha = 0.934927
                beta = 0.209711
            elif w == 0.1:
                alpha = 0.911836
                beta = 0.132639
            elif w == 0.01:
                alpha=0.816639
                beta=0.0508982

            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 
                
        elif n_p == 20:
            if w == 1:
                alpha = 0.917
                beta = 0.45
            elif w == 0.5:
                alpha=0.904329
                beta = 0.319537 
            elif w == 0.28:
                alpha=0.891909
                beta=0.244056
            elif w == 0.1:
                alpha = 0.845526
                beta = 0.157512
            elif w == 0.01:
                alpha=0.8
                beta=0.05

            else:
                print "No saved parameters for (n_p, w) = ", n_p, " ", w 
                
        
        else:
            print "\n\nNo saved parameters for n_p=", n_p , "\n"
    
        return alpha, beta, 0    
    
    else:        
        print "\n\nUnknown type ", system, "\n"
        return None


def consistencyCheck(cmlArgs):
  
    #If no system is specified, we select Qdots
    if (cmlArgs[cmlMAPg['system']] != "def"):
        system = cmlArgs[cmlMAPg['system']] 
    else:
        system = SYSTEMS["QDots"]

    #if no number of particles specified -> n_p=2        
    if (cmlArgs[cmlMAPg['n_p']] != "def"):
        n_p = int(cmlArgs[cmlMAPg['n_p']])
    else:
        n_p = 2
    
    #if no systemConstant is selected, we set it to 1 for qdots, and 2 for atoms etc.
    if (cmlArgs[cmlMAPg['systemConstant']] != "def"):
        systemConstant = float(cmlArgs[cmlMAPg['systemConstant']])
    else:
        if system in [SYSTEMS["QDots"], SYSTEMS["QDots3D"], SYSTEMS["DoubleWell"]]:
            systemConstant = 1.0
        elif system in [SYSTEMS["Atoms"], SYSTEMS["Diatom"]]:
            systemConstant = n_p
            cmlArgs[cmlMAPg['systemConstant']] = str(n_p)
        else:
            raise ValueError("unknown system %s" % system)
            sys.exit(1)
  
    #Get the COM distance
    if system in ["Diatom", "DoubleWell"]:
        if cmlArgs[cmlMAPg['R']] == "def":
            R = 3
        else:
            R = float(cmlArgs[cmlMAPg['R']])
    else:
        R = None
  
  
    alpha, beta, R = varParameterMap(n_p, systemConstant, system, R)

    #No col -> no jast and alpha=1
    if (cmlArgs[cmlMAPg['use_coulomb']] == "0"):
        cmlArgs[cmlMAPg['use_jastrow']] = "0"
    
        if (cmlArgs[cmlMAPvp['alpha']] == "def" or cmlArgs[cmlMAPvp['alpha']] == "0"):
            alpha=1
            cmlArgs[cmlMAPvp['alpha']] = "1"
                 
         
  
    if cmlArgs[cmlMAPvp['alpha']] == "def":
        cmlArgs[cmlMAPvp['alpha']] = str(alpha)
    if cmlArgs[cmlMAPvp['beta']] == "def":
            cmlArgs[cmlMAPvp['beta']] = str(beta)
    if cmlArgs[cmlMAPg['R']] == "def":
        cmlArgs[cmlMAPg['R']] = str(R)
   
    #In case of minimization, we set the inital params equal to the saved if any.
    if cmlArgs[cmlMAPg['doMIN']]=="1":
        if cmlArgs[cmlMAPm['alpha']] == 'def' and cmlArgs[cmlMAPm['beta']] == 'def':
            if not alpha == 0 and not beta == 0:
                cmlArgs[cmlMAPm['alpha']] = str(alpha)
                cmlArgs[cmlMAPm['beta']] = str(beta)
        
    
            
        
     

    return cmlArgs
            
        


def convertToCMLargs(arglist):
    nInputs = len(cmlMAPg) + len(cmlMAPv) \
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
            if name == "sampling" and key == "g":
                val = SAMPLING[val]
            elif name == "system" and key == "g":
                val = SYSTEMS[val]
                
            index = eval("cmlMAP" + key + "[name]")
            val = valConvert(val)
            cmlArgs[index] = val

    cmlArgs = consistencyCheck(cmlArgs)
    cmlArgs = " ".join(cmlArgs)
    
    return cmlArgs
    
def sendVersion(superDir):
    os.system("git branch -v > " + pjoin(superDir, "version.txt"))

def createDistFolder(dirPath):
    
    if dirPath == "def":
        print "Something's wrong... path should never be default."
        return
        
    if not os.path.exists(pjoin(dirPath, "walker_positions")):
        #print "making folder %s" % pjoin(dirPath, "walker_positions")
        initializeDir(dirPath, "walker_positions", date=False)
         
def generateJobScript(Args, path, n_cores):

        with open(pjoin(paths.CODE, 'jobScriptTemplate.slurm'), 'r') as f:
            rawJob = f.read()
        f.close()
        
        superDir, subdir = os.path.split(path)
        superName = os.path.split(superDir)[1]   
        subdirNew = "$SCRATCH/" + superName + "/" + subdir
        
        rawJob = rawJob.replace("__dirName__", subdir)
        rawJob = rawJob.replace("__codeName__", subdir)
        rawJob = rawJob.replace("__nCpus__", str(n_cores))
        rawJob = rawJob.replace("__superDir__", superDir)
        rawJob = rawJob.replace("__superName__", superName)
        rawJob = rawJob.replace("__code__", pjoin(paths.programPath, misc.QMC2programName))
        rawJob = rawJob.replace("__exec__", misc.QMC2programName)
        rawJob = rawJob.replace("__subDir__", subdirNew)     
#        rawJob = rawJob.replace("__subdirOld__", path)
        
        Args = Args.replace(path + "/", subdirNew + "/")
        rawJob = rawJob.replace("__args__", Args)
        rawJob = rawJob.replace("__homeScratch__", paths.scratchPath)
        
        time = raw_input("Expected time? [hh:mm:ss]: ")
        mem = raw_input("Memory use pr CPU? max abel=3900 [Mb]: ")
        mem += "M"
        rawJob = rawJob.replace("__T__", time)     
        rawJob = rawJob.replace("__MEM__", mem)
        
        
        with open(pjoin(paths.CODE, "jobScripts", subdir + ".slurm"), 'w') as f:
            f.write(rawJob)
        f.close()

        
         
         
def initRuns(CMLargs, dirs, superDir, stdoutToFile, mpiFlag, openGUI, n_cores, makeJobScript):

    if openGUI:
        if superDir is None:
            jobDir = paths.scratchPath
        else:            
            jobDir = superDir
        subprocess.Popen(["python", pjoin(paths.toolsPath, 'qmcGUI.py'), jobDir, "> " + pjoin(jobDir, "GUI_out.txt")])

    i = 0
    for CMLarg in CMLargs:
        print "Running job ", dirs[i]
        
        createDistFolder(CMLarg.split()[cmlMAPg["runpath"]])               
        if makeJobScript and superDir is not None:
            
            generateJobScript(CMLarg, dirs[i], n_cores)            
            
        else:
            stdout = (" > %s/stdout.txt" % dirs[i])*stdoutToFile
            MPIrun = ("mpiexec -n %d " % n_cores)*mpiFlag        
        
            os.system(MPIrun + pjoin(paths.programPath, misc.QMC2programName) \
                        + " " + CMLarg + stdout)
    
        i+=1

    print "All queued jobs finished."        
        
    if superDir:
        sendVersion(superDir)
        
        if stdoutToFile and not makeJobScript:
            subprocess.call(["python", pjoin(paths.toolsPath, "output2tex.py"), superDir])
            os.system('cat %s' % pjoin(superDir, 'texTable.tex'))
        #shutil.copy(pjoin(paths.toolsPath, "output2tex.py"), superDir)
        
    
def getTupleString(pre, suff):
        if pre.startswith('var'):
    	      pre = 'vp'
        elif pre == "mP":
            pre = 'g'
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
    
    mainCPPFile = open(pjoin(paths.CODE, 'apps', 'default', 'defaultmain.cpp'), 'r')
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
      
def getCodename(argv):
    
    codename = "QMCrun"
    for arg in argv:
        p = re.findall("codename\s*\=\s*(.+)", arg)
        
        if p:
            codename = p[0]
            argv.remove(arg)            
            
    return codename
    
def getJobFlag(argv):

    for arg in argv:
        if arg == "-JOB":
            argv.remove("-JOB")
            return True            

    return False
    

def main():
    
    if not consistentMap():
        print "The map is inconsistent with the C++ reader."
        ans = raw_input("Proceed? [y]/n")
        if ans in ["n", "no"]:
            sys.exit(1)

    else:
        print "Map consistent"
    
    stdoutToFile, mpiFlag, openGUI, n_cores = parseCML(sys.argv)
    codename = getCodename(sys.argv)
    makeJobScript = getJobFlag(sys.argv)

    print "MPI nodes: ", n_cores
    
    if len(sys.argv) == 1:
        CMLargs, dirs, superDir = parseFiles(openGUI, codename)
    else:
        dirs = [pjoin(paths.scratchPath, "QMC_SCRATCH")]
        setOutPath(sys.argv, dirs[0])
        CMLargs = [convertToCMLargs(sys.argv[1:])]
        
        superDir = None
    
    initRuns(CMLargs, 
             dirs, 
             superDir, 
             stdoutToFile, 
             mpiFlag,
             openGUI, 
             n_cores,
             makeJobScript)
  
  
if __name__ == "__main__":  
    main()
    
