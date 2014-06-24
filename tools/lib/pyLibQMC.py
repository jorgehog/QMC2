# -*- coding: utf-8 -*-

import os, time, re, sys, shutil
from inspect import getfile, currentframe
from os.path import join, split, dirname, abspath



class paths:
    lib = dirname(abspath(getfile(currentframe())))
    toolsPath = split(lib)[0]
    CODE = split(toolsPath)[0]
    src = join(CODE, "src")
    scratchPath = join(CODE ,"scratch")
    iniFilePath = join(CODE ,"iniFiles")
    programPath = "/home/jorgen/code/build-QMC2-Desktop_Qt_5_2_1_GCC_64bit-Release/apps/default"
    
    

class misc:
    QMC2programName = "QMC2_default"
    
class verifyPaths:

    l1_paths = [paths.CODE, paths.src, paths.toolsPath]    
    l2_paths = [paths.iniFilePath, paths.scratchPath, paths.programPath, join(paths.scratchPath, "QMC_SCRATCH"), join(paths.scratchPath, "QMC_SCRATCH", "walker_positions")]

    
    def l1(self, path):
        print "\nNecessary path: %s doesn't exist on the file system.\n" % path
        sys.exit(1)
    
    def l2(self, path):
        print "\nNecessary path: %s doesn't exist on the file system." % path
        a = raw_input("Do you want to create it[y]?")
        if a in ["", "Y", "y", "yes", "Yes"]:
            os.mkdir(path)


    def __init__(self):
        
        for path in self.l1_paths:
            if not os.path.exists(path):
                self.l1(path)
        for path in self.l2_paths:
            if not os.path.exists(path):
                self.l2(path)
                
        
                
verifyPaths()
            

def add_date(filename):

    #in case the file has e.g. a .txt, we want to sandwitch the date and not
    #append it directly
    fileEnding = ""
    if len(filename.split('.')) > 1:   
        fileEnding += "." + ".".join(filename.split('.')[1:])
    
    #asctime converts 'raw local time' to a more refined format.
    date = time.asctime(time.localtime()).replace(' ','_')

    originalFilename = filename.split('.')[0]


    return originalFilename + date + fileEnding;
    
def parseCML(argv):
    
    #Checking flags
    if "stdoutToFile" in argv:
        stdoutToFile = True
        argv.remove("stdoutToFile")
    else:
        stdoutToFile = False
        
    if "noMPI" in argv:
        mpiFlag = False
        argv.remove("noMPI")
    else:
        mpiFlag = True
        
    if "noGUI" in argv:
        openGUI = False
        argv.remove("noGUI")
    else:
        openGUI = True
        
    if re.findall("\-n \d+", " ".join(argv)):
        n_cores = int(re.findall("\-n (\d+)", " ".join(argv))[0])
        argv.remove("-n")
        argv.remove(str(n_cores))
    else:
        n_cores = 4

    if mpiFlag is False:
        n_cores = 1
    
    return stdoutToFile, mpiFlag, openGUI, n_cores

def main():
    spacing = 20
    print """
#===========================================
#Paths:
#===========================================
"""   
    print "CODE:".ljust(spacing) + "%s" % paths.CODE
    print "src:".ljust(spacing) + "%s" % paths.src
    print "toolsPath:".ljust(spacing) + "%s" % paths.toolsPath
    print "scratchPath:".ljust(spacing) + "%s" % paths.scratchPath
    print "iniFiles:".ljust(spacing) + "%s" % paths.iniFilePath
    print """
#==========================================
# Misc:
#==========================================
"""
    print "QMC2programName:".ljust(spacing) + "%s" % misc.QMC2programName
    

if __name__ == "__main__":               
    main()
