# -*- coding: utf-8 -*-

import os, time, re
from os.path import join as pjoin

class paths:
    HOME = os.path.expanduser('~')
    CODE = pjoin(HOME ,"MASTER/QMC2")
    IDEPath = pjoin(HOME ,"NetBeansProjects/nbQMC2")
    src = pjoin(CODE, "src")
    toolsPath = pjoin(CODE ,"tools")
    scratchPath = pjoin(HOME ,"scratch")
    iniFilePath = pjoin(CODE ,"iniFiles")
    programPath =pjoin(HOME ,"NetBeansProjects/nbQMC2/dist/Debug/MPI-Linux-x86")
    buildPath = pjoin(HOME, "tmp/testMake/build")    
    #programPath = HOME + "/NetBeansProjects/nbQMC2/dist/Debug/GNU-Linux-x86"
    

class misc:
    QMC2programName = "nbqmc2"


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
    print "HOME:".ljust(spacing) + "%s" % paths.HOME 
    print "CODE:".ljust(spacing) + "%s" % paths.CODE
    print "IDEPath:".ljust(spacing) + "%s" % paths.IDEPath
    print "src:".ljust(spacing) + "%s" % paths.src
    print "toolsPath:".ljust(spacing) + "%s" % paths.toolsPath
    print "scratchPath:".ljust(spacing) + "%s" % paths.scratchPath
    print "iniFiles:".ljust(spacing) + "%s" % paths.iniFilePath
    print "ProgramPath:".ljust(spacing) + "%s" % paths.programPath
    
    print """
#==========================================
# Misc:
#==========================================
"""
    print "QMC2programName:".ljust(spacing) + "%s" % misc.QMC2programName
    

if __name__ == "__main__":               
    main()
