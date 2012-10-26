# -*- coding: utf-8 -*-

import os, time

class paths:
    HOME = os.path.expanduser('~')
    CODE = HOME +  "/scratch/tryfixQMC"
    IDEPath = HOME + "/NetBeansProjects/nbQMC_debug/nbQMC2_debug"
    toolsPath = CODE + "/tools"
    scratchPath = HOME + "/scratch/debug"
    iniFilePath = CODE + "/iniFiles"
    programPath = HOME + "/NetBeansProjects/nbQMC_debug/nbQMC2_debug/dist/Debug/GNU-Linux-x86"
    

class variables:
    QMC2programName = "nbqmc2_debug"


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
    

def main():
    spacing = 100
    print """
#===========================================
#Paths:
#===========================================
"""   
    print "HOME:\t\t%s" % paths.HOME.ljust(spacing) 
    print "CODE:\t\t%s" % paths.CODE.ljust(spacing)
    print "IDEPath:\t\t%s" % paths.IDEPath.ljust(spacing)
    print "toolsPath:\t\t%s" % paths.toolsPath.ljust(spacing)
    print "scratchPath:\t%s" % paths.scratchPath.ljust(spacing)
    print "iniFiles:\t\t%s" % paths.iniFilePath.ljust(spacing)
    print "ProgramPath:\t%s" % paths.programPath.ljust(spacing)
    
    print """
#==========================================
# Variables:
#==========================================
"""
    print "QMC2programName:\t%s" % variables.QMC2programName.ljust(spacing)
    

if __name__ == "__main__":
    import sys    

    if len(sys.argv) == 2:
        if sys.argv[1] == "-set_tool_path":
            sys.path.append(paths.toolsPath)
            print "path set successfully"
            
    main()
