# -*- coding: utf-8 -*-

import sys, os

class paths:
    HOME = os.path.expanduser('~')
    CODE = HOME +  "/MASTER/QMC2"
    IDEPath = HOME + "/NetBeansProjects/nbQMC2"
    toolsPath = CODE + "/tools"
    scratchPath = HOME + "/scratch"
    iniFilePath = CODE + "/iniFiles"
    programPath = HOME + "/NetBeansProjects/nbQMC2/dist/Debug/GNU-Linux-x86"
    

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

    if len(sys.argv) == 2:
        if sys.argv[1] == "-set_tool_path":
            sys.path.append(paths.toolsPath)
            print "path set successfully"
            
    main()
