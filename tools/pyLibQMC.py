# -*- coding: utf-8 -*-

import os

class paths:
    IDEPath = os.path.expanduser('~') + "/NetBeansProjects/nbQMC2"
    toolsPath = os.path.expanduser('~') + "/MASTER/QMC2/tools"
    
    allPaths = {"toolsPath" : toolsPath, "IDEPath": IDEPath}





def main():
    spacing = 50
    print """
#===========================================
#Paths:
#===========================================
"""    
    print "IDEPath:%s" % paths.IDEPath.rjust(spacing)
    print "toolsPath:%s" % paths.toolsPath.rjust(spacing)
    

if __name__ == "__main__":
    import sys    

    if len(sys.argv) == 2:
        if sys.argv[1] == "-set_tool_path":
            sys.path.append(paths.toolsPath)
            print "path set successfully"
            
    main()
