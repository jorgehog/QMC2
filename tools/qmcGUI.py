# -*- coding: utf-8 -*-

import sys
from pyLibQMC import paths

try:
    import DCVizGUI
except:
    print "No DCViz installation found: GUI not supported."
    sys.exit(0)
    
if __name__ == "__main__":
    path = paths.scratchPath
    if len(sys.argv) > 1:
        path = sys.argv[1]
    
    sys.exit(DCVizGUI.main(path))
