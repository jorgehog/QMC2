# -*- coding: utf-8 -*-

import sys, os
from pyLibQMC import paths

sys.path.append(os.path.join(paths.toolsPath, "DCViz", "GUI"))

import DCVizGUI

if __name__ == "__main__":
    path = paths.scratchPath
    if len(sys.argv) > 1:
        path = sys.argv[1]
    
    DCVizGUI.main(path)
