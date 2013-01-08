# -*- coding: utf-8 -*-

import sys, os
from pyLibQMC import paths

sys.path.append(os.path.join(paths.toolsPath, "DCV", "GUI"))

import DCVizGUI

if __name__ == "__main__":
	DCVizGUI.main(paths.scratchPath)
