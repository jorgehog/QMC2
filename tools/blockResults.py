
import sys, os, re

from PySide.QtCore import *
from PySide.QtGui import *

from pyLibQMC import paths

sys.path.append(os.path.join(paths.toolsPath, "DCViz", "GUI"))

import DCVizGUI

openGUI = True
stdoutToFile = False
mpiFlag = True
n_cores = 4

def selectFile():

    dialog = QFileDialog()
    dialog.setDirectory(os.getcwd())
    dialog.setNameFilter("All arma files (*.arma)")
    dialog.setFileMode(QFileDialog.ExistingFile)
        
    if dialog.exec_():          
        return dialog.selectedFiles()

def main():
    
    global stdoutToFile, mpiFlag, openGUI, n_cores
    
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
    
    app = QApplication(sys.argv)
    
    active = True
    while active:
        blockFile = selectFile()
        
        if not blockFile:
            print "No file selected. Breaking..."
            break
        
        
            




    app.exit()        


if __name__ == "__main__":
    main()