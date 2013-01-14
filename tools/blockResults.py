
import sys, os, subprocess, re
from os.path import join as pjoin

from PySide.QtCore import *
from PySide.QtGui import *

from pyLibQMC import paths, parseCML, misc

sys.path.append(pjoin(paths.toolsPath, "DCViz", "GUI"))

import DCVizGUI

def selectFile():

    dialog = QFileDialog()
    dialog.setDirectory(pjoin(paths.scratchPath, "QMC_SCRATCH"))
    dialog.setNameFilter("All arma files (*.arma)")
    dialog.setFileMode(QFileDialog.ExistingFile)
        
    if dialog.exec_():          
        return dialog.selectedFiles()


class paramForm(QDialog):
    
    def __init__(self, blockFile, n_cores, parent=None):
        super(paramForm, self).__init__(parent)
        
        self.n_cores = n_cores        
        
        #Sniff number of samples
        with open(blockFile, 'r') as f:
            f.readline()
            self.size = f.readline().split()[1]
        
        
        self.headLabel = QLabel("Number of samples: %g" % int(self.size), self)        
        
        self.runButton = QPushButton("Run")
            
        self.labelNb = QLabel("<font>Nblocks:</font>", self)
        self.labelMax = QLabel("<font>Max size:</font>", self)
        self.labelMin = QLabel("<font>Min size:</font>", self)
        
        self.editNb = QLineEdit(str(int(self.size)/1000), self)
        self.labelNb.setBuddy(self.editNb)
        
        self.editMax = QLineEdit(str(int(self.size)/2), self)
        self.labelMax.setBuddy(self.editMax)

        self.editMin = QLineEdit(str(n_cores), self)
        self.labelMin.setBuddy(self.editMin)        
        
        
        layout = QVBoxLayout()

        layout.addWidget(self.headLabel)     
        
        layout.addWidget(self.labelNb)
        layout.addWidget(self.editNb)
        
        layout.addWidget(self.labelMin)       
        layout.addWidget(self.editMin)
       
        layout.addWidget(self.labelMax)
        layout.addWidget(self.editMax)
        
        layout.addStretch()
        
        layout.addWidget(self.runButton)
         
        
        self.setLayout(layout)
        
        
        self.runButton.clicked.connect(self.applyAndRun)
        
        
        
        
    def applyAndRun(self):
    
        isOk = True
        for label, editLine in [(self.labelNb, self.editNb), 
                                (self.labelMax, self.editMax), 
                                (self.labelMin, self.editMin)]:
            try:
                a = int(editLine.text())

                if int(a) != a or a < 1:
                    int("going to exception")
         
                entry = re.findall(r">(.+)<", label.text())[0]
                html = "<font color=green>%s</font>" % entry
                label.setText(html)
                
                
            except:
   
                entry = re.findall(r">(.+)<", label.text())[0]
                html = "<font color=red>%s</font>" % entry
                label.setText(html)
                isOk = False
                
        if not isOk:
            return
         
        self.Nb = int(re.findall("\d+", self.editNb.text())[0])
        self.Max = int(re.findall("\d+", self.editMax.text())[0])
        self.Min = int(re.findall("\d+", self.editMin.text())[0])
        
        
  
        warnedMinOrMaxed = False
        if self.Max > int(self.size)/2:
            self.warning("Max cannot be greater than half the sample size.")
            entry = re.findall(r">(.+)<", self.labelMax.text())[0]
            html = "<font color=red>%s</font>" % entry
            self.labelMax.setText(html)
            warnedMinOrMaxed = True
            
        if self.Min/self.n_cores < 1:
            self.warning("Min cannot be lower than n_cores=%d" % self.n_cores)
            entry = re.findall(r">(.+)<", self.labelMin.text())[0]
            html = "<font color=red>%s</font>" % entry
            self.labelMin.setText(html)
            warnedMinOrMaxed = True
            
        if self.Nb > self.Max/self.n_cores - self.Min/self.n_cores and not warnedMinOrMaxed:
            self.warning("Number of blocks cannot proceed %d" % \
                          (self.Max/self.n_cores - self.Min/self.n_cores))
            entry = re.findall(r">(.+)<", self.labelNb.text())[0]
            html = "<font color=red>%s</font>" % entry
            self.labelNb.setText(html)
            return

        elif self.Nb < 2:
            self.warning("Number of blocks cannot be lower than 2")
            entry = re.findall(r">(.+)<", self.labelNb.text())[0]
            html = "<font color=red>%s</font>" % entry
            self.labelNb.setText(html)
            return
            
            
        if warnedMinOrMaxed:
            return
            
            
        
        self.close()
        
    def warning(self, m):
        QMessageBox.warning(None, "Blocking parameters", m, QMessageBox.Ok)


def main():
    
    stdoutToFile, mpiFlag, openGUI, n_cores = parseCML(sys.argv)

    print "MPI nodes: ", n_cores
    
    app = QApplication(sys.argv)
    
    active = True
    while active:
        blockFile = selectFile()
        
        if not blockFile:
            print "No file selected. Breaking..."
            break
        
        blockFile = blockFile[0]
        
        path, fileName = os.path.split(blockFile)
        
        form = paramForm(blockFile, n_cores)
        form.show()
        
        #Has returncode 0 = Success
        if not form.exec_():
            try:
                nBlocks, maxBlockSize, minBlockSize = form.Nb, form.Max, form.Min
            except:
                print "No sizes set. Breaking..."   
                break                 
        
        mpirun = ""
        if mpiFlag:
            mpirun = "mpiexec -n %d" % n_cores
            
        
        args = mpirun.split() + [pjoin(paths.programPath, misc.QMC2programName), 
                'reblock', fileName.strip("_RAWDATA.arma"), 
                path + "/", str(nBlocks), str(maxBlockSize), str(minBlockSize)]
        
   
        
        if openGUI:
            subprocess.Popen(["python", pjoin(paths.toolsPath, 'qmcGUI.py'),\
            path, ">", pjoin(path, "GUI_out.txt")])
        
        subprocess.call(args)
        active = False
        



    app.exit()        


if __name__ == "__main__":
    main()