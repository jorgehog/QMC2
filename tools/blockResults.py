
import sys, os, subprocess, re
from os.path import join as pjoin

from PySide.QtCore import *
from PySide.QtGui import *

from pyLibQMC import paths, parseCML, misc

sys.path.append(pjoin(paths.toolsPath, "DCViz", "GUI"))

import DCVizGUI

def selectFile(mainDir):

    dialog = QFileDialog()
    dialog.setDirectory(mainDir)
    dialog.setNameFilter("All arma files (*.arma)")
    dialog.setFileMode(QFileDialog.ExistingFile)
        
    if dialog.exec_():          
        return dialog.selectedFiles()


class runAgainForm(QDialog):
    def __init__(self, Q, parent=None):
        super(runAgainForm, self).__init__(parent)

        self.setWindowTitle(Q)
        self.answeredYes = None
        
        self.yes = QPushButton("Yes", self)
        self.yes.clicked.connect(lambda : self.answer(True))
        self.no = QPushButton("No", self)
        self.no.clicked.connect(lambda : self.answer(False))
        
               
        
        layout = QHBoxLayout()
        layout.addWidget(self.yes)
        layout.addWidget(self.no)
        self.setLayout(layout)
        self.resize(300, 50)
        self.show()
        
    def closeEvent(self, event):
        
        if self.answeredYes is None:
            self.answeredYes = False
        
        QDialog.closeEvent(self, event)
        
    def answer(self, ans):
        self.answeredYes = ans
        self.close()
        

                
        

class paramForm(QDialog):
    
    def __init__(self, blockFile, n_cores, lastParams, parent=None):
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
        
        
        if lastParams is None:
            self.editNb = QLineEdit(str(int(self.size)/1000), self)
            self.editMax = QLineEdit(str(int(self.size)/2), self)
            self.editMin = QLineEdit(str(n_cores), self)
        else:
            self.editNb = QLineEdit(str(lastParams[0]), self)
            self.editMax = QLineEdit(str(lastParams[1]), self)
            self.editMin = QLineEdit(str(lastParams[2]), self)
        
        
        self.labelNb.setBuddy(self.editNb)
        self.labelMax.setBuddy(self.editMax)
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
    
    if len(sys.argv) > 1:
        mainDir = sys.argv[1]
    else:
        mainDir = pjoin(paths.scratchPath, "QMC_SCRATCH")
    
    app = QApplication(sys.argv)
    
    active = True
    GUIopened = False
    while active:
        
        blockFile = selectFile(mainDir)
        
        if not blockFile:
            print "No file selected. Breaking..."
            break
        
        blockFile = blockFile[0]
        
        path, fileName = os.path.split(blockFile)
        
        if openGUI and not GUIopened:
            subprocess.Popen(["python", pjoin(paths.toolsPath, 'qmcGUI.py'),\
                path, ">", pjoin(path, "GUI_out.txt")]) 
            GUIopened = True
        
        localSatisfaction = False
        lastParams = None
    
        while not localSatisfaction:
            
            form = paramForm(blockFile, n_cores, lastParams)
            form.show()
            
            #Has returncode 0 = Success
            if not form.exec_():
                try:
                    nBlocks, maxBlockSize, minBlockSize = form.Nb, form.Max, form.Min
                    lastParams = [nBlocks, maxBlockSize, minBlockSize]
                except:
                    print "No sizes set. Breaking..."   
                    break                 
            
            mpirun = ""
            if mpiFlag:
                mpirun = "mpiexec -n %d" % n_cores
            
            args = mpirun.split() + [pjoin(paths.programPath, misc.QMC2programName), 
                    'reblock', fileName.strip("_RAWDATA.arma"), 
                    path + "/", str(nBlocks), str(maxBlockSize), str(minBlockSize)]
            
            subprocess.call(args)
            
            dialog = runAgainForm("Satisfied with result?")
        
            if not dialog.exec_():
                localSatisfaction = dialog.answeredYes
            
        
        dialog = runAgainForm("View more files?")
        if not dialog.exec_():
            active = dialog.answeredYes
        



    app.exit()        


if __name__ == "__main__":
    main()