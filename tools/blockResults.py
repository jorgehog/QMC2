
import sys, os, subprocess, re
from os.path import join as pjoin


from pyLibQMC import paths, parseCML, misc

try:
    from PySide.QtCore import *
    from PySide.QtGui import *
    sys.path.append(pjoin(paths.toolsPath, "DCViz", "GUI"))
    
    forceTerminal = False
       
except:
    
    print "pyside not supported. Terminal usage enabled."
    forceTerminal = True
    
    #Evil haxx
    QDialog = list

try:
    sys.path.append(pjoin(paths.toolsPath, "DCViz", "src"))

    from DCViz_classes import Blocking, MIN_OUT, DMC_OUT
except:
    print "Warning: Displaying results not aviable.."
    

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
        
        self.size = sniffSize(blockFile)   
        
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

def sniffSize(blockFile):
    
    with open(blockFile, 'r') as f:
        f.readline()
        size = f.readline().split()[1]
    f.close()
    
    return size


def selectFileGUI(mainDir):

    dialog = QFileDialog()
    dialog.setDirectory(mainDir)
    dialog.setNameFilter("All arma files (*.arma)")
    dialog.setFileMode(QFileDialog.ExistingFile)
        
    if dialog.exec_():
        
        if dialog.selectedFiles() is not None:
            return dialog.selectedFiles()[0]
            
        return dialog.selectedFiles()



def selectFileFromList(path):
    
    i = 0
    selections = []
    
    print
    for cont in os.listdir(path):
        if os.path.isdir(pjoin(path, cont)):
            print "[%d] %s" % (i, cont)
            selections.append([pjoin(path, cont), "dir"])
            i += 1
            
    for cont in os.listdir(path):
        if re.findall("blocking.+\.arma", cont):
            print "[%d] %s" % (i, cont)
            selections.append([pjoin(path, cont), "file"])
            i += 1
    
    j = raw_input("\nSelect file/folder ('h' for help) #")
    
    try:
        
        if j == "h":

            s = 20         
            
            print "\n%s : up one level" % "'..'".ljust(s)
            print "%s : quit mainloop" % "'q'".ljust(s)
            print "%s : display blocking result n (if any)" % "'display [n]'".ljust(s)
            print "%s : display ASGD result n (if any)" % "'display min [n]'".ljust(s)
            print "%s : display DMC result n (if any)" % "'display dmc [n]'".ljust(s)
                        
            
            selected = selectFileFromList(path)
            
        elif j == "..":
            return [os.path.split(path)[0], "dir"]
            
        elif j == "q":
            return None
            
        elif j.startswith("display"):
            
            if j.split()[1] == "min":
                n = int(j.split()[2])
                
                blockFile = selections[n][0]
      
                if os.path.isdir(blockFile):
                    print "Cannot display directory..."
                else:
                    displayMinTerminal(blockFile)
                    
            elif j.split()[1] == "dmc":
                n = int(j.split()[2])
                
                blockFile = selections[n][0]
      
                if os.path.isdir(blockFile):
                    print "Cannot display directory..."
                else:
                    displayDMCTerminal(blockFile)
                    
            else:
                n = int(j.split()[1])
                blockFile = selections[n][0]
      
                if os.path.isdir(blockFile):
                    print "Cannot display directory..."
                else:
                    displayTerminal(blockFile)
            
            selected = selectFileFromList(path)
                        
        else:
            selected = selections[int(j)]
            
    except:
        print "error in selecting element"
        selected = selectFileFromList(path)

    return selected

def selectFileTerminal(thisDir):

    selected = selectFileFromList(thisDir)

    if selected is None:
         fileName = None
    elif selected[1] == "file":
         fileName = selected[0]
    elif selected[1] == "dir":
        fileName = selectFileTerminal(selected[0])
    else:
        print "Something went wrong.."
  
    return fileName
    


def getParamsGUI(blockFile, lastParams, n_cores):
    
    form = paramForm(blockFile, n_cores, lastParams)
    form.show()
    
    #Has returncode 0 = Success
    if not form.exec_():
        try:
            nBlocks, maxBlockSize, minBlockSize = form.Nb, form.Max, form.Min
            lastParams = [nBlocks, maxBlockSize, minBlockSize]
        except:
            print "No sizes set. Breaking..."   
            return

    return lastParams
    
def getParamsTerminal(blockFile, lastParams, n_cores):
    
    size = int(sniffSize(blockFile))
    
   
    MaxSet = False
    MinSet = False
    NbSet = False
     
   
    while not (MaxSet & MinSet & NbSet):
        
        try:
            
            Max = raw_input("Max block size (< %d) = " % (size/2))
            
            if Max == "q":
                return None
                
            Max = int(Max)
            MaxSet = True
            
        except:
            
            print "Answer could not be cast to integer"
            MaxSet = False
        
        try:
            
            Min =  raw_input("Min block size (> %s) = " % n_cores)
            
            if Min == "q":
                return None
            
            Min = int(Min)
            MinSet = True
            
        except:
            
            print "Answer could not be cast to integer"
            MinSet = False
            
            
        if MaxSet and MinSet:
            
            try:
                
                lim = (Max-Min)/n_cores
                nb = raw_input("nBlocks (2 < nB < %d) = " % lim)
                
                if nb == "q":
                    return None                                
                
                nb = int(nb)
                NbSet = True
                
            except:
                
                print "Answer could not be cast to integer"
                NbSet = False
        
    return [nb, Max, Min]


def getYesNoGUI(Q):
    
    dialog = runAgainForm(Q)
        
    if not dialog.exec_():
        A = dialog.answeredYes
    
    return A

def getYesNoTerminal(Q):
    
    ans = raw_input("\n" + Q + "  ([y]/n): ")
    print
    
    if ans in ["y", "yes", "Y", "YES", "1", "True", "true", ""]:
        return True
        
    elif ans in ["n", "no", "N", "NO", "0", "False", "false"]:
        return False
        
    else:
        print "Answer not understood."
        return getYesNoTerminal(Q)


def displayTerminal(blockFile):
    
    path = blockFile.strip("_RAWDATA.arma") + ".dat"
    
    if not os.path.exists(path):
        print "Spesified file has not been blocked. Display canceled."
        return
    
    displayTool = Blocking(path, dynamic=False)
    displayTool.mainloop()

def displayMinTerminal(blockFile):
    
    path = pjoin(os.path.split(blockFile)[0], "ASGD_out.dat")
    
    if not os.path.exists(path):
        print "Spesified min file does not exist. Display canceled."
        return
    
    displayTool = MIN_OUT(path, dynamic=False)
    displayTool.mainloop()
    
def displayDMCTerminal(blockFile):
    
    path = pjoin(os.path.split(blockFile)[0], "DMC_out.dat")
    
    if not os.path.exists(path):
        print "Spesified dmc file does not exist. Display canceled."
        return
    
    displayTool = DMC_OUT(path, dynamic=False)
    displayTool.mainloop()
    
def getJobFlag(argv):

    for arg in argv:
        if arg == "-JOB":
            argv.remove("-JOB")
            return True            

    return False 

def dumpJobScript(nBlocks, maxBlockSize, minBlockSize, 
                              n_cores, 
                              path,
                              fileName):
                                  
    with open(pjoin(paths.CODE, 'jobScriptTemplate.slurm'), 'r') as f:
        rawJob = f.read()
    f.close()  
    
    rawFile = fileName.strip("_RAWDATA.arma")
    mainDir = os.path.split(path)[1]
   
    args = "reblock %s $SCRATCH/ %d %d %d" % (rawFile,
                                              nBlocks, 
                                              maxBlockSize,
                                              minBlockSize)
    
    rawJob = rawJob.replace("__codeName__", "blocking")
    rawJob = rawJob.replace("__T__", "01:00:00")
    rawJob = rawJob.replace("__MEM__", "3900")
    rawJob = rawJob.replace("__nCpus__", str(n_cores))
    rawJob = rawJob.replace("__superDir__", pjoin(path, fileName))
    rawJob = rawJob.replace("__superName__", rawFile + ".dat")
    rawJob = rawJob.replace("__homeScratch__", path)
    rawJob = rawJob.replace("__code__", pjoin(paths.programPath, misc.QMC2programName))
    rawJob = rawJob.replace("__args__", args)
    rawJob = rawJob.replace(" > __subDir__/stdout.txt", "")
    rawJob = rawJob.replace("__exec__", misc.QMC2programName)
    rawJob = rawJob.replace("__dirName__", mainDir)
    
    with open(pjoin(
                paths.CODE, 
                "jobScripts", 
                "%s_%s.slurm" % (rawFile, mainDir)), 
                        'w') as f:
            f.write(rawJob)
    f.close()
    
def initializeCoumbRuns(mainDir, paramsDict, 
                        n_cores, mpiFlag,
                        makeJobScript):
    
    for root, dir_, files in os.walk(mainDir):
        
        for file_ in files:
            abspath = pjoin(root, file_)            
        
            if re.findall("blocking_.+?\.arma", file_):
                if "DMC" in file_:
                    params = paramsDict["dmc"]
                elif "VMC" in file_:
                    params = paramsDict["vmc"]
                else:
                    continue
                
                initRun(n_cores, mpiFlag,
                        file_, root,
                        abspath,
                        params,
                        makeJobScript)

        
def initRun(n_cores, mpiFlag,
            fileName, path,
            blockFile,
            lastParams,
            makeJobScript,
            display=None):
                
    nBlocks, maxBlockSize, minBlockSize = lastParams
    
    print "Initializing blocking job at "
    print blockFile
    print
    print "with params [n_b, maxSize, minSize] = "
    print lastParams
    print    

    if not makeJobScript:
    
        mpirun = ""
        if mpiFlag:
            mpirun = "mpiexec -n %d" % n_cores
        
        args = mpirun.split() + [pjoin(paths.programPath, misc.QMC2programName), 
                'reblock', fileName.strip("_RAWDATA.arma"), 
                path + "/", str(nBlocks), str(maxBlockSize), str(minBlockSize)]                
        
        subprocess.call(args)
        if display is not None:
            display(blockFile)
        
    else:
        dumpJobScript(nBlocks, maxBlockSize, minBlockSize, 
                      n_cores, 
                      path,
                      fileName)

def main():

    stdoutToFile, mpiFlag, openGUI, n_cores = parseCML(sys.argv)
    
    openGUI = openGUI&(not forceTerminal)

    makeJobScript = getJobFlag(sys.argv)

    if forceTerminal: 
        display = displayTerminal
        selectFile = selectFileTerminal
        getYesNo = getYesNoTerminal
        getParams = getParamsTerminal
        
    else:
        display = lambda a : None
        selectFile = selectFileGUI
        getYesNo = getYesNoGUI
        getParams = getParamsGUI

    print "MPI nodes: ", n_cores
    
    if len(sys.argv) > 1:
        mainDir = sys.argv[1]
        if not os.path.exists(mainDir):
            print "First CML arg must be an existing path."
            sys.exit(1)
            
        fixedParams = {"dmc": [0,0,0], "vmc": [0,0,0]}
        
        #If we supply the correct CML args, we can
        #draw all runs under one comb, saving time.
        if len(sys.argv) > 2:
            for arg in sys.argv[2:]:
                key, val = arg.split("=")
                try:
                    fixedParams[key] = eval(val)
                except:
                    print "Invalid CML arg: ", arg
                    sys.exit(1)
            
            print "successfully initialized blocking runs in path:"
            print mainDir
            print 
            print "with parameters:"
            print fixedParams
            print
            proceed = raw_input("Continue? [y]")
            
            if proceed in ["", "y", "Y", "yes"]:
                initializeCoumbRuns(mainDir, fixedParams, 
                                    n_cores, mpiFlag,
                                    makeJobScript)
                sys.exit(0)
                
            else:
                print "goodbye!"
                sys.exit(1)

            
    else:
        mainDir = pjoin(paths.scratchPath, "QMC_SCRATCH")
    
    if not forceTerminal:
        try:
            app = QApplication(sys.argv)
        except RuntimeError:
            app = QCoreApplication.instance() 
    
    active = True
    GUIopened = False
    while active:
        
        blockFile = selectFile(mainDir)
        
        if not blockFile:
            print "No file selected. Breaking..."
            break
        
        path, fileName = os.path.split(blockFile)
        
        if openGUI and not GUIopened:
            subprocess.Popen(["python", pjoin(paths.toolsPath, 'qmcGUI.py'),\
                path, ">", pjoin(path, "GUI_out.txt")]) 
            GUIopened = True
        
        localSatisfaction = False
        lastParams = None
    
        while not localSatisfaction:
             
            lastParams = getParams(blockFile, lastParams, n_cores)
                
            if lastParams is None:
                break

            initRun(n_cores, mpiFlag,
                    fileName, path,
                    blockFile,
                    lastParams,
                    makeJobScript,
                    display)
                                        
            
            if makeJobScript:
                localSatisfaction = True
            else:
                localSatisfaction = getYesNo("Satisfied with result?")
            
        
        active = getYesNo("View more files?")
    
    if not makeJobScript:
        if getYesNo("Save figs to file?"):
            
            plotTool = Blocking(toFile=True)        
            
            for root, dirs, files in os.walk(mainDir):
                
                if root == mainDir:
                    continue
                
                for outfile in files:
                    if re.findall(plotTool.nametag, outfile):
                        plotTool.filepath = pjoin(root, outfile)
                        plotTool.mainloop()
    
    if not forceTerminal:
        app.exit()


if __name__ == "__main__":
    main()