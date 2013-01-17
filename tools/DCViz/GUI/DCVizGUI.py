# -*- coding: utf-8 -*-

import sys, os, re, time, inspect

try:
    from PySide.QtCore import *
    from PySide.QtGui import *
except:
    print "\n"
    print "You need pyside in order to run the GUI!"
    print "sudo apt-get install pyside-tools"
    print "\n"
    sys.exit(1)
try:
    import matplotlib
except:
    print "\n"
    print "You need matplotlib in order to run this library!"
    print "sudo apt-get install python-matplotlib"
    print "\n" 
    sys.exit(1)

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

#Adding the srcDir to the local pythonpath in order to avoid global pythonpath sets
thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
srcDir = os.path.join(os.path.dirname(thisDir), "src")
sys.path.append(srcDir)
from DCViz_classes import *



class ThreadComm(QObject):
    stopSignal = Signal()
    plotSignal = Signal()
    setFigureSignal = Signal()
    beginSignal = Signal()

class jobThread(QThread):
    def __init__(self, mode, parent=None):
        QThread.__init__(self, parent)
        self.mode = mode
        self.mode.parent = parent
        
        self.comm = ThreadComm()
        self.comm.stopSignal.connect(self.stop)
        self.comm.beginSignal.connect(self.begin)
        
    def run(self):
        self.mode.canStart = False
        self.mode.stopped = False
        self.mode.mainloop()
            
    def stop(self):
        self.mode.stopped = True
    def begin(self):
        self.mode.canStart = True
    


class DCVizPlotWindow(QMainWindow):
    
    def __init__(self, activeMode, parent = None):
        super(DCVizPlotWindow, self).__init__()
        
        self.comm = ThreadComm()
        self.comm.plotSignal.connect(self.displayPlots)
        self.comm.setFigureSignal.connect(self.setFigures)
        
        self.parent = parent
        self.activeMode = activeMode
        self.callStop = True
        self.setAttribute(Qt.WA_DeleteOnClose)

    def displayPlots(self):
  
        for k in range(len(self.activeMode.figures)):
            fig = self.activeMode.figures[k]
            canvas = FigureCanvas(fig[0])
            dockWidget = QDockWidget(self)
            dockWidget.setWidget(canvas)
         
            self.addDockWidget(Qt.DockWidgetArea(1), dockWidget)

        self.show()

    def setFigures(self):
        
        self.activeMode.set_figures()
        self.parent.job.comm.beginSignal.emit()
        
        
    
    def closeEvent(self, event):
   
        if self.activeMode.dynamic and self.parent is not None:
            if self.callStop:
                self.parent.stop()
        
        self.activeMode.close()
        self.parent.plotWinClosed = True
            
        QMainWindow.closeEvent(self, event)


class DCVizGUI(QMainWindow):

    def __init__(self, masterDir):
        
        super(DCVizGUI, self).__init__()
        
        self.masterDir = masterDir
        
        self.comm = ThreadComm()
        
        self.loadDefaults()

        self.setupUI()
        
    def loadDefaults(self):
         
        self.path = os.getcwd()
        self.loadImages()
        
        self.activeMode = None
        self.plotWin = None
        self.job = None
        
        self.dynamic = False
        self.started = False
        self.plotWinClosed = True
        
        #Initial terminal output flag
        self.hideSource = False
        
        self.terminalSilence = False
        self.warningSilence = False

        self.loadExtern()
       
    def setupUI(self):
        
        menubar = self.menuBar()

        #File menu
        filemenu = menubar.addMenu('&File')

        openAction = QAction('&Open', self)
        openAction.setShortcut('Ctrl+O')
        openAction.setStatusTip('Open single file')
        openAction.triggered.connect(self.openfile)
        filemenu.addAction(openAction)

        pathAction = QAction('&Open path', self)
        pathAction.setShortcut('Ctrl+S')
        pathAction.setStatusTip('Open folder')
        pathAction.triggered.connect(self.setpath)
        filemenu.addAction(pathAction)

        filemenu.addSeparator()

        cleanAction = QAction("&Reset", self)
        cleanAction.triggered.connect(self.resetModeSelector)
        filemenu.addAction(cleanAction)
        
        filemenu.addSeparator()
        
        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Esc')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        filemenu.addAction(exitAction)

        
        #Options menu
        optmenu = menubar.addMenu('&Options')

        showConfigAction = QAction('&Display Config', self)
        showConfigAction.setShortcut('Ctrl+C')
        showConfigAction.triggered.connect(lambda : self.showInfo("Config", self.config))
        optmenu.addAction(showConfigAction)

        reloadConfigAction = QAction('&Reload Config', self)
        reloadConfigAction.setShortcut('Ctrl+R')
        reloadConfigAction.triggered.connect(lambda : self.loadExtern(True))
        optmenu.addAction(reloadConfigAction)
        #::::::::::::::::::::::::::::::::::::::::::::::::::

        #  StartButton ::::::::::::::::::::::::::::::::::::
        self.startStopButton = QPushButton(self)
        self.startStopButton.setIcon(self.img["play"])
        self.startStopButton.clicked.connect(self.startOrStop)
        self.startStopButton.setFlat(True)
        #::::::::::::::::::::::::::::::::::::::::::::::::::
    
        #  DynamicCheckBox ::::::::::::::::::::::::::::::::
        self.dynamicCheckBox = QCheckBox("Real-time", self)
        self.dynamicCheckBox.stateChanged.connect(self.flipDynamic)
        #::::::::::::::::::::::::::::::::::::::::::::::::::



        #  Mode selection dropdown menu :::::::::::::::::::
        self.modeSelector = QComboBox(self)
        self.resetModeSelector()
        #::::::::::::::::::::::::::::::::::::::::::::::::::


        #  Refresh timer ::::::::::::::::::::::::::::::::::
        self.dtSlider = QSlider(Qt.Horizontal, self)
        self.dtSlider.setMinimum(0) 
        self.dtSlider.setMaximum(100)
        self.dtSlider.valueChanged.connect(lambda : self.updateDt(False))
        self.dtSlider.sliderReleased.connect(lambda : self.updateDt(True))
       
        
        self.sliderUnit = QLabel("<font size=3>s</font>", self)
        
        self.sliderDisplay = QLineEdit(self)
        self.sliderDisplay.setReadOnly(True)

        #Init values
        self.dtSlider.setValue(self.refreshDtConfig*10)
        if (int(self.refreshDtConfig) - self.refreshDtConfig) == 0:
            initVal = str(int(self.refreshDtConfig))
        else:
            initVal = str(self.refreshDtConfig)
        self.sliderDisplay.clear()
        self.sliderDisplay.insert(str(initVal))

        #init as disabled
        self.dtSlider.setEnabled(False)
        self.sliderDisplay.setEnabled(False)
        self.sliderUnit.setEnabled(False)
        #::::::::::::::::::::::::::::::::::::::::::::::::::


        #  GUI Setup ::::::::::::::::::::::::::::::::::::::
        initPos = (300, 300)
        size = (300, 150)
        buttonSize = QSize(50,50)
        
        self.setGeometry(*(initPos+size))
 
        

        self.startStopButton.resize(buttonSize)
        self.startStopButton.setIconSize(buttonSize)
        self.startStopButton.move(30,40)


        self.modeSelector.resize(QSize(165,25))
        self.modeSelector.move(110,53)


        self.dynamicCheckBox.move(20, size[1] - 50)


        self.dtSlider.resize(QSize(size[0]-55, 20))
        self.dtSlider.move(0, size[1]-20)
        
        self.sliderDisplay.resize(QSize(30, 20))
        self.sliderDisplay.move(size[0]-52, size[1]-20)
        
        self.sliderUnit.resize(QSize(20, 15))
        self.sliderUnit.move(size[0]-20, size[1]-18)

        self.show()
        #::::::::::::::::::::::::::::::::::::::::::::::::::   
        
    def startOrStop(self):
        
        if self.started:
            self.stop()
        else:
            self.start()
            
    def resetPlotWindow(self):
        
        if not self.plotWinClosed:
            self.plotWin.close()
        

        self.plotWin = DCVizPlotWindow(self.activeMode, self)
        self.plotWin.setWindowTitle("DCViz plots")
        self.plotWinClosed = False
        
    
        
    def start(self):
            
        try:
            self.activeMode = self.modeMap[self.modeSelector.currentText()]
        except:
            self.raiseWarning("Select a data set.")
            return
        
        
        if self.dynamic:
            self.startStopButton.setIcon(self.img["stop"])
            self.started = True
        
        self.activeMode.dynamic = self.dynamic
        self.activeMode.plotted = False
        
        self.resetPlotWindow()
        self.job = jobThread(self.activeMode, self.plotWin)
            
        if not self.job.isRunning():
            self.terminalTracker("Job", "Starting.")
            self.job.start()
           
        else:
            print "BUG THREAD: This should never happen..."
    
        
    def stop(self):
        
        self.startStopButton.setIcon(self.img["play"])
        self.started = False
        
        if self.job.isRunning():
        
            self.terminalTracker("Job", "Stoping... ", hold="on")
            self.job.comm.stopSignal.emit()
            self.plotWin.callStop = False
            
            i=0
            while self.job.isRunning():
                time.sleep(0.01)
                i+=1
                if i > 500:
                    print "TIMEOUT: Job didn't exit."
                    sys.exit(1)
                
            self.terminalTracker("Job", "Stopped.")

        else:
            print "BUG THREAD: This should never happen..."


    def openfile(self):
        
        filename, _ = QFileDialog.getOpenFileName(self, \
                                    'Choose file to display', self.masterDir)

        if not filename:
            return
        
        self.detectModetype(filename)
        self.updateModeSelector()
        
    def setpath(self):

        self.path = QFileDialog.getExistingDirectory(self, \
                                    'Choose path...', self.masterDir)
        if not self.path:
            return
        
        for content in [os.path.join(self.path, filename) for filename in os.listdir(self.path)]:
            if os.path.isfile(content):
                self.detectModetype(content)

        self.updateModeSelector()
        
    def detectModetype(self, filename):

        s = 20

        for mode in self.uniqueModes:
            if re.findall(mode.nametag, filename):
                
                self.terminalTracker("Detector", "matched [%s] with [%s]" %  \
                          (os.path.split(filename)[-1].center(s), \
                            self.uniqueModesNames[self.uniqueModes.index(mode)].center(s)))
                            
                args = [filename]
                mode = mode(*args, useGUI=True)
                
                if self.checkConsistency(mode):
                    self.raiseWarning("Similar dataset previously selected: " + str(mode))
                    return
                
                self.modeMap[str(mode)] = mode
                return
                
        self.terminalTracker("Detector", "Unable to fetch data from %s" % os.path.split(filename)[-1])

    def checkConsistency(self, mode):
        return str(mode) in self.modeMap.keys()

    def loadExtern(self, onlyConfig=False):
        
        config = open(os.path.join(thisDir,"config.txt"), 'r')
        self.config = config.read()
        config.close()

        self.refreshedDt = re.findall(r"dynamic refresh interval \[seconds\]\s*=\s*(\d+\.?\d*)", self.config)[0]
        self.refreshDtConfig = float(self.refreshedDt)
        self.terminalSilence = not bool(int(re.findall("terminal tracker\s*=\s*([01])", self.config)[0]))
        self.warningSilence = not bool(int(re.findall("warnings on\s*=\s*([01])", self.config)[0]))

        if not onlyConfig: 
            self.autodetectModes()
            
        #Static overriding
        for mode in self.uniqueModes:
            mode.delay = self.refreshDtConfig
        
        
    def flipDynamic(self):
        self.dynamic = not self.dynamic
        self.dtSlider.setEnabled(self.dynamic)
        self.sliderDisplay.setEnabled(self.dynamic)
        self.sliderUnit.setEnabled(self.dynamic)  


    def updateDt(self, updateConfig):

        newDt = self.dtSlider.value()/10.
        self.sliderDisplay.clear()
        
        if (int(newDt) - newDt) == 0:
            newDt = int(newDt)
            
        self.sliderDisplay.insert(str(newDt))
  
        
        if not updateConfig:
            return
        
        config = open(os.path.join(thisDir,"config.txt"), 'w')
        for line in self.config.split("\n"):
            if line.startswith("dynamic refresh interval [seconds]"):
                config.write(self.config.replace(line, line.replace(self.refreshedDt, str(newDt))))
                config.close()
                self.loadExtern(True)
                return
        

    def autodetectModes(self):
        classfile = open(os.path.join(srcDir, 'DCViz_classes.py'), 'r')
        raw = classfile.read()
        classfile.close()

        self.uniqueModesNames = re.findall('^class (\w+)\(DCVizPlotter\):', raw, re.MULTILINE)
        self.uniqueModes = [eval(subclass) for subclass in self.uniqueModesNames]
        
        self.terminalTracker("Detector", "Found subclasses %s" \
                                % str(self.uniqueModesNames).strip("]").strip("["))
        
        if not self.uniqueModesNames:
            self.raiseWarning("No subclass implementations found.")

        for mode in self.uniqueModes:
            instance = mode()
            try:
                nametag = instance.nametag
            except:
                self.raiseWarning("Subclass %s has no attribute 'nametag' (output filename identifier)." % \
                                  self.uniqueModesNames[self.uniqueModes.index(mode)])
        
    
    def updateModeSelector(self, silent=False):

        if not self.modeMap:
            return
        
        self.modeSelector.clear()
        self.modeSelector.addItems(["Select data..."] + self.modeMap.keys())
        self.modeSelector.model().setData(self.modeSelector.model().index(0,0), Qt.NoItemFlags, Qt.UserRole - 1)
        self.modeSelector.setEnabled(True)
        

    def resetModeSelector(self):
        self.modeMap = {}   
        self.modeSelector.clear()
        self.terminalTracker("GUI", "Data reset.")
        self.modeSelector.addItem("No data")
        self.modeSelector.setEnabled(False)


    def raiseWarning(self, s):

        self.terminalTracker("Warning", s)
        
        if self.warningSilence:
            return
        
        QMessageBox.warning(None, "DCViz GUI", s, QMessageBox.Ok)

    def showInfo(self, h, s):
        self.terminalTracker("Info(%s)" % h, "\n" + s + "\n")
        QMessageBox.information(None, h, s, QMessageBox.Ok)
        
        
    def terminalTracker(self, source, message, hold='off'):
        if self.terminalSilence:
            return

        s = 10

        if hold=='on':
            self.hideSource = True
            print "[%s]:  %s" % (source.center(s), message),
        elif hold=='off':
            if not self.hideSource:
                print "[%s]: " % source.center(s),
            print message
            self.hideSource = False
        else:
            print "Error: hold=%s is not a leagal value." % hold
        
    def loadImages(self):
        self.img = {}
        self.img['play'] = QIcon(os.path.join(thisDir, 'Images', 'play.gif'))
        self.img['stop'] = QIcon(os.path.join(thisDir, 'Images', 'stop.gif'))

def main(masterDir = None):
    
    if masterDir is None:
        masterDir = os.getcwd()
    
    app = QApplication(sys.argv)
    app.setStyle(QStyleFactory.create("Cleanlooks"))
    #app.setWindowIcon(QIcon("Images/QMC.ico"))
                            
    win = DCVizGUI(masterDir)
    win.setWindowTitle('DCViz GUI')
    
    return app.exec_()




if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except:
        path = None
    sys.exit(main(path))
