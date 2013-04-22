# -*- coding: utf-8 -*-

import re, numpy, time, signal, os
import matplotlib.pylab as plab
from os.path import join as pjoin


class dataGenerator:
    def __init__(self, data):
        self.data = data
        
        if len(data.shape) == 2:
            self.n, self.m = data.shape
            self.getD = self.get2Ddata
        else:
            self.m = data.shape
            self.n = 1
            self.getD = self.get1Ddata
        
        self.shape = data.shape
        self.fullshape = (self.n, self.m)
        self.size = data.size

    def get1Ddata(self, i):
        if i == slice(0, 9223372036854775807, None):
            return self.data
        elif i == 0:
            return self.data
        elif i == slice(0, 1, None):
            return self.data
        else:
            raise IndexError("Index out of bounds.")

    def get2Ddata(self, i):
        return self.data[:, i]

    def __iter__(self):
        for i in range(self.m):
            yield self.getD(i)

    def __len__(self):
        return self.n

    def __getitem__(self, i):
        return self.getD(i)

    def __str__(self):
        return str(self.data)

    def __add__(self, other):
            
        if isinstance(other, dataGenerator):
            return self.data + other.data
        else:        
            raise NotImplementedError("add/sub only works for two datasets. Use 'sum([sum(d) for d in data])'")

    def __sub__(self, other):    
            return self.data + (-other.data)


class DCVizPlotter:
    
    figMap = {}

    armaBin = False
    fileBin = False
    Ncols = None
    
    isFamilyMember = False
    familyName = "unnamed"
    familyFileNames = []

    skippedRows = []
    skipCols = 0
    skipRows = None
    
    nFig = 0

    delay = 3

    parent = None
    
    stack = "V"
    
    canStart = False
    
    def __init__(self, filepath=None, dynamic=False, useGUI=False, toFile=False):
        self.dynamic = dynamic
        self.useGUI = useGUI
        
        self.plotted = False
        self.stopped = False
        self.SIGINT_CAPTURED = False
            
        self.figures = []
    
        self.filepath = filepath
        self.file = None
        
        self.toFile = toFile
        
        signal.signal(signal.SIGINT, self.signal_handler)
        
        if self.Ncols is None and self.fileBin:
            self.Error("You need to specify the number of cols 'Ncols' in order to read a binary file.")
            
    
    
    def signal_handler(self, signal, frame):
        print "[%s] Ending session..." % "DCViz".center(10)
        self.SIGINT_CAPTURED = True
    
    def __str__(self):
        
        if self.isFamilyMember:
            return self.familyName + " (family)"
        return ".".join(os.path.split(self.filepath)[-1].split(".")[0:-1])

        
    def get_data(self, setUpFamily):

        if setUpFamily:

            familyHome, _ = os.path.split(self.filepath)
            familyNames = [name for name in os.listdir(familyHome)\
                            if re.findall(self.nametag, name) and "tmp" not in name]
           
            familyMembers = sorted([pjoin(familyHome, name) for name in familyNames])
            N = len(familyMembers)
            
            self.familySkippedRows = [0]*N
            data = [0]*N
            self.familyFileNames = [0]*len(data)     
            
            for i in range(N):
             
                self.file = None
                self.filepath = familyMembers[i]
            
                self.familyFileNames[i] = os.path.basename(familyMembers[i])
                
                data[i] = self.get_data(setUpFamily=False)
                self.familySkippedRows.append(self.skippedRows)
                
            self.file.close()
            return data
            
        data = []
        #attempt to reload file untill data is found.
        t0 = time.time()
        while len(data) == 0:
            self.reload()
            
            if self.armaBin:
                data = self.unpackArmaMatBin(self.file)      
            elif self.fileBin:
                data = self.unpackBinFile(self.file)
            else:
                data = numpy.array(self.rx.findall(self.file.read()), numpy.float)
                
            if time.time() - t0 > 10.0:
               self.Error("TIMEOUT: File was empty for too long...")
               self.file.close()
               return
               
        self.file.close()
        
        return dataGenerator(data)
    
    def unpackBinFile(self, binFile):
 
        data = numpy.fromfile(binFile, dtype=numpy.float64)
  
        m = self.Ncols
        n = int(data.size/m)
        
        data.resize(n, m)
        
        return data
            
    
    def unpackArmaMatBin(self, armaFile):   

        armaFormat =  armaFile.readline()
        
        n, m = armaFile.readline().strip().split()
        n = int(n)       
        m = int(m)
        
        data = numpy.fromfile(armaFile, dtype=numpy.float64)
        data.resize((n,m))
        
        return data
    
    def set_figures(self):
        
        if self.stack not in ["H", "V"]:
            self.Error("Invalid stack argument %s. (Choose either H or V)" % self.stack)
        
        s = ""
        i = 0
        self.figures = []
        for fig in self.figMap.keys():
            s += "self.%s = plab.figure(); " % fig
            s += "self.i%s = self.add_figure(self.%s); " % (fig, fig)
        
            subFigs = self.figMap[fig]
            nFigs = len(subFigs)
            
            for j in range(nFigs):
                
                if self.stack == "V":
                    s += "self.%s = self.%s.add_subplot(%d, 1, %d); " % (subFigs[j], fig, nFigs, j+1)
                else:
                    s += "self.%s = self.%s.add_subplot(1, %d, %d); " % (subFigs[j], fig, nFigs, j+1)
                s += "self.add_subfigure(self.%s, self.i%s); " % (subFigs[j], fig)
            
            i += 1
        exec(s)
      
      
    def manageFigures(self):
        if not self.useGUI:
            self.set_figures()
        else:
            
            self.parent.comm.setFigureSignal.emit()
            
            i = 0
            while not self.canStart:
                time.sleep(0.01)
                i+=1
                if i > 500:
                    self.Error("TIMEOUT: Figures wasn't set...")
                    return          


    def Error(self, s):
        if self.useGUI:
            self.parent.parent.terminalTracker("DCViz", s)
        else:
            print "[%s] %s" % ("DCViz".center(10), s)


    def waitForGreenLight(self):

        if self.dynamic:
            light = 'red'
            while light == 'red':
                
                if self.stopped or self.SIGINT_CAPTURED:
                    return True
                    
                light = self.load_sample()
                time.sleep(1)
        else:
            light = self.load_sample()
            if light == 'red':
                return True
        
        return False
                
    def shouldReplot(self):
        return not self.stopped and not self.SIGINT_CAPTURED
    def shouldBreak(self):
        return self.stopped or self.SIGINT_CAPTURED

    def showFigures(self):
        
        if not self.useGUI:
            if not self.plotted:
                self.show(drawOnly = self.toFile)
            else:
                self.show(drawOnly=True)
        else:
            try:
                self.show(drawOnly=True)
            except:
                #Exception needed in order for 
                #the thread to survive being closed by GUI
                pass
            
            if not self.plotted:
                self.parent.comm.plotSignal.emit()

    def sleep(self):
        for i in range(int(self.delay)):
            time.sleep(1)
            if self.shouldBreak():
                break

        time.sleep(self.delay - int(self.delay))

    def mainloop(self):
        
        if not self.armaBin:
            breakMe = self.waitForGreenLight()
            
            if breakMe:
                return

        self.manageFigures()

        while (self.shouldReplot()):
            self.clear()
           
            data = self.get_data(setUpFamily = self.isFamilyMember)
  
            self.plot(data)  
            self.showFigures()
            self.plotted = True

                
            if self.dynamic:
                self.sleep()     
            else:
                if not self.useGUI and not self.toFile:
                    raw_input("[%s] Press any key to exit" % "DCViz".center(10))
                break
                
        if not self.useGUI:
            if self.toFile:
                self.saveFigs()
            self.close()
        
            
            
    def load_sample(self):
        
        self.reload()
  
        sample = self.file.read()
        
        if not sample:
            self.Error("No data in file...")
            return "red"
            
        if self.fileBin:
            return "green"
            
        skipRows, self.Ncols = self.sniffer(sample)

        self.file.seek(0)
        
        self.skippedRows = []
        for i in range(skipRows):
            self.skippedRows.append(self.file.readline().strip())
       
        anyNumber = r'[\+\-]?\d+\.?\d*[eE]?[\+\-]?\d*'
        self.rx = re.compile(r'\s*.+?'*self.skipCols + (r'\s*(%s)\s+' % anyNumber)*(self.Ncols - self.skipCols-1) + r'(%s)\s*[\n$]' % anyNumber)    
        
        return "green"
        
    
    def sniffer(self, sample):
        
        sampleList = [row.split() for row in sample.split("\n")]

        #If user has specified the number of rows to skip
        if self.skipRows is not None:
            return self.skipRows, len(sampleList[self.skipRows])

        nLast = len(sampleList[-1])
        i = len(sampleList) - 1

        if nLast != len(sampleList[-2]):
            nLast = len(sampleList[-2])
            i -= 1
        
        while len(sampleList[i]) == nLast and i >= 0:
            i -= 1

        #(nRows to skip = i+1, nCols = nLast)
        return i+1, nLast    
    
    def reload(self):
    
        if self.file:
            if not self.file.closed: 
                self.file.close()
        
        self.file = open(self.filepath , "r")
   
        self.skippedRows = []
        if not self.armaBin and self.skipRows is not None:        
            for i in range(self.skipRows):
                self.skippedRows.append(self.file.readline().strip())
              
    def saveFigs(self):
        path, fname = os.path.split(self.filepath)
        
        dirname = "DCViz_out"
        dirpath = pjoin(path, dirname)
        
        if not os.path.exists(dirpath):
            os.mkdir(dirpath)
        
        i = 0
        for fig in zip(*self.figures)[0]:
            
            figname = ".".join(fname.split(".")[0:-1]) + "_" + str(i) + ".png"
            figpath = pjoin(dirpath, figname)
            fig.savefig(figpath)

            i += 1
            
        
            
    def add_figure(self, fig):
        self.figures.append([fig])
        self.nFig += 1
        return len(self.figures) - 1
        
    def add_subfigure(self, subfig, i):
        self.figures[i].append(subfig)
        return len(self.figures[i]) - 1
        
    def show(self, drawOnly=False):
        for fig in self.figures:
            fig[0].canvas.draw()
            if not drawOnly:
                fig[0].show()
        
    def clear(self):
        for fig in self.figures:
            for subfig in fig[1:]:
                subfig.clear()
                subfig.axes.clear()
                subfig.axes.legend_ = None
        
    def close(self):
        self.clear()
        
        for fig in self.figures:
            plab.close(fig[0])
        
        if self.file is not None:
            if not self.file.closed:
                self.file.close()
        
    def plot(self, data):
        "I am virtual"
    
    
