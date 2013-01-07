# -*- coding: utf-8 -*-

import re, numpy, time, sys, signal, os
import mayavi.mlab as mlab
import matplotlib.pylab as plab


class DCVizPlotter:
    
    figMap = {}
    delay = 3
    parent = None
    
    def __init__(self, filepath=None, dynamic=False, useGUI=False):
        self.dynamic = dynamic
        self.useGUI = useGUI
        
        self.plotted = False
        self.stopped = False
        self.SIGINT_CAPTURED = False
            
        self.filepath = filepath
        self.file = None
        
        signal.signal(signal.SIGINT, self.signal_handler)
    
    
    def signal_handler(self, signal, frame):
        print "Ending session..."
        self.SIGINT_CAPTURED = True
        
        
    def set_delay(self, t):
        if int(t) < 1:
            print "invalid delay time. Delay set to 1s."
            t = 1
        
        self.delay = t
    
    def __str__(self):
        return ".".join(os.path.split(self.filepath)[-1].split(".")[0:-1])
        
    def get_data(self):
    
        self.reload()
    
        data = numpy.array(self.rx.findall(self.file.read()), numpy.float)

        self.file.close()
        
        output = []
        for i in range(data.shape[1]):
            output.append(data[:,i])
        
        return tuple(output)
        
    def set_figures(self):
        
        s = ""
        i = 0
        self.figures = []
        for fig in self.figMap.keys():
            s += "self.%s = plab.figure(%d); " % (fig, i)
            s += "self.i%s = self.add_figure(self.%s); " % (fig, fig)
        
            subFigs = self.figMap[fig]
            nFigs = len(subFigs)
            
            for j in range(nFigs):
                s += "self.%s = self.%s.add_subplot(%d, 1, %d); " % (subFigs[j], fig, nFigs, j+1)
                s += "self.add_subfigure(self.%s, self.i%s); " % (subFigs[j], fig)
            
            i += 1
        exec(s)
        #print s.replace("; ", "\n")                

    def mainloop(self):
        
        if self.dynamic:
            light = 'red'
            while light == 'red':
                light = self.load_sample()
                time.sleep(1)
        else:
            light = self.load_sample()
            if light == 'red':
                sys.exit(1)


        self.set_figures()
        
        while (not self.stopped and not self.SIGINT_CAPTURED):

            if self.plotted:
                self.clear()
                print "Replotting..."
           
            data = self.get_data()
    
            self.plot(data)
            
            if not self.useGUI:
                self.show()
            else:
                if self.plotted and self.dynamic:
                    try:
                        self.show(drawOnly=True)
                    except:
                        #Exception needed in order for the thread to survive being closed by GUI
                        pass
                if not self.plotted:
                    self.parent.comm.plotSignal.emit()
                 
                    
            self.plotted = True
                
            if self.dynamic:
                for i in range(int(self.delay)):
                    time.sleep(1)
                    if self.stopped or self.SIGINT_CAPTURED:
                        break

                time.sleep(self.delay - int(self.delay))
                
            else:
                if not self.useGUI:
                    raw_input("Press any key to exit")
                break
                
        
        if not self.useGUI:
            self.close()
            
            
    def load_sample(self):
        
        self.reload()
        sample = self.file.readline()
        if not sample:
            print "No data in file"
            return "red"
        
        self.N = len(sample.split())
        anyNumber = r'[\+\-]?\d+\.?\d*[eE]?[\+\-]?\d*'
        self.rx = re.compile((r'(%s)\s+' % anyNumber)*(self.N-1) + r'(%s)[\n$]' % anyNumber)    
        
        return "green"
    
    def reload(self):
    
        if self.file:
            if not self.file.closed: 
                self.file.close()
            
        self.file = open(self.filepath , "r")
            
    def add_figure(self, fig):
        self.figures.append([fig])
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
                subfig.axes.legend_ = None
        
    def close(self):
        self.clear()
        
        for fig in self.figures:
            plab.close(fig[0])
        
        if not self.file.closed:
            self.file.close()
        
    def plot(self, data):
        "I am virtual"
    
    
