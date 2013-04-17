
# -*- coding: utf-8 -*-

import sys, re, os, inspect

try:
    import numpy
except:
    print "\n"
    print "You need numpy in order to run this library!"
    print "sudo apt-get install python-numpy-dev"
    print "\n" 
    sys.exit(1)

from matplotlib import rc, pylab, colors, ticker

#~ Paths include
classes_thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(classes_thisDir)

from DCViz_sup import DCVizPlotter, dataGenerator

try:
    rc('text', usetex=True)
    rc('font', family='serif')
    pylab.title("$\LaTeX$")
    pylab.draw()
    pylab.clf()
    
except:
    print "Neccessary latex packages not installed. Disabling latex support."
    rc('text', usetex=False)   



class myTestClass(DCVizPlotter):
    nametag =  'testcase\d\.dat' #filename with regex support
    
    #1 figure with 1 subfigure
    figMap = {'fig1': ['subfig1']}
    
    #skip first row. (the function __str__ is printed here)
    skipRows = 1    
    
    def plot(self, data):
        column1 = data[0]

        self.subfig1.set_title('I have $\LaTeX$ support!')
              
        self.subfig1.set_ylim([-1,1])
          
        self.subfig1.plot(column1)
  

class myTestClassFamily(DCVizPlotter):
    nametag =  'testcaseFamily\d\.dat' #filename with regex support
    
    #1 figure with 3 subfigures
    figMap = {'fig1': ['subfig1', 'subfig2', 'subfig3']}
    
    #skip first row. (the function __str__ is printed here)
    skipRows = 1    

    #Using this flag will read all the files matching the nametag
    #(in the same folder.) and make them aviable in the data arg    
    isFamilyMember = True
    familyName = "testcase"
    
    def plot(self, data):
        
        #figures[0] is 'fig1' figures. the 0'th element is the
        #self.fig1 itself. Subfigures are always index [1:]
        mainFig = self.figures[0][0]  
        mainFig.suptitle('I have $\LaTeX$ support!')        
        subfigs = self.figures[0][1:]
    
        #Notice we plot fileData.data and not fileData alone.
        #The dataGenerator class is used to speed up file reading;
        #looping over family members and directly plotting means
        #we send a dataGenerator instance to matplotlib.
        #in order to get the numpy object, we send the data.
        #Alternatively, we could send data[:]
        for subfig, fileData in zip(subfigs, data):
            subfig.plot(fileData.data)
            subfig.set_ylim([-1,1])
        

class EnergyTrail(DCVizPlotter):
    
    nametag = "\w+.trailingE.arma"
    
    armaBin=True
    isFamilyMember=True
    
    figMap = {"fig1": ['eFig']}
    c = ["#C0C0C0", '#008000', '#008000'] 
    l = ["-", "--", "-"]
    
    def plot(self, data):

        eDMC = 65.7
        
        for i, data_ in enumerate(data):
            name = self.familyFileNames[i].split(".")[0].replace("_", " ")[1:]
            self.eFig.plot(numpy.cumsum(data_.data)/numpy.linspace(1,data_.m,data_.m),
                           self.l[i], c=self.c[i], label=name)
                           
        self.eFig.plot([1,data_.m], [eDMC, eDMC], 'k-.', label="DMC Energy")
        self.eFig.legend()
        
        self.eFig.set_xlabel('Cycle')
        self.eFig.set_ylabel(r'$\langle E\rangle $', rotation=0)
        
        self.eFig.axes.get_yaxis().get_label().set_fontsize(20)
        self.eFig.axes.get_xaxis().get_label().set_fontsize(20)

        formatter = ticker.ScalarFormatter(useOffset=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-3, 2))
        self.eFig.axes.get_xaxis().set_major_formatter(formatter)


class Blocking(DCVizPlotter):
    
    nametag = "blocking_\w+_out\d*\.dat"
    figMap = {"Fig": ["blockFig"]}
    
    nameMap = {"0": r"$\alpha$", "1": r"$\beta$", "": ""}
    
    def plot(self, data):
        
        Fig, blockFig = self.Fig, self.blockFig
        blockSize, error = data
        
#        fileName = os.path.basename(self.filepath)
#        title = "Blocking data %s" % (fileName.split("_")[1] + " %s" % \
#              (self.nameMap[re.findall("blocking_\w+_out(\d*)\.dat", self.filepath)[0]]))
        blockFig.plot(blockSize, error, 
                      '*', color='#008000', markersize=10)  
        
        formatter = ticker.ScalarFormatter(useOffset=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-3, 4))
        blockFig.axes.get_yaxis().set_major_formatter(formatter)
        blockFig.axes.get_yaxis().get_label().set_fontsize(30)
        blockFig.axes.get_yaxis().offsetText.set_size(20)
        blockFig.axes.get_xaxis().get_label().set_fontsize(20)

        
#        blockFig.set_title(title)
        blockFig.set_xlabel(r'Block size')
        blockFig.set_ylabel(r'$\sigma$', rotation=0)
        

class DMC_OUT(DCVizPlotter):
    
    nametag = "DMC_out\.dat"
    figMap = {"Fig": ["N_plot", "E_plot"]}
    dt = 0.001
        
        
    def plot(self, data):
        
        E, Eavg, N, Navg, ET = data
        
        t = numpy.linspace(0, self.dt*(len(E) - 1), len(E))
        N_plot, E_plot = self.N_plot, self.E_plot
        
        lw=2
        
        # E PLOTS
        E_plot.plot(t, E, 'k--', label="dmc E", linewidth=lw)
        E_plot.plot(t, ET, '#008000', label="trial E", linewidth=lw)
        
        E_plot.legend()
        E_plot.set_title('Energy convergeance')
        E_plot.set_xlabel(r'$\tau = n*\delta \tau$ [s]')
        E_plot.set_ylabel(r'E [Ha]')
        E_plot.ticklabel_format(useOffset=False, axis='y')
        
        E_plot.axes.get_yaxis().get_label().set_fontsize(20)
        E_plot.axes.get_xaxis().get_label().set_fontsize(20)
        
        # N PLOTS
        N_plot.plot(t, N, '#008000', linewidth=lw)
        N_plot.plot(t, Navg, 'k--', linewidth=lw)
        
        N_plot.set_ylabel(r"N_W($\tau$)")
        N_plot.axes.get_xaxis().set_visible(False)
        N_plot.set_title('Walker population')
        N_plot.ticklabel_format(useOffset=False, axis='y')
        
        N_plot.axes.get_yaxis().get_label().set_fontsize(20)

class radial_out(DCVizPlotter):
    
    nametag = "radial_out.+\.arma"
    figMap = {"fig1": ["radialFig"]}
    
    armaBin = True
    isFamilyMember=True
    
    def plot(self, data):
        
        cut = data[0].n/20     
        cut= 0
        print data[0].n
        
        max_edge = 0        
        for i in range(len(data)):

            edge = float(re.findall("edge(\d+\.?\d*)\.arma", self.familyFileNames[i])[0])
            if edge > max_edge:
                max_edge = edge;
            
            if "dmc" in self.familyFileNames[i]:
                method = "dmc"
                dmc = data[i].data
            elif "vmc" in self.familyFileNames[i]:
                method = "vmc"
                vmc = data[i].data
            
            r = numpy.linspace(0, edge, data[i].n)
            self.radialFig.plot(r, data[i][0], label=method);
        self.radialFig.plot(r, 2*dmc - vmc, "k-")
        self.radialFig.legend()    
        self.radialFig.axes.set_xlim(r[cut], max_edge)
        self.radialFig.set_xlabel('r')
        self.radialFig.set_ylabel(r'$|\psi(r)|^2$')
        
        
        

class dist_out(DCVizPlotter):
    
    nametag = "dist_out.+\.arma"
    figMap = {"fig1": ["subfigHist2d"]}

    armaBin = True
        
    def plot(self, data):
    
        print self.filepath
        print re.findall("edge(\d+\.?\d*)\.arma", self.filepath)
        edge = float(re.findall("edge(\d+\.?\d*)\.arma", self.filepath)[0])
        
        extent = [-edge, edge, -edge, edge]
        im = self.subfigHist2d.imshow(data.data, extent=extent)
        self.subfigHist2d.set_ylabel(r'y')
        self.subfigHist2d.set_xlabel(r'x')
        self.fig1.colorbar(im)

#        self.subfigDist1d.set_ylabel(r'$|P(r)|^2$')
        
class testBinFile(DCVizPlotter):
    
    nametag = "testBin.+\.arma"
    figMap = {"fig": ["subfig"]}

    fileBin = True
    Ncols = 2
    
    skipRows = 2    

    def plot(self, data):
        x, y = data
        self.subfig.plot(x, y, '*')
        
class MIN_OUT(DCVizPlotter):
    
    nametag = "ASGD_out\.dat"
    figMap = {"E_fig"    : ["E_plot"], 
              "step_fig" : ["step_plot"],
              "param_fig": ["param_plot", "grad_plot"]}
    
    indexmap = {0: r"\alpha", 1: r"\beta"}
    l = ["-", "--"]
    c = ['#008000', 'k']    
    
    def plot(self, data):

        n_params = (self.Ncols - 3)/2

        E, Eavg, step = dataGenerator(data[:3])
        
        E_plot, step_plot, param_plot, grad_plot = self.E_plot, self.step_plot, \
                                                    self.param_plot, self.grad_plot
        
        #~ E PLOTS
        E_plot.plot(E, self.l[0], c=self.c[0], label="E")
        E_plot.plot(Eavg, self.l[1], c=self.c[1], label="average E")
        
        E_plot.legend()
        #E_plot.set_title('Energy convergeance')
        E_plot.set_xlabel(r'Cycle')
        E_plot.set_ylabel(r'E [Ha]')
        E_plot.axes.get_yaxis().get_label().set_fontsize(30)
        E_plot.axes.get_xaxis().get_label().set_fontsize(30)
        E_plot.ticklabel_format(useOffset=False, axis='y')
        
        #~ Step plot
        step_plot.plot(step, self.l[0], c=self.c[0])
        #step_plot.set_title('Step length')
        step_plot.set_xlabel('Cycle')
        step_plot.set_ylabel('Step')
        step_plot.axes.get_yaxis().get_label().set_fontsize(30)
        step_plot.axes.get_xaxis().get_label().set_fontsize(30)

        lw = 2
        #~ Param plots
        for i in range(0, 2*n_params, 2):
            param_plot.plot(data[3 + i], self.l[i/2], c=self.c[i/2],label=r'$%s$' % self.indexmap[i/2], linewidth=lw)    
            grad_plot.plot(data[4 + i], self.l[i/2], c=self.c[i/2], linewidth=lw)
        
        param_plot.set_ylabel(r'$\alpha_i$')
        #param_plot.set_title('Variational parameters')
        param_plot.axes.get_xaxis().set_visible(False)
        
        if n_params > 1:
            param_plot.legend()
        param_plot.axes.get_yaxis().get_label().set_fontsize(30)
        param_plot.axes.get_xaxis().get_label().set_fontsize(30)
        param_plot.axes.get_yaxis().labelpad = 20
        
        #grad_plot.set_title('Energy derivatives')
        grad_plot.set_ylabel(r'$\frac{\partial E}{\partial \alpha_i}$', rotation=0)
        grad_plot.set_xlabel('Cycle')
        grad_plot.axes.get_yaxis().get_label().set_fontsize(30)
        grad_plot.axes.get_xaxis().get_label().set_fontsize(30)
            
            
        
def testbedJorgen(dynamic):
    
    from pyLibQMC import paths
    path = os.path.join(paths.scratchPath, "QMC_SCRATCH")

    setting = "MIN"

    if setting == "MIN" or setting == "ALL":
        filepath = os.path.join(path, "ASGD_out.dat")
        plot_tool = MIN_OUT(filepath, dynamic)
    
        if setting == "ALL":
            print "MIN OUTPUT STARTING"
            plot_tool.mainloop()

    if setting == "DMC" or setting == "ALL":
        dt = 0.001
        filepath = os.path.join(path, "DMC_out.dat")
        plot_tool = DMC_OUT(filepath, dynamic)
        plot_tool.dt = dt
    
        if setting == "ALL":
            print "DMC OUTPUT STARTING"
            plot_tool.mainloop()

    if setting != "ALL":
        plot_tool.mainloop()
        
if __name__ == "__main__":
    print "Dynamic = FALSE"
    testbedJorgen(False)
    print "Dynamic = TRUE"
    testbedJorgen(True)
