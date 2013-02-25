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

from matplotlib import rc, pylab

#~ Paths include
classes_thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(classes_thisDir)

from DCViz_sup import DCVizPlotter


useTex = True
if useTex:    
    #~ Enabling latex support in matplotlib
    rc('text', usetex=True)
    rc('font', family='serif')


class myTestClass(DCVizPlotter):
    nametag =  'testcase\d\.dat' #filename with regex support
    
    #1 figure with 1 subfigure
    figMap = {'fig1': ['subfig1']}
    
    #skip first row. (the function __str__ is printed here)
    skipRows = 1    
    
    def plot(self, data):
   
        column1 = data
        

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
    
        for subfig, fileData in zip(subfigs, data):
            subfig.plot(fileData)
            subfig.set_ylim([-1,1])
        

class Blocking(DCVizPlotter):
    
    nametag = "blocking_\w+_out\d*\.dat"
    figMap = {"Fig": ["blockFig"]}
    
    nameMap = {"0": r"$\alpha$", "1": r"$\beta$", "": ""}
    
    def plot(self, data):
        
        Fig, blockFig = self.Fig, self.blockFig
        blockSize, error = data
        
        fileName = os.path.basename(self.filepath)
        title = "Blocking data %s" % (fileName.split("_")[1] + " %s" % \
              (self.nameMap[re.findall("blocking_\w+_out(\d*)\.dat", self.filepath)[0]]))
        
        blockFig.plot(blockSize, error, 'b+')  
        
        blockFig.set_title(title)
        blockFig.set_xlabel(r'Block size')
        blockFig.set_ylabel(r'$\sigma$')
        

class DMC_OUT(DCVizPlotter):
    
    nametag = "DMC_out\.dat"
    figMap = {"Fig": ["N_plot", "E_plot"]}
    dt = 0.001
        
        
    def plot(self, data):
        
        E, Eavg, N, Navg, ET = data
        
        t = numpy.linspace(0, self.dt*(len(E) - 1), len(E))
        N_plot, E_plot = self.N_plot, self.E_plot
        
        
        # E PLOTS
        E_plot.plot(t, E, 'b', label="dmc E")
        E_plot.plot(t, ET, 'r', label="trial E")
        
        E_plot.legend(loc=4)
        E_plot.set_title('Energy convergeance')
        E_plot.set_xlabel(r'$\tau = n*\delta \tau$ [s]')
        E_plot.set_ylabel(r'E [Ha]')
        E_plot.ticklabel_format(useOffset=False, axis='y')
        
        # N PLOTS
        N_plot.plot(t, N, 'b')
        N_plot.plot(t, Navg, 'r')
        
        N_plot.set_ylabel(r"N_W($\tau$)")
        N_plot.axes.get_xaxis().set_visible(False)
        N_plot.set_title('Walker population')
        N_plot.ticklabel_format(useOffset=False, axis='y')

class dist_out(DCVizPlotter):
    
    nametag = "dist_out.+\.arma"
    figMap = {"fig1" : ["subfigRadii"], "fig2": ["subfigHist2d"], "fig3": ["subfigDist1d"]}

    armaBin = True

    isFamilyMember = True
    familyName = "Dist"
  
  
    def plot(self, data):
        
        n_p = len(data[0][0])  
        dim = len(data[0])
        n = len(data)

        cores = re.findall("dist_out(\d+)_\d+\.arma", " ".join(self.familyFileNames[::1000]))     
        nCores = max([int(core) for core in cores]) + 1     
        N = n/nCores
    
        I = [0 for i in range(nCores)]
        
        
        xyz = numpy.empty([nCores, N, n_p, dim])
        rMean = numpy.zeros([nCores, n_p])
        R = numpy.zeros([nCores, n_p, N])
        
        for i in xrange(N*nCores):
            
            #find which core this walker belongs to
            c = int(re.findall("dist_out(\d+)_\d+\.arma", 
                               self.familyFileNames[i])[0])
            
            for j in xrange(dim):
                xyz[c, I[c], :, j] = data[i][j]
                
            I[c]+=1
            
        for c in xrange(nCores):
            for i in xrange(n_p):
                R[c, i, :] = numpy.sqrt((xyz[c]**2)[:, i, :].sum(1))
                rMean[c, i] = R[c, i, :].sum()
      

        indexes = rMean.argsort(1)
        
        rMean.sort(1)
        rMean = rMean.sum(0)/(N*nCores)
        
        
        for i in range(nCores):
            R[i, :, :] = R[i, indexes[i], :]
        
        RSorted = numpy.zeros([n_p, N*nCores])
        R.resize([n_p*nCores, N])
        
        for i in range(n_p):
            RSorted[i, :] = R[i::n_p].reshape(N*nCores)
            
        R = RSorted
        
        nBins = 100
        
        xyz.resize(N*nCores*n_p, dim)
        
        self.subfigHist2d.hexbin(xyz[:,0], xyz[:,1])
#        self.subfigHist2d.set_xlim([-0.5, 0.5])
#        self.subfigHist2d.set_ylim([-0.5, 0.5])
                
        self.subfigHist2d.set_xlabel(r'x')
        self.subfigHist2d.set_ylabel(r'y')
        
        for i in range(n_p):
            self.subfigDist1d.hist(R[i, :], 
                                   nBins, 
                                   histtype='step')
                                   
        self.subfigDist1d.set_xlabel(r'$r = \sqrt{x^2 + y^2' + ' +z^2'*(dim==3) + '}$')

        self.subfigRadii.scatter(range(1, n_p+1), rMean, s=10, marker='o', c='b')
        self.subfigRadii.set_ylim(0, rMean[n_p/2]*2)
        self.subfigRadii.axes.get_xaxis().set_major_locator(pylab.MaxNLocator(integer=True))
        self.subfigRadii.set_xlabel(r'Particle')
        self.subfigRadii.set_ylabel(r'$\langle r \rangle$')
        
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
    c = ['b', 'g']    
    
    def plot(self, data):

        n_params = (self.Ncols - 3)/2
        
        E, Eavg, step = data[:3]
        
        E_plot, step_plot, param_plot, grad_plot = self.E_plot, self.step_plot, \
                                                    self.param_plot, self.grad_plot
        
        #~ E PLOTS
        E_plot.plot(E, 'b', label="E")
        E_plot.plot(Eavg, 'r', label="average E")
        
        E_plot.legend()
        E_plot.set_title('Energy convergeance')
        E_plot.set_xlabel(r'cycle')
        E_plot.set_ylabel(r'E [Ha]')
        E_plot.ticklabel_format(useOffset=False, axis='y')
        
        #~ Step plot
        step_plot.plot(step, 'b')
        step_plot.set_title('step length')
        step_plot.set_xlabel('cycle')
        step_plot.set_ylabel('step')
        E_plot.ticklabel_format(useOffset=False, axis='y')
        
        #~ Param plots
        for i in range(0, 2*n_params, 2):
            param_plot.plot(data[3 + i], self.c[i/2],label=r'$%s$' % self.indexmap[i/2])    
            grad_plot.plot(data[4 + i], self.c[i/2])
        
        param_plot.set_ylabel(r'$\alpha_i$')
        param_plot.set_title('Variational parameters')
        param_plot.axes.get_xaxis().set_visible(False)
        param_plot.legend()
        
        grad_plot.set_title('Energy derivatives')
        grad_plot.set_ylabel(r'$\frac{\partial E}{\partial \alpha_i}$ [Ha]')
        grad_plot.set_xlabel('cycle')
            
            
        
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
