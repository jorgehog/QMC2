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

from random import choice

from matplotlib import rc

#~ Paths include
classes_thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(classes_thisDir)

from DCViz_sup import DCVizPlotter

#~ Enabling latex support in matplotlib
rc('text', usetex=True)
rc('font', family='serif')


class myTestClass(DCVizPlotter):
	nametag =  'testcase\d\.dat' #filename with regex support
	
	#1 figure with 2 subfigures and another figure with only one.
	figMap = {'fig1': ['subfig1', 'subfig2'], 'fig2': ['subfig3']}
	
	def plot(self, data):
		
          column1, column2, column3 = data
          fig1, subfig1, subfig2, fig2, subfig3 = \
				self.fig1, self.subfig1, self.subfig2, self.fig2, self.subfig3
		
          c = ['r', 'g', 'b']
          subfig1.plot(column1, choice(c))
          subfig2.plot(column2, choice(c))
          subfig3.plot(column3, choice(c))
          subfig3.set_title('I have $\LaTeX$ support!')
  
          subfig1.set_xlim([0,12000])
          subfig1.set_ylim([-1,1])
          subfig2.set_xlim([0,12000])
          subfig2.set_ylim([-1,1])
          subfig3.set_xlim([0,12000])
 
		

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
    
    nametag = "dist.+\.arma"
    figMap = {"fig": ["subfig"]}

    armaBin = True
    
    isFamilyMember = True
    familyName = "Dist"

#    def plot(self, data):
#        
#        n_p = len(data[0][0])
#        n = len(data)
#        print "n_p=", n_p
#        print "n=", n
#        X = numpy.zeros(n_p)
#        Y = numpy.zeros(n_p)
#        for rMat in data:
#            x, y = rMat;
#            self.subfig.plot(x,y, 'r*')
#            X += x;
#            Y += y;
#        
#        X /= n
#        Y /= n
#
#        self.subfig.plot(X, Y, '*')     
    def plot(self, data):
        
        
        n_p = len(data[0][0])        
        n = len(data)
        scale = n/100


        r = numpy.zeros(n)
        
        for i in range(n):
            x, y = data[i];
            r[i] = numpy.sqrt(x**2 + y**2).sum()
        r /= n_p
        Y, X = numpy.histogram(r, bins=numpy.linspace(0, r.max(), n/scale))

        self.subfig.plot(X[:-1], Y, '*')
       
    
        
class MIN_OUT(DCVizPlotter):
    
    nametag = "ASGD_out\.dat"
    figMap = {"E_fig"    : ["E_plot"], 
              "step_fig" : ["step_plot"],
              "param_fig": ["param_plot", "grad_plot"]}
    
    indexmap = {0: r"\alpha", 1: r"\beta"}
    c = ['b', 'g']    
    
    def plot(self, data):

        n_params = (self.N - 3)/2
        
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

    setting = "ALL"

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
