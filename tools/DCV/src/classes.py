# -*- coding: utf-8 -*-

try:
    import numpy
except:
    print "\n"
    print "You need numpy in order to run this library!"
    print "sudo apt-get install python-numpy-dev"
    print "\n" 
    sys.exit(1)

import time, sys, re, os, inspect


from random import choice

import matplotlib.pylab as plab
from matplotlib import rc

#~ Paths include
classes_thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(classes_thisDir)

from dcv_sup import DCVizPlotter

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
        
        param_plot.set_ylabel(r'\alpha_i')
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
