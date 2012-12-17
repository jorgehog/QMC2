import numpy
import time, sys, re

import mayavi.mlab as mlab
import matplotlib.pylab as plab
from matplotlib import rc

from pyLibQMC import paths, plot_tools

rc('text', usetex=True)
rc('font', family='serif')

#~ 
#~  
#~ 


path = paths.scratchPath+ "/QMC_SCRATCH/"

class DMC_OUT(plot_tools):
	def __init__(self, filename, path, dt = 0.001, dynamic=False):
		plot_tools.__init__(self, filename, path, dynamic)

		self.dt = dt
	
	def set_figures(self):
		
		self.fig = plab.figure()
		self.iFig = self.add_figure(self.fig)
		
		self.E_plot = self.fig.add_subplot(2,1,2)
		self.N_plot = self.fig.add_subplot(2,1,1)
		
		self.add_subfigure(self.E_plot, self.iFig)
		self.add_subfigure(self.N_plot, self.iFig)
		
		
	def plot(self, data):
		
		E, Eavg, N, Navg, ET = data
		
		t = numpy.linspace(0, dt*(len(E) - 1), len(E))
		fig, N_plot, E_plot = self.fig, self.N_plot, self.E_plot
		
		
		#~ E PLOTS
		E_plot.plot(t, E, 'b', label="dmc E")
		E_plot.plot(t, ET, 'r', label="trial E")
		
		E_plot.legend(loc=4)
		E_plot.set_title('Energy convergeance')
		E_plot.set_xlabel(r'$\tau = n*\delta \tau$ [s]')
		E_plot.set_ylabel(r'E [Ha]')
		E_plot.ticklabel_format(useOffset=False, axis='y')
		
		#~ N PLOTS
		N_plot.plot(t, N, 'b')
		N_plot.plot(t, Navg, 'r')
		
		N_plot.set_ylabel(r"N_W($\tau$)")
		N_plot.axes.get_xaxis().set_visible(False)
		N_plot.set_title('Walker population')
		N_plot.ticklabel_format(useOffset=False, axis='y')
		
		
class MIN_OUT(plot_tools):
	def __init__(self, filename, path, dynamic):
		plot_tools.__init__(self, filename, path, dynamic)
		
		self.indexmap = {0: r"\alpha", 1: r"\beta"}
		self.c = ['b', 'g']
		
	def set_figures(self):
		
		#Adding the energy plot
		self.E_fig = plab.figure(1)
		self.iE_fig = self.add_figure(self.E_fig)
		self.E_plot = self.E_fig.add_subplot(1,1,1)
		self.add_subfigure(self.E_plot, self.iE_fig)
		
		#Adding the step plot
		self.step_fig = plab.figure(2)
		self.istep_fig = self.add_figure(self.step_fig)
		self.step_plot = self.step_fig.add_subplot(1,1,1)
		self.add_subfigure(self.step_plot, self.istep_fig)
		
		#Adding the parameter plot
		self.param_fig = plab.figure(3)
		self.iparam_fig = self.add_figure(self.param_fig)
		self.param_plot = self.param_fig.add_subplot(2,1,1)
		self.add_subfigure(self.param_plot, self.iparam_fig)
		
		#Addinf the gradient subplot
		self.grad_plot = self.param_fig.add_subplot(2,1,2)
		self.add_subfigure(self.grad_plot, self.iparam_fig)
		
		
		self.n_params = (self.N - 3)/2
	
		
		
	def plot(self, data):
		
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
		for i in range(0, 2*self.n_params, 2):
			param_plot.plot(data[3 + i], self.c[i/2],label=r'$%s$' % self.indexmap[i/2])	
			grad_plot.plot(data[4 + i], self.c[i/2])
		
		param_plot.set_ylabel(r'\alpha_i')
		param_plot.set_title('Variational parameters')
		param_plot.axes.get_xaxis().set_visible(False)
		param_plot.legend()
		
		grad_plot.set_title('Energy derivatives')
		grad_plot.set_ylabel(r'$\frac{\partial E}{\partial \alpha_i}$ [Ha]')
		grad_plot.set_xlabel('cycle')
			
			
		

dynamic = True
#~ setting = "MIN"
setting = "ALL"

if setting == "DMC" or setting == "ALL":
	dt = 0.001
	fileName = "DMC_out.dat"
	plot_tool = DMC_OUT(fileName, path, dt, dynamic)
	
	if setting == "ALL":
		plot_tool.mainloop()
	
if setting == "MIN" or setting == "ALL":
	fileName = "ASGD_out.dat"
	plot_tool = MIN_OUT(fileName, path, dynamic)
	
	if setting == "ALL":
		plot_tool.mainloop()

if setting != "ALL":
	plot_tool.mainloop()



"""
t = numpy.linspace(0, 2*numpy.pi, 1000)
x = numpy.cos(t)
y = numpy.sin(t)
z = numpy.linspace(0, 10, 1000)
	
mlab.plot3d(x,y,z)
	
		
mlab.show()
"""

"""
rerunMode = False

cmd = ""

fig = plab.figure()	
E_plot = fig.add_subplot(1,1,1)
#~ N_plot = fig.add_subplot(2,2,1)




while (cmd != "q"):
	
	E, Eavg, N, Navg, ET = plot_tool.get_data()
	
	E *= numpy.random.uniform(0.99999,1.00001, E.shape)

	if len(E) == len(Eavg) == len(N) == len(Navg) == len(ET):

		dt_array = numpy.linspace(0,dt*len(E), len(E) +1)
		print dt_array[:5]
		
		E_plot.plot(E)
		#~ E_plot.plot(
		E_plot.ticklabel_format(useOffset=False, axis='y')
		fig.canvas.draw()
		fig.show()
		
		
		
		if rerunMode:
			try:
				print "replotting in 5 sec.."
				time.sleep(5)
				E_plot.clear()
			except KeyboardInterrupt:
				fig.close()
				break
		else:
			raw_input("Press any key to exit")
			cmd = "q"
	else:
		print "Unequal data set loaded. Loaded during filewrite? Reloading in 1 sec..."
		try:
			time.sleep(1)
		except KeyboardInterrupt:
			print "Exiting.."
			break

	
	
print "closing"
plot_tool.close()
"""
