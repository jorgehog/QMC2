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

fileName = "DMC_out.dat"
path = paths.scratchPath+ "/QMC_SCRATCH/"

class DMC_OUT(plot_tools):
	def __init__(self, filename, path, dt = 0.001, dynamic=False):
		plot_tools.__init__(self, filename, path, dynamic)
		
		self.fig = plab.figure()
		self.iFig = self.add_figure(self.fig)
		
		self.E_plot = self.fig.add_subplot(2,1,2)
		self.N_plot = self.fig.add_subplot(2,1,1)
		
		self.add_subfigure(self.E_plot, self.iFig)
		self.add_subfigure(self.N_plot, self.iFig)
		
		self.dt = dt
	
	def plot(self, data):
		
		E, Eavg, N, Navg, ET = data
		
		#~ data manip
		#~ ET = numpy.cumsum(ET)/numpy.cumsum(numpy.ones(ET.shape))
		#~ E *= numpy.random.uniform(0.9999,1.0001, E.shape)
		#~ N *= numpy.random.uniform(0.99,1.01, N.shape)
		#~ 
		
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
		
		
		
dt = 0.001
dynamic = True

plot_tool = DMC_OUT(fileName, path, dt, dynamic)
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
