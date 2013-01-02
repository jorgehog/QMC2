# -*- coding: utf-8 -*-

import re, numpy, time, sys
import mayavi.mlab as mlab
import matplotlib.pylab as plab

class dcv_plotter:
	
	figMap = {}
	
	def __init__(self, filepath, dynamic=False):
		self.dynamic = dynamic
		self.filepath = filepath 
		self.file = open(self.filepath , "r")
		
		self.set_delay(5)
		
		self.figures = []
	
		
	def set_delay(self, t):
		if int(t) < 1:
			print "invalid delay time. Delay set to 1s."
			t = 1
		if isinstance(t, float):
			print "delay must be in seconds. ", t, " changed to ", int(t)
		self.delay = int(t)
		
		
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
		for fig in self.figMap.keys():
			s += "self.%s = plab.figure(%d); " % (fig, self.i)
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
		
		while (True):
	
			data = self.get_data()
	
			self.plot(data)
			try:
				self.show()
			except:
				print "Windows closed prematurely. Ending..."
				sys.exit(1)
		
			if self.dynamic:
				
					#~ print "replotting in %d s..." % self.delay
					#~ for i in range(int(self.delay)):
						#~ time.sleep(1)

				ans = raw_input("press return to plot again. (q to break)")
						
				self.clear()
					
				if ans == "q":
					print "Ending session..."
					self.close()
					break
			else:
				raw_input("Press any key to exit")
				
				self.close()
			
				break
			
			
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
		if not self.dynamic:
			self.file.seek(0)
			return None
		if not self.file.closed: 
			self.file.close()
		self.file = open(self.filepath , "r")
			
	def add_figure(self, fig):
		self.figures.append([fig])
		return len(self.figures) - 1
		
	def add_subfigure(self, subfig, i):
		self.figures[i].append(subfig)
		return len(self.figures[i]) - 1
		
	def show(self):
		for fig in self.figures:
			fig[0].canvas.draw()
			fig[0].show()
		
	def clear(self):
		for fig in self.figures:
			for subfig in fig[1:]:
				subfig.clear()
		
	def close(self):
		
		self.clear()
		
		for fig in self.figures:
			plab.close(fig[0])
		
		if not self.file.closed:
			self.file.close()
		
	def plot(self, data):
		"I am virtual"
	
	
