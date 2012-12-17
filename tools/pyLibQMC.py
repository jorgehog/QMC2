# -*- coding: utf-8 -*-

import os, time, sys, re, numpy

class paths:
    HOME = os.path.expanduser('~')
    CODE = HOME +  "/MASTER/QMC2"
    IDEPath = HOME + "/NetBeansProjects/nbQMC2"
    toolsPath = CODE + "/tools"
    scratchPath = HOME + "/scratch"
    iniFilePath = CODE + "/iniFiles"
    programPath = HOME + "/NetBeansProjects/nbQMC2/dist/Debug/GNU-Linux-x86"
    

class misc:
    QMC2programName = "nbqmc2"
    

class plot_tools:
	
	def __init__(self, filename, path, dynamic=False):
		self.dynamic = dynamic
		self.filepath = path + filename 
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
				
		while (True):
	
			data = self.get_data()
	
			self.plot(data)
			self.show()
		
			if self.dynamic:
				try:
					print "replotting in %d s..." % self.delay
					for i in range(int(self.delay)):
						time.sleep(1)
						
					self.clear()
				except KeyboardInterrupt:
					self.clear()
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
		self.rx = re.compile(r'(\d+\.?\d*)\s+'*(self.N-1) + r'(\d+\.?\d*)[\n$]')	
		
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
		
		if not self.file.closed:
			self.file.close()
		
	def plot(self, data):
		"I am virtual"
		


def add_date(filename):

    #in case the file has e.g. a .txt, we want to sandwitch the date and not
    #append it directly
    fileEnding = ""
    if len(filename.split('.')) > 1:   
        fileEnding += "." + ".".join(filename.split('.')[1:])
    
    #asctime converts 'raw local time' to a more refined format.
    date = time.asctime(time.localtime()).replace(' ','_')

    originalFilename = filename.split('.')[0]


    return originalFilename + date + fileEnding;
    

def main():
    spacing = 20
    print """
#===========================================
#Paths:
#===========================================
"""   
    print "HOME:".ljust(spacing) + "%s" % paths.HOME 
    print "CODE:".ljust(spacing) + "%s" % paths.CODE
    print "IDEPath:".ljust(spacing) + "%s" % paths.IDEPath
    print "toolsPath:".ljust(spacing) + "%s" % paths.toolsPath
    print "scratchPath:".ljust(spacing) + "%s" % paths.scratchPath
    print "iniFiles:".ljust(spacing) + "%s" % paths.iniFilePath
    print "ProgramPath:".ljust(spacing) + "%s" % paths.programPath
    
    print """
#==========================================
# Misc:
#==========================================
"""
    print "QMC2programName:".ljust(spacing) + "%s" % misc.QMC2programName
    

if __name__ == "__main__":   

    if len(sys.argv) == 2:
        if sys.argv[1] == "-set_tool_path":
            sys.path.append(paths.toolsPath)
            print "path set successfully"
            
    main()
