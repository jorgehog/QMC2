
import os, sys, threading, time, numpy

from math import sin, cos, log, pi

GUIpath = os.path.join(os.path.split(os.getcwd())[-2], 'GUI')
sys.path.append(GUIpath)

import DCVizGUI

class jobThread(threading.Thread):
    j = 0
    stopEvent = threading.Event()    
    funcs = [sin, cos, lambda x: sin(log(abs(x) + 1))*cos(x)]
    family = ""    
    
    def run(self):
        
        F = self.funcs[self.j]
        
        k = 0
        x = numpy.linspace(0, 2*numpy.pi, 100)
            
        dt = 0.1
        P = 10
        w = 2*numpy.pi/P

        while not self.stopEvent.isSet() and k < 600/dt:
            f = open('testcase%s%d.dat' % (self.family, self.j), 'w')
            f.write(str(F)) 
      
            for xi in x:
                f.write("%g\n" % F(xi + k*dt*w))
            f.close()
            time.sleep(dt)
            
            k += 1
        

    

nThreads = 3
jobs = []
family = ""
if len(sys.argv) > 1:
    if sys.argv[1] == "Family":
        family = sys.argv[1]

for i in range(nThreads):
    thread = jobThread()
    thread.j = i
    thread.family = family
    thread.start()
    jobs.append(thread)


success = DCVizGUI.main(os.getcwd())

print "closing job threads...",

for job in jobs:
    job.stopEvent.set()
    
print "closed."

sys.exit(success)

