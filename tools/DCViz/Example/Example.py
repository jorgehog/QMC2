
import os, sys, threading, time

from math import sin, cos, log, pi

GUIpath = os.path.join(os.path.split(os.getcwd())[-2], 'GUI')
sys.path.append(GUIpath)

import DCVizGUI

class jobThread(threading.Thread):
    i = 0
    stopEvent = threading.Event()    
    
    
    def run(self):
        f = open('testcase%d.dat' % self.i, 'w')
        i = 1
        while i < 12000 and not self.stopEvent.isSet():
            x = pi/1500*(self.i + 1)*i
            f.write("%g %g %g\n" % (sin(x), cos(x), log(x)))
            
            time.sleep(0.01)
            i += 1
                
        f.close()

    

nThreads = 3
jobs = []

for i in range(nThreads):
    thread = jobThread()
    thread.i = i
    thread.start()
    jobs.append(thread)


success = DCVizGUI.main(os.getcwd())

print "closing job threads...",

for job in jobs:
    job.stopEvent.set()
    
print "closed."

sys.exit(success)

