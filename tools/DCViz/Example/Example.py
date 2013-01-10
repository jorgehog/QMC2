
import os, sys, threading

GUIpath = os.path.join(os.path.split(os.getcwd())[-2], 'GUI')
sys.path.append(GUIpath)

import DCVizGUI

class jobThread(threading.Thread):
    i = 0
    
    def run(self):
        os.system("python %s %d" % (os.path.join(os.getcwd(), 'fileGenerator.py'), self.i))

    

nThreads = 3
jobs = []

print "\n\n*****press ctrl+c after exiting GUI to end threads*****\n\n"
for i in range(nThreads):
    thread = jobThread()
    thread.i = i
    thread.start()
    jobs.append(thread)
    

DCVizGUI.main(os.getcwd())


