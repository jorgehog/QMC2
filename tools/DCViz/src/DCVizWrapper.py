# -*- coding: utf-8 -*-

import os, sys, inspect, re

#Adding the srcDir to the local pythonpath in order to avoid global pythonpath sets
srcDir= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(srcDir)
from DCViz_classes import *

def terminalTracker(head, msg):
    print  "[%s] %s" % (head.center(10), msg)

def autodetectModes():
    classfile = open(os.path.join(srcDir, 'DCViz_classes.py'), 'r')
    raw = classfile.read()
    classfile.close()

    uniqueModesNames = re.findall('^class (\w+)\(DCVizPlotter\):', raw, re.MULTILINE)
    uniqueModes = [eval(subclass) for subclass in uniqueModesNames]
    
    terminalTracker("Detector".center(10), "Found subclasses %s" \
                            % str(uniqueModesNames).strip("]").strip("["))
    
    if not uniqueModesNames:
        terminalTracker("Warning", "No subclass implementations found.")

    for mode in uniqueModes:
        instance = mode()
        try:
            instance.nametag
        except:
            terminalTracker("Warning", "Subclass %s has no attribute 'nametag' (output filename identifier)." % \
                              uniqueModesNames[uniqueModes.index(mode)])
                              
    return uniqueModes
    
    
def main(path, dynamic):
    
    modes = autodetectModes()
    root, name = os.path.split(path)    
    
    matchedMode = None
    for mode in modes:
        if re.findall(mode.nametag, name):
            matchedMode = mode
            break
    
    if matchedMode is None:
        terminalTracker("Warning", "found no matching nametags for specified filename")
        sys.exit(1)
        
    if dynamic:
        terminalTracker("DCViz", "Interrupt dynamic mode with CTRL+C")
        
    instance = matchedMode(path, dynamic=dynamic)
    instance.mainloop()
    
if __name__ == "__main__":
    dynamic = False
    
    try:
        if "-d" in sys.argv:
            dynamic = True
            sys.argv.remove("-d")
        path = sys.argv[1]
            
    except:
        print "please supply a path as cmdline arg"
        sys.exit(1)
        
    main(path, dynamic)