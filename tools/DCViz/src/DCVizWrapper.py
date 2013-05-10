# -*- coding: utf-8 -*-

import os, sys, inspect, re
from os.path import join as pjoin

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
    
def matchMode(modes, path, silent=False):
    root, name = os.path.split(path)    
    
    matchedMode = None
    for mode in modes:
        if re.findall(mode.nametag, name):
            matchedMode = mode
            break
    
    if matchedMode is None:
        if not silent:
            terminalTracker("Warning", "Found no matching nametags for specified filename")
        return
    
    terminalTracker("DCViz", "Matched [%s] with [%s]" % (name, str(matchedMode)))
    return matchedMode

def main(path, dynamic):
    
    if not os.path.exists(path):
        terminalTracker("Warning", "No such file...")
        return
        
    modes = autodetectModes()
    matchedMode = matchMode(modes, path)
    
    if not matchedMode:
        return;
        
    if dynamic:
        terminalTracker("DCViz", "Interrupt dynamic mode with CTRL+C")
        
    instance = matchedMode(path, dynamic=dynamic)
    instance.mainloop()

def mainToFile(path):
    
    if not os.path.isdir(path):
        raise Exception("Supplied path must be to a directory.")
    
    modes = autodetectModes()
    for root, dirs, files in os.walk(path):
    
        init = True    
    
        for outfile in files:
            matchedMode = matchMode(modes, pjoin(root, outfile), silent=True)

            if matchedMode is None:
                continue
            
            if matchedMode.isFamilyMember:
                if init:
                    thisTag = matchedMode.nametag
                else:
                    if re.findall(thisTag, outfile):
                        continue
            
            
            
            init = False
            plotTool = matchedMode(pjoin(root, outfile), toFile=True)
            plotTool.mainloop()
            

if __name__ == "__main__":
    dynamic = False
    toFile = False    
    
    try:
        if "-d" in sys.argv:
            dynamic = True
            sys.argv.remove("-d")
        elif "-f" in sys.argv:
            toFile = True
            sys.argv.remove("-f")
            
        path = sys.argv[1]
            
    except:
        print "Please supply a path as cmdline arg"
        sys.exit(1)
            
    if toFile:
        mainToFile(path)
    else:
        main(path, dynamic)