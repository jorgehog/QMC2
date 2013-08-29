# -*- coding: utf-8 -*-

import sys, os, re
from os.path import join as pjoin

from pyLibQMC import paths

def stripPath(s0, s1):
    
    i = 0
    while i < len(s0):
        
        if s0[i] != s1[i]:
            break

        i+=1
    
    return s1[i:]
        
        

def getClassFiles():
    classes = []
    ignore = ["HartreeFock.cpp", "Bosons.cpp"]
    print "src", paths.src
    for dirPath, dirname, filenames in os.walk(paths.src):
        for filename in filenames:
            p = re.findall("(.+\.cpp)", filename)
            
            if p and filename not in ignore:
                classes.append(pjoin(dirPath, filename))
    
    return classes
  
def getCIncludes(classes):

    CINCLUDES = pjoin(paths.src, "QMCheaders.h ")    
    for f in classes:
        
        CINCLUDES += pjoin(f) + " "

    return CINCLUDES

if __name__ == "__main__":
    classes = getClassFiles()

    CC = "mpicxx.openmpi"
    CFLAGS = "-O3"
    
    LFLAGS = "-larmadillo -llapack -lblas"
    
    EXEC = "qmc2"
    DISTPATH = paths.buildPath
    print pjoin(DISTPATH, EXEC)
    
    CINCLUDE = getCIncludes(classes)

    compString = "%s %s %s -o %s %s" % (CC,
                                        CFLAGS,
                                        CINCLUDE,
                                        pjoin(DISTPATH, EXEC),
                                        LFLAGS)
    
    
    
    print "compiling: "
    print "%s %s" % (CC, CFLAGS)
    print
    for inc in CINCLUDE.split():
        print inc
    print
    print "-o %s %s" % (pjoin(DISTPATH, EXEC), LFLAGS)

    if len(sys.argv) > 1:
    
        if sys.argv[1] == "compile":
            os.system(compString)
  
    else:
        print "Skipped compilation."
        
