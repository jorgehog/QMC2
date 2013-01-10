# -*- coding: utf-8 -*-

import os, time, sys
from math import sin, cos, log, pi
from random import random



f = open('testcase%d.dat' % int(sys.argv[1]), 'w')
i = 1
while i < 12000:
    x = pi/1500*(int(sys.argv[1]) + 1)*i
    f.write("%g %g %g\n" % (sin(x), cos(x), log(x)))
    
    time.sleep(0.01)
    i += 1
        
f.close()


