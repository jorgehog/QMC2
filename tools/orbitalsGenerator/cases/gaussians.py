# -*- coding: utf-8 -*-

from sympy import printing, Rational, sqrt, exp, latex, Symbol
                   
import sys, os
from os.path import join as pjoin
from math import log

sys.path.append(os.getcwd())

from orbitalsGenerator_super import (orbitalGenerator,
                                     x, y, z, r2, r2_3d)

a = Symbol('a', real=True, positive=True)

class gaussians(orbitalGenerator):
    
    dim = 3
    
    def __init__(self, M):
        
        super(gaussians, self).__init__(M, "gaussians")
        
        self.xp = []        
        self.yp = []
        self.zp = []        
       
        for i in range(M/2):
            self.xp.append(x**i)
            self.yp.append(y**i)
            self.zp.append(z**i)            

          
        self.expFactor = exp(-a*r2_3d)
        
        self.figsPrPage = 5
        
       
    def simplifyLocal(self, expr, qNums, subs=True):
        
        expFac = self.genericFactor(qNums, basic=True)

        if expFac not in expr:
            return expr

        expr = (expr.collect(expFac)/expFac).factor(x, y, z).collect(a).simplify().collect(x).collect(y).collect(z)
        
        expr = expr*expFac

        if subs:
            expr = expr.subs(r2_3d, r2)
        
        return expr
        
    def getNameFromIndex(self, i):
        return "".join([str(k) for k in self.stateMap[i]])
    
    def makeStateMap(self):

        states = []
        for px in range(self.maxImplemented):
            for py in range(self.maxImplemented):
                for pz in range(self.maxImplemented):
                    states.append([px, py, pz])
        
        states = sorted(states, key=lambda x: sum(x))
        
        for i, p in enumerate(states):
#            print p
            if i == self.maxImplemented/2:
                break 
            
            self.stateMap[i] = p
#            if p[1] == p[2] == 0:
#                print
#        print self.stateMap
#        for p in zip(*states):
#            for q in p:
#                print q, "& ",
#            print
            
    def setupOrbitals(self):
      
        for i, stateMap in self.stateMap.items():
            xp, yp, zp = stateMap
            self.orbitals[i] = self.xp[xp]*self.yp[yp]*self.zp[zp]*self.expFactor
   

    def genericFactor(self, qNums, basic=False):

        if basic:
            return exp(-a*r2_3d)
        
        return exp(-a*r2)


    def initCPPbasis(self):
        
        self.cppBasis.dim = 3
        
        self.cppBasis.setName(self.name)
        
        self.cppBasis.setConstVars('double a', 
                                   'double* exp_factor')
                                   
        self.cppBasis.setMembers('double a',
                                 'double* exp_factor',
                                 'double P',
                                 'double x',
                                 'double y',
                                 'double z',
                                 'double x2',
                                 'double y2',
                                 'double z2')
        
    def getCReturn(self, expr, i):
         return "P*(*exp_factor);" 
         
    def getCPre(self, expr, i):
         return "    P = %s;" % printing.ccode(expr/self.genericFactor(i))
         
         
    def extraToFile(self, path):
        
        s = "int px = powers(0);\nint py = powers(1);\nint pz = powers(2);\n\n"
        for i in range(self.maxImplemented/2):
            px,py,pz = self.stateMap[i]
            s += """if ((px == %d) && (py == %d) && (pz == %d)) {
phi = new gaussians___N__(expFactor, alpha);            
del_phi_x = new dell_gaussians___N___x(expFactor, alpha);
del_phi_y = new dell_gaussians___N___y(expFactor, alpha);
del_phi_z = new dell_gaussians___N___z(expFactor, alpha);
lapl_phi = new lapl_gaussians___N__(expFactor, alpha);
}
""".replace("__N__", "%d%d%d" % (px, py, pz)) % (px, py, pz)
            
            
        
        with open(pjoin(path, "%s_powerselector.cpp" % self.name), 'w') as f:
            f.write(s)
            f.close()
            