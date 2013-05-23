# -*- coding: utf-8 -*-

from sympy import printing, Rational, hermite, sqrt, exp, latex
                   
import sys, os
from os.path import join as pjoin

sys.path.append(os.getcwd())

from orbitalsGenerator_super import (orbitalGenerator,
                                     x, y, z, k, r2, r2_3d)

class HOOrbitals3D(orbitalGenerator):
    
    dim = 3
    
    def __init__(self, M):
        
       
        super(HOOrbitals3D, self).__init__(M, "HarmonicOscillator3D")
        
        self.Hx = []        
        self.Hy = []
        self.Hz = []        
       
        nShells = 1 + (M > 2) + (M > 8) + (M > 20) + (M > 40) 
        
        self.nShells = nShells
        print nShells

        for i in range(self.nShells):
            self.Hx.append(hermite(i, k*x))
            self.Hy.append(hermite(i, k*y))
            self.Hz.append(hermite(i, k*z))            

          
        self.expFactor = exp(-Rational(1,2)*k**2*(x**2 + y**2 + z**2))
        
  
       
    def simplifyLocal(self, expr, qNums, subs=True):
        
        expFac = self.genericFactor(qNums, basic=True)

        if expFac not in expr:
            return expr
        
        expr = (expr.collect(expFac)/expFac).expand().collect(k)

        expr = expr.factor(k)*expFac

        if subs:
            expr = expr.subs(x**2 + y**2 + z**2, r2)
        
        return expr
    
    def makeStateMap(self):

        states = []
        for nx in range(self.nShells):
            for ny in range(self.nShells):
                for nz in range(self.nShells):
                    qn = [nx, ny, nz]
                    
                    if sum(qn) >= self.nShells:
                        continue
                    
                    states.append(qn)
                    
        
        states = sorted(states, key=lambda x: sum(x))
        
        for i, p in enumerate(states):
#            print p
            self.stateMap[i] = p
#            if p[1] == p[2] == 0:
#                print
        
        for p in zip(*states):
            for q in p:
                print q, "& ",
            print
            
    def setupOrbitals(self):
      
        for i, stateMap in self.stateMap.items():
            nx, ny, nz = stateMap
            
            self.orbitals[i] = self.Hx[nx]*self.Hy[ny]*self.Hz[nz]*self.expFactor
   
        
    def texOrbitalEq(self):
        return r"""
Orbitals are constructed in the following fashion:
\begin{equation*}
\phi(\vec r)_{n_x, n_y} = H_{n_x}(kx)H_{n_y}(ky)H_{n_z}(kz)e^{-\frac{1}{2}k^2r^2}
\end{equation*}   

where $k = \sqrt{\omega\alpha}$, with $\omega$ being the oscillator frequency and $\alpha$ being the variational parameter.  
"""
        
    def __str__(self):
        
        s = self.beginTable()
        for i, h in enumerate(self.Hx):
            s += "$H_{%d}(kx)$ & $%s$\\\\\n" % (i, latex(h))
            
        s += "\hline\n"
        
        for i, h in enumerate(self.Hy):
            s+= "$H_{%d}(ky)$ & $%s$\\\\\n" % (i, latex(h))
        
        s += "\hline\n"        
        
        for i, h in enumerate(self.Hz):
            s+= "$H_{%d}(kz)$ & $%s$\\\\\n" % (i, latex(h))
            
        s += self.endTable(caption="Hermite polynomials used to construct orbital functions")
        s += "\\clearpage\n"
        
        return s + orbitalGenerator.__str__(self)

    def genericFactor(self, qNums, basic=False):

        if basic:
            return exp(-Rational(1,2)*k**2*r2_3d)
        
        return exp(-Rational(1,2)*k**2*r2)


    def initCPPbasis(self):
        
        self.cppBasis.dim = 3
        
        self.cppBasis.setName(self.name)
        
        self.cppBasis.setConstVars('double* k', 
                                   'double* k2', 
                                   'double* exp_factor')
                                   
        self.cppBasis.setMembers('double* k', 
                                 'double* k2', 
                                 'double* exp_factor',
                                 'double H',
                                 'double x',
                                 'double y',
                                 'double z',
                                 'double x2',
                                 'double y2',
                                 'double z2')
        
    def getCReturn(self, expr, i):
         return "H*(*exp_factor);" 
         
    def getCPre(self, expr, i):
         return "    H = %s;" % printing.ccode(expr/self.genericFactor(i))
         
         
    def extraToFile(self, path):
        
        s = ""
        for i in range(self.maxImplemented/2):
            nx,ny,nz = self.stateMap[i]
            s += "qnums(%d, 0) = %d; qnums(%d, 1) = %d; qnums(%d, 2) = %d;\n" % (i, nx,
                                                                                 i, ny,
                                                                                 i, nz)
        
        with open(pjoin(path, "%s_qnums.cpp" % self.name), 'w') as f:
            f.write(s)
            f.close()
            