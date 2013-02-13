from sympy import printing, Rational, hermite, sqrt, exp, latex
                   
import sys, os

sys.path.append(os.getcwd())

from orbitalsGenerator_super import (orbitalGenerator,
                                     x, y, k, r2, r2_2d)

class HOOrbitals(orbitalGenerator):
    
    dim = 2
    
    def __init__(self, M, doInit=True, toCPP=False):
        
        self.name = "HarmonicOscillator"        
        
        self.setMax(M)
        self.Hx = []        
        self.Hy = []        
        
        self.nShells = int(0.5*(sqrt(1 + 4*self.maxImplemented) - 1))
        for i in range(self.nShells):
            self.Hx.append(hermite(i, k*x))
            self.Hy.append(hermite(i, k*y))
          
        self.expFactor = exp(-Rational(1,2)*k**2*(x**2 + y**2))
        
        super(HOOrbitals, self).__init__(doInit, toCPP)
       
    def simplifyLocal(self, expr, qNums, subs=True):
        expFac = self.genericFactor(qNums, basic=True)
        expr = (expr.collect(expFac)/expFac).expand().collect(k)
        expr = expr.factor(k)*expFac
        
        if subs:
            expr = expr.subs(x**2 + y**2, r2)
        
        return expr
    
    def makeStateMap(self):

        i = 0
        for shell in range(1, self.nShells +1):
            for nx in range(0,shell):
                ny = shell - 1 - nx
                self.stateMap[i] = [nx, ny]
                i += 1
    
    def setupOrbitals(self):
      
        for i, stateMap in self.stateMap.items():
    
            nx, ny = stateMap
            
            self.orbitals[i] = self.Hx[nx]*self.Hy[ny]*self.expFactor
        
    def texOrbitalEq(self):
        return r"""
Orbitals are constructed in the following fashion:
\begin{equation*}
\phi(\vec r)_{n_x, n_y} = H_{n_x}(kx)H_{n_y}(ky)e^{-\frac{1}{2}k^2r^2}
\end{equation*}   

where $k = \sqrt(\omega\alpha)$, with $\omega$ being the oscillator frequency and $\alpha$ being the variational parameter.  
"""
        
    def __str__(self):
        
        s = self.beginTable()
        for i, h in enumerate(self.Hx):
            s += "$H_{%d}(kx)$ & $%s$\\\\\n" % (i, latex(h))
            
        s += "\hline\n"
        
        for i, h in enumerate(self.Hy):
            s+= "$H_{%d}(ky)$ & $%s$\\\\\n" % (i, latex(h))
            
        s += self.endTable(caption="Hermite polynomials used to construct orbital functions")
        s += "\\clearpage\n"
        
        return s + orbitalGenerator.__str__(self)

    def genericFactor(self, qNums, basic=False):

        if basic:
            return exp(-Rational(1,2)*k**2*r2_2d)
        
        return exp(-Rational(1,2)*k**2*r2)


    def initCPPbasis(self):
        
        self.cppBasis.dim = 2
        
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
                                 'double x2',
                                 'double y2')
        
    def getCReturn(self, expr, i):
         return "H*(*exp_factor);" 
         
    def getCPre(self, expr, i):
         return "    H = %s;" % printing.ccode(expr/self.genericFactor(i))
