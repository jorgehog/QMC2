from sympy import (Ylm, 
                   sin, 
                   cos, 
                   exp, 
                   sympify, 
                   ratsimp, 
                   pi, 
                   cancel, 
                   sqf, 
                   printing)
                   
from sympy.physics.hydrogen import R_nl
import sys, os

sys.path.append(os.getcwd())

from orbitalsGenerator_super import (orbitalGenerator,
                                     x, y, z, 
                                     x2, y2, z2,
                                     theta, phi, 
                                     r2d, r3d, r, r_2d,
                                     k)

class hydrogenicOrbitals(orbitalGenerator):
    
    dim = 3
    figsPrPage = 3
    
    def __init__(self, M, doInit=True, toCPP=False):

        self.name = "hydrogenic"        
        
        self.setMax(M)
        
        nShells = 0
        while nShells*(nShells+1)*(2*nShells+1)/6 < self.maxImplemented/2:
            nShells += 1
        self.nShells = nShells

        self.R = {}
        self.S = {}        
        for n in range(1, nShells+1):
            self.R[n] = {}
            for l in range(n):
                self.R[n][l] = self.getRadialFunc(n, l)
        
        for l in range(nShells):
            self.S[l] = {}
            for m in range(-l, l+1):
                self.S[l][m] = self.getSphericalFunc(l, m)
        
        self.expFactor = exp(-k*r3d)

        
        super(hydrogenicOrbitals, self).__init__(doInit, toCPP)

    def texOrbitalEq(self):
        return r"""
Orbitals are constructed in the following fashion:
\begin{equation*}
\phi(\vec r)_{n, l, m} = L_{n - l - 1}^{2l + 1}\Big(\frac{2r}{n}k\Big)S_{l}^{m}(\vec r)e^{-\frac{r}{n}k}
\end{equation*}   

where $n$ is the principal quantum number, $k = \alpha Z$ with $Z$ being the nucleus charge and 
$\alpha$ being the variational parameter.
$$l = 0,\, 1,\, ...,\, (n-1)$$ 
$$m = -l,\, (-l + 1),\, ...,\, (l-1),\, l$$
  
\newpage
"""      
        
    
    def getRadialFunc(self, n, l):
        return R_nl(n, l, r3d, Z=k)
    
    def sphere2Cart(self, func):

        func = func.subs(sin(2*phi), sympify(2)*z*r2d/r3d**2)
        func = func.subs(cos(2*phi), (z**2-r2d**2)/r3d**2)

        func = func.subs(cos(theta), x/r2d)
        func = func.subs(sin(theta), y/r2d)
        func = func.subs(cos(phi), z/r3d)
        func = func.subs(sin(phi), r2d/(r3d))
        
        return func
    
    def getSphericalFunc(self, l, m):
        
        Y = Ylm(l, m, theta, phi)    
        
        if m >= 0:
            S = self.get_real(Y)
        else:
            S = self.get_imag(Y)
            
        return self.sphere2Cart(S)
        

    def simplifyLocal(self, expr, qNums, subs=True):
   
        generic = self.genericFactor(qNums)
        expr = expr.collect(generic)/generic
        
        expr = expr.subs(r2d, r_2d).subs(r3d, r)
       
       
        expr = ratsimp(expr)
        expr = expr.factor(pi) 
        
    
        numer, denom = expr.as_numer_denom()
        
       
        numer = numer.subs(x*x, x2).subs(y*y, y2).subs(z*z, z2)


        numFac = numer.as_independent(x2, y2, z2, x, y, z, r, r_2d)[0]     
        if numFac != 0:
            numer = cancel(numer/numFac)
        else:
            numFac = 1

     
        numer = sqf(numer)                
        numer = self.factorRadiiTerms(numer)

    
        if r_2d in numer:

            numer = numer.subs(r**2, r_2d**2 + z2).collect(r_2d)
            numer = numer.subs(z2, r**2 - r_2d**2)
        
            numer = numer.expand()
            terms = numer.as_ordered_terms()
            for term in terms:
                tmp = numer
                numer = numer.subs(term, term.subs(r**2, r_2d**2 + z2)).expand().collect(r_2d)
              
                if len(str(numer)) > len(str(tmp)):
                    numer = tmp
           
            tmp = numer
            numer = numer.subs(z2, r**2 - r_2d**2)
                
            if len(str(numer)) > len(str(tmp)):
                numer = tmp

            numer = numer.factor(r_2d)
  
        
        else:
            numer = numer.factor(r)


        expr = numFac*numer/denom*self.genericFactor(qNums)  
        
        expr = expr.subs(x2, x*x).subs(y2, y*y).subs(z2, z*z)
        
        if not subs:
            expr = expr.subs(r_2d, r2d).subs(r, r3d)
            
        
        return expr    
        
    
    def makeStateMap(self):
                
        states = {}
        for n in range(1, self.nShells+1):
            states[n] = [sorted(range(-l, l+1), key=lambda x : abs(x)) for l in range(n)]
         
            for i in range(len(states[n])):
                states[n][i] = [-x for x in states[n][i]]
        
        i = 0
        for n in sorted(states.keys()):
            for l, mList in enumerate(states[n]):
                 for m in mList:
                     self.stateMap[i] = [n, l, m]
                     i += 1
                    
        
    def setupOrbitals(self):
        
        self.orbitals = [0 for i in range(self.maxImplemented/2)]

        for i, nlm in sorted(self.stateMap.items(), key=lambda x: x[0]):
            n, l, m = nlm
            self.orbitals[i] = self.R[n][l]*self.S[l][m]
        
    def genericFactor(self, qNums, basic=False):
        
        if basic:
            return exp(-k*r3d/qNums[0])
        
        return exp(-k*r/qNums[0])

    
    def initCPPbasis(self):
    
        self.cppBasis.dim = 3
        
        self.cppBasis.setName(self.name)
        
        self.cppBasis.setConstVars('double* k', 
                                   'double* k2', 
                                   'double* r22d',
                                   'double* r2d',
                                   'double* exp_factor')
                                   
        self.cppBasis.setMembers('double* k', 
                                 'double* k2', 
                                 'double* exp_factor',
                                 'double psi',
                                 'double x',
                                 'double y',
                                 'double z',
                                 'double x2',
                                 'double y2',
                                 'double z2',
                                 'double* r22d',
                                 'double* r2d')
        
    def getCReturn(self, expr, i):
         return "psi*(*exp_factor);" 
         
    def getCPre(self, expr, i):
         qNums = self.stateMap[i]
         return "    psi = %s;" % printing.ccode(expr/self.genericFactor(qNums))

    def makeOrbConstArgs(self, args, i):
        n = self.stateMap[i][0]
        args = args.replace('exp_factor', 'exp_factor_n%d' % n)
        return args
    
   