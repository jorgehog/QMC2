from sympy import (diff, 
                   trigsimp,
                   assoc_legendre,
                   sin, 
                   cos, 
                   exp, 
                   sympify, 
                   printing,
                   Symbol)
                   
from sympy.physics.hydrogen import R_nl
import sys, os
from os.path import join as pjoin

sys.path.append(os.getcwd())

from orbitalsGenerator_super import (orbitalGenerator,
                                     x, y, z, 
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
        print nShells
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
        func = func.subs(cos(2*phi), 1 - r2d**2/r3d**2)

        func = func.subs(cos(theta), x/r2d)
        func = func.subs(sin(theta), y/r2d)
        func = func.subs(cos(phi), z/r3d)
        func = func.subs(sin(phi), r2d/(r3d))
        
        return func
    
#    def getSphericalFunc(self, l, m):
#        
#        Y = Ylm(l, m, theta, phi)    
#        
#        if m >= 0:
#            S = self.get_real(Y)
#        else:
#            S = self.get_imag(Y)
#            
#        return self.sphere2Cart(S)
        
    def getSphericalFunc(self, l, m):
   
        if m < 0:
            m = abs(m)
            fac = sin(m*theta)
        else:
            fac = cos(m*theta)
            

        P = assoc_legendre(l, m, cos(phi))
        
        res = fac*P
        res = self.sphere2Cart(trigsimp(res))
        
        #Line which takes care of the stupid Abs when the argument is obviously real..
        res = res.subs(r2d, r_2d).subs(r3d, r).subs(r, r3d).subs(r_2d, r2d)

        return res;        
        
        

    def simplifyLocal(self, expr, qNums, subs=True):
            
        expr = expr.factor().collect(k).subs(r3d, r).subs(x**2 + y**2 + z**2, r*r).factor()

        if not subs:
            expr = expr.subs(r, r3d)
        
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
    
    def extraToFile(self, path):

        raw = """double %sOrbitals::get_dell_alpha_phi(const Walker* walker, int qnum, int i){
    
    double dphi;
    
    __code__
    
    return dphi;
    
}""" % self.name

        shell = """if (qnum == _q_) {

    __necessities__
    
        //__simple__
        
        dphi = __expr__
        
    } else """
    
        

        Z = Symbol('Z', positive=True, real=True)

        code = "    "

        for i in range(self.maxImplemented/2):

            psi = self.orbitals[i]
            qNums = self.stateMap[i] 
            genFac = self.genericFactor(qNums, basic=False)
           
            kdiff = diff(psi, k).factor(genFac)/psi*Z
            
            simple = self.makeReadable(str(kdiff))

            nec, necList = self.getNecessities(kdiff)
            expr = printing.ccode(kdiff) + ";"
            expr = self.replaceCCode(expr, necList)

            #hack to get the right indent
            nec = "\n".join([" "*4  + nec_i for nec_i in nec.split("\n")])

            subCode = shell
            subCode = subCode.replace("\n\n    __necessities__", nec)\
                             .replace("__expr__", expr)\
                             .replace("__simple__", simple)\
                             .replace("_q_", str(i))
            
            code += subCode
   
        code = code.strip("else ")
       
        ccode = raw.replace("__code__", code)
      
        with open(pjoin(path, "%sOrbitalsAlphaDeriv.cpp" % self.name), 'w') as f:
            f.write(ccode)
            f.close()
        
        
        
        
        
        
        
        
        