from sympy import hermite, Symbol, sqrt, diff, exp, Ylm, laguerre_l as laguerre, re, im, collect
from sympy import cos, sin, latex, I, pi, sympify
from math import ceil
import re as regxp
import sympy.galgebra.latex_ex as tex


x = Symbol('x', real=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)

theta = Symbol('theta', real=True)
phi = Symbol('phi', real=True)

r2_2d = x**2 + y**2
r2d = sqrt(r2_2d)

r2_3d = r2_2d + z**2
r3d = sqrt(r2_3d)

r = Symbol('r', real=True)
r2 = Symbol('r2', real=True)
r_2d = Symbol('r_2d', real=True)

k = Symbol('k', real=True)

class orbitalGenerator(object):
    
    def __init__(self):
    
        self.stateMap = {}
        self.orbitals = [0]*int(ceil(self.maxImplemented/2.))
        
        self.makeStateMap()
        self.setupOrbitals()
        self.getGradientsAndLaplacians()
        self.simplify()
        
    def simplify(self):
        for i in range(len(self.orbitals)):
                   
            self.orbitals[i] = self.simplifyLocal(self.orbitals[i])
            
            for j in range(self.dim): 
                self.gradients[i][j] = self.simplifyLocal(self.gradients[i][j])
       
            self.Laplacians[i] = self.simplifyLocal(self.Laplacians[i])
           

    def getGradientsAndLaplacians(self):
        self.Laplacians = []
        self.gradients = []
        
        for i, orbital in enumerate(self.orbitals):
            
            dx = diff(orbital, x)                        
            ddx = diff(dx, x)
            
            dy = diff(orbital, y)
            ddy = diff(dy, y)

            gradient = [dx, dy]            
            Laplace = ddx + ddy
            
            if self.dim == 3:
                dz = diff(orbital, z)
                ddz = diff(dz, z)
                gradient.append(dz)
                Laplace += ddz
     
            self.gradients.append(gradient)
            self.Laplacians.append(Laplace)
            
    def setMax(self, M):
        self.maxImplemented = M

    def simplifyLocal(self, expr):
        return expr

    def __str__(self):
        
        s = ""
        xi = ['x', 'y', 'z']
        
        for key in sorted(self.stateMap.keys()):
            s += "END"
            s += "%d : " % key
            for n in self.stateMap[key]:
                s += "%d " % n
            s += "\n"
            s += "Orbital = %s\n" % latex(self.orbitals[key])
            s += "Gradient:\n"
         
            for i in range(self.dim):
                s+= "%s: %s\n" % (xi[i], latex(self.gradients[key][i]))
           
                
            s += "Laplace = %s\n" % latex(self.Laplacians[key])
    
            
        return s

class HOOrbitals(orbitalGenerator):
    dim = 2
    
    def __init__(self):
        
        self.setMax(30)
        self.Hx = []        
        self.Hy = []        
        
        self.nShells = int(0.5*(sqrt(1 + 4*self.maxImplemented) - 1))
        for i in range(self.nShells):
            self.Hx.append(hermite(i, k*x))
            self.Hy.append(hermite(i, k*y))
            
        
        super(HOOrbitals, self).__init__()
        
    def simplifyLocal(self, expr):
        expr = (expr.collect(self.expFactor)/self.expFactor).expand().collect(k)
        expr = expr.factor(k).subs(1.0*x, x).subs(1.0*y, y)
        return (expr*self.expFactor).subs(x**2 + y**2, r2)
    
    def makeStateMap(self):

        i = 0
        for shell in range(1, self.nShells +1):
            for nx in range(0,shell):
                ny = shell - 1 - nx
                self.stateMap[i] = [nx, ny]
                i += 1
    
    def setupOrbitals(self):
  
        self.expFactor = exp(-0.5*k**2*(x**2 + y**2))
        
        for i, stateMap in self.stateMap.items():
    
            nx, ny = stateMap
            
            self.orbitals[i] = self.Hx[nx]*self.Hy[ny]*self.expFactor
        
            
    def __str__(self):
        
        s = ""
        for i, h in enumerate(self.Hx):
            s += "H%d(kx) = %s\n" % (i, latex(h))
        
        for i, h in enumerate(self.Hy):
            s+= "H%d(ky) = %s\n" % (i, latex(h))
    
        
        return s + orbitalGenerator.__str__(self)


class hydrogenicOrbitals(orbitalGenerator):
    dim = 3
    
    def __init__(self):
        self.setMax(10)
        
        nShells = 0
        while nShells*(nShells+1)*(2*nShells+1)/6 < self.maxImplemented/2:
            nShells += 1
        self.nShells = nShells

        self.R = {}
        self.Y = {}        
        for n in range(1, nShells+1):
            self.R[n] = {}
            for l in range(n):
                self.R[n][l] = self.getRadialFunc(n, l)
        
        for l in range(nShells):
            self.Y[l] = {}
            for m in range(-l, l+1):
                self.Y[l][m] = self.getSphericalFunc(l, m)
   
        super(hydrogenicOrbitals, self).__init__()
        
    
    def getRadialFunc(self, n, l):
        return laguerre(n - l - 1, 2*l + 1, 2*r3d/n)*exp(-r3d/n)
    
    def sphere2Cart(self, func):
        func = func.subs(cos(theta), x/r2d).subs(sin(theta), y/r2d)
        func = func.subs(cos(phi), z/r3d).subs(sin(phi), r2d/r3d)
        return func
    
    def getSphericalFunc(self, l, m):
        return self.sphere2Cart(self.get_real(Ylm(l, m, theta, phi)))
        
    def simplifyLocal(self, expr):
        return expr.subs(str(r2d), r_2d).subs("(x**2 + y**2 + z**2)**(1/2)", r).subs(r2_3d, r2)
    
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
        
    def get_real(self, expr):
        
        s = str(expr.factor(exp(x)))


        print s
        if not "I" in s:
            print "nothing"
            return expr
        else:
            if regxp.findall('exp\((.+)\)', s):
                arguments = regxp.findall('exp\((.+?)\)', s)
                for argument in arguments:
                    if not "I" in argument:
                        continue
                    print argument
                    s = s.replace('exp(%s)' % (argument), 'cos(%s)' % argument.replace("I*", "")).replace("I", "")
        print "after:"
        print sympify(s.replace('cos()', 'cos(1)'))
        return sympify(s.replace('cos()', 'cos(1)'))
                    
        
    def setupOrbitals(self):
        
        self.orbitals = [0 for i in range(self.maxImplemented/2)]

        for i, nlm in sorted(self.stateMap.items(), key=lambda x: x[0]):
            n, l, m = nlm
            self.orbitals[i] = (self.R[n][l]*self.Y[l][m]).simplify()
        


def latexInit():
    s = r"""\documentclass[a4paper,10pt]{article}
\usepackage[utf8]{inputenc}
\usepackage{amsmath}

\title{}
\author{}
\date{}

\begin{document}

\maketitle

\begin{equation*}
\begin{split}
"""
    return s
    
def latexEnd():
    s = r"""\end{split}
\end{equation*}


\end{document}
"""
    return s

def newEq():
    return "\end{split}\n\end{equation*}\n\\begin{equation*}\n\\begin{split}\n"

def main():
    obj = HOOrbitals()
#    obj = hydrogenicOrbitals()
    
    with open('/home/jorgmeister/scratch/orbitals.tex', 'w') as f:        
        f.write(latexInit())
        f.write(str(obj).replace("\n", r"\\" + "\n").replace("END", newEq()))
        f.write(latexEnd())
        f.close()
    

       
    
    #atoms = hydrogenicOrbitals()
    #print atoms

    

if __name__ == "__main__":
    main()
