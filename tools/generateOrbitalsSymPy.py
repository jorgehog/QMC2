from sympy import *
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

r = Symbol('r', real=True, positive = True)
r2 = Symbol('r^2', real=True, positive = True)
r_2d = Symbol('r_2d', real=True, positive = True)

k = Symbol('k', real=True, positive = True)

class orbitalGenerator(object):
    
    genericFactor = sympify(1)    
    
    def __init__(self):
    
        self.stateMap = {}
        self.orbitals = [0]*int(ceil(self.maxImplemented/2.))
        
        self.makeStateMap()
        self.setupOrbitals()
        self.getGradientsAndLaplacians()
        self.simplify()
        
    def simplify(self):
        for i in range(len(self.orbitals)):
                   
            self.orbitals[i] = self.simplifyLocal(self.orbitals[i], self.stateMap[i])
            
            for j in range(self.dim): 
                self.gradients[i][j] = self.simplifyLocal(self.gradients[i][j], self.stateMap[i])
       
            self.Laplacians[i] = self.simplifyLocal(self.Laplacians[i], self.stateMap[i])
           

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

    def texOrbitalEq(self):
        return ""

    def beginTable(self):
        return r"""
\begin{table}
\begin{center}
\begin{tabular}{c|l}
"""

    def endTable(self, caption = ""):
        return r"""\end{tabular}
\caption{%s}
\end{center}
\end{table}

""" % caption

    def __str__(self):
        
        s = ""
        xi = [r'\vec i', r'\vec j', r'\vec k']
        phi = r"\phi(\vec r)"
        
        for key in sorted(self.stateMap.keys()):
            s += self.beginTable()

            qNums = str(self.stateMap[key]).strip("]").strip("[")
            
            genericFactor = self.genericFactor(self.stateMap[key])            
            
            s += "$\phi_{%d} \\rightarrow \phi_{%s}$\\\\\n" % (key, qNums)
            s += "\hline\n"
            s += "$%s$ & $%s$\\\\\n" % (phi, latex(self.orbitals[key]/genericFactor))
            s += "\hline\n"
         
            for i in range(self.dim):
                s+= "$%s\cdot \\nabla %s$ & $%s$\\\\\n" % \
                        (xi[i], phi, latex(self.gradients[key][i]/genericFactor))
           
            s += "\hline\n"
            s += "$\\nabla^2 %s$ & $%s$\\\\\n" % (phi, latex(self.Laplacians[key]/genericFactor))
            s += self.endTable(caption="Orbital expressions %s : %s. Factor $%s$ is omitted." % \
                                (self.__class__.__name__, qNums, latex(genericFactor)))
            
            if (key+1)%5 == 0:
                s += "\\clearpage\n"
            
        return s
    
    def genericFactor(self, qNums):
        return sympify(1)
        
    def get_real(self, expr):
        
        
#        print "\n::::REAL CALL\n----Before:::----"
#        print expr
        expr = expr.factor(exp(x))
        if not I in expr:
#            print "-----nothing----"
            return expr
        else:
            
            s = str(expr)
            if regxp.findall('exp\((.+)\)', s):
                arguments = regxp.findall('exp\((.+?)\)', s)
                for argument in arguments:
                    if not "I" in argument:
                        continue
#                    print "found imaginary exp function: ", argument
                    newArg = argument.replace("I*", "")                    
                    s = s.replace('exp(%s)' % (argument), 'cos(%s)' % argument.replace("I*", ""))
                    expr = expr.subs(exp(eval(argument)), cos(eval(newArg)))
#        print "----after:::----"
#        print expr
#        print "END OF FUNC\n"
        return expr
        
    def get_imag(self, expr):
#        print "\n::::IMAGINARY CALL\n----Before:::----"
#        print expr
        expr = expr.factor(exp(x))
        if not I in expr:
#            print "-----THIS EXPR IS REAL RETURN ZERO----"
            return sympify(0)
        else:
            

            s = str(expr)
            if regxp.findall('exp\((.+?)\)', s):
                arguments = regxp.findall('exp\((.+?)\)', s)
                for argument in arguments:
                    if not "I" in argument:
                        continue
#                    print "found imaginary exp function: ", argument
                    newArg = argument.replace("I*", "")       
                    s = s.replace('exp(%s)' % (argument), 'sin(%s)' % argument.replace("I*", ""))
                    expr =  expr.subs(exp(eval(argument)), sin(eval(newArg)))
#        print "----after:::----"
#        print expr
#        print "END OF FUNC\n"
        return expr

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
          
        self.expFactor = exp(-Rational(1,2)*k**2*(x**2 + y**2))
        
        super(HOOrbitals, self).__init__()
        
    def simplifyLocal(self, expr, qNums):
        expr = (expr.collect(self.expFactor)/self.expFactor).expand().collect(k)
        expr = expr.factor(k)
        return (expr*self.expFactor).subs(x**2 + y**2, r2)
    
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

where $k = \omega\alpha$ is the scaled oscillator frequency.  
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

    def genericFactor(self, qNums):
        return exp(-Rational(1,2)*k**2*r2)

class hydrogenicOrbitals(orbitalGenerator):
    dim = 3
    
    def __init__(self):
        self.setMax(10)
        
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
        
        self.expFactor = exp(-r3d)

        
        super(hydrogenicOrbitals, self).__init__()
        
        
    
    def getRadialFunc(self, n, l):
        return laguerre_l(n - l - 1, 2*l + 1, 2*r3d/n)*exp(-r3d/n)
    
    def sphere2Cart(self, func):

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
            
        return self.sphere2Cart(S)*r2_3d**l
 
#   def simplifyLocal(self, expr):
#        expr = (expr.collect(self.expFactor)/self.expFactor).expand().collect(k)
#        expr = expr.factor(k).subs(1.0*x, x).subs(1.0*y, y)
#        return (expr*self.expFactor).subs(x**2 + y**2, r2)

    def factorCore(self, hit, dim, s, newS, xi):
        if hit:
            factors = []
            lowest = 100
            for i in range(dim):
                if int(hit[i][1]) < lowest:
                    lowest = int(hit[i][1])
                print hit[i], 
             
            print
            print lowest
            
            for i in range(dim):
                orig = hit[i][0] + "*" +xi[i] + "**" + hit[i][1] + hit[i][2]
                s = s.replace(orig, "")
                factors.append(hit[i][0] + (int(hit[i][1]) - lowest)*("*" + xi[i]) + hit[i][2])
            print "---factors---"
            print factors
            print "-------------"
            if len(set(factors)) == 1:
                newS += " + " + factors[0] + "*(x**2 + y**2" + " + z**2"*(dim==3) + ")**(%d/2)" % lowest
    
            else:
                print "MISMATCH IN FACTORS"
                
        return s, newS

    def factorRadiiFunc(self, expr, dim):

        s = str(expr).replace("- ", "-")
        xi = ['x', 'y', 'z']
        
        hits = []
        for i in range(dim):
            hits.append(regxp.findall("([^\s]+?)\*%s\*\*(\d+)(\*?[^\s]*)" % xi[i], s))
        
        print "\n------%dD------GOT EXPRESSION---------------" % dim
        print expr

        print "---hits---"
        print hits
        print "----------"
        newS = ""
        for hit in zip(*hits):
            s, newS = self.factorCore(hit, dim, s, newS, xi)

        if len(hits) == 3:
            if len(hits[0]) == 1:
                s, newS = self.factorCore(hits[0], dim, s, newS)
        
        if newS:
            print "----------CONCATINATED RADII TERMS----------"
            s = "+".join(s.replace("+  +", "").split())
            print newS
            
            expr = sympify(r"%s" % (newS + " + " + s))
            expr = eval(str(expr))
            print expr
        
        print "-------------------END FACTORIZATION-----------\n"        
        return expr
     
                
    def factorRadiiTerms(self, expr):
        expr = self.factorRadiiFunc(expr, dim=3)

        expr = expr.subs(r2d, r_2d)
        expr = expr.subs(r3d, r)
        expr = expr.subs(r2_3d, r*r)
        expr = expr.subs(r2_2d, r_2d*r_2d)

        expr = self.factorRadiiFunc(expr, dim=2)
#        
        expr = expr.subs(r2d, r_2d)
        expr = expr.subs(r3d, r)
        expr = expr.subs(r2_3d, r*r) 
        expr = expr.subs(r2_2d, r_2d*r_2d)
        
        return expr
        

    def simplifyLocal(self, expr, qNums):
        generic = exp(-r3d/qNums[0])
        expr = expr.collect(generic)/generic
        
        expr = expr.subs(r2d, r_2d)
        expr = expr.subs(r3d, r)
        expr = expr.subs(r2_3d, r*r)
        expr = expr.subs(r2_2d, r_2d*r_2d)
        
        expr = ratsimp(expr)
        
        numer, denom = expr.as_numer_denom()

        print "---------START---%s----" % str(qNums).strip("[").strip("]")
        print numer
        

        recursionDepth = 2
        for i in range(recursionDepth):
            numer = self.factorRadiiTerms(numer)
        
        numer = numer.factor(pi)
        numer = self.factorRadiiTerms(numer)
        
        if r_2d in numer:
            numer = numer.factor(r_2d)
        else:
            numer = numer.factor(r)
        print "-------------------------"
        print numer
        print "-----------END-----------"        
        
        expr = numer/denom
        
#        expr = expr.factor([x, y, z])
        

        
        return expr*self.genericFactor(qNums)
    
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
        
    def genericFactor(self, qNums):
        return exp(-r/qNums[0])

def texThis(thing):

    s = r"""\documentclass[a4paper,10pt]{article}
\usepackage[utf8]{inputenc}
\usepackage{amsmath}

\title{SymPy generated orbital functions}
\author{}
\date{}

\begin{document}

\maketitle

%s

%s

\end{document}

""" % (thing.texOrbitalEq(), thing)
    return s
    
def main():
#    orbitalSet = HOOrbitals()
    orbitalSet = hydrogenicOrbitals()
    
    with open('/home/jorgmeister/scratch/orbitals.tex', 'w') as f:        
        f.write(texThis(orbitalSet))
        f.close()

    #atoms = hydrogenicOrbitals()
    #print atoms

    

if __name__ == "__main__":
    main()
