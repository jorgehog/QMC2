from sympy import *
from math import ceil
import re as regxp
import sympy.galgebra.latex_ex as tex

cycle = lambda a: [a[-1]] + a[:-1]
def copyList2(a):
    
    b = [[]]*len(a)
    for i in range(len(a)):
        b[i] = a[i][:]
        
    return b

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

x2, y2, z2 = symbols('x2 y2 z2', real=True, positive=True)


class orbitalGenerator(object):
    
    genericFactor = sympify(1)    
    xi = ['x', 'y', 'z']
    xi2 = ['x2', 'y2', 'z2']
    
    def __init__(self, doInit=True):
        if not doInit:
            return
    
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
    
    def __init__(self, M, doInit=True):
        
        self.setMax(M)
        self.Hx = []        
        self.Hy = []        
        
        self.nShells = int(0.5*(sqrt(1 + 4*self.maxImplemented) - 1))
        for i in range(self.nShells):
            self.Hx.append(hermite(i, k*x))
            self.Hy.append(hermite(i, k*y))
          
        self.expFactor = exp(-Rational(1,2)*k**2*(x**2 + y**2))
        
        super(HOOrbitals, self).__init__(doInit)
        
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
    
    def __init__(self, M, doInit=True):
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
        
        self.expFactor = exp(-r3d)

        
        super(hydrogenicOrbitals, self).__init__(doInit)

    def texOrbitalEq(self):
        return r"""
Orbitals are constructed in the following fashion:
\begin{equation*}
\phi(\vec r)_{n, l, m} = L_{n - l - 1}^{2l + 1}(2r/n)S_{l}^{m}(\vec r)e^{-\frac{r}{n}}
\end{equation*}   

where $n$ is the principal quantum numer. $l = 0, 1, ..., (n-1)$. $m = -l, -l + 1, ..., l$.  
"""      
        
    
    def getRadialFunc(self, n, l):
        return laguerre_l(n - l - 1, 2*l + 1, 2*r3d/n)*exp(-r3d/n)
    
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
            
        return self.sphere2Cart(S)*r2_3d**l

     
 
    def factorRadiiTerms(self, expr):
   
        print "\n------------GOT EXPRESSION---------------"
        print expr
        
 
        raw = str(expr).replace("- ", "-").replace("+ ", "+")
        raw = regxp.sub("\((\w)", "(+\g<1>", raw)
        raw = regxp.sub("^(\w)", "+\g<1>", raw)
        s = raw        
        print s
  
        hits, totalScore = self.getOptimizedHits(raw)

        #Nothing can be done        
        if totalScore == 0:
            print "Nothing to do."
            return expr

        for hit in zip(*hits):
            s = self.factorCore(hit, s)
#        return expr

        if s != raw:
            print "----------CONCATINATED RADII TERMS----------"
  
            print s
            
            expr = sympify(s)
            expr = eval(str(expr))
            print expr
     
        
        print "-------------------END FACTORIZATION-----------\n"        
        return expr

    def getExp(self, hit, i):
        "Input I: xi^I match. No match means a 1 exponent."

        thisExp = hit[i][1]

        if not thisExp:
            thisExp = 1
        else:
            thisExp = int(thisExp)
            
        return thisExp
        

    def getLowestExp(self, hit):
        """Finds the lowest exponent in the match, and splits up expressions
        so that ax^4 + ax^2y^2 -> ax^2x^2 + ax^2y^2 -> ax^2(x^2 + y^2)
        """
        
        lowest = 100
        for i in range(self.dim):
            
            #Placein None to optimize factorization
            if hit[i] is None:
                continue            

            thisExp = self.getExp(hit, i)
            
            if thisExp < lowest:
                lowest = thisExp
                
        return lowest

    def getFactors(self, hit):
        
        factors = [0,0,0]
        for i in range(self.dim):
                
            #Placein None to optimize factorization
            if hit[i] is None:
                factors[i] = None
                continue                
            
            #First: Collects the exponent from the expression.
            #Compresses it so it's FACTOR*xi^2 only (which means x2^1).
            thisExp = self.getExp(hit, i)
            pre = hit[i][0]

            if not pre and type(pre) is str:
                pre = "+"
            
            xiTerms = self.xi2[i] + ("**%d" % (thisExp - 1))*(thisExp - 1)
            suff = hit[i][2]
            factors[i] = "%s*%s"
            
        return factors

    def markReplacements(self, hit, s):
        
        orig = ["SHOULD NEVER APPEAR", 
                "SHOULD NEVER APPEAR", 
                "SHOULD NEVER APPEAR"]
        
        for i in range(self.dim):
            
            #Placein None to optimize factorization
            if hit[i] is None:
                continue             
            
            #The original expression
            orig[i] = hit[i][0] + "*"*(hit[i][0] not in ["+", "-"]) + (self.xi2[i]  + "**" + hit[i][1]).strip("*") + hit[i][2]
        
            #Set a trigger for where to insert the new expression
            #Adding fake ( or lineshift to not disturb other factors
            s = s.replace("(" + orig[i], "(__%s__" % self.xi2[i])
            s = s.replace(" " + orig[i], " __%s__" % self.xi2[i])

        return s, orig
            
    def checkFactors(self, factors):
       
       strippedPlus = []
       for factor in factors:
           if factor:
               if factor[0] not in ["+", "-"]:
                   strippedPlus.append("+" + factor)
               else:
                   strippedPlus.append(factor)
           else:
               #This will match an empty string
               if type(factor) is str:
                   factor = "+"
               strippedPlus.append(factor)
       
       return len(set(strippedPlus)) == 1
       

    def factorCore(self, hit, s):

#        lowest = self.getLowestExp(hit)
           
        factors = self.getFactors(hit)
        
        s, orig = self.markReplacements(hit, s)
            
        print "-------factors-------"
        for fac in factors:
            print fac
        print "---------------------"
            
    
        print
        print lowest

        
        #If all the factors match it means we got a x^2 y^2 z^2 term
        if self.checkFactors(factors):
            newS = factors[1] + "*"*(factors[1] not in ["+", "-"]) + "r**%d" % (2*lowest)
            print factors[1]
            print newS
            print s
            print orig
         
            #Insert the new factor over the prev x expression
            s = s.replace("__x2__", newS)
                
            #Remove y and z expressions
            for i in range(1, self.dim):
                s = s.replace("__%s__" % self.xi2[i], "")
            
        #as above, only a match for 2d only
        elif self.checkFactors(factors[:2]) and factors.count(None) != 2:
            newS = factors[1] + "*r_2d**%d" % (2*lowest)
            
            print newS
            print s
            print orig
            
            #overwrite the x-expression
            s = s.replace("__x2__", newS)
            
            #delete the y-expression
            s = s.replace("__y2__", "")
            
            #Reinsert the z-expression
            s = s.replace("__z2__", orig[2])
            
        
            
        else:
            for i in range(self.dim):
                s = s.replace("__%s__" % self.xi2[i], orig[i])
                
                
        print "NEW S: ", s
        print "---------------------------------------------------------"
                
                
        return s

        
    def getFactorScore(self, *xyzFactors):

        score = 0        

        lowest = self.getLowestExp(xyzFactors)
                
        factors = self.getFactors(xyzFactors, lowest)

           
        #If all the factors match it means we got a x^2 y^2 z^2 term 
        if self.checkFactors(factors): 
            score += 3
        
        #as above, only a match for 2d only
        elif self.checkFactors(factors[:2]):
            score += 2
        
        return score        
        
    def getOptimizedHits(self, raw):
        
        hits = []
        maxLen = 0
        for i in range(self.dim):
            hits.append(regxp.findall("([^\s\(]*?)\*?%s\*?\*?(\d*)(\*?[^\s\)]*)" % self.xi2[i], raw))
            if len(hits[-1]) > maxLen:
                maxLen = len(hits[-1])
                
        print "---hits---"
        

        for hit in hits:
            print hit

        print "----------"
        
                
        print "GOT MAXLEN: ", maxLen
        
        for i in range(self.dim):
            if len(hits[i]) != maxLen:
                for k in range(maxLen - len(hits[i])):
                    hits[i].append(None)
           


        goodHits = []
        totalScore = 0
        for i, xFac in enumerate(hits[0]):
            
            maxScore = 0
            localOptimal = []
            jM = i
            kM = i
            
            for j, yFac in enumerate(hits[1]):
                for k, zFac in enumerate(hits[2]):
                    
                    score = self.getFactorScore(xFac, yFac, zFac)
                    
                    if score > maxScore:

                        localOptimal = [xFac, yFac, zFac]
                        maxScore = score
                        
                        jM = j
                        kM = k
                        
            
            totalScore += maxScore
            if localOptimal:
                
                print maxScore
                print localOptimal
                print "------------------------_"
                
                for hit in hits:
                    print hit
        
                print "----------"            
                
                raw_input()
                
                
            
                        
                if maxScore == 2:
                    hits[1].pop(jM)
                    goodHits.append(localOptimal[:2] + [None])
                    
                elif maxScore == 3:
                    hits[1].pop(jM)
                    hits[2].pop(kM)
                    goodHits.append(localOptimal[:])
            
        #Transpose
        goodHits = zip(*goodHits)
        
        print "GOT MAXSCORE ", totalScore
        print "AT CONFIG "
        for hit in goodHits:
            print hit
            
        if goodHits:
            raw_input()
        
        return goodHits, totalScore
        

        

    def simplifyLocal(self, expr, qNums):
        if qNums != [2, 1, 0]:
            return expr
        print expr
        generic = exp(-r3d/qNums[0])
        expr = expr.collect(generic)/generic
        
        expr = expr.subs(r2d, r_2d)
        expr = expr.subs(r3d, r)
        expr = expr.subs(r2_3d, r*r)
        expr = expr.subs(r2_2d, r_2d*r_2d)
        
        expr = ratsimp(expr)
        expr = expr.factor(pi)   
        
        numer, denom = expr.as_numer_denom()
        
        numer = numer.subs(x*x, x2)
        numer = numer.subs(y*y, y2)
        numer = numer.subs(z*z, z2)

        print "\n\n\n\n---------START---%s----" % str(qNums).strip("[").strip("]")
        print numer
        print "-----------------"
        print denom

        numer = numer.factor(r)
        ###TEST####
        numer = x2**4 + y2**4
        ####
        numer = self.factorRadiiTerms(numer)     
        import sys        
        sys.exit(1)
        if r_2d in numer:
            numer = numer.factor(r_2d)
        else:
            numer = numer.factor(r)        
            
        numer = self.factorRadiiTerms(numer)
        
        #Magic line to eliminate z from '2d' expressions
        if r_2d in numer:
            numer = numer.subs(r*r, r_2d*r_2d + z2).subs(z2, r**2 - r_2d**2).factor(r_2d)

        numerTest = self.factorRadiiTerms(numer)
        if numerTest != numer:
            print "THIS MADE SENCE HERE:"
            print numer
            print numerTest
            raw_input()
        while numerTest != numer:    
            numer = numerTest
            print "lolol-------------------------------\n\n\n\----LOOLO.---\n\n"
            numerTest = self.factorRadiiTerms(numer) 
        
        
        print "--------WA THIS------------"
        print str(numerTest).replace("**", "^")
        print "---------------------------"
        expr = numerTest/denom
        
        print expr
        
        print "---------------END---------------"

        
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
    orbitalSet = hydrogenicOrbitals(10)
    
    with open('/home/jorgmeister/scratch/orbitals.tex', 'w') as f:        
        f.write(texThis(orbitalSet))
        f.close()

    #atoms = hydrogenicOrbitals()
    #print atoms

    

if __name__ == "__main__":
    main()
