from sympy import (sympify, 
                   sqrt, 
                   symbols, 
                   cos, 
                   sin, 
                   exp, 
                   printing, 
                   latex, 
                   diff, 
                   I)
                   
import re as regxp
from math import ceil
from subprocess import call, PIPE

import sys, os
from os.path import join as pjoin

sys.path.append(pjoin(os.getcwd()))
from CPPBasis import CPPBasis

x, y, z, theta, phi = symbols('x y z theta phi', real=True)

r2_2d = x**2 + y**2
r2_3d = r2_2d + z**2

r2d = sqrt(r2_2d)
r3d = sqrt(r2_3d)

r, r2, r_2d, k, x2, y2, z2 = symbols('r r^2 r_2d k x2 y2 z2', real=True, positive = True)

class orbitalGenerator(object):
    
    genericFactor = sympify(1)    
    xi = ['x', 'y', 'z']
    xi2 = ['x2', 'y2', 'z2']
    
    figsPrPage = 5
    
    cppBasis = CPPBasis()
    
    def __init__(self, M, name = "someOrbitals"):
        
        self.name = name
        self.setMax(M)
        
        self.stateMap = {}
        self.orbitals = [0]*int(ceil(self.maxImplemented/2.))
    
    def closedFormify(self):
        
        self.makeStateMap()

        self.setupOrbitals()
        
        self.clearOrbitalNumericFactors()
        self.getGradientsAndLaplacians()
        self.simplify()
   
        self.initCPPbasis()
        self.cppBasis.make(self)

    def getNameFromIndex(self, i):
        return str(i)

    def extraToFile(self, path):
        return;        
        
    def initCPPbasis(self):
        raise NotImplementedError("initCPPbasis not in class " + self.__class__.__name__)        
        
    def simplify(self):
        for i in range(len(self.orbitals)):
        
            self.orbitals[i] = self.simplifyLocal(self.orbitals[i], self.stateMap[i])
            
            for j in range(self.dim): 
                self.gradients[i][j] = self.simplifyLocal(self.gradients[i][j], self.stateMap[i])
                
            self.Laplacians[i] = self.simplifyLocal(self.Laplacians[i], self.stateMap[i])

    def factorRadiiTerms(self, expr):
 
        raw = str(expr).replace("- ", "-").replace("+ ", "+")
        raw = regxp.sub("\((\w)", "(+\g<1>", raw)
        raw = regxp.sub("^(\w)", "+\g<1>", raw)
        raw = regxp.sub("([\+\-])", "\g<1>1*", raw)
        s = raw        
  
        hits = self.getOptimizedHits(raw)

        #Nothing can be done        
        if hits is None:
            return expr

        for hit in zip(*hits):
            s = self.factorCore(hit, s)

        if s != raw:  
            expr = sympify(s)
            expr = eval(str(expr))
  
        return expr

    def getExp(self, hit, i):
        "Input I: xi^I match. No match means a 1 exponent."

        thisExp = hit[i][1]

        if not thisExp:
            thisExp = 1
        else:
            thisExp = int(thisExp)
            
        return thisExp
        
    def getFactors(self, hit):
        
        factors = [0,0,0]
        for i in range(self.dim):
                
            #If we get a None, it's a dummy factor fill
            if hit[i] is None:
                factors[i] = None
                continue                
            
            #First: Collects the exponent from the expression.
            #Compresses it so it's FACTOR*xi^2 only (which means x2^1).
            thisExp = self.getExp(hit, i)
        
            pre = hit[i][0]
            
            xiTerms = (self.xi2[i] + "**%d" % (thisExp - 1))*bool(thisExp - 1)
            suff = hit[i][2]
            factors[i] = pre + ("*" + xiTerms.replace("**1", ""))*bool(xiTerms.replace("**1", "")) + ("*" + suff)*bool(suff)
            
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
            
            sPre = s            
            #Set a trigger for where to insert the new expression
            #Adding fake ( or lineshift to not disturb other factors
            s = s.replace("(" + orig[i], "(__%s__" % self.xi2[i])
            s = s.replace(" " + orig[i], " __%s__" % self.xi2[i])
            
            #This means the term is at the start of the line
            if s == sPre:
                s = s.replace(orig[i] + " ", "__%s__ " % self.xi2[i])

        return s, orig
            
    def checkFactors(self, factors):
       
       if factors[0] is None:
           return False
               
       return len(set(factors)) == 1
       

    def factorCore(self, hit, s):
           
        factors = self.getFactors(hit)
        
        s, orig = self.markReplacements(hit, s)
        
        #If all the factors match it means we got a x^2 y^2 z^2 term
        if self.checkFactors(factors):
            newS = factors[1] + "*r**2"
   
            #Insert the new factor over the prev x expression
            s = s.replace("__x2__", newS)
                
            #Remove y and z expressions
            for i in range(1, self.dim):
                s = s.replace("__%s__" % self.xi2[i], "")
            
        #as above, only a match for 2d only
        elif self.checkFactors(factors[:2]) and factors.count(None) != 2:
            newS = factors[1] + "*r_2d**2"
            
            #overwrite the x-expression
            s = s.replace("__x2__", newS)
            
            #delete the y-expression
            s = s.replace("__y2__", "")
            
            #Reinsert the z-expression
            s = s.replace("__z2__", orig[2])
    
        else:
            for i in range(self.dim):
                s = s.replace("__%s__" % self.xi2[i], orig[i])

        return s

        
    def getFactorScore(self, *xyzFactors):

        score = 0        
                
        factors = self.getFactors(xyzFactors)
           
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

        #Fill the hits with None elements to make them all of same length
        for i in range(self.dim):
            if len(hits[i]) != maxLen:
                for k in range(maxLen - len(hits[i])):
                    hits[i].append(None)
           

        #Test all combinations of factors to get the most compact factorization
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
                
                if maxScore == 2:
                    hits[1].pop(jM)
                    goodHits.append(localOptimal[:2] + [None])
                    
                elif maxScore == 3:
                    hits[1].pop(jM)
                    hits[2].pop(kM)
                    goodHits.append(localOptimal[:])
            
        #Transpose
        goodHits = zip(*goodHits)
        
        if totalScore == 0:
            return None
            
        
        return goodHits     
      
    def clearOrbitalNumericFactors(self):
        "Cleans off any numerical factors. NB: Only use on phi orbitals."        
        
        for i in range(self.maxImplemented/2):
            qNums = self.stateMap[i]
 
            simple = self.simplifyLocal(self.orbitals[i], qNums, subs=False)
            
      
            
            
            numer, denom = simple.as_numer_denom()
#            numer = self.simplifyLocal(numer, qNums, subs=False)
            
#            if i != 4:
#                continue
        
            if numer != 1:
      
                numFac = numer.as_independent(x, y, z)[0]
                singleFactor = len(numer.as_ordered_terms()) == 1 and numFac == numer

                   
#                if numFac != 1 and (numFac not in numer.as_ordered_terms() or singleFactor):     
                numer = numer/abs(numFac)
               
            if denom != 1:
                
                numFac = denom.as_independent(x, y, z)[0]
                singleFactor = len(denom.as_ordered_terms()) == 1 and numFac == denom
               
#                if numFac != 1 and (numFac not in denom.as_ordered_terms() or singleFactor):
                denom = denom/abs(numFac)
            
            numer = self.simplifyLocal(numer, qNums, subs=False)
            self.orbitals[i] = numer/denom
      
                
    
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
            
            if (key+1)%self.figsPrPage == 0:
                s += "\\clearpage\n"
            
        return s
    
    def genericFactor(self, qNums, basic=False):
        return sympify(1)
        
    def get_real(self, expr):
        
        expr = expr.factor(exp(x))
        if not I in expr:
            return expr
        else:
            
            s = str(expr)
            if regxp.findall('exp\((.+)\)', s):
                arguments = regxp.findall('exp\((.+?)\)', s)
                for argument in arguments:
                    if not "I" in argument:
                        continue

                    newArg = argument.replace("I*", "")                    
                    s = s.replace('exp(%s)' % (argument), 'cos(%s)' % argument.replace("I*", ""))
                    expr = expr.subs(exp(eval(argument)), cos(eval(newArg)))

        return expr
        
    def get_imag(self, expr):

        expr = expr.factor(exp(x))
        if not I in expr:
            return sympify(0)
        else:
            
            s = str(expr)
            if regxp.findall('exp\((.+?)\)', s):
                arguments = regxp.findall('exp\((.+?)\)', s)
                for argument in arguments:
                    if not "I" in argument:
                        continue
                    
                    newArg = argument.replace("I*", "")       
                    s = s.replace('exp(%s)' % (argument), 'sin(%s)' % argument.replace("I*", ""))
                    expr =  expr.subs(exp(eval(argument)), sin(eval(newArg)))

        return expr        

    def getTeX(self):

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

""" % (self.texOrbitalEq(), str(self))

        return s

    def TeXToFile(self, path):
        with open(pjoin(path, '%s.tex' % self.name), 'w') as f:
             f.write(self.getTeX())
             f.close()
        
        cwd = os.getcwd()
        os.chdir(path)
        call(["pdflatex", pjoin(path, '%s.tex' % self.name)], stdout=PIPE)
        os.chdir(cwd)

    def makeCPP(self):
        self.setupSuperclass()
        self.setupSubclasses()

            
    def replaceCCode(self, code, necList):
 
        hits = regxp.findall("pow\((\w+), (\d+)\)", code)
        for hit in hits:
            
            d = int(hit[1])
            
            replc = ""
            orig = "pow(%s, %d)" % (hit[0], d)
            
            #k-factors
            if hit[0] == "k":
                replc = "k2*"*(d/2) + "k*"*(d%2 != 0)
                
            #r factors
            elif hit[0] == 'r':
                replc = 'r^2*'*(d/2) + "r*"*(d%2 != 0)
            
            #r_2d factors
            elif hit[0] == 'r_2d':
                replc = 'r_2d^2*'*(d/2) + 'r_2d*'*(d%2 != 0)                
                
            #x factors
            for x in self.xi:
                x2 = "%s2" % x
                if hit[0] == x:
                    if x2 in necList:
                        replc = ("%s*" % x2)*(d/2) + \
                                ("%s" % x)*(d%2 != 0)
                    else:
                        replc =  ("%s*" % x)*d


            if replc:
                code = code.replace(orig, replc.strip("*"))
 
        #replacing r^2, r_2d^2, r_2d is trivial
        code = code.replace("r^2", "walker->get_r_i2(i)")
        code = code.replace("r_2d^2", 'r22d')
        code = code.replace("r_2d", 'r2d')
        
        #replacing r however is not so much, since walke[r] will match...
        code = regxp.sub("([^\w])(r)([^\w])", "\g<1>walker->get_r_i(i)\g<3>", code)
   
        #inserting pointer syntax
        code = code.replace("r2d", "(*r2d)").replace("r22d", "(*r22d)").replace("k2", "(*k2)")
        code = regxp.sub("([^\w*]?)k([^2\w])", "\g<1>(*k)\g<2>", code)
        
        return code
        
     

    def getCCode(self, expr, i):

        exprS = str(expr)

        nec, necList = self.getNecessities(expr)
 
        simple = self.makeReadable(exprS)

        pre = self.getCPre(expr, i)      
        pre = self.replaceCCode(pre, necList)
        
        ret = self.getCReturn(expr, i)         
        ret = self.replaceCCode(ret, necList)        
    
        return nec, simple, pre, ret   
               
    def makeReadable(self, s):
        s = s.replace("**", "^")
        return s
    
    def getNecessities(self, expr):
        
        s = printing.ccode(expr)        
        
        nec = []
        necS = ""
        necS2 = ""

        for i, x in enumerate(self.xi):
            l = len(regxp.findall("pow\(%s\, \d+\)" % x, s))
       
            if regxp.findall("[^\w]?%s[^\w]" % x, s):
                nec.append(x)
                necS += "    %s = walker->r(i, %d);\n" % (x, i)

            if l > 0:
                x2 = x + "2"
                nec.append(x2)
                necS2 += "    %s = %s*%s;\n" % (x2, x, x)

                
        necS = ("%s\n%s" % (necS, necS2)).strip("\n")
        
        #Manually add the lineshifts so that if no nec, then no lineshift        
        if necS:
            necS = "\n\n" + necS
        
        return necS, nec
        
    def getCReturn(self, expr, i):
         raise NotImplementedError("Function not implemented in subclass.")
         
    def getCPre(self, expr, i):
         raise NotImplementedError("Function not implemented in subclass.")
    
    def CPPToFile(self, path):
        
        with open(pjoin(path, '%s.cpp' % self.name), 'w') as f:        
             f.write(self.cppBasis.getCPP())
             f.close()
        
        with open(pjoin(path, '%s.h' % self.name), 'w') as f:        
             f.write(self.cppBasis.getH())
             f.close()
            
        with open(pjoin(path, '%s_orbConst.cpp' % self.name), 'w') as f:        
             f.write(self.cppBasis.getOrbitalConstructorArgs())
             f.close()
             
    def makeOrbConstArgs(self, args, i):
        return args