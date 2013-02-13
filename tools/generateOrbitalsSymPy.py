from sympy import *
from math import ceil
import re as regxp
from sympy.physics.hydrogen import R_nl

from pyLibQMC import paths
from os.path import join as pjoin

cycle = lambda a: [a[-1]] + a[:-1]
def copyList2(a):
    
    b = [[]]*len(a)
    for i in range(len(a)):
        b[i] = a[i][:]
        
    return b

x, y, z, theta, phi = symbols('x y z theta phi', real=True)

r2_2d = x**2 + y**2
r2_3d = r2_2d + z**2

r2d = sqrt(r2_2d)
r3d = sqrt(r2_3d)

r, r2, r_2d, k, x2, y2, z2 = symbols('r r^2 r_2d k x2 y2 z2', real=True, positive = True)


class CPPBasis:
    
    def __init__(self):

        self.dim = 3        

        self.xi = ['x', 'y', 'z']
        
        self.includes = ["../../QMCheaders.h"]        
        
        self.name = 'undefined'        
        
        self.superClassInfo = {'memberVars'   : {'double*' : [] ,
                                                 'double'  : [] ,
                                                 'int*'    : [] ,
                                                 'int'     : []},
    
                               'constArgs'    : {'double*' : [] ,
                                                 'double'  : [] ,
                                                 'int*'    : [] ,
                                                 'int'     : []}}
      
        self.supClassShell = '''
class __name__ : public BasisFunctions {
protected:

__members__

public:

    __name__(__constArgs__);

};
'''
                                        
        self.subClassShell = '''
class __name__ : public __superName__ {
public:

    __name__(__constArgs__);
    virtual double eval(const Walker* walker, int i);

};
'''

        self.supClassConstShell = '''__name__::__name__(__constArgs__) {
__statements__
}
'''

        self.subClassConstShell = '''
__name__::__name__(__constArgs__)
: __superName__(__constArgsRaw__) {
__statements__
}
'''
    
        self.evalShell = '''
double __name__::eval(const Walker* walker, int i) {

__necessities__
    
    //__simpleExpr__
    
__preCalc__
    return __return__
    
}
'''
     
     
        self.rawCPP = """
__includes__

//Superclass Constructor
__superClassConst__


/*
    Subclass Constructors
*/

__subclassConst__


/*
    Subclass Eval functions
*/

__subclassEval__
"""        

        self.rawH = """
        
#ifndef __nameUpper___H
#define __nameUpper___H 

//Superclass 
__superClass__

/*
    Subclasses
*/

__subClass__

#endif /* __nameUpper___H */
"""
        
    
    def getSeparator(self, size=20, header=""):
        l = (int(size) - len(header))/2
        if l < 0:
            l = 0
        
        return "\n" + "/*\n    " + "-"*l + header + "-"*l + "\n*/" + "\n"

    def setupClasses(self):
        
        superName = self.name
        
        self.supClassShell = self.supClassShell.replace("__name__", superName)
        self.subClassShell = self.subClassShell.replace("__superName__", superName) 
        
        self.supClassConstShell = self.supClassConstShell.replace("__name__", superName)
        self.subClassConstShell = self.subClassConstShell.replace("__superName__", superName)
        
        args = ""
        argsRaw = ""
        for ctype in self.superClassInfo['constArgs'].keys():
            
            names = self.superClassInfo['constArgs'][ctype]
            
            if not names:
                continue
            
            for name in names:
                args += "%s %s, " % (ctype, name)
                argsRaw += "%s, " % name
        
        args = args.strip(", ")
        argsRaw = argsRaw.strip(", ") 
        self.argsRaw = argsRaw
        
        self.supClassConstShell = \
            self.supClassConstShell.replace("__constArgs__", args)
        self.subClassConstShell = \
            self.subClassConstShell.replace("__constArgs__", args)
        self.subClassConstShell = \
            self.subClassConstShell.replace("__constArgsRaw__", argsRaw)
            
        self.supClassShell = self.supClassShell.replace("__constArgs__", args)
        self.subClassShell = self.subClassShell.replace("__constArgs__", args)
        
        supConstStatements = ""
        for arg in argsRaw.split(","):
            supConstStatements += "    this->%s = %s;\n" % (arg.strip(), arg.strip())
        supConstStatements = supConstStatements.strip("\n")
        
        self.supClassConstShell = \
            self.supClassConstShell.replace("__statements__", supConstStatements)
   
        members = ""
        for ctype in self.superClassInfo['memberVars'].keys():
            
            names = self.superClassInfo['memberVars'][ctype]
            
            if not names:
                continue
            
            for name in names:
                members += "    %s %s;\n" % (ctype, name)
         
            members += "\n"
        
        members = members.strip("\n")
        
        self.supClassShell = self.supClassShell.replace("__members__", members)
        
        self.rawH = self.rawH.replace("__superClass__", self.supClassShell)
        self.rawH = self.rawH.replace("__nameUpper__", superName.upper())

        includes = ""        
        for include in self.includes:
            includes += "#include \"%s\"\n" % include

        self.rawCPP = self.rawCPP.replace("__includes__", includes)
        self.rawCPP = self.rawCPP.replace("__superClassConst__", self.supClassConstShell)
        

    
    def makeSubclasses(self, name, cppConstStatements = ""):
        
        EVAL = {}        
        
        #H part
        H    = self.subClassShell.replace("__name__", name)
        CPP  = self.subClassConstShell.replace("__name__", name)
        EVAL["phi"] = self.evalShell.replace("__name__", name)
        
        nameDell = "dell_" + name + "_"    
        
        
        for i in range(self.dim):
            xi = self.xi[i]
            dimName = nameDell + xi
            H    += self.subClassShell.replace("__name__", dimName)
            CPP  += self.subClassConstShell.replace("__name__", dimName)
            EVAL["d" + xi] = self.evalShell.replace("__name__", dimName)
            
        laplName = "lapl_" + name
        
        H   += self.subClassShell.replace("__name__", laplName)
        CPP += self.subClassConstShell.replace("__name__", laplName)   
        EVAL["lapl"] = self.evalShell.replace("__name__", laplName)
        
        CPP = CPP.replace("__statements__", cppConstStatements)
        
        return H, CPP, EVAL
    
    def setConstVars(self, *args):
        for arg in args:
            ctype, name = arg.split()
            self.superClassInfo['constArgs'][ctype].append(name)
            
    def setMembers(self, *args):
        for arg in args:
            ctype, name = arg.split()
            self.superClassInfo['memberVars'][ctype].append(name)
    
    def setName(self, name):
        self.name = name

    def updateOrbConstState(self, orbConstState, orbitals, i):
        
        bf = 'basis_functions'
        bfi = "_#_%s_d_[%d] = new _#_%s_%d_dim_(%s);\n" % (bf, 
                                                           i, 
                                                           self.name, 
                                                           i, 
                                                           orbitals.makeOrbConstArgs(self.argsRaw,
                                                                                     i))
                                                           
        orbConstState[0] += bfi.replace("_#_", "").replace("_dim_", "").replace("_d_", "")
        
        for j in range(self.dim):
            orbConstState[1] += \
            bfi.replace("_#_", "dell_").replace("_dim_", "_%s" % self.xi[j]).replace("_d_", "[%s]" % j)

        orbConstState[2] += bfi.replace("_#_", "lapl_").replace("_dim_", "").replace("_d_", "")
    
        return orbConstState
    
    def updateExpressions(self, constCPP, evalCPP, constH, orbitals, i):
        
        size = 60        
        sep = self.getSeparator(size=size, header="  END %d  " % i)    
        
        name = "%s_%d" % (self.name, i)        
        H, CPP, EVAL = self.makeSubclasses(name)
        
        constH += H + sep
        constCPP += CPP + sep
   
        evalCPP += self.getEvalFunc(EVAL, orbitals, i) + sep
        
        
        return constCPP, evalCPP, constH
    
    def getEvalExpr(self, orbitals, raw, expr, i):
        
        nec, simple, preCalc, ret = orbitals.getCCode(expr, i)
        
        return raw.replace("__necessities__", nec).\
                replace("__preCalc__", preCalc).\
                replace("__return__", ret).\
                replace("__simpleExpr__", simple)
    
    def getEvalFunc(self, EVAL, orbitals, i):
           
        s = ""
    
        #Phi part 
        s += self.getEvalExpr(orbitals, 
                              EVAL['phi'], 
                              orbitals.orbitals[i], 
                              i)
        
        #Dell part
        for j, key in enumerate(['d' + k for k in self.xi[:self.dim]]):
            s += self.getEvalExpr(orbitals, 
                                  EVAL[key],
                                  orbitals.gradients[i][j],
                                  i)
            
        #Laplace part
        s += self.getEvalExpr(orbitals, 
                              EVAL['lapl'], 
                              orbitals.Laplacians[i],
                              i)
        
        return s
        
    def make(self, orbitals):
        
        self.setupClasses()
        
        orbConstState = ["", "", ""]
        
        constCPP = ""
        evalCPP = ""
        constH = ""
        for i in range(orbitals.maxImplemented/2):
        
            orbConstState = self.updateOrbConstState(orbConstState, orbitals, i)
            
            constCPP, evalCPP, constH =\
                self.updateExpressions(constCPP, 
                                       evalCPP, 
                                       constH, 
                                       orbitals, 
                                       i)
            
        
        
        self.constStatesOrbitals = "\n".join(orbConstState)
        
        self.rawH = self.rawH.replace("__subClass__", constH)
        self.rawCPP = self.rawCPP.replace("__subclassConst__", constCPP)
        self.rawCPP = self.rawCPP.replace("__subclassEval__", evalCPP)
    
    def getH(self):
        return self.rawH
    
    def getCPP(self):
        return self.rawCPP
    
    def getOrbitalConstructorArgs(self):
        return self.constStatesOrbitals
        
 
        
class orbitalGenerator(object):
    
    name = "orbitals"    
    
    genericFactor = sympify(1)    
    xi = ['x', 'y', 'z']
    xi2 = ['x2', 'y2', 'z2']
    
    figsPrPage = 5
    
    cppBasis = CPPBasis()
    
    def __init__(self, doInit, toCPP):
        if not doInit:
            return
    
        self.stateMap = {}
        self.orbitals = [0]*int(ceil(self.maxImplemented/2.))
        
        self.makeStateMap()
        self.setupOrbitals()
        self.clearOrbitalNumericFactors()
        self.getGradientsAndLaplacians()
        self.simplify()
        
        if toCPP:
            self.initCPPbasis()
            self.cppBasis.make(self)
        
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
            numer = numer.factor()
            
            if numer != 1:
      
                numFac = numer.as_independent(x, y, z)[0]
                singleFactor = len(numer.as_ordered_terms()) == 1 and numFac == numer
               
                if numFac != 1 and (numFac not in numer.as_ordered_terms() or singleFactor):     
                    numer = numer/numFac

                
            if denom != 1:
                
                numFac = denom.as_independent(x, y, z)[0]
                singleFactor = len(denom.as_ordered_terms()) == 1 and numFac == denom
               
                if numFac != 1 and (numFac not in denom.as_ordered_terms() or singleFactor):
                    denom = denom/numFac

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
                replc = "k2*"*(d/2) + "*k*"*(d%2 != 0)
                
            #r factors
            elif hit[0] == 'r':
                replc = 'r^2*'*(d/2) + "r"*(d%2 != 0)
            
            #r_2d factors
            elif hit[0] == 'r_2d':
                replc = 'r_2d^2'*(d/2) + 'r_2d'*(d%2 != 0)                
                
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

                
        necS = ("%s\n%s" % (necS, necS2))
        
        return necS.strip("\n"), nec
        
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
        
    
class HOOrbitals(orbitalGenerator):
    
    dim = 2
    
    def __init__(self, M, doInit=True, toCPP=False):
        
        self.name = "alphaHO"        
        
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

where $k = \omega\alpha$, with $\omega$ being the oscillator frequency and $\alpha$ being the variational parameter.  
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
    
    

def main():
    
#     orbitalSet = HOOrbitals(42, toCPP=True)
#
#      orbitalSet.TeXToFile(paths.scratchPath)
#      orbitalSet.CPPToFile(paths.scratchPath)

#   
#
     orbitalSet = hydrogenicOrbitals(10, toCPP=True)
#     
     orbitalSet.TeXToFile(paths.scratchPath)
     orbitalSet.CPPToFile(paths.scratchPath)
#
#

    

if __name__ == "__main__":
    main()
