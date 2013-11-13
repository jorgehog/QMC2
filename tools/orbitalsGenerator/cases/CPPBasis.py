


class CPPBasis:
    
    def __init__(self):

        self.dim = 3        

        self.xi = ['x', 'y', 'z']
        
        self.includes = ["../../Walker/Walker.h"]        
        
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

        self.supClassConstShell = '''__name__::__name__(__constArgs__) {
__statements__
}
'''


                                        
        self.subClassShell = '''
class __name__ : public __superName__ {
public:

    __name__(__constArgs__) : __superName__(__constArgsRaw__) {
        __statements__
    }
        
    virtual double eval(const Walker* walker, int i);

};
'''

#        self.subClassConstShell = '''
#__name__::__name__(__constArgs__)
#: __superName__(__constArgsRaw__) {
#__statements__
#}
#'''
    
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
    Subclass Eval functions
*/

__subclassEval__
"""        

        self.rawH = """#ifndef __nameUpper___H
#define __nameUpper___H 

#include "../BasisFunctions.h"

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
#        self.subClassConstShell = self.subClassConstShell.replace("__superName__", superName)
        
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
#        self.subClassConstShell = \
#            self.subClassConstShell.replace("__constArgs__", args)
#        self.subClassConstShell = \
#            self.subClassConstShell.replace("__constArgsRaw__", argsRaw)
            
        self.supClassShell = self.supClassShell.replace("__constArgs__", args)
        self.subClassShell = self.subClassShell.replace("__constArgs__", args).replace("__constArgsRaw__", argsRaw)
        
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
        thisHeader = superName + ".h"
        if thisHeader not in self.includes:
            self.includes = [thisHeader] + self.includes
        
        for include in self.includes:
            includes += "#include \"%s\"\n" % include

        self.rawCPP = self.rawCPP.replace("__includes__", includes)
        self.rawCPP = self.rawCPP.replace("__superClassConst__", self.supClassConstShell)
        

    
    def makeSubclasses(self, name, cppConstStatements = ""):
        
        EVAL = {}        
        
        #H part
        H    = self.subClassShell.replace("__name__", name).replace("__statements__", cppConstStatements)
#        CPP  = self.subClassConstShell.replace("__name__", name)
        EVAL["phi"] = self.evalShell.replace("__name__", name)
        
        nameDell = "dell_" + name + "_"    
        
        for i in range(self.dim):
            xi = self.xi[i]
            dimName = nameDell + xi
            H    += self.subClassShell.replace("__name__", dimName).replace("__statements__", cppConstStatements)
#            CPP  += self.subClassConstShell.replace("__name__", dimName)
            EVAL["d" + xi] = self.evalShell.replace("__name__", dimName)
            
        laplName = "lapl_" + name
        
        H   += self.subClassShell.replace("__name__", laplName).replace("__statements__", cppConstStatements)
#        CPP += self.subClassConstShell.replace("__name__", laplName)   
        EVAL["lapl"] = self.evalShell.replace("__name__", laplName)
        
#        CPP = CPP.replace("__statements__", cppConstStatements)
        
#        return H, CPP, EVAL
        return H, EVAL
    
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
    
    def updateExpressions(self, evalCPP, constH, orbitals, i):
        
        size = 60        
        sep = self.getSeparator(size=size, header="  END %d  " % i)    
        
        name = "%s_%s" % (self.name, orbitals.getNameFromIndex(i))        
#        H, CPP, EVAL = self.makeSubclasses(name)
        H, EVAL = self.makeSubclasses(name)
        
        constH += H + sep
#        constCPP += CPP + sep
   
        evalCPP += self.getEvalFunc(EVAL, orbitals, i) + sep
        
        
#        return constCPP, evalCPP, constH
        return evalCPP, constH
    
    def getEvalExpr(self, orbitals, raw, expr, i):
        
        nec, simple, preCalc, ret = orbitals.getCCode(expr, i)
        
        walkerUnused = True
        iUnused = True
        for e in [nec, simple, preCalc, ret]:
            if "walker" in e:
                walkerUnused = False
            if "(i" in e or "i)" in e or ", i" in e or ",i" in e:
                iUnused = False

        if (walkerUnused or iUnused) and not nec:
            nec = "\n\n" + nec                
        if walkerUnused:
            nec += "    (void) walker;\n"
        if iUnused:
            nec += "    (void) i;\n"
        
        return raw.replace("\n\n__necessities__", nec).\
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
        
        evalCPP = ""
        constH = ""
        for i in range(orbitals.maxImplemented/2):
        
            orbConstState = self.updateOrbConstState(orbConstState, orbitals, i)
            
            evalCPP, constH =\
                self.updateExpressions(evalCPP, 
                                       constH, 
                                       orbitals, 
                                       i)
            
        
        
        self.constStatesOrbitals = "\n".join(orbConstState)
        
        self.rawH = self.rawH.replace("__subClass__", constH)
#        self.rawCPP = self.rawCPP.replace("__subclassConst__", constCPP)
        self.rawCPP = self.rawCPP.replace("__subclassEval__", evalCPP)
    
    def getH(self):
        return self.rawH
    
    def getCPP(self):
        return self.rawCPP
    
    def getOrbitalConstructorArgs(self):
        return self.constStatesOrbitals
        