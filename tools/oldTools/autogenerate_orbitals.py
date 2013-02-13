import re

rawH = """class __name__ : public alphaHO {
public:

    __name__(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};"""

rawCPP_const = """__name__::__name__(double* k, double* k2, double* exp_factor)
: alphaHO(k, k2, exp_factor) {

}"""

rawCPP_eval = """double __name__::eval(const Walker* walker, int i) {
    __NECESSITIES__
    
    //__cheader__
    H = __expr__;
    return H * (*exp_factor);
    
}"""

def setNecessities(string, expr):
    
    getX2 = "\n    x2 = walker->r(i, 0)*walker->r(i, 0);"
    getY2 = "\n    y2 = walker->r(i, 1)*walker->r(i, 1);"

    x2p = "x2"
    y2p = "y2"
    
    necessities = ""

    if len(re.findall(x2p, expr)) > 1:
        necessities += getX2
    elif len(re.findall(x2p, expr)) == 1:
        string = string.replace("x2", "x*x")
    if len(re.findall(y2p, expr)) > 1:
        necessities += getY2
    elif len(re.findall(y2p, expr)) == 1:
        string = string.replace("y2", "y*y")
        
    
    return string.replace("__NECESSITIES__", necessities)

def translateNames(string):
    
    string = string.replace("x", "walker->r(i, 0)")   
    string = string.replace("walker->r(i, 0)2", "x2")
    string = string.replace("ewalker->r(i, 0)p", "exp")   
    string = string.replace("_walker->r(i, 0)", "_x")
    
    string = string.replace("y", "walker->r(i, 1)")   
    string = string.replace("walker->r(i, 1)2", "y2") 
    string = string.replace("_walker->r(i, 1)", "_y")
    
    string = string.replace("r2", "walker->get_r_i2(i)")
    
    string = string.replace("k2*", "(*k2)*").replace("k*", "(*k)*")
    
    return string

maxImp = 15;
names = ["alphaHO_%d" % i for i in range(0, maxImp)]
lapl_name = ["lapl_" + name for name in names]
dell_name = ["dell_" + name for name in names]

expr = ["1",
        "2*k*y",
        "2*k*x",
        "4*k2*y2 - 2",
        "4*k2*x*y",
        "4*k2*x2 - 2",
        "8*k2*k*y2*y - 12*k*y",
        "2*k*x*(4*k2*y2 - 2)",
        "2*k*y*(4*k2*x2 - 2)",
        "8*k2*k*x2*x - 12*k*x",
        "16*k2*k2*y2*y2 - 48*k2*y2 + 12",
      #  "2*k*x*(8*k2*k*y2*y - 12*k*y)",
        "2*k*x*(k*y*(8*k2*y2 -12))",
        "4*(2*k2*x2 - 1)*(2*k2*y2 - 1)",
         "2*k*y*(k*x*(8*k2*x2 -12))",
      #  "2*k*y*(8*k2*k*x2*x - 12*k*x)",
        "16*k2*k2*x2*x2 - 48*k2*x2 + 12"]

lapl_expr = ["k2*(k2*r2 - 2)", 
            "2*k2*k*y*(k2*r2-4)", 
            "2*k2*k*x*(k2*r2-4)", 
            "2*k2*(2*k2*y2-1)*(k2*r2-6)", 
            "4*k2*k2*x*y*(k2*r2-6)", 
            "2*k2*(2*k2*x*x-1)*(k2*r2-6)", 
            "4*k2*k*y*(2*k2*y2-3)*(k2*r2-8)", 
            "4*k2*k*x*(2*k2*y2-1)*(k2*r2-8)", 
            "4*k2*k*y*(2*k2*x2-1)*(k2*r2-8)", 
            "4*k2*k*x*(2*k2*x2-3)*(k2*r2-8)", 
            "4*k2*(4*k2*k2*y2*y2 - 12*k2*y2 + 3)*(k2*r2-10)",
            "8*k2*k2*x*y*(2*k2*y2-3)*(k2*r2-10)", 
            "4*k2*(2*k2*x2-1)*(2*k2*y2-1)*(k2*r2-10)",
            "8*k2*k2*x*y*(2*k2*x2-3)*(k2*r2-10)",
            "4*k2*(4*k2*k2*x2*x2 - 12*k2*x2 + 3)*(k2*r2-10)"]

dell_expr0 = ["-k2*x",
                "-2*k2*k*x*y", 
                "-2*k*(k2*x2-1)", 
                "x*(k2*2-4*k2*k2*y2)", 
                "-4*k2*y*(k2*x2-1)", 
                "-2*k2*x*(2*k2*x2-5)", 
                "4*k2*k*x*y*(3-2*k2*y2)", 
                "-4*k*(k2*x2-1)*(2*k2*y2-1)",
                "4*k2*k*x*y*(5-2*k2*x2)", 
                "-4*k*(2*k2*k2*x2*x2 - 9*k2*x2 +3)", 
                "-4*k2*x*(4*k2*k2*y2*y2 - 12*k2*y2 + 3)", 
                "-8*k2*y*(k2*x2-1)*(2*k2*y2-3)",
                "-4*k2*x*(2*k2*x2-5)*(2*k2*y2-1)", 
                "-8*k2*y*(2*k2*k2*x2*x2 - 9*k2*x2 + 3)", 
                "-4*k2*x*(4*k2*k2*x2*x2 - 28*k2*x2 + 27)"]

dell_expr1 = ["-k2*y", 
                "-2*k*(k2*y2-1)", 
                "-2*k2*k*x*y", 
                "-2*k2*y*(2*k2*y2-5)", 
                "-4*k2*x*(k2*y2-1)", 
                "y*(k2*2 - 4*k2*k2*x2)", 
                "-4*k*(2*k2*k2*y2*y2 - 9*k2*y2 + 3)", 
                "4*k2*k*x*y*(5-2*k2*y2)",
                "-4*k*(k2*y2-1)*(2*k2*x2-1)",
                "4*k2*k*x*y*(3-2*k2*x2)", 
                "-4*k2*y*(4*k2*k2*y2*y2 - 28*k2*y2 + 27)",
                "-8*k2*x*(2*k2*k2*y2*y2 - 9*k2*y2 + 3)", 
                "-4*k2*y*(2*k2*x2-1)*(2*k2*y2-5)", 
                "-8*k2*x*(k2*y2-1)*(2*k2*x2-3)", 
                "-4*k2*y*(4*k2*k2*x2*x2 - 12*k2*x2 + 3)"]

hRawStr = ""
cppRawStrConst = ""
cppRawStrEval = ""
for i in range(maxImp):
    
    hRawStr += rawH.replace("__name__", names[i]) + "\n\n"
    hRawStr += rawH.replace("__name__", dell_name[i] + "_x") + "\n\n"
    hRawStr += rawH.replace("__name__", dell_name[i] + "_y") + "\n\n"
    hRawStr += rawH.replace("__name__", lapl_name[i]) + "\n\n//-----------------------------------//\n\n"
   
    cppRawStrConst += rawCPP_const.replace("__name__", names[i]) + "\n\n"
    cppRawStrConst += rawCPP_const.replace("__name__",  dell_name[i] + "_x")\
                      + "\n\n"
    cppRawStrConst += rawCPP_const.replace("__name__", dell_name[i] + "_y")\
                      + "\n\n"
    cppRawStrConst += rawCPP_const.replace("__name__", lapl_name[i])\
                      + "\n\n//-----------------------------------//\n\n"
  
    tmpEval = rawCPP_eval.replace("__name__", names[i])
    tmpEval = tmpEval.replace("__expr__", expr[i])
    tmpEval = setNecessities(tmpEval, expr[i])
    tmpEval = translateNames(tmpEval)
    tmpEval = tmpEval.replace("__cheader__", expr[i])
    cppRawStrEval += tmpEval + "\n\n"       
         
    tmpEval = rawCPP_eval.replace("__name__", dell_name[i] + "_x")
    tmpEval = tmpEval.replace("__expr__", dell_expr0[i])
    tmpEval = setNecessities(tmpEval, dell_expr0[i])
    tmpEval = translateNames(tmpEval)
    tmpEval = tmpEval.replace("__cheader__", dell_expr0[i])
    cppRawStrEval += tmpEval + "\n\n"

    tmpEval = rawCPP_eval.replace("__name__", dell_name[i] + "_y")
    tmpEval = tmpEval.replace("__expr__", dell_expr1[i])
    tmpEval = setNecessities(tmpEval, dell_expr1[i])
    tmpEval = translateNames(tmpEval)
    tmpEval = tmpEval.replace("__cheader__",dell_expr1[i])
    cppRawStrEval += tmpEval + "\n\n"

    tmpEval = rawCPP_eval.replace("__name__", lapl_name[i])
    tmpEval = tmpEval.replace("__expr__", lapl_expr[i])
    tmpEval = setNecessities(tmpEval, lapl_expr[i])
    tmpEval = translateNames(tmpEval)
    tmpEval = tmpEval.replace("__cheader__", lapl_expr[i])
    cppRawStrEval += tmpEval + "\n\n//-----------------------------------//\n\n"




#WRITING THE HFILE
hFile = open("../oldFiles/alphaHO.h", 'r')

preH = ""
start = True
suffH = ""
end = False
for line in hFile:
    if line.startswith("//START OF GENERATED FUNCTIONS"):
        preH += str(line+"\n");
        start = False;
    elif line.startswith("//END OF GENERATED FUNCTIONS"):
        end = True;
        
    if start:
        preH += str(line);
    elif end:
        suffH += str(line);
        
    
hFile.close()

H = preH + hRawStr + suffH;

hFile = open("../src/BasisFunctions/alphaHO/alphaHO.h", "w")
hFile.write(H)
hFile.close()


#print H


#WRITING THE CPPFILE
cppFile = open("../oldFiles/alphaHO.cpp", 'r')


preCPP = ""
start = True
suffCPP = ""
end = False

for line in cppFile:
    if line.startswith("//START OF CONSTRUCTORS"):
        preCPP+=str(line + "\n")
        start = False
    elif line.startswith("//END OF EVAL FUNCTIONS"):
        end = True
    
    if start:
        preCPP += str(line)
    elif end:
        suffCPP += str(line)


cppFile.close()

CPP = preCPP + cppRawStrConst + "//END OF CONSTRUCTORS\n//START OF EVAL FUNCTIONS\n\n" + \
      cppRawStrEval + suffCPP
     
#print CPP
cppFile = open("../src/BasisFunctions/alphaHO/alphaHO.cpp", "w")
cppFile.write(CPP)
cppFile.close()
    
