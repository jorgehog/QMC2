rawH = """class __name__ : public alphaHO {
public:

    __name__(double* k, double* k2);
    virtual double eval(const Walker* walker, int i) const;

};"""

rawCPP_const = """__name__::__name__(double* k, double* k2)
: alphaHO(k, k2) {

}"""

rawCPP_eval = """double __name__::eval(const Walker* walker, int i) const {

    double r2 = walker->get_r_i2(i);
__NECESSITIES__

    double H = __expr__;

    return H * exp(-0.5 * (*k2) * r2);
    
}"""

def setPosParam(string, expr, i):
    
    getX =  "    double x = walker->r(i, 0);\n"
    getY =  "    double y = walker->r(i, 1);\n"
    getX2 = "    double x2 = x*x;\n"
    getY2 = "    double y2 = y*y;\n"

    necessities = ""
    if expr[i][1][0]: #use x
        necessities += getX
        
    if expr[i][1][1]: #use y
        necessities += getY
        
    if expr[i][1][2]: #use x2
        if (getX not in necessities):
           necessities += getX 
        necessities += getX2
        
    if expr[i][1][3]: #use y2
        if (getY not in necessities):
           necessities += getY
        necessities += getY2
        
    
    return string.replace("__NECESSITIES__", necessities)



maxImp = 15;
names = ["alphaHO_%d" % i for i in range(0, maxImp)]
lapl_name = ["lapl_" + name for name in names]
dell_name = ["dell_" + name for name in names]

expr = [["1",[0,0,0,0]],
        ["2*k*y",[0,1,0,0]],
        ["2*k*x",[1,0,0,0]],
        ["4*k2*y2 - 2",[0,0,0,1]],
        ["4*k2*x*y",[1,1,0,0]],
        ["4*k2*x2 - 2",[0,0,1,0]],
        ["8*k2*k*y2*y - 12*k*y",[0,1,0,1]],
        ["2*k*x*(4*k2*y2 - 2)",[1,0,0,1]],
        ["2*k*y*(4*k2*x2 - 2)",[0,1,1,0]],
        ["8*k2*k*x2*x - 12*k*x",[1,0,1,0]],
        ["16*k2*k2*y2*y2 - 48*k2*y2 + 12",[0,0,0,1]],
        ["2*k*x*(8*k2*k*y2*y - 12*k*y)",[1,1,0,1]],
        ["4*(2*k2*x2 - 1)*(2*k2*y2 - 1)",[0,0,1,1]],
        ["2*k*y*(8*k2*k*x2*x - 12*k*x)",[1,1,1,0]],
        ["16*k2*k2*x2*x2 - 48*k2*x2 + 12",[0,0,1,0]]]

lapl_expr = [["k2*(k2*r2 - 2)", [0,0,0,0]],
            ["2*k2*k*y*(k2*r2-4)", [0,1,0,0]],
            ["2*k2*k*x*(k2*r2-4)", [1,0,1,1]],
            ["2*k2*(2*k2*y2-1)*(k2*r2-6)", [0,0,0,1]],
            ["4*k2*k2*x*y*(k2*r2-6)", [1,1,0,0]],
            ["2*k2*(2*k2*x*x-1)*(k2*r2-6)", [1,0,0,0]],
            ["4*k2*k*y*(2*k2*y2-3)*(k2*r2-8)", [0,1,0,1]],
            ["4*k2*k*x*(2*k2*y2-1)*(k2*r2-8)", [1,0,0,1]],
            ["4*k2*k*y*(2*k2*x2-1)*(k2*r2-8)", [0,1,1,0]],
            ["4*k2*k*x*(2*k2*x2-3)*(k2*r2-8)", [1,0,1,0]],
            ["4*k2*(4*k2*k2*y2*y2 - 12*k2*y2 + 3)*(k2*r2-10)", [0,0,0,1]],
            ["8*k2*k2*x*y*(2*k2*y2-3)*(k2*r2-10)", [1,1,0,1]],
            ["4*k2*(2*k2*x2-1)*(2*k2*y2-1)*(k2*r2-10)", [0,0,1,1]],
            ["8*k2*k2*x*y*(2*k2*x2-3)*(k2*r2-10)", [1,1,1,0]],
            ["4*k2*(4*k2*k2*x2*x2 - 12*k2*x2 + 3)*(k2*r2-10)", [0,0,1,0]]]

dell_expr0 = [["-k2*x", [1,0,0,0]],
                ["-2*k2*k*x*y", [1,1,0,0]],
                ["-2*k*(k2*x2-1)", [0,0,1,0]],
                ["x*(k2*2-4*k2*k2*y2)", [1,0,0,1]],
                ["-4*k2*y*(k2*x2-1)", [0,1,1,0]],
                ["-2*k2*x*(2*k2*x2-5)", [1,0,1,0]],
                ["4*k2*k*x*y*(3-2*k2*y2)", [1,1,0,1]],
                ["-4*k*(k2*x2-1)*(2*k2*y2-1)", [0,0,1,1]],
                ["4*k2*k*x*y*(5-2*k2*x2)", [1,1,1,0]],
                ["-4*k*(2*k2*k2*x2*x2 - 9*k2*x2 +3)", [0,0,1,0]],
                ["-4*k2*x*(4*k2*k2*y2*y2 - 12*k2*y2 + 3)", [1,0,0,1]],
                ["-8*k2*y*(k2*x2-1)*(2*k2*y2-3)", [0,1,1,1]],
                ["-4*k2*x*(2*k2*x2-5)*(2*k2*y2-1)", [1,0,1,1]],
                ["-8*k2*y*(2*k2*k2*x2*x2 - 9*k2*x2 + 3)", [0,1,1,0]],
                ["-4*k2*x*(4*k2*k2*x2*x2 - 28*k2*x2 + 27)", [1,0,1,0]]]

dell_expr1 = [["-k2*y", [0,1,0,0]],
                ["-2*k*(k2*y2-1)", [0,0,0,1]],
                ["-2*k2*k*x*y", [1,1,0,0]],
                ["-2*k2*y*(2*k2*y2-5)", [0,1,0,1]],
                ["-4*k2*x*(k2*y2-1)", [1,0,0,1]],
                ["y*(k2*2 - 4*k2*k2*x2)", [0,1,1,0]],
                ["-4*k*(2*k2*k2*y2*y2 - 9*k2*y2 + 3)", [0,0,0,1]],
                ["4*k2*k*x*y*(5-2*k2*y2)", [1,1,0,1]],
                ["-4*k*(k2*y2-1)*(2*k2*x2-1)", [0,0,1,1]],
                ["4*k2*k*x*y*(3-2*k2*x2)", [1,1,1,0]],
                ["-4*k2*y*(4*k2*k2*y2*y2 - 28*k2*y2 + 27)", [0,1,0,1]],
                ["-8*k2*x*(2*k2*k2*y2*y2 - 9*k2*y2 + 3)", [1,0,0,1]],
                ["-4*k2*y*(2*k2*x2-1)*(2*k2*y2-5)", [0,1,1,1]],
                ["-8*k2*x*(k2*y2-1)*(2*k2*x2-3)", [1,0,1,1]],
                ["-4*k2*y*(4*k2*k2*x2*x2 - 12*k2*x2 + 3)", [0,1,1,0]]]

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
    tmpEval = tmpEval.replace("__expr__", expr[i][0].replace("k2*", "(*k2)*").replace("k*", "(*k)*"))
    tmpEval = setPosParam(tmpEval, expr, i)
    cppRawStrEval += tmpEval + "\n\n"       
         
    tmpEval = rawCPP_eval.replace("__name__", dell_name[i] + "_x")
    tmpEval = tmpEval.replace("__expr__", dell_expr0[i][0].replace("k2*", "(*k2)*").replace("k*", "(*k)*"))
    tmpEval = setPosParam(tmpEval, dell_expr0, i)
    cppRawStrEval += tmpEval + "\n\n"

    tmpEval = rawCPP_eval.replace("__name__", dell_name[i] + "_y")
    tmpEval = tmpEval.replace("__expr__", dell_expr1[i][0].replace("k2*", "(*k2)*").replace("k*", "(*k)*"))
    tmpEval = setPosParam(tmpEval, dell_expr1, i)
    cppRawStrEval += tmpEval + "\n\n"

    tmpEval = rawCPP_eval.replace("__name__", lapl_name[i])
    tmpEval = tmpEval.replace("__expr__", lapl_expr[i][0].replace("k2*", "(*k2)*").replace("k*", "(*k)*"))
    tmpEval = setPosParam(tmpEval, lapl_expr, i)
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
      
cppFile = open("../src/BasisFunctions/alphaHO/alphaHO.cpp", "w")
cppFile.write(CPP)
cppFile.close()
    
