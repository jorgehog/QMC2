rawH = """class __name__ : public function {
protected:
    double alpha;
    double w;

public:

    __name__(double alpha, double w);
    
    virtual double eval(const Walker* walker, int qnum, int i, int d = 0);

};"""

rawCPP_const = """__name__::__name__(double alpha, double w) {
    this->w = w;
    this->alpha = alpha;
}"""

rawCPP_eval = """double __name__::eval(const Walker* walker, int qnum, int i, int d) {
    double r2 = walker->get_r_i2(i);

    double x = walker->r(i, 0);
    double y = walker->r(i, 1);

    double H = __expr__;

    return H * exp(-0.5 * alpha * w * r2);
}"""

lapl_name = ["lapl_HO_%d" % i for i in range(1,16)]
dell_name = ["dell_HO_%d" % i for i in range(1,16)]

lapl_expr = ["w * alpha * (w * alpha * r2 - dim)", "2 * x * w * alpha * (w * alpha * r2 - 4)", "2 * y * w * alpha * (w * alpha * r2 - 4)", "2 * (4 + w * alpha * (2 - 12 * x * x + w * alpha * r2 * (2 * x * x - 1)))", "2 * (4 + w * alpha * (2 - 12 * y * y + w * alpha * r2 * (2 * y * y - 1)))", "4 * alpha * w * x * y * (alpha * w * r2 - 6)", "4 * x * (aw2 * r2 * (2 * y * y - 1) + 4 * alpha * w * (1 - 4 * y * y) + 4)", "4 * y * (aw2 * r2 * (2 * x * x - 1) + 4 * alpha * w * (1 - 4 * x * x) + 4)", "4 * y * (aw2 * r2 * (2 * y * y - 3) + 4 * alpha * w * (3 - 4 * y * y) + 12)", "4 * x * (aw2 * r2 * (2 * x * x - 3) + 4 * alpha * w * (3 - 4 * x * x) + 12)", "8 * x * y * (aw2 * r2 * (2 * x * x - 3) + 2 * alpha * w * (9 - 10 * x * x) + 12)", "8 * x * y * (aw2 * r2 * (2 * y * y - 3) + 2 * alpha * w * (9 - 10 * y * y) + 12)", "4 * (aw2 * r2 * (3 + 4 * y * y * (y * y - 3)) + 2 * alpha * w * (36 * y * y - 3 - 20 * y * y * y * y) + 24 * (y * y - 1))", "4 * (aw2 * r2 * (3 + 4 * x * x * (x * x - 3)) + 2 * alpha * w * (36 * x * x - 3 - 20 * x * x * x * x) + 24 * (x * x - 1))", "4 * (aw2 * r2 * (4 * x * x * y * y + 1 + 2 * r2) - 2 * alpha * w * (1 + 6 * r2 + 20 * x * x * y * y) + 8 * (r2 + 1))"]

dell_expr0 = ["-w * alpha * walker->r(particle,d)", "2 * (1 - alpha * w * x * x)", "-2 * alpha * w * y*x", "2 * x * (4 + alpha * w - 2 * x * x * alpha * w)", "(4 * y * y - 2)*(-alpha * w * x)", "4 * y - 4 * x * x * y * w*alpha", "-4 * (-1 + alpha * w * x * x)*(-1 + 2 * y * y)", "-4 * x * y * (alpha * w * (2 * x * x - 1) - 4)", "4 * alpha * w * y * (3 - 2 * y * y) * x", "4 * (x * x * (alpha * w * (3 - 2 * x * x) + 6) - 3)", "-8 * y * (alpha * w * x * x * (2 * x * x - 3) - 3 * (2 * x * x - 1))", "-8 * y * (2 * y * y - 3)*(alpha * w * x * x - 1)", "-4 * alpha * w * x * (4 * y * y * y * y - 12 * y * y + 3)", "-4 * x * (3 * (alpha * w + 8) + 4 * alpha * w * x * x * x * x - 4 * x * x * (3 * alpha * w + 4))", "-4 * x * (2 * y * y - 1)*(alpha * w * (2 * x * x - 1) - 4)"]

dell_expr1 = ["-w * alpha * walker->r(particle,d)", "-2 * alpha * w * x*y", "2 * (1 - alpha * w * y * y)", "(4 * x * x - 2)*(-alpha * w * y)", "2 * y * (4 + alpha * w - 2 * y * y * alpha * w)", "4 * x - 4 * x * y * y * w*alpha", "-4 * x * y * (alpha * w * (2 * y * y - 1) - 4)", "-4 * (-1 + alpha * w * y * y)*(-1 + 2 * x * x)", "4 * (y * y * (alpha * w * (3 - 2 * y * y) + 6) - 3)", "4 * alpha * w * x * (3 - 2 * x * x) * y", "-8 * x * (2 * x * x - 3)*(alpha * w * y * y - 1)", "-8 * x * (alpha * w * y * y * (2 * y * y - 3) - 3 * (2 * y * y - 1))", "-4 * y * (3 * (alpha * w + 8) + 4 * alpha * w * y * y * y * y - 4 * y * y * (3 * alpha * w + 4))", "-4 * alpha * w * y * (4 * x * x * x * x - 12 * x * x + 3)", "-4 * y * (2 * x * x - 1)*(alpha * w * (2 * y * y - 1) - 4)"]

hRawStr = ""
cppRawStrConst = ""
cppRawStrEval = ""
for i in range(15):
    hRawStr += rawH.replace("__name__", dell_name[i] + "_0") + "\n\n"
    hRawStr += rawH.replace("__name__", dell_name[i] + "_1") + "\n\n"
    hRawStr += rawH.replace("__name__", lapl_name[i]) + "\n\n"
   
    cppRawStrConst += rawCPP_const.replace("__name__",  dell_name[i] + "_0")\
                      + "\n\n"
    cppRawStrConst += rawCPP_const.replace("__name__", dell_name[i] + "_1")\
                      + "\n\n"
    cppRawStrConst += rawCPP_const.replace("__name__", lapl_name[i])\
                      + "\n\n"
    tmpEval = rawCPP_eval.replace("__name__", dell_name[i] + "_0")
    tmpEval = tmpEval.replace("__expr__", dell_expr0[i])
    cppRawStrEval += tmpEval + "\n\n"

    tmpEval = rawCPP_eval.replace("__name__", dell_name[i] + "_1")
    tmpEval = tmpEval.replace("__expr__", dell_expr1[i])
    cppRawStrEval += tmpEval + "\n\n"

    tmpEval = rawCPP_eval.replace("__name__", lapl_name[i])
    tmpEval = tmpEval.replace("__expr__", lapl_expr[i])
    cppRawStrEval += tmpEval + "\n\n"




#WRITING THE HFILE
hFile = open("function.h", 'r')

preH = ""
for line in hFile:
    preH += str(line);
hFile.close()

hEnding = preH[-25:];
preH = preH[:-25]

H = preH + hRawStr + hEnding;

hFile = open("function.h", "w")
hFile.write(H)
hFile.close()





#WRITING THE CPPFILE
cppFile = open("functions.cpp", 'r')

preCPPconst = ""
preCPPeval = ""

toConst = True
for line in cppFile:
    if str(line).startswith("//END OF CONSTRUCTORS"):
        toConst = False
    else:
        if toConst:
            preCPPconst += str(line);
        else:
            preCPPeval += str(line);
cppFile.close()

CPP = preCPPconst + cppRawStrConst + "//END OF CONSTRUCTORS\n\n" + \
      preCPPeval + cppRawStrEval

cppFile = open("functions.cpp", "w")
cppFile.write(CPP)
cppFile.close()
    
