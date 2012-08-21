import sys
from scipy.special import hermite

fileName = "HO_functions" + ".cpp"

#basisFile = open(fileName, 'w')

fileStr = "";
funcStr = """class __NAME__ : public function {
protected:
    __FUNCTION_VARIABLES__

public:
    
    __NAME__(__CONSTRUCTOR_ARGS__){
        __CONSTRUCTOR_TASKS__
    }
    
    virtual double eval(const Walker* walker, int qnum, int i, int d=0){
        int x = walker->r(i,0);
        int y = walker->r(i,1);
        return __EXPRESSION__;
    }
};"""


def analyzeFormula(rawFormula):
    constructorArgs = []
    terms = rawFormula.split("+")

    for term in terms:
        factors = term.split("*")
        for factor in factors:
            if factor.startswith("ARG"):
                constructorArgs.append([factor.split("_")[1],\
                                         factor.split("_")[2]])
            
            
    return constructorArgs

def setFuncVarAndConstructor(args, funcStr):
    
    funcVar = ""
    constructorArgs = ""
    constructorTasks = ""
    
    for arg in args:
        argtype, argname = arg

        funcVar += argtype.lower() + " " + argname + ";\n    " 
        constructorArgs += argtype.lower() + " " + argname + ", "
        constructorTasks += "this->" +argname + " = " + argname + ";\n        "
        

    if len(args) != 0:
        constructorArgs = constructorArgs[:-2] #stripping last comma
        constructorTasks = constructorTasks[:-9]
        funcVar = funcVar[:-5]

    funcStr = funcStr.replace("__FUNCTION_VARIABLES__", funcVar)
    funcStr = funcStr.replace("__CONSTRUCTOR_ARGS__", constructorArgs)
    funcStr = funcStr.replace("__CONSTRUCTOR_TASKS__", constructorTasks)

    return funcStr

def generateExpressions(rawFormula, qnums):
    terms = rawFormula.split("+")
    nTerms = len(terms)
    fermiLevel = len(qnums);
                 
    expressions = [""]*fermiLevel;

    levelDependentTerms = [0]*nTerms;
    for i in range(len(levelDependentTerms)):
        levelDependentTerms[i] = []

    levelIndependentTerms = [""]*nTerms
    
    for n in range(nTerms):
        term = terms[n]
        for factor in term.split("*"):
            factor = replaceWords(factor)

            
            if factor.startswith("LEVEL"):
                levelDependentTerms[n].append(factor.split("_")[1:])
            else:
                if factor.startswith("ARG"):
                    factor = factor.split("_")[-1]
                    
                levelIndependentTerms[n] += factor + "*";
        if levelIndependentTerms[n] != "":
            levelIndependentTerms[n] = levelIndependentTerms[n][:-1]

  
    for n in range(nTerms):
        func_type = ""
        for i in range(len(levelDependentTerms[n])):
            func_type += levelDependentTerms[n][i][0]

        #2D hermite based function
        if func_type == "HH":
            for I in range(len(qnums)):
                qnum = qnums[I]
                thisStrX = ""
                thisStrY = ""

                Hermite_x = hermite(qnum[0]);
                Hermite_y = hermite(qnum[1]);

                #BUGGY CBA
                maxPower = qnum[0];
                for i in Hermite_x:
                    if abs(i) > 0.001:
                        thisStrX += str(int(round(i))) + "*" + "x*"*maxPower;
                    
                        if len(thisStrX) != 0:
                            thisStrX = thisStrX[:-1] + "+"

                    maxPower -= 1
                        
                maxPower = qnum[1];
                for i in Hermite_y:
                    if abs(i) > 0.001:
                        thisStrY += str(int(round(i))) + "*" +"y*"*maxPower;
                    
                        if len(thisStrY) != 0:
                            thisStrY = thisStrY[:-1]+ "+"

                    maxPower -= 1

                
            expressions[I] += "(" + thisStrX + ")" + "*" + "(" + thisStrY + ")"\
                             + "*" + levelIndependentTerms[n]

    print expressions
    
                 
def replaceWords(word):
    word = word.replace("R2", "walker->get_r_i2(i)")
    #Add more here




    return word
                
                                           

def getQnumsHO(fermiShell):
    qnums = []
    fermiLevel = fermiShell*(fermiShell + 1);
    for shell in range(1, fermiShell +1):
        for nx in range(0,shell):
            ny = shell - 1 - nx
            qnums.append([nx, ny])
    print len(qnums)
    return qnums


print funcStr

fermiShell = 5;

qnums = getQnumsHO(fermiShell)

generatingFormula = "LEVEL_H_nx_x*LEVEL_H_ny_y*exp(-0.5*ARG_DOUBLE_w*ARG_DOUBLE_alpha*R2)"

constructorArgs =  analyzeFormula(generatingFormula)

funcStr = setFuncVarAndConstructor(constructorArgs, funcStr);

generateExpressions(generatingFormula, qnums)

print funcStr




