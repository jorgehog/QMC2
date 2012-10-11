# -*- coding: utf-8 -*-

#from pyLibQMC import paths

import os, sys, re



def getVmcE(path):
    stdout = open(path + "/stdout.txt", 'r')
    stdoutRaw = "\n".join(stdout.readlines())
    stdout.close()
    
    pattern = "VMC energy: (\d+\.?\d*)"
    
    return float(re.findall(pattern, stdoutRaw)[0])

def getDmcE(path):
    stdout = open(path + "/stdout.txt", 'r')
    stdoutRaw = "\n".join(stdout.readlines())
    stdout.close()
    
    pattern = "dmcE:\s*(\d+\.?\d*)\s*\|\s*Nw:\s*\d+\|\s*100\.?[0]*%"
    
    #[0] = therm; [1] = production
    return float(re.findall(pattern, stdoutRaw)[1])

def getAlpha(path):
    stdout = open(path + "/stdout.txt", 'r')
    stdoutRaw = "\n".join(stdout.readlines())
    stdout.close()
    
    pattern = "Finished minimizing. Final parameters:"+\
    "\s*Alpha:\s*(\d+\.?\d*)\s*Beta:\s*\d+\.?\d*"

    return float(re.findall(pattern, stdoutRaw)[0])

def getBeta(path):
    stdout = open(path + "/stdout.txt", 'r')
    stdoutRaw = "\n".join(stdout.readlines())
    stdout.close()
    
    pattern = "Finished minimizing. Final parameters:"+\
    "\s*Alpha:\s*\d+\.?\d*\s*Beta:\s*(\d+\.?\d*)"

    return float(re.findall(pattern, stdoutRaw)[0])

def setGenPar(path, param):
    
    for filename in os.listdir(path):
        if filename.split(".")[1] == "ini":
            ininame = filename
    
    ini = open(path + "/" + ininame, 'r')
    iniRaw = "\n".join(ini.readlines())
    ini.close()

    npp = "n_p\s*=\s*(\d+)"
    wp = "w\s*=\s*(\d+\.?\d*)"
    if re.findall(npp, iniRaw):
        np = int(re.findall(npp, iniRaw)[0])
    if re.findall(wp, iniRaw):
        w = float(re.findall(wp, iniRaw)[0])
    
    param['n_p'] = np
    param['w'] = w

def dumpData(runparams):
    print "Found data:"
    for runparam in runparams:
        for key in runparam.keys():
            print key, "\t", runparam[key]
        print "------------"
   
def dumpDataList(tableParams):
    for i in range(len(tableParams)):
        for j in range(len(tableParams[0])):
            print tableParams[i][j], "\t",
        print
    
def compressData(data):
    for j in [0, 1]:
        for i in range(1, len(data)):
            if (data[i][j] == data[i-1][j]):
                data[i][j] = ""
            elif data[i-1][j] == "":
                k=i-2
                while data[k][j] == "":
                    k-=1
                if data[i][j] == data[k][j]:
                    data[i][j] = ""

        
def zeros(x,y):
    matrix = []
    for i in range(x):
        tmp = []
        for j in range(y):
            tmp.append(0)
        matrix.append(tmp)
    return matrix
    
def generateTex(table, mapping):
    raw = r"""\begin{table}
\begin{center}
\label{}
\begin{tabular}{__placement__}
__HEADER__
\hline
__TABLE__
\end{tabular}
\caption{}
\end{center}
\end{table}
"""
    
    spacing = 8    
    
    placement = ""
    for i in range(len(mapping)):
        placement += "c"
    
    header = ""
    for key in mapping:
        key = key.replace("alpha", r"$\alpha$")
        key = key.replace("beta", r"$\beta$")
        key = key.replace("n_p", "N")
        key = key.replace("w", r"$\omega$")
        
        header += " %s &" % key.center(spacing)
    header = (header[:-1] + r"\\")


    textable=""
    for line in table:
        texline = ""
        for val in line:
            texline += " %s &" % str(val).center(spacing)
        textable += texline[:-1] + r"\\" + "\n"
    textable = textable.strip("\n")
    texcode = raw.replace("__placement__", placement)
    texcode = texcode.replace("__HEADER__", header)
    texcode = texcode.replace("__TABLE__", textable)
    
    return texcode
    
    
    
    
def main():
    
    dirCont = [os.getcwd() + "/" + i for i in os.listdir(os.getcwd())] 

    runpaths = []
    for thing in dirCont:
        if os.path.isdir(thing):
            endings = [name.split(".")[1] for name in os.listdir(thing)]
            if "ini" in endings:
                runpaths.append(thing)
        
        
    if not runpaths:
       print "Found no runs"
       sys.exit(1) 
                
    runparams = []
    for i in range(len(runpaths)):
        runparams.append({})
                
    for runparam in runparams:
        runparam['n_p'] = 2
        runparam['w'] = 1.


    for i in range(len(runpaths)):
        runpath = runpaths[i]
        runparam = runparams[i]
        
        setGenPar(runpath, runparam)
        runparam['E_{VMC}'] = getVmcE(runpath)
        runparam['E_{DMC}'] = getDmcE(runpath)
        runparam['alpha'] = getAlpha(runpath)
        runparam['beta'] = getBeta(runpath)
  
    
    mapping = ['n_p', 'w', 'E_{VMC}', 'E_{DMC}', 'alpha', 'beta']
    tableParams = zeros(len(runparams), len(mapping))

    for i in range(len(runparams)):
        for key in runparams[i].keys():
            tableParams[i][mapping.index(key)] = runparams[i][key]

    tableParams.sort()
 
    compressData(tableParams)

    texString = generateTex(tableParams, mapping)
    
    ofile = open(os.getcwd() + "/texTable.tex", 'w')
    ofile.write(texString)
    ofile.close()
            
    

        
if __name__ == "__main__":
    main()