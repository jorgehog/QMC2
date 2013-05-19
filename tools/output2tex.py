# -*- coding: utf-8 -*-

import os, sys, re
from os.path import join as pjoin



def getR(path):
    for filename in os.listdir(path):
        if not os.path.isdir(pjoin(path, filename)):
            if filename.split(".")[1] == "ini":
                ininame = filename
    
    ini = open(pjoin(path, ininame), 'r')
    iniRaw = ini.read()
    ini.close()
    
    R = None
    if re.findall("R\s*=\s*(.+)[\n$]", iniRaw):
        R = re.findall("R\s*=\s*(.+)[\n$]", iniRaw)[0]
        
    return R

def getVmcE(path):
    stdout = open(pjoin(path, "stdout.txt"), 'r')
    stdoutRaw = stdout.read()
    stdout.close()
    
    pattern = "VMC energy: (\-?\d+\.?\d*)"
    
    r = re.findall(pattern, stdoutRaw)
    
    if r:
        return float(r[0])
    else:
        return "N/A"


def getDmcE(path):
    stdout = open(pjoin(path, "stdout.txt"), 'r')
    stdoutRaw = stdout.read()
    stdout.close()
    
    pattern = "dmcE:\s*(\-?\d+\.?\d*).+\|\s*100\.?[0]*%"
    
    r = re.findall(pattern, stdoutRaw)
    
    if r:
        #r[0] = therm; r[1] = production

	#Both thermalization and main cycles succeeded.
	if len(r) ==2:
        	return float(r[1])
	else:
	#Run aborted after thermalization.
		return "~" + r[0]
    else:
	#Run aborted.
        return "N/A"


def getAlpha(path):
    stdout = open(pjoin(path, "stdout.txt"), 'r')
    stdoutRaw = stdout.read()
    stdout.close()
    
    pattern = "Finished minimizing. Final parameters:"+\
    "\s*Alpha:\s*(\d+\.?\d*)(\s*Beta:\s*\d+\.?\d*)?"
    
    r = re.findall(pattern, stdoutRaw)

    if r:
        return float(r[0][0])
    else:
        return "N/A"


def getBeta(path):
    stdout = open(pjoin(path, "stdout.txt"), 'r')
    stdoutRaw = stdout.read()
    stdout.close()
    
    pattern = "Finished minimizing. Final parameters:"+\
    "\s*Alpha:\s*\d+\.?\d*\s*Beta:\s*(\d+\.?\d*)"

    r = re.findall(pattern, stdoutRaw)

    if r:
        return float(r[0])
    else:
        return "N/A"


def getSystem(path):
    for filename in os.listdir(path):
        if not os.path.isdir(pjoin(path, filename)):
            if filename.split(".")[1] == "ini":
                ininame = filename
    
    ini = open(pjoin(path, ininame), 'r')
    iniRaw = ini.read()
    ini.close()
    
    system = "QDots"
    if re.findall("system\s*=\s*(\w+)[\n$]", iniRaw):
        system = re.findall("system\s*=\s*(\w+)[\n$]", iniRaw)[0]
        
    return system
    
    

def setGenPar(path, param):
    
    for filename in os.listdir(path):
        if not os.path.isdir(pjoin(path, filename)):
            if filename.split(".")[1] == "ini":
                ininame = filename
    
    ini = open(pjoin(path, ininame), 'r')
    iniRaw = ini.read()
    ini.close()

    npp = "n_p\s*=\s*(\d+)"
    scp = "systemConstant\s*=\s*(\d+\.?\d*)"
    if re.findall(npp, iniRaw):
        np = int(re.findall(npp, iniRaw)[0])
        param['n_p'] = np
    if re.findall(scp, iniRaw):
        sc = float(re.findall(scp, iniRaw)[0])
        param['sc'] = sc
    

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
    
def generateTex(table, mapping, system):
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
        placement += "c" + "|"*(i==1)
    
    header = ""
    for key in mapping:
        key = key.replace("alpha", r"$\alpha$")
        key = key.replace("beta", r"$\beta$")
        key = key.replace("n_p", "N")
        if system == "QDots":
            key = key.replace("sc", r"$\omega$")
        elif system == "Atoms":
            key = key.replace("sc", r"$\mathrm{Z}$")
        key = key.replace("E_DMC", "$\mathrm{E_{DMC}}$")
        key = key.replace("E_VMC", "$\mathrm{E_{VMC}}$")
        
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
    
    
def getMap(runpath):

    mapVarPar = False
    mapBeta = False
    mapVMC = False
    mapDMC = False
    
    mapVarPattern = "doMIN\s*=\s*1"
    notMapBetaPattern = "use_coulomb\s*=\s*0|use_jastrow\s*=\s*0"
    mapVMCPattern = "doVMC\s*=\s*1"
    mapDMCPattern = "doDMC\s*=\s*1"    
    
    mapping = ['n_p', 'sc']
   
         
        
    for filename in os.listdir(runpath):
        if not os.path.isdir(pjoin(runpath, filename)):
            if filename.split(".")[1] == "ini":
                ininame = filename
                ini = open(pjoin(runpath, ininame), 'r')
                iniRaw = ini.read()
                ini.close()
                
                if re.findall("system\s*=\s*Diatom\s*", iniRaw) or re.findall("system\s*=\s*DoubleWell\s*", iniRaw):
                    mapR = True
                if re.findall(mapVarPattern, iniRaw):
                    mapVarPar = True
                if re.findall(mapVMCPattern, iniRaw):
                    mapVMC = True
                if re.findall(mapDMCPattern, iniRaw):
                    mapDMC = True
                if not re.findall(notMapBetaPattern, iniRaw):
                    mapBeta = True

                if mapR:
                    mapping.append('R')
                if mapVMC:
                    mapping.append('E_VMC')
                if mapDMC:
                    mapping.append('E_DMC')
                if mapVarPar:
                    mapping.append('alpha')
                    if mapBeta:
                        mapping.append('beta')   
                        
    
    return mapping

    
    
def main():
    
    if len(sys.argv) == 1:
        print "Path must be set as cmlarg"
        sys.exit(1)
    else:
        mainDir = sys.argv[1]
    
    dirCont = [pjoin(mainDir, dir_) for dir_ in os.listdir(mainDir) \
                if os.path.isdir(pjoin(mainDir, dir_))] 

    runpaths = []
    for thing in dirCont:
  
        endings = [name.split(".")[1] for name in os.listdir(thing) \
                    if not os.path.isdir(pjoin(thing, name))]

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
        runparam['sc'] = 1.


    for i in range(len(runpaths)):
        runpath = runpaths[i]
        runparam = runparams[i]
        
        setGenPar(runpath, runparam)
        runparam['E_VMC'] = getVmcE(runpath)
        runparam['E_DMC'] = getDmcE(runpath)
        runparam['alpha'] = getAlpha(runpath)
        runparam['beta'] = getBeta(runpath)
        runparam['system'] = getSystem(runpath)
        runparam['R'] = getR(runpath)
        runparam['mapping'] = getMap(runpath) 
  
    
    systems = {}
    for runparam in runparams:
        if runparam["system"] not in systems.keys():
            systems[runparam["system"]] = [runparam]
        else:
            systems[runparam["system"]].append(runparam)
    
    texFile = ""
    
    #Initially loop and make a table for each system type
    for system in systems.keys():

        #Extract all different 'run-styles' from each system type
        maps = {}
        for runparam in systems[system]:
          
            if str(runparam['mapping']) not in maps.keys():
                maps[str(runparam['mapping'])] = [runparam]
            else:
                maps[str(runparam['mapping'])].append(runparam)

        #Loop and make table for each system types run style        
        for mapping in maps.keys():
            tableParams = zeros(len(maps[mapping]), len(eval(mapping)))

            for i in range(len(maps[mapping])):
                for key in eval(mapping):
                    tableParams[i][eval(mapping).index(key)] = maps[mapping][i][key]

            tableParams.sort()
 
            compressData(tableParams)
        

            texString = generateTex(tableParams, eval(mapping), system)
        
            texFile += "%s : %s\n\n" % (system, mapping) + texString + "\n\n"
    
    ofile = open(pjoin(mainDir, "texTable.tex"), 'w')
    ofile.write(texFile)
    ofile.close()
            
    

        
if __name__ == "__main__":
    main()
