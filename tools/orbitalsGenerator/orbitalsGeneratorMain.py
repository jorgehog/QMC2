import sys, os
from os.path import join as pjoin

sys.path.append(pjoin(os.getcwd(), "cases"))

import HO, hydrogenic



def main():
    
    outPath = pjoin(os.getcwd(), "output")    
    
#    orbitalSet = HO.HOOrbitals(56, toCPP=True)
    orbitalSet = hydrogenic.hydrogenicOrbitals(18, toCPP=True)
#     
    orbitalSet.TeXToFile(outPath)
    orbitalSet.CPPToFile(outPath)
    orbitalSet.extraToFile(outPath)
    
    

if __name__ == "__main__":
    main()