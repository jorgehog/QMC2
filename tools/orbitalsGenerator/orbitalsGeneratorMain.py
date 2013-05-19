import sys, os
from os.path import join as pjoin

sys.path.append(pjoin(os.getcwd(), "cases"))

import HO, hydrogenic, HO_3D



def main():
    
    outPath = pjoin(os.getcwd(), "output")    
    
#    orbitalSet = HO.HOOrbitals(56, toCPP=True)
#    orbitalSet = hydrogenic.hydrogenicOrbitals(36, toCPP=True)
    orbitalSet = HO_3D.HOOrbitals3D(40, toCPP=True)
#     
    orbitalSet.TeXToFile(outPath)
    orbitalSet.CPPToFile(outPath)
    orbitalSet.extraToFile(outPath)
    
    

if __name__ == "__main__":
    main()
