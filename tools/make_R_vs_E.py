import os, re, sys
from pyLibQMC import paths

sys.append(os.path.join(paths.toolspath, "DCViz", "src"))

from DCViz_classes import R_vs_E

cwd = sys.argv[1]

R = []
pot = []
i = 0
k = 0
s = ""

for dir in os.listdir(cwd):
    
    path = os.path.join(cwd, dir)
    if not os.path.isdir(path):
        continue
    file = os.path.join(path, "stdout.txt")
    k += 1
    try:
       f = open(file, "r")
       F = f.read()
       f.close()
       
       core = re.findall("DiAtomCore (-?\d+\.?\d*)", F)[0]
       col = re.findall("Coulomb (-?\d+\.?\d*)", F)[0]
       kin = re.findall("Kinetic (-?\d+\.?\d*)", F)[0]

    except (IOError, IndexError):
        print "no such file!", path
        continue
    
    i += 1
    

    
    R = re.findall("R(\d+)", path)[-1]

    s += "%g %g %g %g\n" % (float(R[0] + "." + R[1:]), float(core),   float(col), float(kin))

print i/float(k)*100, "%"

with open(os.path.join(cwd, "R_vs_E.dat"), 'w') as F:
    F.write(s)

R_vs_E(os.path.join(cwd, "R_vs_E.dat"))
R_vs_E.mainloop()