import sys, os, numpy as np

from math import sqrt
from mayavi import mlab


def parseCML():
    if len(sys.argv) == 1:
        raise Exception("Path to 3D file must be given as first cml arg.")
    else:
        path = sys.argv[1]
        
        if not os.path.exists(path):
            raise Exception("%s does not exist on your file system." % path)
    
    if len(sys.argv) == 4:
        vmin, vmax = [float(x) for x in sys.argv[2:]]
    else:
        vmin = 0
        vmax = 1
        
    return path, vmin, vmax       
        
def loadArmaCube(path):
    """reading a armadillo binary cube representing a 3D histogram"""
    
    with open(path, 'rb') as binFile:
        binFile.readline()
        n, m, l = [int(n_i) for n_i in binFile.readline().strip().split()]
        
        data = np.fromfile(binFile, dtype=np.float64)
        
    print ("Loaded %d data points" % data.shape),
    
    data.resize(n, m, l)            
    
    print "reshaped to a %s array. " % str(data.shape)

    return data




def earthSpherify(data):
    """Creates a spherical representation of the data with a slice to the center"""
    
    n, m, l = data.shape    
    
    f = lambda i, n: (i*2.0/(n-1) - 1)**2
    f = np.vectorize(f)

    D_i = f(xrange(0, n), n)
    D_j = f(xrange(0, m), m)
    D_k = f(xrange(0, l), l)
    
    #Create the sphere
    for i , d_i in enumerate(D_i):
        for j, d_j in enumerate(D_j):
            for k, d_k in enumerate(D_k):
                
                if sqrt(d_i + d_j + d_k) > 1:
                    data[i, j, k] = 0

    #Create the slice
    data[n/2:, m/2:, l/2:] = 0;
    
#    data[np.where(data < 0.1*data.mean())] = 0;
    
    return data;

def sliceXZ(data):
    n, m, l = data.shape

    data[:, l/2:, :] = 0;

    return data

def main():

    path, vmin, vmax = parseCML()

    data = loadArmaCube(path)
    
    try:
        path2 = path.replace("dmc", "vmc")
        data2 = loadArmaCube(path2)
        data = 2*data - data2
        print "Pure success!"
    except:
        pass
    
    print data.sum()*(2*49.9667/199.)**2
    if "Atoms" in path or "QDots3D" in path:
        data = earthSpherify(data)
    elif "Diatom" in path:
        data = sliceXZ(data)
    
    data[np.where(data < 0)] = 0
    mlab.pipeline.volume(mlab.pipeline.scalar_field(data), vmin=vmin, vmax=vmax)
    mlab.show()
    
    
if __name__ == "__main__":
    main()
