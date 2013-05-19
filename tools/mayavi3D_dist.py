import sys, os, numpy as np

from math import sqrt
from mayavi import mlab


def parseCML():
    if len(sys.argv[1]) < 2:
        raise Exception("Path to 3D file must be given as first cml arg.")
    else:
        path = sys.argv[1]
        
        if not os.path.exists(path):
            raise Exception("%s does not exist on your file system." % path)
        
    return path        
        
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

    path = parseCML()

    data = loadArmaCube(path)
    
    if "Atoms" in path or "QDots3D" in path:
        data = earthSpherify(data)
    elif "Diatom" in path:
        data = sliceXZ(data)
    
    mlab.pipeline.volume(mlab.pipeline.scalar_field(data), vmin=0, vmax=1)
    mlab.show()
    
    
if __name__ == "__main__":
    main()
