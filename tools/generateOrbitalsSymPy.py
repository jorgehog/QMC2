from sympy import hermite, Symbol, sqrt, diff, exp, Ylm, laguerre_l as laguerre, re, im
from sympy import cos, sin
from math import ceil

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

theta = Symbol('theta')
phi = Symbol('phi')

r2_2d = x**2 + y**2
r2d = sqrt(r2_2d)

r2_3d = r2_2d + z**2
r3d = sqrt(r2_3d)

r = Symbol('r')
r2 = Symbol('r2')

k = Symbol('k')

class orbitalGenerator(object):
    
    def __init__(self):
    
        self.stateMap = {}
        self.orbitals = [0]*int(ceil(self.maxImplemented/2.))
        
        self.makeStateMap()
        self.setupOrbitals()
        self.getGradientsAndLaplacians()
        self.simplify()
        
    def simplify(self, stop = False):
        for i in range(len(self.orbitals)):
            if self.dim == 3:
                self.orbitals[i] = self.orbitals[i].subs('(x**2 + y**2 + z**2)**(1/2)', r)
                self.orbitals[i] = self.orbitals[i].subs('x**2 + y**2 + z**2', r2)
                self.orbitals[i] = self.orbitals[i].factor(r2_2d)
            else:
                self.orbitals[i] = self.orbitals[i].subs('(x**2 + y**2)**(1/2)', r)
                self.orbitals[i] = self.orbitals[i].subs('1.0*k**2*x**2 + 1.0*k**2*y**2', k**2*r2)
                self.orbitals[i] = self.orbitals[i].factor(r2_3d)
                
        for i in range(len(self.gradients)):
            for j in range(self.dim):
                
                if self.dim == 3:
                    self.gradients[i][j] = self.gradients[i][j].subs('(x**2 + y**2 + z**2)**(1/2)', r)
                    self.gradients[i][j] = self.gradients[i][j].subs('x**2 + y**2 + z**2', r2)
                    self.gradients[i][j] = self.gradients[i][j].factor(r2_2d)
                else:
                    self.gradients[i][j] = self.gradients[i][j].subs('(x**2 + y**2)**(1/2)', r)
                    self.gradients[i][j] = self.gradients[i][j].subs('1.0*k**2*x**2 + 1.0*k**2*y**2', k**2*r2)
                    self.gradients[i][j] = self.gradients[i][j].factor(r2_3d)
        
        for i in range(len(self.Laplacians)):
             
            if self.dim == 3:
                self.Laplacians[i] = self.Laplacians[i].subs('(x**2 + y**2 + z**2)**(1/2)', r)
                self.Laplacians[i] = self.Laplacians[i].subs('x**2 + y**2 + z**2', r2)
                self.Laplacians[i] = self.Laplacians[i].factor(r2_2d)
            else:
                self.Laplacians[i] = self.Laplacians[i].subs('(x**2 + y**2)**(1/2)', r)
                self.Laplacians[i] = self.Laplacians[i].subs('1.0*k**2*x**2 + 1.0*k**2*y**2', k**2*r2)
                self.Laplacians[i] = self.Laplacians[i].factor(r2_3d)

        if not stop:
            self.simplify(True)

    def getGradientsAndLaplacians(self):
        self.Laplacians = []
        self.gradients = []
        
        for i, orbital in enumerate(self.orbitals):
            
            dx = diff(orbital, x)                        
            ddx = diff(dx, x)
            
            dy = diff(orbital, y)
            ddy = diff(dy, y)

            gradient = [dx, dy]            
            Laplace = ddx + ddy
            
            if self.dim == 3:
                dz = diff(orbital, z)
                ddz = diff(dz, z)
                gradient.append(dz)
                Laplace += ddz
     
            self.gradients.append([dx.simplify() for dx in gradient])
            self.Laplacians.append(Laplace.simplify())
            
    def setMax(self, M):
        self.maxImplemented = M

    def __str__(self):
        
        s = ""
        xi = ['x', 'y', 'z']
        
        for key in sorted(self.stateMap.keys()):
            s += "%2d -> " % key
            for n in self.stateMap[key]:
                s += "%d " % n
            s += "--------------------------\n"
            s += "Orbital = %s\n" % self.orbitals[key]
            s += "Gradient:\n"
         
            for i in range(self.dim):
                s+= "%s: %s\n" % (xi[i], self.gradients[key][i])
           
                
            s += "Laplace = %s\n\n" % (self.Laplacians[key])
    
            
        return s.strip("\n")

class HOOrbitals(orbitalGenerator):
    dim = 2
    
    def __init__(self):
        
        self.setMax(30)
        self.Hx = []        
        self.Hy = []        
        
        self.nShells = int(0.5*(sqrt(1 + 4*self.maxImplemented) - 1))
        for i in range(self.nShells):
            self.Hx.append(hermite(i, k*x))
            self.Hy.append(hermite(i, k*y))
            
        
        super(HOOrbitals, self).__init__()
        
        
    def makeStateMap(self):

        i = 0
        for shell in range(1, self.nShells +1):
            for nx in range(0,shell):
                ny = shell - 1 - nx
                self.stateMap[i] = [nx, ny]
                i += 1
    
    def setupOrbitals(self):
  
        expFactor = exp(-0.5*k*k*r2_2d)
        
        for i, stateMap in self.stateMap.items():
    
            nx, ny = stateMap
            
            self.orbitals[i] = self.Hx[nx]*self.Hy[ny]*expFactor
        
            
    def __str__(self):
        
        s = ""
        for i, h in enumerate(self.Hx):
            s += "H%d(kx) = %s\n" % (i, h)
        s += "\n"
        
        for i, h in enumerate(self.Hy):
            s+= "H%d(ky) = %s\n" % (i, h)
        s += "\n"
    
        
        return s + orbitalGenerator.__str__(self)


class hydrogenicOrbitals(orbitalGenerator):
    dim = 3
    
    def __init__(self):
        self.setMax(10)
        
        nShells = 0
        while nShells*(nShells+1)*(2*nShells+1)/6 < self.maxImplemented/2:
            nShells += 1
        self.nShells = nShells

        self.R = {}
        self.Y = {}        
        for n in range(1, nShells+1):
            self.R[n] = {}
            for l in range(n):
                self.R[n][l] = self.getRadialFunc(n, l)
        
        for l in range(nShells):
            self.Y[l] = {}
            for m in range(-l, l+1):
                self.Y[l][m] = self.getSphericalFunc(l, m)
   
        super(hydrogenicOrbitals, self).__init__()
        
        
    def getRadialFunc(self, n, l):
        return laguerre(n - l - 1, 2*l + 1, 2*r3d/n)*exp(-r3d/n)
    
    def sphere2Cart(self, func):
        func = func.subs(cos(theta), x/r2d).subs(sin(theta), y/r2d)
        func = func.subs(cos(phi), z/r3d).subs(sin(phi), r2d/r3d)
        return func
    
    def getSphericalFunc(self, l, m):
        return self.sphere2Cart(re(Ylm(l, m, theta, phi)))
    
    def makeStateMap(self):
                
        states = {}
        for n in range(1, self.nShells+1):
            states[n] = [sorted(range(-l, l+1), key=lambda x : abs(x)) for l in range(n)]
         
            for i in range(len(states[n])):
                states[n][i] = [-x for x in states[n][i]]
        
        i = 0
        for n in sorted(states.keys()):
            for l, mList in enumerate(states[n]):
                 for m in mList:
                     self.stateMap[i] = [n, l, m]
                     i += 1
        
        
    def setupOrbitals(self):
        
        self.orbitals = [0 for i in range(self.maxImplemented/2)]

        for i, nlm in sorted(self.stateMap.items(), key=lambda x: x[0]):
            n, l, m = nlm
            self.orbitals[i] = (self.R[n][l]*self.Y[l][m]).simplify()
        


def main():        
    qdots = HOOrbitals()
    print qdots
    #atoms = hydrogenicOrbitals()
    #print atoms
    

if __name__ == "__main__":
    main()
