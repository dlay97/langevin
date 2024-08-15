"""
Solves the Langevin equation for the simple test case
    $$ V(x,y) = ax^2 + by^2, m_{ij} \sim V(x,y) $$
and compares against the numerical solution
"""

import sys, os

lgvPath = '..//..//../src/'
if lgvPath not in sys.path:
    sys.path.insert(0,lgvPath)
import langevin as lgv

import numpy as np
import matplotlib.pyplot as plt

import h5py

class QuadraticPES:
    def __init__(self,a=1,b=3):
        self.a = a
        self.b = b
        
    def __call__(self,coords):
        slc = (len(coords.shape)-1)*(slice(None),)
        return self.a*coords[slc+(0,)]**2 + self.b*coords[slc+(1,)]**2
    
    def grad(self,coords):
        slc = (len(coords.shape)-1)*(slice(None),)
        
        ret = np.zeros(coords.shape[:-1]+(2,))
        ret[slc+(0,)] = 2*self.a*coords[slc+(0,)]
        ret[slc+(1,)] = 2*self.b*coords[slc+(1,)]
        return ret
    
class QuadraticInertia:
    def __init__(self,a,b,c):
        self.a = a
        self.b = b
        self.c = c
        
    def __call__(self,coords,returnFlat=False):
        shp = coords.shape[:-1]
        slc = (len(shp))*(slice(None),)
        if returnFlat:
            ret = np.zeros(shp+(3,))
            ret[slc+(0,)] = self.a[0]*coords[slc+(0,)]**2 + self.a[1]*coords[slc+(1,)]**2
            ret[slc+(1,)] = self.b[0]*coords[slc+(0,)]**2 + self.b[1]*coords[slc+(1,)]**2
            ret[slc+(2,)] = self.c[0]*coords[slc+(0,)]**2 + self.c[1]*coords[slc+(1,)]**2
        else:
            ret = np.zeros(shp+(2,2))
            ret[slc+(0,0)] = self.a[0]*coords[slc+(0,)]**2 + self.a[1]*coords[slc+(1,)]**2
            ret[slc+(0,1)] = self.b[0]*coords[slc+(0,)]**2 + self.b[1]*coords[slc+(1,)]**2
            ret[slc+(1,0)] = ret[slc+(0,1)]
            ret[slc+(1,1)] = self.c[0]*coords[slc+(0,)]**2 + self.c[1]*coords[slc+(1,)]**2
        return ret
    
    def grad(self,coords,returnFlat=False):
        shp = coords.shape[:-1]
        slc = (len(shp))*(slice(None),)
        #Arranged as (coords,) + (component,) + (grad,), so 'returnFlat' either
        #gives inertArr as a length 3 vector, or a 2x2 matrix
        if returnFlat:
            ret = np.zeros(shp+(3,2))
            for i in range(2):
                ret[slc+(0,i)] = 2*self.a[i]*coords[slc+(i,)]
                ret[slc+(1,i)] = 2*self.b[i]*coords[slc+(i,)]
                ret[slc+(2,i)] = 2*self.c[i]*coords[slc+(i,)]
        else:
            ret = np.zeros(shp+(2,2,2))
            for i in range(2):
                ret[slc+(0,0,i)] = 2*self.a[i]*coords[slc+(i,)]
                ret[slc+(0,1,i)] = 2*self.b[i]*coords[slc+(i,)]
                ret[slc+(1,0,i)] = ret[slc+(0,1,i)]
                ret[slc+(1,1,i)] = 2*self.c[i]*coords[slc+(i,)]
        return ret
    
pes = QuadraticPES()
invInert = QuadraticInertia([1,2],[0.1,0.2],[3,4])

x = np.arange(-5,5.1,0.1)
y = np.arange(-6,6.1,0.1)

meshArr = np.swapaxes(np.meshgrid(x,y),0,-1)
zz = pes(meshArr)

gamma = np.array([[0.5,0.2],
                  [0.2,2.2]])

startCoords = np.array([2.,4])
startMomenta = np.zeros(2)

y0 = np.hstack((startCoords,startMomenta))

tMax = 1

exactSolver = lgv.DeterministicSolution(gamma,pes,pes.grad,invInert,invInert.grad)
dtSave = 0.1 #See scipy.integrate.solve_ivp for documentation
sol = exactSolver.solve((0,tMax),y0,t_eval=np.arange(0,tMax+dtSave,dtSave))

print('Final exact solution:',sol.y[:2,-1])

# sys.exit()

#%% Solving with Langevin
nIndependentRuns = 20
nRuns = 1000

pythonSeed = 10 #for reproducibility purposes
np.random.seed(pythonSeed)
runSeeds = np.random.randint(1,10**5,nIndependentRuns)
print(runSeeds)

start = np.stack(nRuns*(startCoords,))

uniqueCoords = (x,y)
neckVals = np.ones(zz.shape) #expect to not stop early due to scission

meanSols = np.zeros((2,nIndependentRuns))
for i in range(nIndependentRuns):
    print('Starting run ',i)
    lgvObj = lgv.Jax(0,uniqueCoords,gamma,neckVals,zz,invInert(meshArr,returnFlat=True))
    
    #A, energyOnOuterTurningLine, dt, maxTime, scissionNeckVal
    nonArrayParams = [10000,pes(startCoords),10.**(-4),tMax,0.5]
    lgvObj.set_langevin_params(*nonArrayParams)
    
    start = np.stack(nRuns*(startCoords,))
    
    res = lgvObj.run(start,runSeeds[i])
    meanSols[:,i] = res[0].mean(axis=0)
    print('Langevin Mean:',res[0].mean(axis=0))

#%% Plotting

fig, ax = plt.subplots()
ax.contour(x,y,zz.T,colors='gray',levels=np.arange(0,100,10))

ax.scatter(sol.y[0],sol.y[1],label='Deterministic (Numerical)',zorder=100)

ax.scatter(*meanSols,marker='x',label='Langevin')

ax.set(xlim=(0,3),ylim=(-1,5),
       xlabel='x',ylabel='y',title=r'Langevin vs Exact, $V(x)=x^2$, $M_{ij}=\delta_{ij}$')
ax.legend()

os.makedirs('plots',exist_ok=True)
fig.savefig('plots/'+str(nRuns)+'.pdf',bbox_inches='tight')
