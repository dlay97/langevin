"""
Solves the Langevin equation for the simple test case
    $$ V(x,y) = ax^2 + by^2, m_{ij} = \delta_{ij} $$
and compares against the numerical solution

Currently (Sep 19 2023) only goes through the Fortran implementation
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
    
class IdentityInertia:
    #Assumes 2D
    def __call__(self,coords,returnFlat=False):
        shp = coords.shape[:-1]
        slc = (len(shp))*(slice(None),)
        if returnFlat:
            ret = np.zeros(shp+(3,))
            ret[slc+(0,)] = 1
            ret[slc+(2,)] = 1
        else:
            ret = np.zeros(shp+(2,2))
            ret[slc+(0,0)] = 1
            ret[slc+(1,1)] = 1
        return 2*ret
    
    def grad(self,coords,returnFlat=False):
        shp = coords.shape[:-1]
        #Arranged as (coords,) + (inertArr,) + (grad,), so 'returnFlat' either
        #gives inertArr as a length 3 vector, or a 2x2 matrix
        if returnFlat:
            ret = np.zeros(shp+(3,2))
        else:
            ret = np.zeros(shp+(2,2,2))
        return ret

pes = QuadraticPES()
invInert = IdentityInertia()

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

#%% Solving with Langevin
nIndependentRuns = 20
nRuns = 100

pythonSeed = 10 #for reproducibility purposes
np.random.seed(pythonSeed)

#NOTE: this should be run outside of any IDE. The redirecting of the output
#likes to fight with IDEs (at least, it doesn't agree well with Spyder)
for pid in range(nIndependentRuns):
    seed = np.random.randint(10**5,size=33)

    start = np.stack(nRuns*(startCoords,))
    
    uniqueCoords = (x,y)
    neckVals = np.ones(zz.shape) #expect to not stop early due to scission
    with lgv.stdout_redirected(to=str(pid)+'.out'):
        lgvObj = lgv.Fortran(pid,uniqueCoords,gamma,neckVals,zz,
                             invInert(meshArr,returnFlat=True))
    
        #A, energyOnOuterTurningLine, dt, maxTime, scissionNeckVal, savePaths, seed
        nonArrayParams = [252,pes(startCoords),10.**(-4),1.,0.5,True,seed]
        lgvObj.set_langevin_params(*nonArrayParams)
    
        allCoords, allMomenta, lastIter, finishedSuccessfully = lgvObj.run(start,True)

#%% Plotting
fig, ax = plt.subplots()
ax.contour(x,y,zz.T,colors='gray',levels=np.arange(0,100,10))

ax.scatter(sol.y[0],sol.y[1],label='Deterministic (Numerical)',zorder=100)

for pid in range(nIndependentRuns):
    with h5py.File('logs-'+str(nRuns)+'/'+str(pid).zfill(6)+'.lgv','r') as h5File:
        allCoords = np.array(h5File['allCoords'])
    
    meanTrajectory = allCoords.mean(axis=0)

    if pid == 0:
        ax.plot(*meanTrajectory.T,label='Langevin')
    else:
        ax.plot(*meanTrajectory.T)

ax.set(xlim=(0,3),ylim=(-1,5),
       xlabel='x',ylabel='y',title=r'Langevin vs Exact, $V(x)=x^2$, $M_{ij}=\delta_{ij}$')
ax.legend()

os.makedirs('plots',exist_ok=True)
fig.savefig('plots/'+str(nRuns)+'.pdf',bbox_inches='tight')
