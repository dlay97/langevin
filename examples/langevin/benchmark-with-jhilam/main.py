import sys, os

lgvPath = '..//..//../src/'
if lgvPath not in sys.path:
    sys.path.insert(0,lgvPath)
import langevin as lgv

import numpy as np
import matplotlib.pyplot as plt

import h5py
import pandas as pd

import jax.numpy as jnp
from jax import jit
from jax.scipy.ndimage import map_coordinates
from functools import partial

def trimmed_matrix_inversion(df,cols,eigClipVal=10**(-4)):
    """
    Takes eigenvalue decomposition, trims eigenvalues, then rebuilds original
    matrix. Then, inverts that matrix.
    """
    #Inverse of n(n+1)/2 unique values in an n x n symmetric matrix
    nCoords = int(0.5*(np.sqrt(1+8*len(cols))-1))
    
    mat = np.zeros((len(df),nCoords,nCoords))
    idx = np.array(np.triu_indices(nCoords)).T
    
    for (cIter,c) in enumerate(cols):
        mat[(slice(None),)+tuple(idx[cIter])] = df[c]
        mat[(slice(None),)+tuple(idx[cIter][::-1])] = df[c]
    
    #Manipulations verified by plotting different components
    vals, vecs = np.linalg.eigh(mat)
    vals = vals.clip(eigClipVal)
    
    fullVals = np.zeros(mat.shape)
    for i in range(nCoords):
        fullVals[:,i,i] = vals[:,i]
    
    mat = np.linalg.inv(vecs) @ fullVals @ vecs
    matInv = np.linalg.inv(mat)
    
    return np.array([matInv[(slice(None),)+tuple(c)] for c in idx]).T

pesDf = pd.read_csv('jhilams-run/pot2030.in',sep='\s+',header=None,
                    names=['Q20','Q30','EHFB'])

uniqueCoords = [np.unique(pesDf[c]) for c in ['Q20','Q30']]
shp = tuple([len(u) for u in uniqueCoords])
zz = pesDf['EHFB'].to_numpy().reshape(shp)

fig, ax = plt.subplots()
cf = ax.contourf(*uniqueCoords,zz.T.clip(-5,30),cmap='Spectral_r',levels=30)
plt.colorbar(cf,ax=ax)

inertDf = pd.read_csv('jhilams-run/mas2030.in',sep='\s+',header=None,
                      names=['Q20','Q30','M22','M23','M33'])
invInertArr = trimmed_matrix_inversion(inertDf,['M22','M23','M33'])
invInertArr = invInertArr.reshape(shp+(-1,))

gamma = np.array([[0.5,0.2],
                  [0.2,2.2]])

initialCoordsDf = pd.read_csv('jhilams-run/prob.in',sep='\s+',header=None,
                              names=['Q20','Q30','prob'])
startCoords = initialCoordsDf[['Q20','Q30']].to_numpy()[0].reshape((1,2))

startMomenta = np.zeros((len(startCoords),2))

y0 = np.hstack((startCoords,startMomenta))

tMax = 10000

# exactSolver = lgv.DeterministicSolution(gamma,pes,pes.grad,invInert,invInert.grad)
# dtSave = 0.1 #See scipy.integrate.solve_ivp for documentation
# sol = exactSolver.solve((0,tMax),y0,t_eval=np.arange(0,tMax+dtSave,dtSave))

# print('Final exact solution:',sol.y[:2,-1])

#%% Solving with Langevin
nRuns = 1000

pythonSeed = 10 #for reproducibility purposes
np.random.seed(pythonSeed)
runSeed = np.random.randint(1,10**5)

start = np.stack(nRuns*(startCoords,)).reshape((-1,2))

neckDf = pd.read_csv('jhilams-run/scission.in',sep='\s+',header=None,
                     names=['Q20','Q30','neck','Z','A'])
neckArr = neckDf['neck'].to_numpy().reshape(shp)

lgvObj = lgv.Jax(0,uniqueCoords,gamma,neckArr,zz,invInertArr)
    
#A, energyOnOuterTurningLine, dt, maxTime, scissionNeckVal, seed
nonArrayParams = [252,0.97,10.**(-3),tMax,0.5,runSeed]
# nonArrayParams = [252,0.97,10.**(-1),tMax,0.5,runSeed]
lgvObj.set_langevin_params(*nonArrayParams)
    
res = lgvObj.run(start)

# sys.exit()
#%% Plotting

ax.contour(*uniqueCoords,zz.T,colors='white',levels=[0.97,])

ax.scatter(*startCoords.T,color='green',marker='.',zorder=100)

with h5py.File('logs/000000.lgv','r') as h5File:
    finalCoords = np.array(h5File['allCoords'])
    finishedSuccessfully = np.array(h5File['finishedSuccessfully'])
    lastIter = np.array(h5File['lastIter'])
    
print(lastIter)
print(finishedSuccessfully.sum())
pesBounds = [(u.min(),u.max()) for u in uniqueCoords]
for dimIter in range(2):
    isOobOnAxis = np.logical_or(finalCoords[:,dimIter]<pesBounds[dimIter][0],
                                finalCoords[:,dimIter]>pesBounds[dimIter][1])
    finishedSuccessfully = np.logical_and(finishedSuccessfully,~isOobOnAxis)
print(finishedSuccessfully.sum())
ax.scatter(*finalCoords[finishedSuccessfully].T,marker='.',color='black',
           zorder=100)

isInBounds = np.ones(finalCoords.shape[0],dtype=bool)
for dimIter in range(2):
    isOobOnAxis = np.logical_or(finalCoords[:,dimIter]<pesBounds[dimIter][0],
                                finalCoords[:,dimIter]>pesBounds[dimIter][1])
    isInBounds = np.logical_and(isInBounds,~isOobOnAxis)
print(isInBounds.sum())

ax.scatter(*finalCoords[isInBounds].T,marker='x',color='green',
           zorder=80)
    
ax.contour(*uniqueCoords,neckArr.T,levels=[0.5,],colors=['magenta',])
# ax.scatter(sol.y[0],sol.y[1],label='Deterministic (Numerical)',zorder=100)

# ax.scatter(*meanSols,marker='x',label='Langevin')

# ax.set(xlim=(0,3),ylim=(-1,5),
#        xlabel='x',ylabel='y',title=r'Langevin vs Exact, $V(x)=x^2$, $M_{ij}=\delta_{ij}$')
# ax.legend()

os.makedirs('plots',exist_ok=True)
fig.savefig('plots/'+str(nRuns)+'.pdf',bbox_inches='tight')

