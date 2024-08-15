import sys, os

lgvPath = '..//..//../src/'
if lgvPath not in sys.path:
    sys.path.insert(0,lgvPath)
import langevin as lgv

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import multiprocess as mp
import time

def invert_mass_tensor(df):
    #Handling non-invertible mass tensors, due to previous interpolation error or
    #non-converged DFT calculations. Can be done more cleverly
    singularMassMin = 10**(-6)
    
    massTensor = np.zeros((len(df),2,2))
    massTensor[:,0,0] = df['M22'].clip(singularMassMin)
    massTensor[:,1,0] = df['M32']
    massTensor[:,0,1] = df['M32']
    massTensor[:,1,1] = df['M33'].clip(singularMassMin)
    
    singularInds = np.where(np.linalg.det(massTensor) < singularMassMin)
    
    massTensor[singularInds,0,1] = 0
    massTensor[singularInds,1,0] = 0
    
    det = massTensor[:,0,0]*massTensor[:,1,1] - massTensor[:,0,1]**2
    
    invDf = df.copy()
    invDf['M22_inv'] = massTensor[:,1,1]/det
    invDf['M32_inv'] = -massTensor[:,0,1]/det
    invDf['M33_inv'] = massTensor[:,0,0]/det
    
    return invDf

def fortran_wrapper(maxTime):
    #taken from np.random.randint(10**5,size=33)
    seed1 = np.array([29457, 28683, 35644, 20981, 62472, 55480, 78568, 24366, 26984,
                      98338, 42705, 89382, 40858, 44949, 59073, 77907, 31910, 44698,
                      89181, 73492, 90923, 15733, 24527, 89399, 13051, 83380, 69470,
                      83347, 83110, 35746, 80080, 88737, 85734])

    lgvObj = lgv.Fortran(0,uniqueCoords,gamma,neckVals,pesArr,invMassArr,dx=[0.5,0.5])

    #A, energyOnOuterTurningLine, dt, maxTime, scissionNeckVal, useSeed, seed
    nonArrayParams = [252,0.97,10.**(-4),maxTime,0.5,True,seed1]
    lgvObj.set_langevin_params(*nonArrayParams)

    allCoords, allMomenta, lastIter, finishedSuccessfully = lgvObj.run(start,False)
    
    return 

def python_wrapper(maxTime):
    lgvObj = lgv.Python(uniqueCoords,gamma,neckVals,pesArr,invMassArr,dx=[0.5,0.5])
    
    lgvObj.run(start,252,0.97,10.**(-4),maxTime,0.5,True,5,
               0)
    return

pesDf = pd.read_csv('pot2030.in',header=None,sep='\s+',names=['Q20','Q30','EHFB'])
inertDf = pd.read_csv('mas2030.in',header=None,sep='\s+',names=['Q20','Q30','M22','M32','M33'])
neckDf = pd.read_csv('scission.in',header=None,sep='\s+',names=['Q20','Q30','qn','Z','A'])
initPtsDf = pd.read_csv('prob.in',header=None,sep='\s+',names=['Q20','Q30','prob'])

inertDf = invert_mass_tensor(inertDf)

uniqueCoords = [np.unique(pesDf[col]) for col in ['Q20','Q30']]
shp = tuple([len(u) for u in uniqueCoords])

pesArr = pesDf['EHFB'].to_numpy().reshape(shp)

nEventsPerCoord = 10
start = initPtsDf[['Q20','Q30']].to_numpy()
start = np.tile(start.T,nEventsPerCoord).T

neckVals = neckDf['qn'].to_numpy().reshape(shp)

inertArrays = [inertDf[m].to_numpy().reshape(shp) for m in ['M22_inv','M32_inv','M33_inv']]
invMassArr = np.stack(inertArrays,axis=-1)

gamma = np.array([[0.5,0.2],
                  [0.2,2.2]])

maxTimes = [0.01,0.1,1,10]
maxTests = 10

fortranTimes = []
for maxTime in maxTimes:
    runtimeArr = np.zeros(maxTests)
    
    for i in range(maxTests):
        t0 = time.time()
        fortran_wrapper(maxTime)
        t1 = time.time()
        
        runtimeArr[i] = t1 - t0
        
    fortranTimes.append(runtimeArr)
    
# print('Starting Python version')
# pythonTimes = []
# for maxTime in maxTimes:
#     print(maxTime)
#     runtimeArr = np.zeros(maxTests)
    
#     for i in range(maxTests):
#         print(i)
#         t0 = time.time()
#         python_wrapper(maxTime)
#         t1 = time.time()
        
#         runtimeArr[i] = t1 - t0
        
#     pythonTimes.append(runtimeArr)

fortranTimes = np.array(fortranTimes)
np.savetxt('fortranTimes.dat',fortranTimes)

# pythonTimes = np.array(pythonTimes)
# np.savetxt('pythonTimes.dat',pythonTimes)

#%%
fig, ax = plt.subplots()

maxTimes = np.array(maxTimes)
nIters = maxTimes/(10**(-4))

keys = ['Fortran','Python','JAX-1','JAX-2']
for (fIter,f) in enumerate(['fortranTimes.dat','pythonTimes.dat','jax-inner.dat','jax-loop.dat']):
    arr = np.loadtxt(f)
    mean = np.mean(arr,axis=1)
    std = np.std(arr,axis=1)

    ax.errorbar(nIters[:len(mean)],mean,yerr=std,marker='.',capsize=5,label=keys[fIter])

ax.set(xscale='log',yscale='log',
       xlabel='Number of Iterations',ylabel='Run Time (s)',
       title='Langevin Benchmark (Preliminary)')
ax.legend()

fig.savefig('runtime.pdf',bbox_inches='tight')

