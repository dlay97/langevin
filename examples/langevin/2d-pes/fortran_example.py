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

#multiprocessing appears to just... work
def parallel_wrapper(args):
    pid, seed = args
    with lgv.stdout_redirected(to=str(pid)+'.out'):
        lgvObj = lgv.Fortran(pid,uniqueCoords,gamma,neckVals,pesArr,invMassArr,dx=[0.5,0.5])
    
        #A, energyOnOuterTurningLine, dt, maxTime, scissionNeckVal, useSeed, seed
        nonArrayParams = [252,0.97,10.**(-4),10.,0.5,True,seed]
        lgvObj.set_langevin_params(*nonArrayParams)
        t0 = time.time()
        allCoords, allMomenta, lastIter, finishedSuccessfully = lgvObj.run(start,False)
        t1 = time.time()
        print('Run time: %.3e s'%(t1-t0))
    
    return #out, err

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

neckInterp = lgv.FastMultilinearInterp(uniqueCoords,neckVals)

# print(neckInterp(start))

inertArrays = [inertDf[m].to_numpy().reshape(shp) for m in ['M22_inv','M32_inv','M33_inv']]
invMassArr = np.stack(inertArrays,axis=-1)

gamma = np.array([[0.5,0.2],
                  [0.2,2.2]])

#taken from np.random.randint(10**5,size=33)
seed1 = np.array([29457, 28683, 35644, 20981, 62472, 55480, 78568, 24366, 26984,
                  98338, 42705, 89382, 40858, 44949, 59073, 77907, 31910, 44698,
                  89181, 73492, 90923, 15733, 24527, 89399, 13051, 83380, 69470,
                  83347, 83110, 35746, 80080, 88737, 85734])
seed2 = np.array([39257, 21411, 90327, 95806, 70788, 43235, 27203, 41369, 16239,
                  63337, 89717,    79, 46165, 88125, 63125, 50678, 80579, 45702,
                  20985, 36440, 79398, 95367, 11446,  3354, 36768, 15579, 32551,
                  34840,  9982, 68752, 18751, 87217, 90144])

args = list(zip(range(2),[seed1,seed2]))
with mp.Pool(2) as pool:
    pool.map(parallel_wrapper,args)
# print(res)
# lgvObj = lgv.Langevin(uniqueCoords,'fortran',gamma,neckVals,pesArr,invMassArr,dx=[0.5,0.5])
# #A, energyOnOuterTurningLine, dt, maxTime, scissionNeckVal, savePaths, seed
# nonArrayParams = [252,0.97,10.**(-4),50.,0.5,True,seed1]
# lgvObj.set_langevin_params(nonArrayParams)

# # sys.exit()
# allCoords, allMomenta, lastIter, finishedSuccessfully = lgvObj.run(start,True)

# #%%
# allFinalCoords = np.array([allCoords[i,lastIter[i]-1] for i in range(allCoords.shape[0])])

# fig, ax = plt.subplots()
# cf = ax.contourf(*uniqueCoords,pesArr.T,cmap='Spectral_r',levels=30)
# plt.colorbar(cf,ax=ax)

# ax.scatter(*start.T,color='white')
# ax.contour(*uniqueCoords,neckVals.T,colors=['black'],levels=[0,0.5,1])

# ax.scatter(*allFinalCoords.T,color='red',marker='x')


# print(lgvObj.invInert.shape)
