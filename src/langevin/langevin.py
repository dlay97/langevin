import langevin_module as lgvf #for langevin fortran
import numpy as np
import sys, os
import h5py

from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate._rgi_cython import evaluate_linear_2d

from scipy.linalg import sqrtm

from scipy.integrate import solve_ivp

from contextlib import contextmanager

import warnings
import time

from jax.scipy.ndimage import map_coordinates
from jax import jit, random
import jax.numpy as jnp
import jax.lax as lax

from functools import partial

#Warns when passing arrays to Fortran, which seems an unnecessary warning
warnings.filterwarnings('ignore',message='.*was not an f_contiguous NumPy array.*')

#Warns if trying to store paths larger than maxArrSizeWarn GB; halts if
#trying to store paths larger than maxArrSizeErr GB
maxArrSizeWarn = 1.
maxArrSizeErr = 3.

@contextmanager
def stdout_redirected(to=os.devnull):
    '''
    Redirects *only* stdout, but also from the Fortran code
    From https://stackoverflow.com/a/17954769
    TODO: modify to include stderr
    
    import os

    with stdout_redirected(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    '''
    fd = sys.stdout.fileno()

    ##### assert that Python and C stdio write using the same file descriptor
    ####assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stdout")) == fd == 1

    def _redirect_stdout(to):
        sys.stdout.close() # + implicit flush()
        os.dup2(to.fileno(), fd) # fd writes to 'to' file
        sys.stdout = os.fdopen(fd, 'w') # Python writes to fd

    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        with open(to, 'w') as file:
            _redirect_stdout(to=file)
        try:
            yield # allow code to be run with the redirected stdout
        finally:
            _redirect_stdout(to=old_stdout) # restore stdout.
                                            # buffering and flags such as
                                            # CLOEXEC may be different

class FastMultilinearInterp(RegularGridInterpolator):
    #The implemented interpolator in scipy has a lot of overhead, just checking
    #for out-of-bounds errors everytime it's evaluated. That is totally unnecessary
    #in our case, and removing that checking saves ~30% of the run time. In hindsight,
    #this is some of the overhead that will continue to vanish as we increase
    #the number of points being evaluated - while it will exist in each iteration,
    #as a % of the work done, it will decrease
    def _prepare_xi(self, xi):
        ndim = len(self.grid)
        
        #This reshaping has a negligible impact on the run time of the code,
        #and leaving it in makes this a general-purpose interpolator
        xi_shape = xi.shape
        xi = xi.reshape(-1, xi_shape[-1])
        xi = np.asarray(xi, dtype=float)

        # find nans in input
        # nans = np.any(np.isnan(xi), axis=-1)
        nans = False

        out_of_bounds = None

        return xi, xi_shape, ndim, nans, out_of_bounds
    
    def __call__(self, xi, method=None):
        """
        Ripped straight from the source, for a further 16% speedup over the
        modifications to _prepare_xi
        """
        method = self.method

        xi, xi_shape, ndim, nans, out_of_bounds = self._prepare_xi(xi)
        
        #This overhead is probably necessary for when we look at the 3D case
        if method == "linear":
            indices, norm_distances = self._find_indices(xi.T)
            if (ndim == 2 and hasattr(self.values, 'dtype') and
                    self.values.ndim == 2 and self.values.flags.writeable and
                    self.values.dtype in (np.float64, np.complex128) and
                    self.values.dtype.byteorder == '='):
                # until cython supports const fused types, the fast path
                # cannot support non-writeable values
                # a fast path
                out = np.empty(indices.shape[1], dtype=self.values.dtype)
                result = evaluate_linear_2d(self.values,
                                            indices,
                                            norm_distances,
                                            self.grid,
                                            out)
            else:
                # print('asdf')
                result = self._evaluate_linear(indices, norm_distances)
        elif method == "nearest":
            indices, norm_distances = self._find_indices(xi.T)
            result = self._evaluate_nearest(indices, norm_distances)
        elif method in self._SPLINE_METHODS:
            result = self._evaluate_spline(xi, method)

        return result.reshape(xi_shape[:-1] + self.values.shape[ndim:])


def interp_array_to_different_grid(uniqueCoords,arrOrFunc,dx,returnCoords=False):
    nDims = len(uniqueCoords)
    if len(dx) != nDims:
        raise ValueError('len(dx) != len(uniqueCoords)')
        
    newCoords = []
    
    for dIter in range(nDims):
        diff = np.diff(uniqueCoords[dIter])
        if len(np.unique(diff)) != 1:
            raise ValueError('Initial array should be uniformly space')
        bounds = [uniqueCoords[dIter][0],uniqueCoords[dIter][-1]]
        newCoords.append(np.arange(bounds[0],bounds[1]+dx[dIter],dx[dIter]))
        
    meshArr = np.swapaxes(np.array(np.meshgrid(*newCoords)),0,-1)
    
    if isinstance(arrOrFunc,np.ndarray):
        interp = FastMultilinearInterp(uniqueCoords,arrOrFunc)
        ret = interp(meshArr)
    else:
        ret = arrOrFunc(meshArr)
    
    if returnCoords:
        return ret, newCoords
    else:
        return ret

def setup_grad_and_hess(uniqueCoords,pes,invInert,dx=None):
    nDims = len(uniqueCoords)
    triuInds = np.triu_indices(nDims)
    if dx is None:
        dx = np.array([u[1]-u[0] for u in uniqueCoords])
    assert len(dx) == nDims
    
    slc = nDims*(slice(None),)
    
    pesGrad = np.stack(np.gradient(pes,*dx),axis=-1)
    
    #PES Hessian
    hessList = []
    for dIter in range(nDims):
        hessList.append(np.stack(np.gradient(pesGrad[slc+(dIter,)],*dx,),
                                             axis=-1))
    hessArr = np.stack(hessList,axis=-1)
    pesHess = hessArr[slc+triuInds]
    
    #Interpolating to new grid
    pes, newCoords = interp_array_to_different_grid(uniqueCoords,pes,dx,
                                                    returnCoords=True)
    pesGrad = interp_array_to_different_grid(uniqueCoords,pesGrad,dx)
    pesHess = interp_array_to_different_grid(uniqueCoords,pesHess,dx)
    
    #probably some stuff in here that can be tidied up
    #Inverse inertia gradient
    gradList = []
    for cIter in range(invInert.shape[-1]):
        gradList.append(np.stack(np.gradient(invInert[slc+(cIter,)],*dx,),
                                 axis=-1))
    gradArr = np.stack(gradList,axis=-2)
    invInertGrad = gradArr
    
    #Inverse inertia Hessian    
    invInertHess = np.zeros(invInert.shape+(nDims*(nDims+1)//2,))
    for compIter in range(invInert.shape[-1]):
        hessList = []
        for dIter in range(nDims):
            hessList.append(np.stack(np.gradient(invInertGrad[slc+(compIter,dIter)],
                                                 *dx,
                                                 ),axis=-1))
        hessArr = np.stack(hessList,axis=-1)
        invInertHess[slc+(compIter,)] = hessArr[slc+triuInds]
        
    #Interpolating to new grid
    invInert = interp_array_to_different_grid(uniqueCoords,invInert,dx)
    invInertGrad = interp_array_to_different_grid(uniqueCoords,invInertGrad,dx)
    invInertHess = interp_array_to_different_grid(uniqueCoords,invInertHess,dx)
    
    return newCoords, pes, pesGrad, pesHess, invInert, invInertGrad, invInertHess

class LangevinFileIO:
    def __init__(self,outputDir,fileID):
        os.makedirs(outputDir,exist_ok=True)
        
        self.fileName = os.path.join(outputDir,str(fileID).zfill(6)+'.lgv')
        
    def write_params(self,A,energyOnOuterTurningLine,dt,maxTime,
                     scissionNeckVal,useSeed,seed,method):
        with h5py.File(self.fileName,'w') as h5File:
            h5File.attrs.create('method',method)
            h5File.attrs.create('A',A)
            h5File.attrs.create('energyOnOuterTurningLine',energyOnOuterTurningLine)
            h5File.attrs.create('dt',dt)
            h5File.attrs.create('maxTime',maxTime)
            h5File.attrs.create('scissionNeckVal',scissionNeckVal)
            h5File.attrs.create('useSeed',useSeed)
            h5File.attrs.create('seed',seed)
                
        return None
    
    def write_datasets(self,allCoords,allMomenta,lastIter,finishedSuccessfully,runTime):
        with h5py.File(self.fileName,'a') as h5File:
            h5File.create_dataset('allCoords',data=allCoords)
            h5File.create_dataset('allMomenta',data=allMomenta)
            h5File.create_dataset('lastIter',data=lastIter)
            h5File.create_dataset('finishedSuccessfully',data=finishedSuccessfully)
            h5File.attrs.create('runTime',runTime)
        return None

class Fortran:
    def __init__(self,fileID,uniqueCoords,fric,neckVals,pes,invInert,
                 dx=None,outputDir='logs'):
        self.nDims = len(uniqueCoords)
        if self.nDims != 2:
            raise NotImplementedError
        
        self.fric = fric
        
        self.dx = dx
        if self.dx is None:
            self.dx = np.array([u[1]-u[0] for u in uniqueCoords])
            
        self.uniqueCoords, self.pes, self.pesGrad, self.pesHess, \
            self.invInert, self.invInertGrad, self.invInertHess = \
                setup_grad_and_hess(uniqueCoords,pes,invInert,self.dx)
            
        self.neckVals = interp_array_to_different_grid(uniqueCoords,neckVals,self.dx)
        
        self.ioObj = LangevinFileIO(outputDir,fileID)
        
        lgvf.langevin_2d.setup_arrays(*self.uniqueCoords,self.pes,
                                      self.pesGrad,self.pesHess,
                                      self.invInert,self.invInertGrad,
                                      self.invInertHess,
                                      self.fric,sqrtm(self.fric),
                                      self.neckVals
                                      )
    
    def set_langevin_params(self,A,energyOnOuterTurningLine,dt,maxTime,
                            scissionNeckVal,useSeed,seed):
        lgvf.langevin_2d.setup_params(A,energyOnOuterTurningLine,dt,maxTime,
                                      scissionNeckVal,useSeed,seed)
            
        #Approximately in GB
        self.approxArraySizePerStart = maxTime/dt*4*8*10**(-9)
        
        self.ioObj.write_params(A,energyOnOuterTurningLine,dt,maxTime,
                                scissionNeckVal,useSeed,seed,'fortran')
        return None
    
    def run(self,initialPoints,savePaths):
        if savePaths:
            size = self.approxArraySizePerStart*len(initialPoints)
            if size > maxArrSizeWarn:
                print('Total array size: approx %.3f GB'%size,
                      flush=True)
                warnings.warn('Attempting to allocate more than %.2f GB of memory. Consider disabling "savePaths"'%maxArrSizeWarn)
            if size > maxArrSizeErr:
                raise ValueError('Current configuration attempts to allocate more than %.2f GB of memory. Disable "savePaths" to run Langevin'%maxArrSizeErr)
            
        t0 = time.time()
        allCoords, allMomenta, lastIter, finishedSuccessfully = \
            lgvf.langevin_2d.run(initialPoints,savePaths)
        t1 = time.time()
        
        #Subtract 1 on lastIter b/c Fortran
        self.ioObj.write_datasets(allCoords,allMomenta,lastIter-1,finishedSuccessfully,
                                  t1-t0)
        
        return allCoords, allMomenta, lastIter, finishedSuccessfully
    
class Python:
    def __init__(self,uniqueCoords,fric,neckVals,pes,invInert,
                 dx=None,outputDir='logs'):
        self.nDims = len(uniqueCoords)
        
        self.fric = fric
        
        self.dx = dx
        if self.dx is None:
            self.dx = np.array([u[1]-u[0] for u in uniqueCoords])
            
        neckVals = interp_array_to_different_grid(uniqueCoords,neckVals,dx)
            
        self.uniqueCoords, self.pes, self.pesGrad, self.pesHess, \
            self.invInert, self.invInertGrad, self.invInertHess = \
                setup_grad_and_hess(uniqueCoords,pes,invInert,self.dx)
        
        self.neck_interpolator = FastMultilinearInterp(self.uniqueCoords,neckVals)
        
        self.make_pes_interpolators()
        self.make_inertia_interpolators()
        
        self.fric = fric
        self.randForce = sqrtm(self.fric)
            
        self.outputDir = outputDir
    
    def make_pes_interpolators(self):
        arrInterp = FastMultilinearInterp(self.uniqueCoords,self.pes)
        gradInterps = [FastMultilinearInterp(self.uniqueCoords,g.T) for g in self.pesGrad.T]
        hessInterps = [FastMultilinearInterp(self.uniqueCoords,h.T) for h in self.pesHess.T]
        
        def gradient_interp(points):
            pointsShape = points.shape[:-1]
            
            arrOut = np.zeros(pointsShape+(self.nDims,))
            for i in range(self.nDims):
                arrOut[len(pointsShape)*(slice(None),)+(i,)] = gradInterps[i](points)
                
            return arrOut
        
        triuInds = np.array(np.triu_indices(self.nDims))
        
        def hessian_interp(points):
            pointsShape = points.shape[:-1]
            
            arrOut = np.zeros(pointsShape+(self.nDims,self.nDims))
            for (i,hess) in enumerate(hessInterps):
                hessEval = hess(points)
                #The upper triangular indices
                arrOut[len(pointsShape)*(slice(None),) + tuple(triuInds[:,i])] = hessEval
                #And, the lower ones (that's what the [::-1] does)
                arrOut[len(pointsShape)*(slice(None),) + tuple(triuInds[:,i][::-1])] = hessEval
            
            return arrOut
        
        self.pesFuncs = [arrInterp, gradient_interp, hessian_interp]
        
        return
    
    def make_inertia_interpolators(self):
        #grad should be [gradsOfM11, gradsOfM12, ...]
        #hess should be [hessOfM11, hessOfM12, ...] as the output of hessian defined above
        
        tensorInterps = [FastMultilinearInterp(self.uniqueCoords,f) for f in np.moveaxis(self.invInert,-1,0)]
        
        gradInterps = []
        for i in range(self.nDims*(self.nDims+1)//2):
            componentGradInterps = []
            for j in range(self.nDims):
                componentGradInterps.append(FastMultilinearInterp(self.uniqueCoords,
                                                                  self.invInertGrad[self.nDims*(slice(None),)+(i,j)]))
            gradInterps.append(componentGradInterps)
            
        hessInterps = []
        for i in range(self.nDims*(self.nDims+1)//2):
            componentHessInterps = []
            for j in range(self.nDims*(self.nDims+1)//2):
                componentHessInterps.append(FastMultilinearInterp(self.uniqueCoords,
                                                                  self.invInertHess[self.nDims*(slice(None),)+(i,j)]))
            hessInterps.append(componentHessInterps)
            
        triuInds = np.array(np.triu_indices(self.nDims))
            
        def field_interp(points):
            pointsShape = points.shape[:-1]
            
            arrOut = np.zeros(pointsShape+(self.nDims,self.nDims))
            for (i,interp) in enumerate(tensorInterps):
                tensorArr = interp(points)
                
                #The upper triangular indices
                arrOut[len(pointsShape)*(slice(None),) + tuple(triuInds[:,i])] = tensorArr
                #And, the lower ones (that's what the [::-1] does)
                arrOut[len(pointsShape)*(slice(None),) + tuple(triuInds[:,i][::-1])] = tensorArr
            return arrOut
        
        def gradient_interp(points):
            pointsShape = points.shape[:-1]
            
            #pointsShape, then tensor, then gradient
            arrOut = np.zeros(pointsShape+(self.nDims,self.nDims)+(self.nDims,))
            for (i,componentGradInterps) in enumerate(gradInterps):
                for (j,interp) in enumerate(componentGradInterps):
                    evalArr = interp(points)
                    
                    #The upper triangular indices
                    arrOut[len(pointsShape)*(slice(None),) + tuple(triuInds[:,i])+(j,)] = evalArr
                    #And, the lower ones (that's what the [::-1] does)
                    arrOut[len(pointsShape)*(slice(None),) + tuple(triuInds[:,i][::-1])+(j,)] = evalArr
            
            return arrOut
        
        #For use in hessian_interp
        symmetryIndsDict = {}
        for (i,componentHessInterps) in enumerate(hessInterps):
            for (j,_) in enumerate(componentHessInterps):
                #Upper triangular in tensor and Hessian
                idx1 = tuple(triuInds[:,i])+tuple(triuInds[:,j])
                #Lower in tensor, upper in Hessian
                idx2 = tuple(triuInds[:,i][::-1])+tuple(triuInds[:,j])
                #Upper in tensor, lower in Hessian
                idx3 = tuple(triuInds[:,i])+tuple(triuInds[:,j][::-1])
                #Lower in both
                idx4 = tuple(triuInds[:,i][::-1])+tuple(triuInds[:,j][::-1])
                
                symmetryIndsDict[(i,j)] = [idx1,idx2,idx3,idx4]
        
        # @profile
        def hessian_interp(points):
            pointsShape = points.shape[:-1]
            
            #pointsShape, then tensor, then Hessian
            arrOut = np.zeros(pointsShape+(self.nDims,self.nDims)+(self.nDims,self.nDims))
            for (i,componentHessInterps) in enumerate(hessInterps):
                for (j,interp) in enumerate(componentHessInterps):
                    evalArr = interp(points)
                    
                    for idx in symmetryIndsDict[(i,j)]:
                        arrOut[len(pointsShape)*(slice(None),)+idx] = evalArr
            
            return arrOut
        
        self.invInertFuncs = [field_interp, gradient_interp, hessian_interp]
        
        return 
    
    def make_stopper(self,neckValue):
        pesBounds = [(u.min(),u.max()) for u in self.uniqueCoords]
        
        def stopper(coords):
            needsUpdate = np.ones(coords.shape[0],dtype=bool)
            
            #Out-of-bounds checker
            for dimIter in range(self.nDims):
                isOobOnAxis = np.logical_or(coords[:,dimIter]<=pesBounds[dimIter][0],
                                            coords[:,dimIter]>=pesBounds[dimIter][1])
                needsUpdate = np.logical_and(needsUpdate,~isOobOnAxis)
            
            #If the neck is below neckValue, is done. Also, have to account for
            #out-of-bounds points here, as well
            neckVals = self.neck_interpolator(coords[needsUpdate])
            neckValIsBelowMin = (neckVals<=neckValue)
            needsUpdate[needsUpdate] = np.logical_and(needsUpdate[needsUpdate],~neckValIsBelowMin)
            
            return needsUpdate
        
        return stopper
    
    def run(self,initialPoints,A,E0,dt,maxTime,scissionNeckVal,useSeed,randomSeed,
            fileID):
        #gArr -> self.randForce
        #gammaArr -> self.fric
        
        levelDens = A/10.
        
        ioObj = LangevinFileIO(self.outputDir,fileID)
        ioObj.write_params(A,E0,dt,maxTime,scissionNeckVal,useSeed,randomSeed,'python')
        stopper = self.make_stopper(scissionNeckVal)
        
        #Assume here that g is constant in space
        #Probably dt can be increased - at least 1 trajectory on 10**(-3) seemed fine
        #Tried using np.einsum_path to optimize tensor contractions, but that was actually
        #slower for some reason
        #Using einsumt for threaded summations doesn't make sense (and is actually slower).
        #It makes more sense to spawn multiple single-threaded langevin runs, rather
        #than spawning threads inside the langevin code
        
        if useSeed:
            np.random.seed(randomSeed)
        #for parallel output: https://stackoverflow.com/questions/1501651/log-output-of-multiprocessing-process
        maxIters = int(maxTime/dt)
        
        nRuns, nDims = initialPoints.shape
        needsUpdateMask = np.ones(nRuns,dtype=bool)
        
        coords = initialPoints.copy()
        momentum = np.zeros(coords.shape)
        
        lastIter = np.zeros(nRuns,dtype=int)
        
        t0 = time.time()
        for i in range(maxIters):
            needsUpdateMask = stopper(coords)
            
            if needsUpdateMask.sum() == 0:
                break
            
            coordsToUse = coords[needsUpdateMask]
            momentumToUse = momentum[needsUpdateMask]
            lastIter[needsUpdateMask] += 1
            
            omega = np.random.normal(scale=np.sqrt(2),size=(coordsToUse.shape+(2,)))
            # omega = np.ones(coordsToUse.shape+(2,))
            Gamma1 = dt**0.5 * omega[:,:,0]
            Gamma2 = dt**1.5 * (0.5*omega[:,0] + 1/(2*np.sqrt(3))*omega[:,:,1])
            #Gamma3 I think only is used when considering coordinate-dependent gArr
            Gamma3 = dt**1.5 * (0.5*omega[:,0] - 1/(2*np.sqrt(3))*omega[:,:,1])
            
            pesEval = self.pesFuncs[0](coordsToUse)
            pesGrad = self.pesFuncs[1](coordsToUse)
            pesHess = self.pesFuncs[2](coordsToUse)
            
            mass = self.invInertFuncs[0](coordsToUse)
            massGrad = self.invInertFuncs[1](coordsToUse)
            massHess = self.invInertFuncs[2](coordsToUse)
            
            exint = E0 - pesEval
            exint = exint.clip(0)
            temperature = np.sqrt(exint/levelDens)
            
            g = np.stack(coordsToUse.shape[0]*(self.randForce,))
            g = g*np.sqrt(temperature)[:,None,None]
            
            v = np.einsum("aij,aj->ai",mass,momentumToUse)
            h = -pesGrad -0.5 * np.einsum("ajki,aj,ak->ai",massGrad,momentumToUse,momentumToUse) - \
                np.einsum("ij,aj->ai",self.fric,v)
            
            dvdq = np.einsum("aikj,ak->aij",massGrad,momentumToUse)
            dvdp = mass
            
            dhdq = -pesHess - 0.5 * np.einsum("ak,aklij,al->aij",momentumToUse,massHess,momentumToUse) \
                - np.einsum("ik,akj->aij",self.fric,dvdq)
            dhdp =  -np.einsum("ajki,ak->aij",massGrad,momentumToUse) \
                - np.einsum("ik,akj->aij",self.fric,dvdp)
            
            momentumChange = dt*h + 0.5*(np.einsum("aij,aj->ai",dhdq,v) + np.einsum("aij,aj->ai",dhdp,h))*dt**2 + \
                np.einsum("aij,aj->ai",g,Gamma1) + np.einsum("aij,ajk,ak->ai",dhdp,g,Gamma2)
            positionChange = dt*v + 0.5*(np.einsum("aij,aj->ai",dvdq,v) + np.einsum("aij,aj->ai",dvdp,h))*dt**2 + \
                np.einsum("aij,ajk,ak->ai",dvdp,g,Gamma2)
                    
            coords[needsUpdateMask] += positionChange
            momentum[needsUpdateMask] += momentumChange
            
        t1 = time.time()
        
        ioObj.write_datasets(coords,momentum,lastIter,needsUpdateMask,t1-t0)
        
        return coords
    
class Jax:
    def __init__(self,fileID,uniqueCoords,fric,neckVals,pes,invInert,
                 dx=None,outputDir='logs'):
        self.nDims = len(uniqueCoords)
        if self.nDims != 2:
            raise NotImplementedError
        
        self.fric = fric
        self.randForce = sqrtm(self.fric)
        
        self.dx = dx
        if self.dx is None:
            self.dx = np.array([u[1]-u[0] for u in uniqueCoords])
            
        self.uniqueCoords, self.pes, self.pesGrad, self.pesHess, \
            self.invInert, self.invInertGrad, self.invInertHess = \
                setup_grad_and_hess(uniqueCoords,pes,invInert,self.dx)
            
        self.xIntercept = np.array([u[0] for u in self.uniqueCoords])
            
        self.neckVals = interp_array_to_different_grid(uniqueCoords,neckVals,self.dx)
        
        #Setting up interpolators
        self.pes_interp, self.pes_grad, self.pes_hess = \
            self._make_pes_interpolators()
            
        self.inert_interp, self.inert_grad, self.inert_hess = \
            self._make_inverse_inertia_interpolators()
        
        self.stop_func = self._make_stopper()
        
        self.ioObj = LangevinFileIO(outputDir,fileID)
        
    def set_langevin_params(self,A,energyOnOuterTurningLine,dt,maxTime,
                            scissionNeckVal):
        #Approximately in GB
        self.approxArraySizePerStart = maxTime/dt*4*8*10**(-9)
        
        self.A = A
        self.E0 = energyOnOuterTurningLine
        self.dt = dt
        self.maxTime = maxTime
        self.scissionNeckVal = scissionNeckVal
        # self.ioObj.write_params(A,energyOnOuterTurningLine,dt,maxTime,
        #                         scissionNeckVal,useSeed,seed,'fortran')
        return None
        
    @partial(jit,static_argnames=('self',))
    def rescale_points(self,points):
        return jnp.moveaxis((points - self.xIntercept)/self.dx,-1,0)
    
    def _make_pes_interpolators(self):
        slc = self.nDims*(slice(None),)
        
        @jit
        def arr_interp(points):
            return map_coordinates(self.pes,self.rescale_points(points),1,mode='nearest')
        
        @jit
        def gradient_interp(points):
            pointsShape = points.shape[:-1]
            
            arrOut = jnp.zeros(pointsShape+(self.nDims,))
            for i in range(self.nDims):
                evalArr = map_coordinates(self.pesGrad[slc+(i,)],
                                          self.rescale_points(points),1,mode='nearest')
                arrOut = arrOut.at[len(pointsShape)*(slice(None),)+(i,)].set(evalArr)
                
            return arrOut
        
        triuInds = np.array(np.triu_indices(self.nDims))
        
        @jit
        def hessian_interp(points):
            pointsShape = points.shape[:-1]
            
            arrOut = jnp.zeros(pointsShape+(self.nDims,self.nDims))
            for i in range(self.pesHess.shape[-1]):
                hessEval = map_coordinates(self.pesHess[slc+(i,)],
                                           self.rescale_points(points),1,mode='nearest')
                #The upper triangular indices
                arrOut = arrOut.at[len(pointsShape)*(slice(None),) + tuple(triuInds[:,i])].set(hessEval)
                #And, the lower ones (that's what the [::-1] does)
                arrOut = arrOut.at[len(pointsShape)*(slice(None),) + tuple(triuInds[:,i][::-1])].set(hessEval)
            
            return arrOut
        
        return arr_interp, gradient_interp, hessian_interp
    
    def _make_inverse_inertia_interpolators(self):
        triuInds = np.array(np.triu_indices(self.nDims))
        slc = self.nDims*(slice(None),)
        nInertComps = self.invInert.shape[-1]
        
        @jit
        def field_interp(points):
            pointsShape = points.shape[:-1]
            pointsSlice = len(pointsShape)*(slice(None),)
            
            arrOut = jnp.zeros(pointsShape+(self.nDims,self.nDims))
            
            for i in range(nInertComps):
                tensorArr = map_coordinates(self.invInert[slc+(i,)],self.rescale_points(points),1,
                                            mode='nearest')
                #The upper triangular indices
                arrOut = arrOut.at[pointsSlice + tuple(triuInds[:,i])].set(tensorArr)
                                                                                                 
                #And, the lower ones (that's what the [::-1] does)
                arrOut = arrOut.at[pointsSlice + tuple(triuInds[:,i][::-1])].set(tensorArr)
                
            return arrOut
        
        @jit
        def gradient_interp(points):
            pointsShape = points.shape[:-1]
            pointsSlice = len(pointsShape)*(slice(None),)
            
            #pointsShape, then tensor, then gradient
            arrOut = jnp.zeros(pointsShape+(self.nDims,self.nDims)+(self.nDims,))
            
            for i in range(nInertComps):
                for j in range(self.nDims):
                    evalArr = map_coordinates(self.invInertGrad[slc+(i,j)],
                                              self.rescale_points(points),1,
                                              mode='nearest')
                    
                    #The upper triangular indices
                    arrOut = arrOut.at[pointsSlice + tuple(triuInds[:,i])+(j,)].set(evalArr)
                    
                    #And, the lower ones (that's what the [::-1] does)
                    arrOut = arrOut.at[pointsSlice + tuple(triuInds[:,i][::-1])+(j,)].set(evalArr)
            
            return arrOut
        
        #For use in hessian_interp
        symmetryIndsDict = {}
        for i in range(nInertComps):
            for j in range(nInertComps):
                #Upper triangular in tensor and Hessian
                idx1 = tuple(triuInds[:,i])+tuple(triuInds[:,j])
                #Lower in tensor, upper in Hessian
                idx2 = tuple(triuInds[:,i][::-1])+tuple(triuInds[:,j])
                #Upper in tensor, lower in Hessian
                idx3 = tuple(triuInds[:,i])+tuple(triuInds[:,j][::-1])
                #Lower in both
                idx4 = tuple(triuInds[:,i][::-1])+tuple(triuInds[:,j][::-1])
                
                symmetryIndsDict[(i,j)] = [idx1,idx2,idx3,idx4]
        
        @jit
        def hessian_interp(points):
            pointsShape = points.shape[:-1]
            
            #pointsShape, then tensor, then Hessian
            arrOut = jnp.zeros(pointsShape+(self.nDims,self.nDims)+(self.nDims,self.nDims))
            for i in range(nInertComps):
                for j in range(nInertComps):
                    evalArr = map_coordinates(self.invInertHess[slc+(i,j)],
                                              self.rescale_points(points),1,mode='nearest')
                    
                    for idx in symmetryIndsDict[(i,j)]:
                        arrOut = arrOut.at[len(pointsShape)*(slice(None),)+idx].set(evalArr)
            
            return arrOut
            
        return field_interp, gradient_interp, hessian_interp
    
    def _make_stopper(self):
        pesBounds = [(u.min(),u.max()) for u in self.uniqueCoords]
        
        @jit
        def check_if_run_should_stop(points):
            needsUpdate = jnp.ones(points.shape[0],dtype=bool)
            
            #Out-of-bounds checker
            for dimIter in range(self.nDims):
                isOobOnAxis = jnp.logical_or(points[:,dimIter]<=pesBounds[dimIter][0],
                                             points[:,dimIter]>=pesBounds[dimIter][1])
                needsUpdate = jnp.logical_and(needsUpdate,~isOobOnAxis)
            #If the neck is below neckValue, is done. Also, have to account for
            #out-of-bounds points here, as well
            neckVals = map_coordinates(self.neckVals,self.rescale_points(points),1,
                                       mode='nearest')
            
            neckValIsBelowMin = (neckVals<=self.scissionNeckVal)
            needsUpdate = jnp.logical_and(needsUpdate,~neckValIsBelowMin)
            
            return needsUpdate
        
        return check_if_run_should_stop
    
    def make_loop_stopper(self):
        maxIters = self.maxTime//self.dt
        
        @jit
        def stop(args):
            _, _, _, stepIter = args
            return jnp.where(stepIter>maxIters,False,True)
        return stop
    
    def _make_single_step(self):
        #Define as a wrapper b/c lax.while_loop expects a single argument,
        #so I don't know how to separate 'self' out of the list of arguments
        levelDensity = self.A/10.
        
        
        @jit
        def step(args):
            coords, momentum, key, stepIter = args
            needsUpdateMask = self.stop_func(coords)
            
            key, subkey = random.split(key)
            #Sampling from a normal with $\sigma=\sqrt{2}$
            omega = 2**(0.5)*random.normal(subkey,coords.shape+(2,))
            
            Gamma1 = self.dt**0.5 * omega[:,:,0]
            Gamma2 = self.dt**1.5 * (0.5*omega[:,0] + 1/(2*jnp.sqrt(3))*omega[:,:,1])
            #Gamma3 I think only is used when considering coordinate-dependent gArr
            Gamma3 = self.dt**1.5 * (0.5*omega[:,0] - 1/(2*jnp.sqrt(3))*omega[:,:,1])
            
            pesEval = self.pes_interp(coords)
            pesGradEval = self.pes_grad(coords)
            pesHessEval = self.pes_hess(coords)
            
            massEval = self.inert_interp(coords)
            massGradEval = self.inert_grad(coords)
            massHessEval = self.inert_hess(coords)
            
            exint = self.E0 - pesEval
            exint = exint.clip(0)
            temperature = jnp.sqrt(exint/levelDensity)
            
            g = jnp.stack(coords.shape[0]*(self.randForce,))
            g = g*jnp.sqrt(temperature)[:,None,None]
            
            v = jnp.einsum("aij,aj->ai",massEval,momentum)
            h = -pesGradEval -0.5 * jnp.einsum("ajki,aj,ak->ai",massGradEval,momentum,momentum) - \
                jnp.einsum("ij,aj->ai",self.fric,v)
            
            dvdq = jnp.einsum("aikj,ak->aij",massGradEval,momentum)
            dvdp = massEval
            
            dhdq = -pesHessEval - 0.5 * jnp.einsum("ak,aklij,al->aij",momentum,massHessEval,momentum) - \
                jnp.einsum("ik,akj->aij",self.fric,dvdq)
            dhdp = -jnp.einsum("ajki,ak->aij",massGradEval,momentum) - jnp.einsum("ik,akj->aij",self.fric,dvdp)
            
            momentumChange = self.dt*h + 0.5*(jnp.einsum("aij,aj->ai",dhdq,v) + jnp.einsum("aij,aj->ai",dhdp,h))*self.dt**2 + \
                jnp.einsum("aij,aj->ai",g,Gamma1) + jnp.einsum("aij,ajk,ak->ai",dhdp,g,Gamma2)
            positionChange = self.dt*v + 0.5*(jnp.einsum("aij,aj->ai",dvdq,v) + jnp.einsum("aij,aj->ai",dvdp,h))*self.dt**2 + \
                jnp.einsum("aij,ajk,ak->ai",dvdp,g,Gamma2)
            
            coords = jnp.where(needsUpdateMask[:,None],coords+positionChange,coords)
            momentum = jnp.where(needsUpdateMask[:,None],momentum+momentumChange,momentum)
            
            return [coords, momentum, key, stepIter+1]
        return step
    
    def run(self,initialCoords,seed):
        loop_stopper = self.make_loop_stopper()
        step = self._make_single_step()
        
        key = random.PRNGKey(seed)
        initialGuess = [initialCoords,np.zeros(initialCoords.shape),key,0]
        
        return lax.while_loop(loop_stopper,step,initialGuess)
    
class DeterministicSolution:
    def __init__(self,fric,pes,pesGrad,invInert,invInertGrad):
        self.fric = fric
        self.pes = pes
        self.pesGrad = pesGrad
        self.invInert = invInert
        self.invInertGrad = invInertGrad
        
        self.nDims = self.fric.shape[0]
    
    def func(self,t,y):
        q = y[:self.nDims]
        p = y[self.nDims:]
        
        mInv = self.invInert(q)
        
        ret = np.zeros(len(y))
        
        ret[:self.nDims] = mInv @ p
        
        ret[self.nDims:] = -self.pesGrad(q) - \
            0.5*np.einsum("jki,j,k->i",self.invInertGrad(q),p,p) - \
                self.fric @ mInv @ p
                
        return ret
    
    def solve(self,tSpan,y0,**kwargs):
        return solve_ivp(self.func,tSpan,y0,**kwargs)

"""
Fortran method:
    -Requires grid data to be fed into Fortran code
    -For each value, if array is provided:
        -If derivative(s) not provided, use np.gradient to get them
            -Check this first. Can then loop through all arrays and put them
                on the correct grid size
        -If dx is None, takes raw array
        -If dx != spacing of array, interpolate to desired grid
        


"""

