'''This Python code is an automatically generated wrapper
for Fortran code made by 'fmodpy'. The original documentation
for the Fortran source code follows.

!===================================================================================================
!   Things to feed in from python:
!       -Starting point
!           Include k copies of initial point for k runs
!           Will be serialized here, so probably load balancing is irrelevant
!       -PES/grad/hess
!       -Inverse inertia/grad/hess
!       -Friction and random force (constant)
!       -Random seed
!       -Neck values on grid, for stopping
!       -dt, maxIters, neck criteria
!       -File to output to, and whether to save the paths or not
!
!===================================================================================================
'''

import os
import ctypes
import platform
import numpy

# --------------------------------------------------------------------
#               CONFIGURATION
# 
_verbose = True
_fort_compiler = "gfortran"
_shared_object_name = "langevin_module." + platform.machine() + ".so"
_this_directory = os.path.dirname(os.path.abspath(__file__))
_path_to_lib = os.path.join(_this_directory, _shared_object_name)
_compile_options = ['-fPIC', '-shared', '-O3']
_ordered_dependencies = ['normal.f90', 'langevin-python-module.f90', 'langevin_module_c_wrapper.f90']
_symbol_files = []# 
# --------------------------------------------------------------------
#               AUTO-COMPILING
#
# Try to import the prerequisite symbols for the compiled code.
for _ in _symbol_files:
    _ = ctypes.CDLL(os.path.join(_this_directory, _), mode=ctypes.RTLD_GLOBAL)
# Try to import the existing object. If that fails, recompile and then try.
try:
    # Check to see if the source files have been modified and a recompilation is needed.
    if (max(max([0]+[os.path.getmtime(os.path.realpath(os.path.join(_this_directory,_))) for _ in _symbol_files]),
            max([0]+[os.path.getmtime(os.path.realpath(os.path.join(_this_directory,_))) for _ in _ordered_dependencies]))
        > os.path.getmtime(_path_to_lib)):
        print()
        print("WARNING: Recompiling because the modification time of a source file is newer than the library.", flush=True)
        print()
        if os.path.exists(_path_to_lib):
            os.remove(_path_to_lib)
        raise NotImplementedError(f"The newest library code has not been compiled.")
    # Import the library.
    clib = ctypes.CDLL(_path_to_lib)
except:
    # Remove the shared object if it exists, because it is faulty.
    if os.path.exists(_shared_object_name):
        os.remove(_shared_object_name)
    # Compile a new shared object.
    _command = " ".join([_fort_compiler] + _compile_options + ["-o", _shared_object_name] + _ordered_dependencies)
    if _verbose:
        print("Running system command with arguments")
        print("  ", _command)
    # Run the compilation command.
    import subprocess
    subprocess.run(_command, shell=True, cwd=_this_directory)
    # Import the shared object file as a C library with ctypes.
    clib = ctypes.CDLL(_path_to_lib)
# --------------------------------------------------------------------


class langevin_2d:
    '''!===============================================================================================
!   The solver in 2d
!==============================================================================================='''

    # Declare 'uniquecoords1'
    def get_uniquecoords1(self):
        uniquecoords1_allocated = ctypes.c_bool(False)
        uniquecoords1_dim_1 = ctypes.c_long()
        uniquecoords1 = ctypes.c_void_p()
        clib.langevin_2d_get_uniquecoords1(ctypes.byref(uniquecoords1_allocated), ctypes.byref(uniquecoords1_dim_1), ctypes.byref(uniquecoords1))
        if (not uniquecoords1_allocated.value): return None
        uniquecoords1_size = (uniquecoords1_dim_1.value)
        if (uniquecoords1_size > 0):
            uniquecoords1 = numpy.array(ctypes.cast(uniquecoords1, ctypes.POINTER(ctypes.c_float*uniquecoords1_size)).contents, copy=False)
        else:
            uniquecoords1 = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        return uniquecoords1
    def set_uniquecoords1(self, uniquecoords1):
        if ((not issubclass(type(uniquecoords1), numpy.ndarray)) or
            (not numpy.asarray(uniquecoords1).flags.f_contiguous) or
            (not (uniquecoords1.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'uniquecoords1' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            uniquecoords1 = numpy.asarray(uniquecoords1, dtype=ctypes.c_float, order='F')
        uniquecoords1_dim_1 = ctypes.c_long(uniquecoords1.shape[0])
        clib.langevin_2d_set_uniquecoords1(ctypes.byref(uniquecoords1_dim_1), ctypes.c_void_p(uniquecoords1.ctypes.data))
    uniquecoords1 = property(get_uniquecoords1, set_uniquecoords1)

    # Declare 'uniquecoords2'
    def get_uniquecoords2(self):
        uniquecoords2_allocated = ctypes.c_bool(False)
        uniquecoords2_dim_1 = ctypes.c_long()
        uniquecoords2 = ctypes.c_void_p()
        clib.langevin_2d_get_uniquecoords2(ctypes.byref(uniquecoords2_allocated), ctypes.byref(uniquecoords2_dim_1), ctypes.byref(uniquecoords2))
        if (not uniquecoords2_allocated.value): return None
        uniquecoords2_size = (uniquecoords2_dim_1.value)
        if (uniquecoords2_size > 0):
            uniquecoords2 = numpy.array(ctypes.cast(uniquecoords2, ctypes.POINTER(ctypes.c_float*uniquecoords2_size)).contents, copy=False)
        else:
            uniquecoords2 = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        return uniquecoords2
    def set_uniquecoords2(self, uniquecoords2):
        if ((not issubclass(type(uniquecoords2), numpy.ndarray)) or
            (not numpy.asarray(uniquecoords2).flags.f_contiguous) or
            (not (uniquecoords2.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'uniquecoords2' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            uniquecoords2 = numpy.asarray(uniquecoords2, dtype=ctypes.c_float, order='F')
        uniquecoords2_dim_1 = ctypes.c_long(uniquecoords2.shape[0])
        clib.langevin_2d_set_uniquecoords2(ctypes.byref(uniquecoords2_dim_1), ctypes.c_void_p(uniquecoords2.ctypes.data))
    uniquecoords2 = property(get_uniquecoords2, set_uniquecoords2)

    # Declare 'pes'
    def get_pes(self):
        pes_allocated = ctypes.c_bool(False)
        pes_dim_1 = ctypes.c_long()
        pes_dim_2 = ctypes.c_long()
        pes = ctypes.c_void_p()
        clib.langevin_2d_get_pes(ctypes.byref(pes_allocated), ctypes.byref(pes_dim_1), ctypes.byref(pes_dim_2), ctypes.byref(pes))
        if (not pes_allocated.value): return None
        pes_size = (pes_dim_1.value) * (pes_dim_2.value)
        if (pes_size > 0):
            pes = numpy.array(ctypes.cast(pes, ctypes.POINTER(ctypes.c_float*pes_size)).contents, copy=False)
        else:
            pes = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        pes = pes.reshape(pes_dim_2.value,pes_dim_1.value).T
        return pes
    def set_pes(self, pes):
        if ((not issubclass(type(pes), numpy.ndarray)) or
            (not numpy.asarray(pes).flags.f_contiguous) or
            (not (pes.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'pes' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            pes = numpy.asarray(pes, dtype=ctypes.c_float, order='F')
        pes_dim_1 = ctypes.c_long(pes.shape[0])
        pes_dim_2 = ctypes.c_long(pes.shape[1])
        clib.langevin_2d_set_pes(ctypes.byref(pes_dim_1), ctypes.byref(pes_dim_2), ctypes.c_void_p(pes.ctypes.data))
    pes = property(get_pes, set_pes)

    # Declare 'neckvals'
    def get_neckvals(self):
        neckvals_allocated = ctypes.c_bool(False)
        neckvals_dim_1 = ctypes.c_long()
        neckvals_dim_2 = ctypes.c_long()
        neckvals = ctypes.c_void_p()
        clib.langevin_2d_get_neckvals(ctypes.byref(neckvals_allocated), ctypes.byref(neckvals_dim_1), ctypes.byref(neckvals_dim_2), ctypes.byref(neckvals))
        if (not neckvals_allocated.value): return None
        neckvals_size = (neckvals_dim_1.value) * (neckvals_dim_2.value)
        if (neckvals_size > 0):
            neckvals = numpy.array(ctypes.cast(neckvals, ctypes.POINTER(ctypes.c_float*neckvals_size)).contents, copy=False)
        else:
            neckvals = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        neckvals = neckvals.reshape(neckvals_dim_2.value,neckvals_dim_1.value).T
        return neckvals
    def set_neckvals(self, neckvals):
        if ((not issubclass(type(neckvals), numpy.ndarray)) or
            (not numpy.asarray(neckvals).flags.f_contiguous) or
            (not (neckvals.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'neckvals' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            neckvals = numpy.asarray(neckvals, dtype=ctypes.c_float, order='F')
        neckvals_dim_1 = ctypes.c_long(neckvals.shape[0])
        neckvals_dim_2 = ctypes.c_long(neckvals.shape[1])
        clib.langevin_2d_set_neckvals(ctypes.byref(neckvals_dim_1), ctypes.byref(neckvals_dim_2), ctypes.c_void_p(neckvals.ctypes.data))
    neckvals = property(get_neckvals, set_neckvals)

    # Declare 'pesgrad'
    def get_pesgrad(self):
        pesgrad_allocated = ctypes.c_bool(False)
        pesgrad_dim_1 = ctypes.c_long()
        pesgrad_dim_2 = ctypes.c_long()
        pesgrad_dim_3 = ctypes.c_long()
        pesgrad = ctypes.c_void_p()
        clib.langevin_2d_get_pesgrad(ctypes.byref(pesgrad_allocated), ctypes.byref(pesgrad_dim_1), ctypes.byref(pesgrad_dim_2), ctypes.byref(pesgrad_dim_3), ctypes.byref(pesgrad))
        if (not pesgrad_allocated.value): return None
        pesgrad_size = (pesgrad_dim_1.value) * (pesgrad_dim_2.value) * (pesgrad_dim_3.value)
        if (pesgrad_size > 0):
            pesgrad = numpy.array(ctypes.cast(pesgrad, ctypes.POINTER(ctypes.c_float*pesgrad_size)).contents, copy=False)
        else:
            pesgrad = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        pesgrad = pesgrad.reshape(pesgrad_dim_3.value,pesgrad_dim_2.value,pesgrad_dim_1.value).T
        return pesgrad
    def set_pesgrad(self, pesgrad):
        if ((not issubclass(type(pesgrad), numpy.ndarray)) or
            (not numpy.asarray(pesgrad).flags.f_contiguous) or
            (not (pesgrad.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'pesgrad' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            pesgrad = numpy.asarray(pesgrad, dtype=ctypes.c_float, order='F')
        pesgrad_dim_1 = ctypes.c_long(pesgrad.shape[0])
        pesgrad_dim_2 = ctypes.c_long(pesgrad.shape[1])
        pesgrad_dim_3 = ctypes.c_long(pesgrad.shape[2])
        clib.langevin_2d_set_pesgrad(ctypes.byref(pesgrad_dim_1), ctypes.byref(pesgrad_dim_2), ctypes.byref(pesgrad_dim_3), ctypes.c_void_p(pesgrad.ctypes.data))
    pesgrad = property(get_pesgrad, set_pesgrad)

    # Declare 'peshess'
    def get_peshess(self):
        peshess_allocated = ctypes.c_bool(False)
        peshess_dim_1 = ctypes.c_long()
        peshess_dim_2 = ctypes.c_long()
        peshess_dim_3 = ctypes.c_long()
        peshess = ctypes.c_void_p()
        clib.langevin_2d_get_peshess(ctypes.byref(peshess_allocated), ctypes.byref(peshess_dim_1), ctypes.byref(peshess_dim_2), ctypes.byref(peshess_dim_3), ctypes.byref(peshess))
        if (not peshess_allocated.value): return None
        peshess_size = (peshess_dim_1.value) * (peshess_dim_2.value) * (peshess_dim_3.value)
        if (peshess_size > 0):
            peshess = numpy.array(ctypes.cast(peshess, ctypes.POINTER(ctypes.c_float*peshess_size)).contents, copy=False)
        else:
            peshess = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        peshess = peshess.reshape(peshess_dim_3.value,peshess_dim_2.value,peshess_dim_1.value).T
        return peshess
    def set_peshess(self, peshess):
        if ((not issubclass(type(peshess), numpy.ndarray)) or
            (not numpy.asarray(peshess).flags.f_contiguous) or
            (not (peshess.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'peshess' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            peshess = numpy.asarray(peshess, dtype=ctypes.c_float, order='F')
        peshess_dim_1 = ctypes.c_long(peshess.shape[0])
        peshess_dim_2 = ctypes.c_long(peshess.shape[1])
        peshess_dim_3 = ctypes.c_long(peshess.shape[2])
        clib.langevin_2d_set_peshess(ctypes.byref(peshess_dim_1), ctypes.byref(peshess_dim_2), ctypes.byref(peshess_dim_3), ctypes.c_void_p(peshess.ctypes.data))
    peshess = property(get_peshess, set_peshess)

    # Declare 'inversemetric'
    def get_inversemetric(self):
        inversemetric_allocated = ctypes.c_bool(False)
        inversemetric_dim_1 = ctypes.c_long()
        inversemetric_dim_2 = ctypes.c_long()
        inversemetric_dim_3 = ctypes.c_long()
        inversemetric = ctypes.c_void_p()
        clib.langevin_2d_get_inversemetric(ctypes.byref(inversemetric_allocated), ctypes.byref(inversemetric_dim_1), ctypes.byref(inversemetric_dim_2), ctypes.byref(inversemetric_dim_3), ctypes.byref(inversemetric))
        if (not inversemetric_allocated.value): return None
        inversemetric_size = (inversemetric_dim_1.value) * (inversemetric_dim_2.value) * (inversemetric_dim_3.value)
        if (inversemetric_size > 0):
            inversemetric = numpy.array(ctypes.cast(inversemetric, ctypes.POINTER(ctypes.c_float*inversemetric_size)).contents, copy=False)
        else:
            inversemetric = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        inversemetric = inversemetric.reshape(inversemetric_dim_3.value,inversemetric_dim_2.value,inversemetric_dim_1.value).T
        return inversemetric
    def set_inversemetric(self, inversemetric):
        if ((not issubclass(type(inversemetric), numpy.ndarray)) or
            (not numpy.asarray(inversemetric).flags.f_contiguous) or
            (not (inversemetric.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'inversemetric' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            inversemetric = numpy.asarray(inversemetric, dtype=ctypes.c_float, order='F')
        inversemetric_dim_1 = ctypes.c_long(inversemetric.shape[0])
        inversemetric_dim_2 = ctypes.c_long(inversemetric.shape[1])
        inversemetric_dim_3 = ctypes.c_long(inversemetric.shape[2])
        clib.langevin_2d_set_inversemetric(ctypes.byref(inversemetric_dim_1), ctypes.byref(inversemetric_dim_2), ctypes.byref(inversemetric_dim_3), ctypes.c_void_p(inversemetric.ctypes.data))
    inversemetric = property(get_inversemetric, set_inversemetric)

    # Declare 'inversemetricgrad'
    def get_inversemetricgrad(self):
        inversemetricgrad_allocated = ctypes.c_bool(False)
        inversemetricgrad_dim_1 = ctypes.c_long()
        inversemetricgrad_dim_2 = ctypes.c_long()
        inversemetricgrad_dim_3 = ctypes.c_long()
        inversemetricgrad_dim_4 = ctypes.c_long()
        inversemetricgrad = ctypes.c_void_p()
        clib.langevin_2d_get_inversemetricgrad(ctypes.byref(inversemetricgrad_allocated), ctypes.byref(inversemetricgrad_dim_1), ctypes.byref(inversemetricgrad_dim_2), ctypes.byref(inversemetricgrad_dim_3), ctypes.byref(inversemetricgrad_dim_4), ctypes.byref(inversemetricgrad))
        if (not inversemetricgrad_allocated.value): return None
        inversemetricgrad_size = (inversemetricgrad_dim_1.value) * (inversemetricgrad_dim_2.value) * (inversemetricgrad_dim_3.value) * (inversemetricgrad_dim_4.value)
        if (inversemetricgrad_size > 0):
            inversemetricgrad = numpy.array(ctypes.cast(inversemetricgrad, ctypes.POINTER(ctypes.c_float*inversemetricgrad_size)).contents, copy=False)
        else:
            inversemetricgrad = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        inversemetricgrad = inversemetricgrad.reshape(inversemetricgrad_dim_4.value,inversemetricgrad_dim_3.value,inversemetricgrad_dim_2.value,inversemetricgrad_dim_1.value).T
        return inversemetricgrad
    def set_inversemetricgrad(self, inversemetricgrad):
        if ((not issubclass(type(inversemetricgrad), numpy.ndarray)) or
            (not numpy.asarray(inversemetricgrad).flags.f_contiguous) or
            (not (inversemetricgrad.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'inversemetricgrad' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            inversemetricgrad = numpy.asarray(inversemetricgrad, dtype=ctypes.c_float, order='F')
        inversemetricgrad_dim_1 = ctypes.c_long(inversemetricgrad.shape[0])
        inversemetricgrad_dim_2 = ctypes.c_long(inversemetricgrad.shape[1])
        inversemetricgrad_dim_3 = ctypes.c_long(inversemetricgrad.shape[2])
        inversemetricgrad_dim_4 = ctypes.c_long(inversemetricgrad.shape[3])
        clib.langevin_2d_set_inversemetricgrad(ctypes.byref(inversemetricgrad_dim_1), ctypes.byref(inversemetricgrad_dim_2), ctypes.byref(inversemetricgrad_dim_3), ctypes.byref(inversemetricgrad_dim_4), ctypes.c_void_p(inversemetricgrad.ctypes.data))
    inversemetricgrad = property(get_inversemetricgrad, set_inversemetricgrad)

    # Declare 'inversemetrichess'
    def get_inversemetrichess(self):
        inversemetrichess_allocated = ctypes.c_bool(False)
        inversemetrichess_dim_1 = ctypes.c_long()
        inversemetrichess_dim_2 = ctypes.c_long()
        inversemetrichess_dim_3 = ctypes.c_long()
        inversemetrichess_dim_4 = ctypes.c_long()
        inversemetrichess = ctypes.c_void_p()
        clib.langevin_2d_get_inversemetrichess(ctypes.byref(inversemetrichess_allocated), ctypes.byref(inversemetrichess_dim_1), ctypes.byref(inversemetrichess_dim_2), ctypes.byref(inversemetrichess_dim_3), ctypes.byref(inversemetrichess_dim_4), ctypes.byref(inversemetrichess))
        if (not inversemetrichess_allocated.value): return None
        inversemetrichess_size = (inversemetrichess_dim_1.value) * (inversemetrichess_dim_2.value) * (inversemetrichess_dim_3.value) * (inversemetrichess_dim_4.value)
        if (inversemetrichess_size > 0):
            inversemetrichess = numpy.array(ctypes.cast(inversemetrichess, ctypes.POINTER(ctypes.c_float*inversemetrichess_size)).contents, copy=False)
        else:
            inversemetrichess = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        inversemetrichess = inversemetrichess.reshape(inversemetrichess_dim_4.value,inversemetrichess_dim_3.value,inversemetrichess_dim_2.value,inversemetrichess_dim_1.value).T
        return inversemetrichess
    def set_inversemetrichess(self, inversemetrichess):
        if ((not issubclass(type(inversemetrichess), numpy.ndarray)) or
            (not numpy.asarray(inversemetrichess).flags.f_contiguous) or
            (not (inversemetrichess.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'inversemetrichess' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            inversemetrichess = numpy.asarray(inversemetrichess, dtype=ctypes.c_float, order='F')
        inversemetrichess_dim_1 = ctypes.c_long(inversemetrichess.shape[0])
        inversemetrichess_dim_2 = ctypes.c_long(inversemetrichess.shape[1])
        inversemetrichess_dim_3 = ctypes.c_long(inversemetrichess.shape[2])
        inversemetrichess_dim_4 = ctypes.c_long(inversemetrichess.shape[3])
        clib.langevin_2d_set_inversemetrichess(ctypes.byref(inversemetrichess_dim_1), ctypes.byref(inversemetrichess_dim_2), ctypes.byref(inversemetrichess_dim_3), ctypes.byref(inversemetrichess_dim_4), ctypes.c_void_p(inversemetrichess.ctypes.data))
    inversemetrichess = property(get_inversemetrichess, set_inversemetrichess)

    # Declare 'friction'
    def get_friction(self):
        friction_allocated = ctypes.c_bool(False)
        friction_dim_1 = ctypes.c_long()
        friction_dim_2 = ctypes.c_long()
        friction = ctypes.c_void_p()
        clib.langevin_2d_get_friction(ctypes.byref(friction_allocated), ctypes.byref(friction_dim_1), ctypes.byref(friction_dim_2), ctypes.byref(friction))
        if (not friction_allocated.value): return None
        friction_size = (friction_dim_1.value) * (friction_dim_2.value)
        if (friction_size > 0):
            friction = numpy.array(ctypes.cast(friction, ctypes.POINTER(ctypes.c_float*friction_size)).contents, copy=False)
        else:
            friction = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        friction = friction.reshape(friction_dim_2.value,friction_dim_1.value).T
        return friction
    def set_friction(self, friction):
        if ((not issubclass(type(friction), numpy.ndarray)) or
            (not numpy.asarray(friction).flags.f_contiguous) or
            (not (friction.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'friction' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            friction = numpy.asarray(friction, dtype=ctypes.c_float, order='F')
        friction_dim_1 = ctypes.c_long(friction.shape[0])
        friction_dim_2 = ctypes.c_long(friction.shape[1])
        clib.langevin_2d_set_friction(ctypes.byref(friction_dim_1), ctypes.byref(friction_dim_2), ctypes.c_void_p(friction.ctypes.data))
    friction = property(get_friction, set_friction)

    # Declare 'randomforce'
    def get_randomforce(self):
        randomforce_allocated = ctypes.c_bool(False)
        randomforce_dim_1 = ctypes.c_long()
        randomforce_dim_2 = ctypes.c_long()
        randomforce = ctypes.c_void_p()
        clib.langevin_2d_get_randomforce(ctypes.byref(randomforce_allocated), ctypes.byref(randomforce_dim_1), ctypes.byref(randomforce_dim_2), ctypes.byref(randomforce))
        if (not randomforce_allocated.value): return None
        randomforce_size = (randomforce_dim_1.value) * (randomforce_dim_2.value)
        if (randomforce_size > 0):
            randomforce = numpy.array(ctypes.cast(randomforce, ctypes.POINTER(ctypes.c_float*randomforce_size)).contents, copy=False)
        else:
            randomforce = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        randomforce = randomforce.reshape(randomforce_dim_2.value,randomforce_dim_1.value).T
        return randomforce
    def set_randomforce(self, randomforce):
        if ((not issubclass(type(randomforce), numpy.ndarray)) or
            (not numpy.asarray(randomforce).flags.f_contiguous) or
            (not (randomforce.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'randomforce' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            randomforce = numpy.asarray(randomforce, dtype=ctypes.c_float, order='F')
        randomforce_dim_1 = ctypes.c_long(randomforce.shape[0])
        randomforce_dim_2 = ctypes.c_long(randomforce.shape[1])
        clib.langevin_2d_set_randomforce(ctypes.byref(randomforce_dim_1), ctypes.byref(randomforce_dim_2), ctypes.c_void_p(randomforce.ctypes.data))
    randomforce = property(get_randomforce, set_randomforce)

    # Declare 'ngridpoints'
    def get_ngridpoints(self):
        ngridpoints_dim_1 = ctypes.c_long()
        ngridpoints = ctypes.c_void_p()
        clib.langevin_2d_get_ngridpoints(ctypes.byref(ngridpoints_dim_1), ctypes.byref(ngridpoints))
        ngridpoints_size = (ngridpoints_dim_1.value)
        if (ngridpoints_size > 0):
            ngridpoints = numpy.array(ctypes.cast(ngridpoints, ctypes.POINTER(ctypes.c_int*ngridpoints_size)).contents, copy=False)
        else:
            ngridpoints = numpy.zeros((0,), dtype=ctypes.c_int, order='F')
        return ngridpoints
    def set_ngridpoints(self, ngridpoints):
        if ((not issubclass(type(ngridpoints), numpy.ndarray)) or
            (not numpy.asarray(ngridpoints).flags.f_contiguous) or
            (not (ngridpoints.dtype == numpy.dtype(ctypes.c_int)))):
            import warnings
            warnings.warn("The provided argument 'ngridpoints' was not an f_contiguous NumPy array of type 'ctypes.c_int' (or equivalent). Automatically converting (probably creating a full copy).")
            ngridpoints = numpy.asarray(ngridpoints, dtype=ctypes.c_int, order='F')
        ngridpoints_dim_1 = ctypes.c_long(ngridpoints.shape[0])
        clib.langevin_2d_set_ngridpoints(ctypes.byref(ngridpoints_dim_1), ctypes.c_void_p(ngridpoints.ctypes.data))
    ngridpoints = property(get_ngridpoints, set_ngridpoints)

    # Declare 'dx'
    def get_dx(self):
        dx_dim_1 = ctypes.c_long()
        dx = ctypes.c_void_p()
        clib.langevin_2d_get_dx(ctypes.byref(dx_dim_1), ctypes.byref(dx))
        dx_size = (dx_dim_1.value)
        if (dx_size > 0):
            dx = numpy.array(ctypes.cast(dx, ctypes.POINTER(ctypes.c_float*dx_size)).contents, copy=False)
        else:
            dx = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        return dx
    def set_dx(self, dx):
        if ((not issubclass(type(dx), numpy.ndarray)) or
            (not numpy.asarray(dx).flags.f_contiguous) or
            (not (dx.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'dx' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            dx = numpy.asarray(dx, dtype=ctypes.c_float, order='F')
        dx_dim_1 = ctypes.c_long(dx.shape[0])
        clib.langevin_2d_set_dx(ctypes.byref(dx_dim_1), ctypes.c_void_p(dx.ctypes.data))
    dx = property(get_dx, set_dx)

    # Declare 'massval'
    def get_massval(self):
        massval = ctypes.c_float()
        clib.langevin_2d_get_massval(ctypes.byref(massval))
        return massval.value
    def set_massval(self, massval):
        massval = ctypes.c_float(massval)
        clib.langevin_2d_set_massval(ctypes.byref(massval))
    massval = property(get_massval, set_massval)

    # Declare 'startingeneg'
    def get_startingeneg(self):
        startingeneg = ctypes.c_float()
        clib.langevin_2d_get_startingeneg(ctypes.byref(startingeneg))
        return startingeneg.value
    def set_startingeneg(self, startingeneg):
        startingeneg = ctypes.c_float(startingeneg)
        clib.langevin_2d_set_startingeneg(ctypes.byref(startingeneg))
    startingeneg = property(get_startingeneg, set_startingeneg)

    # Declare 'dt'
    def get_dt(self):
        dt = ctypes.c_float()
        clib.langevin_2d_get_dt(ctypes.byref(dt))
        return dt.value
    def set_dt(self, dt):
        dt = ctypes.c_float(dt)
        clib.langevin_2d_set_dt(ctypes.byref(dt))
    dt = property(get_dt, set_dt)

    # Declare 'maxtime'
    def get_maxtime(self):
        maxtime = ctypes.c_float()
        clib.langevin_2d_get_maxtime(ctypes.byref(maxtime))
        return maxtime.value
    def set_maxtime(self, maxtime):
        maxtime = ctypes.c_float(maxtime)
        clib.langevin_2d_set_maxtime(ctypes.byref(maxtime))
    maxtime = property(get_maxtime, set_maxtime)

    # Declare 'scissionneckval'
    def get_scissionneckval(self):
        scissionneckval = ctypes.c_float()
        clib.langevin_2d_get_scissionneckval(ctypes.byref(scissionneckval))
        return scissionneckval.value
    def set_scissionneckval(self, scissionneckval):
        scissionneckval = ctypes.c_float(scissionneckval)
        clib.langevin_2d_set_scissionneckval(ctypes.byref(scissionneckval))
    scissionneckval = property(get_scissionneckval, set_scissionneckval)

    # Declare 'leveldensity'
    def get_leveldensity(self):
        leveldensity = ctypes.c_float()
        clib.langevin_2d_get_leveldensity(ctypes.byref(leveldensity))
        return leveldensity.value
    def set_leveldensity(self, leveldensity):
        leveldensity = ctypes.c_float(leveldensity)
        clib.langevin_2d_set_leveldensity(ctypes.byref(leveldensity))
    leveldensity = property(get_leveldensity, set_leveldensity)

    # Declare 'maxiters'
    def get_maxiters(self):
        maxiters = ctypes.c_int()
        clib.langevin_2d_get_maxiters(ctypes.byref(maxiters))
        return maxiters.value
    def set_maxiters(self, maxiters):
        maxiters = ctypes.c_int(maxiters)
        clib.langevin_2d_set_maxiters(ctypes.byref(maxiters))
    maxiters = property(get_maxiters, set_maxiters)

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine SETUP_ARRAYS
    
    def setup_arrays(self, uncoords1in, uncoords2in, pesin, pesgradin, peshessin, invmetin, invmetgradin, invmethessin, fricin, randforcein, neckin):
        ''''''
        
        # Setting up "uncoords1in"
        if ((not issubclass(type(uncoords1in), numpy.ndarray)) or
            (not numpy.asarray(uncoords1in).flags.f_contiguous) or
            (not (uncoords1in.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'uncoords1in' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            uncoords1in = numpy.asarray(uncoords1in, dtype=ctypes.c_double, order='F')
        uncoords1in_dim_1 = ctypes.c_long(uncoords1in.shape[0])
        
        # Setting up "uncoords2in"
        if ((not issubclass(type(uncoords2in), numpy.ndarray)) or
            (not numpy.asarray(uncoords2in).flags.f_contiguous) or
            (not (uncoords2in.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'uncoords2in' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            uncoords2in = numpy.asarray(uncoords2in, dtype=ctypes.c_double, order='F')
        uncoords2in_dim_1 = ctypes.c_long(uncoords2in.shape[0])
        
        # Setting up "pesin"
        if ((not issubclass(type(pesin), numpy.ndarray)) or
            (not numpy.asarray(pesin).flags.f_contiguous) or
            (not (pesin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'pesin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            pesin = numpy.asarray(pesin, dtype=ctypes.c_double, order='F')
        pesin_dim_1 = ctypes.c_long(pesin.shape[0])
        pesin_dim_2 = ctypes.c_long(pesin.shape[1])
        
        # Setting up "pesgradin"
        if ((not issubclass(type(pesgradin), numpy.ndarray)) or
            (not numpy.asarray(pesgradin).flags.f_contiguous) or
            (not (pesgradin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'pesgradin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            pesgradin = numpy.asarray(pesgradin, dtype=ctypes.c_double, order='F')
        pesgradin_dim_1 = ctypes.c_long(pesgradin.shape[0])
        pesgradin_dim_2 = ctypes.c_long(pesgradin.shape[1])
        pesgradin_dim_3 = ctypes.c_long(pesgradin.shape[2])
        
        # Setting up "peshessin"
        if ((not issubclass(type(peshessin), numpy.ndarray)) or
            (not numpy.asarray(peshessin).flags.f_contiguous) or
            (not (peshessin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'peshessin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            peshessin = numpy.asarray(peshessin, dtype=ctypes.c_double, order='F')
        peshessin_dim_1 = ctypes.c_long(peshessin.shape[0])
        peshessin_dim_2 = ctypes.c_long(peshessin.shape[1])
        peshessin_dim_3 = ctypes.c_long(peshessin.shape[2])
        
        # Setting up "invmetin"
        if ((not issubclass(type(invmetin), numpy.ndarray)) or
            (not numpy.asarray(invmetin).flags.f_contiguous) or
            (not (invmetin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'invmetin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            invmetin = numpy.asarray(invmetin, dtype=ctypes.c_double, order='F')
        invmetin_dim_1 = ctypes.c_long(invmetin.shape[0])
        invmetin_dim_2 = ctypes.c_long(invmetin.shape[1])
        invmetin_dim_3 = ctypes.c_long(invmetin.shape[2])
        
        # Setting up "invmetgradin"
        if ((not issubclass(type(invmetgradin), numpy.ndarray)) or
            (not numpy.asarray(invmetgradin).flags.f_contiguous) or
            (not (invmetgradin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'invmetgradin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            invmetgradin = numpy.asarray(invmetgradin, dtype=ctypes.c_double, order='F')
        invmetgradin_dim_1 = ctypes.c_long(invmetgradin.shape[0])
        invmetgradin_dim_2 = ctypes.c_long(invmetgradin.shape[1])
        invmetgradin_dim_3 = ctypes.c_long(invmetgradin.shape[2])
        invmetgradin_dim_4 = ctypes.c_long(invmetgradin.shape[3])
        
        # Setting up "invmethessin"
        if ((not issubclass(type(invmethessin), numpy.ndarray)) or
            (not numpy.asarray(invmethessin).flags.f_contiguous) or
            (not (invmethessin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'invmethessin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            invmethessin = numpy.asarray(invmethessin, dtype=ctypes.c_double, order='F')
        invmethessin_dim_1 = ctypes.c_long(invmethessin.shape[0])
        invmethessin_dim_2 = ctypes.c_long(invmethessin.shape[1])
        invmethessin_dim_3 = ctypes.c_long(invmethessin.shape[2])
        invmethessin_dim_4 = ctypes.c_long(invmethessin.shape[3])
        
        # Setting up "fricin"
        if ((not issubclass(type(fricin), numpy.ndarray)) or
            (not numpy.asarray(fricin).flags.f_contiguous) or
            (not (fricin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'fricin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            fricin = numpy.asarray(fricin, dtype=ctypes.c_double, order='F')
        fricin_dim_1 = ctypes.c_long(fricin.shape[0])
        fricin_dim_2 = ctypes.c_long(fricin.shape[1])
        
        # Setting up "randforcein"
        if ((not issubclass(type(randforcein), numpy.ndarray)) or
            (not numpy.asarray(randforcein).flags.f_contiguous) or
            (not (randforcein.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'randforcein' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            randforcein = numpy.asarray(randforcein, dtype=ctypes.c_double, order='F')
        randforcein_dim_1 = ctypes.c_long(randforcein.shape[0])
        randforcein_dim_2 = ctypes.c_long(randforcein.shape[1])
        
        # Setting up "neckin"
        if ((not issubclass(type(neckin), numpy.ndarray)) or
            (not numpy.asarray(neckin).flags.f_contiguous) or
            (not (neckin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'neckin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            neckin = numpy.asarray(neckin, dtype=ctypes.c_double, order='F')
        neckin_dim_1 = ctypes.c_long(neckin.shape[0])
        neckin_dim_2 = ctypes.c_long(neckin.shape[1])
    
        # Call C-accessible Fortran wrapper.
        clib.c_setup_arrays(ctypes.byref(uncoords1in_dim_1), ctypes.c_void_p(uncoords1in.ctypes.data), ctypes.byref(uncoords2in_dim_1), ctypes.c_void_p(uncoords2in.ctypes.data), ctypes.byref(pesin_dim_1), ctypes.byref(pesin_dim_2), ctypes.c_void_p(pesin.ctypes.data), ctypes.byref(pesgradin_dim_1), ctypes.byref(pesgradin_dim_2), ctypes.byref(pesgradin_dim_3), ctypes.c_void_p(pesgradin.ctypes.data), ctypes.byref(peshessin_dim_1), ctypes.byref(peshessin_dim_2), ctypes.byref(peshessin_dim_3), ctypes.c_void_p(peshessin.ctypes.data), ctypes.byref(invmetin_dim_1), ctypes.byref(invmetin_dim_2), ctypes.byref(invmetin_dim_3), ctypes.c_void_p(invmetin.ctypes.data), ctypes.byref(invmetgradin_dim_1), ctypes.byref(invmetgradin_dim_2), ctypes.byref(invmetgradin_dim_3), ctypes.byref(invmetgradin_dim_4), ctypes.c_void_p(invmetgradin.ctypes.data), ctypes.byref(invmethessin_dim_1), ctypes.byref(invmethessin_dim_2), ctypes.byref(invmethessin_dim_3), ctypes.byref(invmethessin_dim_4), ctypes.c_void_p(invmethessin.ctypes.data), ctypes.byref(fricin_dim_1), ctypes.byref(fricin_dim_2), ctypes.c_void_p(fricin.ctypes.data), ctypes.byref(randforcein_dim_1), ctypes.byref(randforcein_dim_2), ctypes.c_void_p(randforcein.ctypes.data), ctypes.byref(neckin_dim_1), ctypes.byref(neckin_dim_2), ctypes.c_void_p(neckin.ctypes.data))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return 

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine SETUP_PARAMS
    
    def setup_params(self, massvalin, startingenegin, dtin, maxtimein, scissionneckvalin, useseed, seed):
        ''''''
        
        # Setting up "massvalin"
        if (type(massvalin) is not ctypes.c_double): massvalin = ctypes.c_double(massvalin)
        
        # Setting up "startingenegin"
        if (type(startingenegin) is not ctypes.c_double): startingenegin = ctypes.c_double(startingenegin)
        
        # Setting up "dtin"
        if (type(dtin) is not ctypes.c_double): dtin = ctypes.c_double(dtin)
        
        # Setting up "maxtimein"
        if (type(maxtimein) is not ctypes.c_double): maxtimein = ctypes.c_double(maxtimein)
        
        # Setting up "scissionneckvalin"
        if (type(scissionneckvalin) is not ctypes.c_double): scissionneckvalin = ctypes.c_double(scissionneckvalin)
        
        # Setting up "useseed"
        if (type(useseed) is not ctypes.c_int): useseed = ctypes.c_int(useseed)
        
        # Setting up "seed"
        if ((not issubclass(type(seed), numpy.ndarray)) or
            (not numpy.asarray(seed).flags.f_contiguous) or
            (not (seed.dtype == numpy.dtype(ctypes.c_int)))):
            import warnings
            warnings.warn("The provided argument 'seed' was not an f_contiguous NumPy array of type 'ctypes.c_int' (or equivalent). Automatically converting (probably creating a full copy).")
            seed = numpy.asarray(seed, dtype=ctypes.c_int, order='F')
        seed_dim_1 = ctypes.c_long(seed.shape[0])
    
        # Call C-accessible Fortran wrapper.
        clib.c_setup_params(ctypes.byref(massvalin), ctypes.byref(startingenegin), ctypes.byref(dtin), ctypes.byref(maxtimein), ctypes.byref(scissionneckvalin), ctypes.byref(useseed), ctypes.byref(seed_dim_1), ctypes.c_void_p(seed.ctypes.data))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return 

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine GET_INTERPOLATION_INDICES
    
    def get_interpolation_indices(self, coords, idx=None):
        ''''''
        
        # Setting up "coords"
        if ((not issubclass(type(coords), numpy.ndarray)) or
            (not numpy.asarray(coords).flags.f_contiguous) or
            (not (coords.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'coords' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            coords = numpy.asarray(coords, dtype=ctypes.c_double, order='F')
        coords_dim_1 = ctypes.c_long(coords.shape[0])
        
        # Setting up "idx"
        if (idx is None):
            idx = numpy.zeros(shape=(2), dtype=ctypes.c_int, order='F')
        elif ((not issubclass(type(idx), numpy.ndarray)) or
              (not numpy.asarray(idx).flags.f_contiguous) or
              (not (idx.dtype == numpy.dtype(ctypes.c_int)))):
            import warnings
            warnings.warn("The provided argument 'idx' was not an f_contiguous NumPy array of type 'ctypes.c_int' (or equivalent). Automatically converting (probably creating a full copy).")
            idx = numpy.asarray(idx, dtype=ctypes.c_int, order='F')
        idx_dim_1 = ctypes.c_long(idx.shape[0])
        
        # Setting up "isoutofbounds"
        isoutofbounds = ctypes.c_int()
    
        # Call C-accessible Fortran wrapper.
        clib.c_get_interpolation_indices(ctypes.byref(coords_dim_1), ctypes.c_void_p(coords.ctypes.data), ctypes.byref(idx_dim_1), ctypes.c_void_p(idx.ctypes.data), ctypes.byref(isoutofbounds))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return coords, idx, isoutofbounds.value

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine BILINEAR_INTERP
    
    def bilinear_interp(self, idx, coords, isoutofbounds, arrin):
        ''''''
        
        # Setting up "idx"
        if ((not issubclass(type(idx), numpy.ndarray)) or
            (not numpy.asarray(idx).flags.f_contiguous) or
            (not (idx.dtype == numpy.dtype(ctypes.c_int)))):
            import warnings
            warnings.warn("The provided argument 'idx' was not an f_contiguous NumPy array of type 'ctypes.c_int' (or equivalent). Automatically converting (probably creating a full copy).")
            idx = numpy.asarray(idx, dtype=ctypes.c_int, order='F')
        idx_dim_1 = ctypes.c_long(idx.shape[0])
        
        # Setting up "coords"
        if ((not issubclass(type(coords), numpy.ndarray)) or
            (not numpy.asarray(coords).flags.f_contiguous) or
            (not (coords.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'coords' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            coords = numpy.asarray(coords, dtype=ctypes.c_double, order='F')
        coords_dim_1 = ctypes.c_long(coords.shape[0])
        
        # Setting up "isoutofbounds"
        if (type(isoutofbounds) is not ctypes.c_int): isoutofbounds = ctypes.c_int(isoutofbounds)
        
        # Setting up "arrin"
        if ((not issubclass(type(arrin), numpy.ndarray)) or
            (not numpy.asarray(arrin).flags.f_contiguous) or
            (not (arrin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'arrin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            arrin = numpy.asarray(arrin, dtype=ctypes.c_double, order='F')
        arrin_dim_1 = ctypes.c_long(arrin.shape[0])
        arrin_dim_2 = ctypes.c_long(arrin.shape[1])
        
        # Setting up "val"
        val = ctypes.c_double()
    
        # Call C-accessible Fortran wrapper.
        clib.c_bilinear_interp(ctypes.byref(idx_dim_1), ctypes.c_void_p(idx.ctypes.data), ctypes.byref(coords_dim_1), ctypes.c_void_p(coords.ctypes.data), ctypes.byref(isoutofbounds), ctypes.byref(arrin_dim_1), ctypes.byref(arrin_dim_2), ctypes.c_void_p(arrin.ctypes.data), ctypes.byref(val))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return idx, coords, isoutofbounds.value, arrin, val.value

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine INTERP_WRAPPER
    
    def interp_wrapper(self, coords, eneggrad=None, eneghess=None, invmet=None, invmetgrad=None, invmethess=None):
        '''!total of 42 evals, matching the previous code, minus (3+6)*2=18 from friction/grad
!and randForce/grad
!trust this wrapper against scipy.interpolate.interpn'''
        
        # Setting up "coords"
        if ((not issubclass(type(coords), numpy.ndarray)) or
            (not numpy.asarray(coords).flags.f_contiguous) or
            (not (coords.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'coords' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            coords = numpy.asarray(coords, dtype=ctypes.c_double, order='F')
        coords_dim_1 = ctypes.c_long(coords.shape[0])
        
        # Setting up "eneg"
        eneg = ctypes.c_double()
        
        # Setting up "eneggrad"
        if (eneggrad is None):
            eneggrad = numpy.zeros(shape=(2), dtype=ctypes.c_double, order='F')
        elif ((not issubclass(type(eneggrad), numpy.ndarray)) or
              (not numpy.asarray(eneggrad).flags.f_contiguous) or
              (not (eneggrad.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'eneggrad' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            eneggrad = numpy.asarray(eneggrad, dtype=ctypes.c_double, order='F')
        eneggrad_dim_1 = ctypes.c_long(eneggrad.shape[0])
        
        # Setting up "eneghess"
        if (eneghess is None):
            eneghess = numpy.zeros(shape=(3), dtype=ctypes.c_double, order='F')
        elif ((not issubclass(type(eneghess), numpy.ndarray)) or
              (not numpy.asarray(eneghess).flags.f_contiguous) or
              (not (eneghess.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'eneghess' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            eneghess = numpy.asarray(eneghess, dtype=ctypes.c_double, order='F')
        eneghess_dim_1 = ctypes.c_long(eneghess.shape[0])
        
        # Setting up "invmet"
        if (invmet is None):
            invmet = numpy.zeros(shape=(3), dtype=ctypes.c_double, order='F')
        elif ((not issubclass(type(invmet), numpy.ndarray)) or
              (not numpy.asarray(invmet).flags.f_contiguous) or
              (not (invmet.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'invmet' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            invmet = numpy.asarray(invmet, dtype=ctypes.c_double, order='F')
        invmet_dim_1 = ctypes.c_long(invmet.shape[0])
        
        # Setting up "invmetgrad"
        if (invmetgrad is None):
            invmetgrad = numpy.zeros(shape=(3, 2), dtype=ctypes.c_double, order='F')
        elif ((not issubclass(type(invmetgrad), numpy.ndarray)) or
              (not numpy.asarray(invmetgrad).flags.f_contiguous) or
              (not (invmetgrad.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'invmetgrad' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            invmetgrad = numpy.asarray(invmetgrad, dtype=ctypes.c_double, order='F')
        invmetgrad_dim_1 = ctypes.c_long(invmetgrad.shape[0])
        invmetgrad_dim_2 = ctypes.c_long(invmetgrad.shape[1])
        
        # Setting up "invmethess"
        if (invmethess is None):
            invmethess = numpy.zeros(shape=(3, 3), dtype=ctypes.c_double, order='F')
        elif ((not issubclass(type(invmethess), numpy.ndarray)) or
              (not numpy.asarray(invmethess).flags.f_contiguous) or
              (not (invmethess.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'invmethess' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            invmethess = numpy.asarray(invmethess, dtype=ctypes.c_double, order='F')
        invmethess_dim_1 = ctypes.c_long(invmethess.shape[0])
        invmethess_dim_2 = ctypes.c_long(invmethess.shape[1])
        
        # Setting up "neckval"
        neckval = ctypes.c_double()
        
        # Setting up "isoutofbounds"
        isoutofbounds = ctypes.c_int()
    
        # Call C-accessible Fortran wrapper.
        clib.c_interp_wrapper(ctypes.byref(coords_dim_1), ctypes.c_void_p(coords.ctypes.data), ctypes.byref(eneg), ctypes.byref(eneggrad_dim_1), ctypes.c_void_p(eneggrad.ctypes.data), ctypes.byref(eneghess_dim_1), ctypes.c_void_p(eneghess.ctypes.data), ctypes.byref(invmet_dim_1), ctypes.c_void_p(invmet.ctypes.data), ctypes.byref(invmetgrad_dim_1), ctypes.byref(invmetgrad_dim_2), ctypes.c_void_p(invmetgrad.ctypes.data), ctypes.byref(invmethess_dim_1), ctypes.byref(invmethess_dim_2), ctypes.c_void_p(invmethess.ctypes.data), ctypes.byref(neckval), ctypes.byref(isoutofbounds))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return eneg.value, eneggrad, eneghess, invmet, invmetgrad, invmethess, neckval.value, isoutofbounds.value

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine TIME_EVOLVE_OLD
    
    def time_evolve_old(self, coordsin, allcoords=None, allmomenta=None):
        '''!==================================================================================
!
!   Time evolves one trajectory. Should take in:
!       coords - the starting coords
!       saveToFile - whether or not to save this run
!       fileName - if save, need this for parallel i/o (may be ID instead)
!       dsetName - if save, need this to specify the dataset number
!
!agrees with Daniel's python code when ran with constants instead of random numbers
!
!=================================================================================='''
        
        # Setting up "coordsin"
        if ((not issubclass(type(coordsin), numpy.ndarray)) or
            (not numpy.asarray(coordsin).flags.f_contiguous) or
            (not (coordsin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'coordsin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            coordsin = numpy.asarray(coordsin, dtype=ctypes.c_double, order='F')
        coordsin_dim_1 = ctypes.c_long(coordsin.shape[0])
        
        # Setting up "allcoords"
        if (allcoords is None):
            allcoords = numpy.zeros(shape=(self.maxiters, 2), dtype=ctypes.c_double, order='F')
        elif ((not issubclass(type(allcoords), numpy.ndarray)) or
              (not numpy.asarray(allcoords).flags.f_contiguous) or
              (not (allcoords.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'allcoords' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            allcoords = numpy.asarray(allcoords, dtype=ctypes.c_double, order='F')
        allcoords_dim_1 = ctypes.c_long(allcoords.shape[0])
        allcoords_dim_2 = ctypes.c_long(allcoords.shape[1])
        
        # Setting up "allmomenta"
        if (allmomenta is None):
            allmomenta = numpy.zeros(shape=(self.maxiters, 2), dtype=ctypes.c_double, order='F')
        elif ((not issubclass(type(allmomenta), numpy.ndarray)) or
              (not numpy.asarray(allmomenta).flags.f_contiguous) or
              (not (allmomenta.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'allmomenta' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            allmomenta = numpy.asarray(allmomenta, dtype=ctypes.c_double, order='F')
        allmomenta_dim_1 = ctypes.c_long(allmomenta.shape[0])
        allmomenta_dim_2 = ctypes.c_long(allmomenta.shape[1])
        
        # Setting up "lastiter"
        lastiter = ctypes.c_int()
        
        # Setting up "finishedsuccessfully"
        finishedsuccessfully = ctypes.c_int()
    
        # Call C-accessible Fortran wrapper.
        clib.c_time_evolve_old(ctypes.byref(coordsin_dim_1), ctypes.c_void_p(coordsin.ctypes.data), ctypes.byref(allcoords_dim_1), ctypes.byref(allcoords_dim_2), ctypes.c_void_p(allcoords.ctypes.data), ctypes.byref(allmomenta_dim_1), ctypes.byref(allmomenta_dim_2), ctypes.c_void_p(allmomenta.ctypes.data), ctypes.byref(lastiter), ctypes.byref(finishedsuccessfully))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return allcoords, allmomenta, lastiter.value, finishedsuccessfully.value

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine TIME_EVOLVE
    
    def time_evolve(self, coordsin, allcoords=None, allmomenta=None):
        '''!==================================================================================
!
!   Time evolves one trajectory. Should take in:
!       coords - the starting coords
!       saveToFile - whether or not to save this run
!       fileName - if save, need this for parallel i/o (may be ID instead)
!       dsetName - if save, need this to specify the dataset number
!
!agrees with Daniel's python code when ran with constants instead of random numbers
!
!=================================================================================='''
        
        # Setting up "coordsin"
        if ((not issubclass(type(coordsin), numpy.ndarray)) or
            (not numpy.asarray(coordsin).flags.f_contiguous) or
            (not (coordsin.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'coordsin' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            coordsin = numpy.asarray(coordsin, dtype=ctypes.c_double, order='F')
        coordsin_dim_1 = ctypes.c_long(coordsin.shape[0])
        
        # Setting up "allcoords"
        if (allcoords is None):
            allcoords = numpy.zeros(shape=(self.maxiters, 2), dtype=ctypes.c_double, order='F')
        elif ((not issubclass(type(allcoords), numpy.ndarray)) or
              (not numpy.asarray(allcoords).flags.f_contiguous) or
              (not (allcoords.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'allcoords' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            allcoords = numpy.asarray(allcoords, dtype=ctypes.c_double, order='F')
        allcoords_dim_1 = ctypes.c_long(allcoords.shape[0])
        allcoords_dim_2 = ctypes.c_long(allcoords.shape[1])
        
        # Setting up "allmomenta"
        if (allmomenta is None):
            allmomenta = numpy.zeros(shape=(self.maxiters, 2), dtype=ctypes.c_double, order='F')
        elif ((not issubclass(type(allmomenta), numpy.ndarray)) or
              (not numpy.asarray(allmomenta).flags.f_contiguous) or
              (not (allmomenta.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'allmomenta' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            allmomenta = numpy.asarray(allmomenta, dtype=ctypes.c_double, order='F')
        allmomenta_dim_1 = ctypes.c_long(allmomenta.shape[0])
        allmomenta_dim_2 = ctypes.c_long(allmomenta.shape[1])
        
        # Setting up "lastiter"
        lastiter = ctypes.c_int()
        
        # Setting up "finishedsuccessfully"
        finishedsuccessfully = ctypes.c_int()
    
        # Call C-accessible Fortran wrapper.
        clib.c_time_evolve(ctypes.byref(coordsin_dim_1), ctypes.c_void_p(coordsin.ctypes.data), ctypes.byref(allcoords_dim_1), ctypes.byref(allcoords_dim_2), ctypes.c_void_p(allcoords.ctypes.data), ctypes.byref(allmomenta_dim_1), ctypes.byref(allmomenta_dim_2), ctypes.c_void_p(allmomenta.ctypes.data), ctypes.byref(lastiter), ctypes.byref(finishedsuccessfully))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return allcoords, allmomenta, lastiter.value, finishedsuccessfully.value

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine RUN
    
    def run(self, initialcoords, savepaths):
        ''''''
        
        # Setting up "initialcoords"
        if ((not issubclass(type(initialcoords), numpy.ndarray)) or
            (not numpy.asarray(initialcoords).flags.f_contiguous) or
            (not (initialcoords.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'initialcoords' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            initialcoords = numpy.asarray(initialcoords, dtype=ctypes.c_double, order='F')
        initialcoords_dim_1 = ctypes.c_long(initialcoords.shape[0])
        initialcoords_dim_2 = ctypes.c_long(initialcoords.shape[1])
        
        # Setting up "savepaths"
        if (type(savepaths) is not ctypes.c_int): savepaths = ctypes.c_int(savepaths)
        
        # Setting up "allcoords"
        allcoords = ctypes.c_void_p()
        allcoords_dim_1 = ctypes.c_long()
        allcoords_dim_2 = ctypes.c_long()
        allcoords_dim_3 = ctypes.c_long()
        
        # Setting up "allmomenta"
        allmomenta = ctypes.c_void_p()
        allmomenta_dim_1 = ctypes.c_long()
        allmomenta_dim_2 = ctypes.c_long()
        allmomenta_dim_3 = ctypes.c_long()
        
        # Setting up "lastiter"
        lastiter = ctypes.c_void_p()
        lastiter_dim_1 = ctypes.c_long()
        
        # Setting up "finishedsuccessfully"
        finishedsuccessfully = ctypes.c_void_p()
        finishedsuccessfully_dim_1 = ctypes.c_long()
    
        # Call C-accessible Fortran wrapper.
        clib.c_run(ctypes.byref(initialcoords_dim_1), ctypes.byref(initialcoords_dim_2), ctypes.c_void_p(initialcoords.ctypes.data), ctypes.byref(savepaths), ctypes.byref(allcoords_dim_1), ctypes.byref(allcoords_dim_2), ctypes.byref(allcoords_dim_3), ctypes.byref(allcoords), ctypes.byref(allmomenta_dim_1), ctypes.byref(allmomenta_dim_2), ctypes.byref(allmomenta_dim_3), ctypes.byref(allmomenta), ctypes.byref(lastiter_dim_1), ctypes.byref(lastiter), ctypes.byref(finishedsuccessfully_dim_1), ctypes.byref(finishedsuccessfully))
    
        # Post-processing "allcoords"
        allcoords_size = (allcoords_dim_1.value) * (allcoords_dim_2.value) * (allcoords_dim_3.value)
        if (allcoords_size > 0):
            allcoords = numpy.array(ctypes.cast(allcoords, ctypes.POINTER(ctypes.c_double*allcoords_size)).contents, copy=False)
            allcoords = allcoords.reshape(allcoords_dim_3.value,allcoords_dim_2.value,allcoords_dim_1.value).T
        elif (allcoords_size == 0):
            allcoords = numpy.zeros(shape=(allcoords_dim_3.value,allcoords_dim_2.value,allcoords_dim_1.value), dtype=ctypes.c_double, order='F')
        else:
            allcoords = None
        
        # Post-processing "allmomenta"
        allmomenta_size = (allmomenta_dim_1.value) * (allmomenta_dim_2.value) * (allmomenta_dim_3.value)
        if (allmomenta_size > 0):
            allmomenta = numpy.array(ctypes.cast(allmomenta, ctypes.POINTER(ctypes.c_double*allmomenta_size)).contents, copy=False)
            allmomenta = allmomenta.reshape(allmomenta_dim_3.value,allmomenta_dim_2.value,allmomenta_dim_1.value).T
        elif (allmomenta_size == 0):
            allmomenta = numpy.zeros(shape=(allmomenta_dim_3.value,allmomenta_dim_2.value,allmomenta_dim_1.value), dtype=ctypes.c_double, order='F')
        else:
            allmomenta = None
        
        # Post-processing "lastiter"
        lastiter_size = (lastiter_dim_1.value)
        if (lastiter_size > 0):
            lastiter = numpy.array(ctypes.cast(lastiter, ctypes.POINTER(ctypes.c_int*lastiter_size)).contents, copy=False)
        elif (lastiter_size == 0):
            lastiter = numpy.zeros(shape=(lastiter_dim_1.value), dtype=ctypes.c_int, order='F')
        else:
            lastiter = None
        
        # Post-processing "finishedsuccessfully"
        finishedsuccessfully_size = (finishedsuccessfully_dim_1.value)
        if (finishedsuccessfully_size > 0):
            finishedsuccessfully = numpy.array(ctypes.cast(finishedsuccessfully, ctypes.POINTER(ctypes.c_int*finishedsuccessfully_size)).contents, copy=False)
        elif (finishedsuccessfully_size == 0):
            finishedsuccessfully = numpy.zeros(shape=(finishedsuccessfully_dim_1.value), dtype=ctypes.c_int, order='F')
        else:
            finishedsuccessfully = None
        
        # Return final results, 'INTENT(OUT)' arguments only.
        return allcoords, allmomenta, lastiter, finishedsuccessfully

langevin_2d = langevin_2d()

