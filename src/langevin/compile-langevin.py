import fmodpy

# kwargs = {'-I':'/usr/include/hdf5/serial'}
# hdf5Dir = '/usr/include/hdf5/serial'
# test = [os.path.join(hdf5Dir,f) for f in os.listdir(hdf5Dir) if f.endswith('.mod')]
fmodpy.configure(rebuild=True)
code = fmodpy.fimport('langevin-python-module.f90',
                      name='langevin_module',
                      dependencies=['normal.f90',],#hdf5Dir+'/hdf5.mod'],#+test,
                      show_warnings=True,
                      # f_compiler_args='-I/usr/include/hdf5/serial',
                      )

# code = fmodpy.fimport('langevin_2d.mod',name='langevin',
#                       f_compiler_args='-I\ /usr/include/hdf5/serial',
#                       dependencies=['normal.f90',],
#                       show_warnings=True,)
# print(code.__dict__.keys())

