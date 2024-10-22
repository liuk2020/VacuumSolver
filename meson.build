project('vacuum', 'c')


py3 = import('python').find_installation(pure: false)
py3_dep = py3.dependency()


add_languages('fortran')

fortran_compiler = meson.get_compiler('fortran')
f90flags = []
f90flags += '-ffree-line-length-none'
f90flags += '-O3'

incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()
inc_np = include_directories(incdir_numpy)

incdir_f2py = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

# Share this object across multiple modules.
fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  dependencies: py3_dep,
  include_directories: [inc_np, inc_f2py],
)
fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [inc_np, inc_f2py],
)


fortran_files = files(
    'vacuum/fortran/pyvariables.f90', 
    'vacuum/fortran/basefunction.f90'
)
translated = custom_target('vacuum_fortran',
  input : fortran_files,
  output: 'vacuum_fortranmodule.c',
  command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'vacuum_fortran'],
  install : true,
  install_dir : py3.get_install_dir()
)

py3.extension_module('vacuum_fortran', 
                     fortran_files, translated,
                     dependencies: fortranobject_dep,
                     install: true,
                     fortran_args: f90flags)


py3.install_sources(
    'vacuum/__init__.py',
    'vacuum/vacuumProblem.py',
    subdir: 'vacuum'
)