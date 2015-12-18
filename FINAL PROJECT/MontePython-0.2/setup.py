#!/usr/bin/env python

import os,sys,shutil,string
from distutils.core import setup, Extension

prefix=sys.prefix
for i in sys.argv:
    if '--prefix' in i:
        prefix=i.split('=')[1]

def list_files(path):
    tmp = os.listdir(path)
    list = []
    for file in tmp:
        list += [string.join([path,file],'/')]
    return list


try:
    import Numeric
except:
    raise "Couldn't find Numeric module. Quiting..."
try:
    import pypar
except:
    raise "Couldn't find pypar module. Quiting..."

name = 'MontePython'    #module name
version = '0.2' #module's version
description = """MontePython: Python Diffusion Monte Carlo
pure estimators"""
author = 'Jon Kristian Nilsen' #me
author_email = 'j.k.nilsen@usit.uio.no'
url = ''

# lib directories:
swig_dir  = 'src/SWIG'
cpp_dir   = 'src/C++'
py_dir    = 'src/Python'
ex_dir    = 'src/example'
mpi_dir   = '/usr/local/include'
mpi_lib   = '/usr/local/lib'
pypar_lib = '/usr/lib/python2.4/site-packages/pypar'
# openmpi:
#mpi_libraries = ['mpi','mpi_cxx','pthread','orte','opal','util','nsl','dl',]
# mpich2:
mpi_libraries = ['mpich','mpichcxx','rt','pthread',]

use_sprng = False

if use_sprng:
    shutil.copy(os.path.join(cpp_dir,'random_sprng.cpp'),
                os.path.join(cpp_dir,'random.cpp'))
    shutil.copy(os.path.join(cpp_dir,'random_sprng.hpp'),
                os.path.join(cpp_dir,'random.hpp'))
    shutil.copy(os.path.join(swig_dir,'%s_cxx_sprng.i'%name),
                os.path.join(swig_dir,'%s_cxx.i'%name))
    sprng_dir = '/usr/include'
    sprng_lib = '/usr/lib'
else:
    shutil.copy(os.path.join(cpp_dir,'random_lecuyer.cpp'),
                os.path.join(cpp_dir,'random.cpp'))
    shutil.copy(os.path.join(cpp_dir,'random_lecuyer.hpp'),
                os.path.join(cpp_dir,'random.hpp'))
    shutil.copy(os.path.join(swig_dir,'%s_cxx_lecuyer.i'%name),
                os.path.join(swig_dir,'%s_cxx.i'%name))
    

# swig to be run manually, so
if use_sprng:
    swig_cmd = 'swig -python -c++ -I../C++ -I%s %s_cxx.i'%(sprng_dir,name)
else:
    swig_cmd = 'swig -python -c++ -I../C++ %s_cxx.i'%(name)
    
os.chdir(swig_dir)
print 'running SWIG:',swig_cmd
if os.path.exists('%s_cxx.py'%name):
    os.remove('%s_cxx.py'%name)
if os.path.exists('%s_cxx_wrap.cxx'%name):
    os.remove('%s_cxx_wrap.cxx'%name)
os.system(swig_cmd)
os.chdir(os.pardir)
os.chdir(os.pardir)

shutil.copy(swig_dir+'/'+name+'_cxx.py','.')

shutil.copy(py_dir+'/'+name+'.py','.')

shutil.copy(py_dir+'/WalkerArray.py','.')


sources = [cpp_dir+'/convert.cpp',cpp_dir+'/functions.cpp',cpp_dir+'/mc.cpp',\
           cpp_dir+'/particle.cpp',cpp_dir+'/QickArray.cpp',\
           cpp_dir+'/random.cpp',cpp_dir+'/walker.cpp',cpp_dir+'/warray.cpp',\
           swig_dir+'/'+name+'_cxx_wrap.cxx',]

if not use_sprng:
    sources += [cpp_dir+'/RngStream.cpp']

include_dirs = [swig_dir,cpp_dir,mpi_dir,]
if use_sprng:
    include_dirs += [sprng_dir]
library_dirs     = [mpi_lib]
libraries    = mpi_libraries+['stdc++','python2.4']

cpp_files  = list_files(cpp_dir)
swig_files = list_files(swig_dir)
ex_files   = list_files(ex_dir)

extra_objects = [os.path.join(pypar_lib,'mpiext.so')]

if use_sprng:
    extra_objects += [os.path.join(sprng_lib,'libsprng.a')]

setup(name = name, version = version, description = description,
      author = author, author_email = author_email, url = url,
      ext_modules = [Extension('_'+name+'_cxx', sources,
                               include_dirs = include_dirs,
                               library_dirs = library_dirs,
			       libraries=libraries,
                               extra_objects=extra_objects)],
      py_modules = ['MontePython','MontePython_cxx','WalkerArray'],
      data_files = [(ex_dir,ex_files),
                    (cpp_dir,cpp_files),
                    (swig_dir,swig_files)]
      )
