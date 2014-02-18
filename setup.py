#!/usr/bin/env python

"""
setup.py file for nlab methods
"""
import numpy
from distutils.core import setup, Extension
import platform as plt


f = plt.system()

if (f =='Darwin'):
	print 'OSX Identitfied'
	ext = ['']
	lib = 0
	swg = ['-c++', '-I./cpp/', '-Wall']
	incd= [numpy.get_include(), './cpp/']
	lbdir = 0
	xtc = ['-O3']
else:
	print 'Assuming Linux'
	ext = 0
	lib = 0
	swg = ['-c++',  '-I./cpp/']
	incd= [numpy.get_include(), './cpp/']
	lbdir = 0
	xtc  = ['-fpermissive', '-O3']

nlab_module = Extension('_nlab',
                           sources=['cpp/Neuron.cxx',
                           			'cpp/Neuron_ML.cxx',
                           			'cpp/Neuron_HH.cxx',
                           			'cpp/Neuron_Osc.cxx',
						   			'cpp/nlab.cxx',
						   			'cpp/ConMat.cxx',
									'nlab.i'],
                           include_dirs = incd,
                           swig_opts=swg,
						   libraries = lib,
                           extra_link_args = ext,
						   extra_compile_args = xtc,
						   library_dirs=lbdir,
                           )


setup (name = 'nlab',
       version = '0.1',
       author      = "FHK",
       description = """SWIG wrapping of nlab""",
       ext_modules = [nlab_module],
       py_modules = ["nlab"],
       )
