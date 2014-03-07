#!/usr/bin/env python
#  
#  For more information, please see: http://software.sci.utah.edu
#  
#  The MIT License
#  
#  Copyright (c) 2009 Scientific Computing and Imaging Institute,
#  University of Utah.
#  
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#  
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#  
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#  
#    Author : Martin Cole

import string
from sys import argv
import sys
import subprocess
import re
import os
import time
import glob
import fileinput
import Utils
import math

if sys.version_info[0] < 3 :
	from thread import start_new_thread
else :
	from _thread import start_new_thread

from math import fabs
from math import ceil
from math import log


print_only = False
use_shell = True
if sys.platform == "win32" :
	use_shell = False

# req: nrrd is unsigned char
def build_label_list(innrrd) :
	hist_prefix = os.path.splitext(innrrd)[0] + ".hist"
	tmp_nrrd = hist_prefix + ".nrrd"
	hist_nrrd =  hist_prefix + ".nhdr"
	# produce an ascii file where each line is the number of 
	unu_cmmd = os.path.join(binary_path,"unu")
	cmmd = unu_cmmd + ' histo -b 256 -i "%s" -o "%s"' % (innrrd,tmp_nrrd)
	Utils.do_system(cmmd,print_only,use_shell)
	cmmd = unu_cmmd + ' save -i "%s" -e ascii -f nrrd -o "%s"' % (tmp_nrrd,hist_nrrd)
	Utils.do_system(cmmd,print_only,use_shell)

	f = open(hist_prefix+".ascii")

	label_list = []
	label_index = 0
	label_occurences = 0
	for line in f:
		try:
			label_occurences = int(line)
		except ValueError:
			label_occurences = 0
		if label_occurences != 0:
			label_list.append(label_index)
		label_index+=1

	f.close()

	# cleanup
	os.remove(hist_prefix+".ascii")
	os.remove(tmp_nrrd)
	os.remove(hist_nrrd)

	return label_list;

def build_mat_list(label_list) :
	mat_list = []
	for label in label_list:
		mat_list.append("mat"+str(label))
	return mat_list

def get_dataset_dims(innrrd) :
	axsizes_re = re.compile("sizes: (\d+) (\d+) (\d+)")
	
	unu_cmmd = os.path.join(binary_path,"unu")
	cmmd = r'"%s" head "%s"' % (unu_cmmd, innrrd)
	
	#print("Cmmd: " + cmmd)
	
	p = subprocess.Popen(cmmd,  shell=use_shell, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	Utils.rec_pid(p.pid)
	stdoutdata = p.communicate()[0]
	
	#print(stdoutdata)

	dim_str = ""

	if (sys.version_info[0] < 3) :
		for l in string.split(stdoutdata, '\n') :
			mo = axsizes_re.match(l.decode())
			if mo != None :
				dim_str = mo.group()
	else :
		for l in stdoutdata.split(bytes('\n','ascii')) :
			mo = axsizes_re.match(l.decode())
			if mo != None :
				dim_str =  mo.group()

	dims = dim_str.split()
	return (dims[1],dims[2],dims[3])

# because this has been tested with our code, we will assume
# that the initial number of bins is 100 in each dimension
def compute_initial_medial_axis_bins(innrrd) :
	return 100

# we will assume that regardless of the number of bins 
# set initially, we want 4 times as many bins in each direction at the most
# refined level of the grid as there are voxels in that direction.
# that is:
# 2^l*init_bins >= 4*max[voxels_x,voxels_y,voxels_z]
# 2^(l-2)       >= max_vox_num / init_bins
# l-2           >= log2 max_vox_num - log2 init_bins
# l             >= 2 + log2 max_vox_num - log2 init_bins
# l             =  ceil(2 + log2 max_vox_num - log2 init_bins)
def compute_refinement_levels(innrrd,init_bins) :
	dims = get_dataset_dims(innrrd)
	dim_max = max(dims)

	refinement_levels = ceil(2.0 + log( float(dim_max),2.0 ) - log( float(init_bins),2.0 ))
	refinement_levels = max(refinement_levels,0)

	return refinement_levels

# by default, we will let the mat radii be a 80% voxel width
# this number is in voxel coordinates and is unaffected by
# spacing information.  Hence you should think of this as
# feature size in voxel space. 
def compute_mat_radii(innrrd) :
	return .8

def compute_max_sizing_field(innrrd) :
	return 3.0

# there doesn't seem to be any great way of determining
# this a priori.  Hence, we default to 100.
def compute_num_particle_iters(innrrd) :
	return 100

# for now, we will simply return the flags required for the
# proper functioning of stage 8.  It's possible that we will
# want to specify defaults for 'a' and 'q' as well.
def compute_tetgen_joined_vol_flags(innrrd) : 
	return "zpAA"

def output_model_config(label_nrrd,
						output_path,
						label_list,
						mat_list,
						mat_radii,
						initial_medial_axis_bins,
						refinement_levels,
						max_sizing_field,
						num_particle_iters,
						tetgen_joined_vol_flags):
	print("# The source label map.")
	label_nrrd = label_nrrd.replace('\\','/')
	print("model_input_file=\""+label_nrrd+"\"")
	print("# The directory where BioMesh3D will write its output.")
	output_path = output_path.replace('\\','/')
	print("model_output_path=\""+output_path+"\"")

	print("# The list of labels that you care about in the label map.")
	mats_str="mats=("
	for i in range(len(label_list)-1) :
		mats_str += "%s, " % label_list[i]
	mats_str += "%s)" % label_list[len(label_list)-1]
	print(mats_str)

	print("# A string identifying the labels in your label map:\n# e.g., 'air', 'skin', 'bone'.")
	mat_names_str="mat_names=("
	for i in range(len(mat_list)-1) :
		mat_names_str += "'%s', " % mat_list[i]
	mat_names_str += "'%s')" % mat_list[len(mat_list)-1]
	print(mat_names_str)

	print("# Minimum feature size in your data in voxel units.")
	print("mat_radii=%f" % mat_radii)
	print("# The number of bins at the lowest resolution level of the medial\n# axis multigrid.")
	print("initial_medial_axis_bins=%d"%initial_medial_axis_bins)
	print("# Depth of the hierarchy for the medial axis computation.  Larger\n# numbers mean longer run times, but better approximation.")
	print("refinement_levels=%d"%refinement_levels)
	print("# Sets a cap on sizing field.  The smaller the number, the more\n# dense the final mesh.")
	print("max_sizing_field=%f"%max_sizing_field)
	print("# Number of passes through the particle placement algorithm.  More\n# iterations mean better feature approximation but longer runtimes.")
	print("num_particle_iters=%d"%num_particle_iters)
	print("# Flags to pass to tetgen for the final pipeline stage.")
	print("tetgen_joined_vol_flags=\""+tetgen_joined_vol_flags+"\"")

if __name__ == "__main__" :

# Check that we got the right arguments......
	if len(argv) < 3 or not( os.path.exists( argv[1] ) ):
		print('Usage: %s [label_nrrd] [model_file]' % argv[0])
		sys.exit(1)

	if len(argv) > 3:
		main_bin_dir = argv[2]
	else:
		scripts_dir = os.path.abspath(os.path.dirname(argv[0]))
		main_bin_dir, dummy = os.path.split(scripts_dir) 
	
	binary_path = main_bin_dir    
	
	if sys.platform == "win32":
		unu_exec = "unu.exe"
	else:
		unu_exec = "unu"
		os.setpgrp()

	test_binary_path = os.path.join(main_bin_dir, unu_exec)
	if not(os.path.exists(test_binary_path)) :
		# are we in a development build directory?
		binary_path = os.path.join(main_bin_dir, "Release")
		test_binary_path = os.path.join(binary_path, unu_exec)
		if not(os.path.exists(test_binary_path)) :
			binary_path = os.path.join(main_bin_dir, "Debug")
			test_binary_path = os.path.join(binary_path, unu_exec)
			if not(os.path.exists(test_binary_path)) :
				binary_path="" 
	 

	label_nrrd = os.path.abspath(argv[1])
	model_file = os.path.abspath(argv[2])

	output_path = "."

	model_path, dummy = os.path.split(model_file)

	#print("label_nrrd " + label_nrrd)

	Utils.output_path = os.path.normpath(os.path.join(model_path,output_path))

	if not os.path.exists(Utils.output_path) :
		os.makedirs(Utils.output_path)

	label_list = build_label_list(label_nrrd)
	mat_list = build_mat_list(label_list)
	mat_radii = compute_mat_radii(label_nrrd)
	initial_medial_axis_bins = compute_initial_medial_axis_bins(label_nrrd)
	refinement_levels = compute_refinement_levels(label_nrrd,initial_medial_axis_bins)
	max_sizing_field = compute_max_sizing_field(label_nrrd)
	num_particle_iters = compute_num_particle_iters(label_nrrd)
	tetgen_joined_vol_flags = compute_tetgen_joined_vol_flags(label_nrrd)

	f = open(model_file,'w')
	
	
	saved_stdout = sys.stdout
	sys.stdout = f
	
	# send relative path to nrrd, so that relocating the data does not break
	# the pipeline.
	rel_label_nrrd = os.path.relpath(label_nrrd,model_path)

	output_model_config(rel_label_nrrd,
						output_path,
						label_list,
						mat_list,
						mat_radii,
						initial_medial_axis_bins,
						refinement_levels,
						max_sizing_field,
						num_particle_iters,
						tetgen_joined_vol_flags)

	sys.stdout = saved_stdout

	f.close()

	sys.exit(0)
