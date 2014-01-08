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

import os
import sys
from sys import argv
import subprocess
import glob
import time
import Utils

print_only = False
use_shell = True
if sys.platform == "win32" :
	use_shell = False

combined_prefix = "particle-union"

def assert_exists(file):
	if not(os.path.exists(file)):
		print("No such file %s. Perhaps you need to run a previous stage."%file)
		sys.exit(-1)

def make_nodes_from_dir(d) :
	filename = os.path.join(d,combined_prefix+".pts")
	if os.path.exists(filename) :
		os.remove(filename)

	# Basically does the same thing as 'cat'
	pts_file = open(filename, 'w')
	for f in glob.glob(os.path.join(d,"particle_params*.pts")):
		f = open(f, 'r')
		pts_file.writelines(f.readlines())
		f.close()
	pts_file.close()

	path_prefix = os.path.join(d,combined_prefix)
	cmmd = r'"%s" %s.pts %s.node' % (os.path.join(binary_path,"pts2node"), path_prefix, path_prefix)
	Utils.do_system(cmmd,print_only,use_shell)

# must be called after combined is created.
def make_surfs(d) :

	assert_exists("medial_axis_param_file.txt")

	path_prefix = os.path.join(d,combined_prefix)
	i = 0
	for nm in mat_names :
		nm = Utils.extract_root_matname(nm)
		path_nm = os.path.join(d,nm)
		cmmd = r'"%s" %s.1 medial_axis_param_file.txt %s.m %d' % (os.path.join(binary_path,"removetets"), path_prefix, path_nm, i)
		#print(cmmd)
		Utils.do_system(cmmd,print_only,use_shell)

		cmmd = r'"%s" %s.m %s.ts.fld' % (os.path.join(binary_path,"MtoTriSurfField"), path_nm, path_nm)
		#print(cmmd)
		Utils.do_system(cmmd,print_only,use_shell)
		i = i + 1


def make_combined_surf(d) :

	assert_exists("medial_axis_param_file.txt")
	path_prefix = os.path.join(d,combined_prefix)
	cmmd = r'"%s" %s.node' % (os.path.join(binary_path,"tetgen"), path_prefix)
	#print(cmmd)
	Utils.do_system(cmmd,print_only,use_shell)
  
	cmmd = r'"%s" %s.1 medial_axis_param_file.txt %s.m 0' % (os.path.join(binary_path,"removetetsall"), path_prefix, path_prefix)
	#print(cmmd)
	Utils.do_system(cmmd,print_only,use_shell)

	cmmd = r'"%s" %s.m %s.ts.fld' % (os.path.join(binary_path,"MtoTriSurfField"), path_prefix, path_prefix)
	#print(cmmd)
	Utils.do_system(cmmd,print_only,use_shell)
	
if __name__ == '__main__' :


# Check that we got the right arguments......
	if len(argv) < 2 or not(os.path.exists(argv[1])):
		print('Usage: %s [model_config] [binary_path]' % argv[0])
		exit(1)

	model_config = argv[1]
	exec(open(model_config).read())
	model_path, dummy = os.path.split(model_config)
	
	if not(os.path.isabs(model_input_file)) :
	  model_input_file = os.path.normpath(os.path.join(model_path, model_input_file))

	if not(os.path.isabs(model_output_path)) :
	  model_output_path = os.path.normpath(os.path.join(model_path, model_output_path))
   
	if len(argv) > 2:
		binary_path = argv[2]
	else:
		binary_path = ""  

# Done checking arguments

	cur_dir = os.getcwd();
	os.chdir(model_output_path);

	Utils.output_path = model_output_path
	Utils.current_stage = 7

	Utils.delete_complete_log()
	Utils.delete_error_log()
	Utils.rec_running_log()
	
	start_time = time.time()

	dir_name = 'junctions'
	make_nodes_from_dir(dir_name) 
	make_combined_surf(dir_name) 
#    make_surfs(dir_name) 

	stop_time = time.time()
	f = open("build-material-interface-mesh-runtime.txt","w")
	f.write("%lf" % (stop_time-start_time))
	f.close()

	Utils.rec_completed_log()
	os.chdir(cur_dir)
