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

import math, os
import time
import sys
from sys import argv
import subprocess
import glob
import Utils


print_only = False
use_shell = True
if sys.platform == "win32" :
	use_shell = False

def assert_exists(file):
	if not(os.path.exists(file)):
		print("No such file %s. Perhaps you need to run a previous stage."%file)
		sys.exit(-1)

if __name__ == '__main__' :


# Check that we got the right arguments......
	if len(argv) < 2 or not(os.path.exists(argv[1])):
		print('Usage: %s [model_config] [binary_path]' % argv[0]);
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
	Utils.current_stage = 5

	Utils.delete_complete_log()
	Utils.delete_error_log()
	Utils.rec_running_log()
	
	p_ptcl_path = '.%sseeds' % (os.sep)
	dummy, trans_name_root = os.path.split(model_input_file)
	trans_name_root = trans_name_root[:-5]
	transform_fname = os.path.join(model_output_path,"%s_pad_transform.tf" % trans_name_root)
	
	if not os.path.exists(p_ptcl_path) :
		os.makedirs(p_ptcl_path)

	gen_seeds_cmd = r'"%s" "%s" "%s"' % (os.path.join(binary_path,"SeedIntersections"), p_ptcl_path, transform_fname)

	start_time = time.time()

	for n in mat_names :
		# convert the crossing vols into NRRDs
		n = Utils.extract_root_matname(n)
		vol_fname = "%s_crossing.nrrd" % n
		assert_exists(vol_fname)

		# add the crossing NRRD to the argument list for 
		# SeedIntersections
		gen_seeds_cmd += " %s_crossing.nrrd" % n 

	print(gen_seeds_cmd)
	Utils.do_system(gen_seeds_cmd,print_only,use_shell)

	os.chdir(p_ptcl_path)

	join_fields_cmmd = r'"%s" -indexdata seeds-all.pc.fld ' % os.path.join(binary_path,"JoinFields")
	for f in glob.glob("*.ptcl"):
		(file_base,file_ext) = os.path.splitext(f)
		field_name = file_base + ".pc.fld"
		ptcl_to_field_cmmd = r'"%s" -nodata %s %s' % (os.path.join(binary_path,"PtclToField"), f, field_name)
		print(ptcl_to_field_cmmd)
		Utils.do_system(ptcl_to_field_cmmd,print_only,use_shell)
		join_fields_cmmd += field_name + " "

	print(join_fields_cmmd)
	Utils.do_system(join_fields_cmmd,print_only,use_shell)

	os.environ['GENERATE_SEEDS_OUTPUT'] = '%s/seeds-all.pc.fld' % os.getcwd()

	os.chdir("..")

	stop_time = time.time()
	f = open("generate-seeds-runtime.txt","w")
	f.write("%lf" % (stop_time-start_time))
	f.close()

	Utils.rec_running_log()
	os.chdir(cur_dir)
