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
import re
import glob
import time
import Utils

combined_prefix = "particle-union"

print_only = False
use_shell = True
if sys.platform == "win32" :
	use_shell = False

def assert_exists(file):
	if not(os.path.exists(file)):
		print("No such file %s. Perhaps you need to run a previous stage."%file)
		sys.exit(-1)

# must be called after combined is created.
def make_vols(d) :

	i = 0
	tetgen_solo_vol_flags = tetgen_joined_vol_flags
	tet_fname = os.path.join(d,"tets-all.fld")
	tet_trans_fname = os.path.join(d,"tets-all_transformed.fld")
	join_fields_cmmd = r'"%s" -indexdata %s' % (os.path.join(binary_path,"JoinFields"), tet_fname)
	join_transformed_fields_cmmd = r'"%s" -indexdata %s' % (os.path.join(binary_path,"JoinFields"), tet_trans_fname)
	for nm in mat_names :
		nm = Utils.extract_root_matname(nm)
		tet_fname = "%s.tets.fld" % os.path.join(d, nm)
		mat_surf = "%s.ts.fld" % os.path.join(d,nm)
		assert_exists(mat_surf)
		cmmd = r'"%s" -cmmd_line %s -main %s -out %s' % (os.path.join(binary_path,"InterfaceWithTetGen"), tetgen_solo_vol_flags, mat_surf, tet_fname)
		#print(cmmd)
		Utils.do_system(cmmd,print_only,use_shell)
		transform_with_transform(tet_fname)  
		join_fields_cmmd += "%s.tets.fld " % os.path.join(d,nm)
		join_transformed_fields_cmmd += "%s.tets_transformed.fld " % os.path.join(d, nm)
		i = i + 1
	Utils.do_system(join_fields_cmmd,print_only,use_shell)
	Utils.do_system(join_transformed_fields_cmmd,print_only,use_shell)

	# go through and label all of the tets 
	labeled_tet_fname = os.path.join(d,"tets-all-labeled.fld")
	cmmd = r'"%s" %s %s %s' % (os.path.join(binary_path,"LabelTets"),os.path.join(model_output_path,"dominant_labelmap.nrrd"),tet_trans_fname,labeled_tet_fname)
	#print("%s\n" % cmmd)
	Utils.do_system(cmmd,print_only,use_shell)

def make_combined_vol(d) :

	combined_surf = "%s.ts.fld" % os.path.join(d,combined_prefix)
	assert_exists(combined_surf)
	tet_fname = "%s.tets.fld" % os.path.join(d, combined_prefix)
	cmmd = r'"%s" -cmmd_line %s -main %s -out %s' % (os.path.join(binary_path,"InterfaceWithTetGen"), tetgen_joined_vol_flags, combined_surf, tet_fname)
	#print(cmmd)
	Utils.do_system(cmmd,print_only,use_shell)

	# go through and label all of the tets 
	labeled_tet_fname = "%s.tets-labeled.fld" % os.path.join(d, combined_prefix)
	cmmd = r'"%s" "%s" %s %s' % (os.path.join(binary_path,"LabelTets"),os.path.join(model_output_path,"dominant_labelmap.nrrd"),tet_fname,labeled_tet_fname)

	Utils.do_system(cmmd,print_only,use_shell)

	transform_with_transform(labeled_tet_fname)

def transform_with_transform(tet_file) :

	assert_exists(tet_file)
	transformed_fname = "%s_transformed.fld" % tet_file[:-4]
	transf = "%s_pad_transform.tf" % model_input_file[:-5]
	dummy,filename = os.path.split(transf)
	transf = os.path.join(model_output_path,filename)
	print("Transforming %s" % tet_file)
	trans_cmmd = r'"%s" -input %s -transform "%s" -output %s' % (os.path.join(binary_path,"TransformFieldWithTransform"), tet_file, transf, transformed_fname)
	Utils.do_system(trans_cmmd,print_only,use_shell)
	print("Done Transforming")

def extract_string_matching_regex(src_str,regex_str) :
	p = re.compile(regex_str)
	m = p.search(src_str)
	str = ""
	if (m != None) :
		str = m.group()
	return str

# In the original BioMesh3D, the user was allowed complete 
# freedom in specifying the tetgen flags.  However, the current
# solution requires that certain flags be specified by the user,
# and if they are not specified, the code fails to work.  Hence,
# this function extracts the degrees of freedom that the user is
# allowed to specify, and ensures that the required flags are 
# present.

def correct_tetgen_flags(tet_gen_flags) :

	fp_regex = '[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?'
	int_regex = '[0-9]+'

	# we have to do a little something extra here as a & q can be specified
	# without arguments
	astr = extract_string_matching_regex(tet_gen_flags,'a('+fp_regex+')?')
	qstr = extract_string_matching_regex(tet_gen_flags,'q('+fp_regex+')?')
	Mstr = extract_string_matching_regex(tet_gen_flags,'M')
	Ystr = extract_string_matching_regex(tet_gen_flags,'Y')
	Sstr = extract_string_matching_regex(tet_gen_flags,'S'+int_regex)
	Tstr = extract_string_matching_regex(tet_gen_flags,'T'+fp_regex)

	# now compose a correct set of parameters    
	fixed_tet_gen_flags = "zpAA" + astr + qstr + Mstr + Ystr + Sstr + Tstr

	return fixed_tet_gen_flags

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
	Utils.current_stage = 8

	Utils.delete_complete_log()
	Utils.delete_error_log()
	Utils.rec_running_log()
	start_time = time.time()    

	dir_name = 'junctions'
	assert_exists(dir_name)

	tetgen_joined_vol_flags = correct_tetgen_flags(tetgen_joined_vol_flags)

	make_combined_vol(dir_name) 
#    make_vols(dir_name) 

	stop_time = time.time()
	f = open("build-volumetric-mesh-runtime.txt","w")
	f.write("%lf" % (stop_time-start_time))
	f.close()

	Utils.rec_completed_log()
	os.chdir(cur_dir)
