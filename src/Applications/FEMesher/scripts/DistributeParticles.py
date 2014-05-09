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
import shutil
import pdb

if sys.version_info[0] < 3 :
	from thread import start_new_thread
else :
	from _thread import start_new_thread

from math import fabs
import Utils


print_only = False
use_shell = True
if sys.platform == "win32" :
	use_shell = False

def dump_log(p) :
	if print_only :
		return

	# Warning: The data read is buffered in memory, so do not use this method if the data size is large or unlimited.
	# See Python docs, section 18.1.2.
	(stdoutdata, stderrdata) = p.communicate()
	f = open("%d.log" % p.pid, "w")
	f.writelines(stdoutdata.decode())
	f.close()

def wait_on_procs(procs) :    

	#once all tight procs are done, make the vols    
	done = print_only
	count = 0;
	while not done :
		time.sleep(1)
		Utils.rec_running_log()
		ndone = 0
		#all tight procs need to be finished before moving on.
		np = 0
		for p in procs :
			if p.poll() != None : #None means process not finished.
				ndone = ndone + 1

		if ndone == len(procs) :
			done = True
		else :
			if count % 30 == 0 :
				#print "%d processes are complete of %d" % (ndone, len(procs))
				print("%d processes are complete of %d pid:%d" % (ndone, len(procs), p.pid))
		count = count + 1

def start_proc(optimize_ps_cmmd,procs,max_procs,n_proc) :
	#print("Start Proc: len(procs:%d Cmmd:%s\n"% (len(procs),optimize_ps_cmmd))
	if (len(procs) >= max_procs) :
		done = False
		while not done :
			time.sleep(1)
			for p in procs :
				if p.poll() != None :
					done = True
					procs.remove(p)
					break
	procstdout = "procs" + str(n_proc) + ".log"
	optimize_ps_cmmd = optimize_ps_cmmd + " > " + procstdout
	print("**** parallel process console output: " + optimize_ps_cmmd)
	#pdb.set_trace()
	proc = subprocess.Popen(optimize_ps_cmmd,  shell=True, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=model_output_path)

	Utils.rec_pid(proc.pid)
	procs.append(proc)
	# the output from processes need to be dumped or the
	# buffer fills up and it sleeps.
	#print("starting thread to dump")
	#print(proc)
	start_new_thread(dump_log, (proc, ))
	

def assert_exists(file):
	if not(os.path.exists(file)):
		print("No such file %s. Perhaps you need to run a previous stage."%file)
		sys.exit(-1)

def gen_isect_file(fname) :
	global g_num_isects
	global g_nquads
	global g_ntrips
	global g_ndoubs
	global g_quad_tri_juncs
	global g_quad_tri_seeds

	#print("test: in gen_isect_file ", fname)

	g_num_isects = 0
	g_quad_tri_juncs = []
	g_quad_tri_seeds = []

	seed_path = "seeds"
	if not os.path.exists(seed_path):
		print("Path %s does not exist.  Perhaps you need to run a previous stage." %seed_path)
		sys.exit(-1)

	input_seed_path = "input_seeds"
	if not os.path.exists(input_seed_path):
		os.mkdir(input_seed_path)

	isect_file_text = ""

	nmats = len(mats)

	# first write out the quad junctions to be optimized
	g_nquads = 0
	quad_text = ""

	for i in range(0,nmats-3) :
		for j in range(i+1,nmats-2) :
			for k in range(j+1,nmats-1) :
				for l in range(k+1,nmats) :
					seed_file = "%s_%s_%s_%s_%s_seed.ptcl" % (os.path.join(seed_path,"q"), i, j, k, l)
					if os.path.exists(seed_file):
						input_seed_file = "%s_%s_%s_%s_%s_seed.ptcl" % (os.path.join(input_seed_path,"q"), i, j, k, l)
						shutil.copyfile(seed_file, input_seed_file)
						g_quad_tri_seeds.append(input_seed_file)
						g_quad_tri_juncs.append(os.path.join("junctions","particle_params.txt_%s.ptcl" % g_nquads))
						quad_text += "%s %s %s %s y %s y\n" % (i, j, k, l, input_seed_file)
						g_nquads+=1

	isect_file_text += "%s\n" % g_nquads + quad_text

	# second write out the triple junctions to be optimized
	g_ntrips = 0
	trip_text = ""

	for i in range(0,nmats-2) :
		for j in range(i+1,nmats-1) :
			for k in range(j+1,nmats) :
				seed_file = "%s_%s_%s_%s_seed.ptcl" % (os.path.join(seed_path,"t"),i,j,k)
				if os.path.exists(seed_file):
					input_seed_file = "%s_%s_%s_%s_seed.ptcl" % (os.path.join(input_seed_path,"t"), i, j, k)
					shutil.copyfile(seed_file, input_seed_file)
					g_quad_tri_seeds.append(input_seed_file)
					g_quad_tri_juncs.append(os.path.join("junctions","particle_params.txt_%s.ptcl" % (g_nquads+g_ntrips)))
					trip_text += "%s %s %s y %s y\n" % (i, j, k, input_seed_file)
					g_ntrips+=1

	isect_file_text += "%s\n" % g_ntrips + trip_text

	# finally write out the double junctions
	g_ndoubs = 0
	doub_text = ""

	for i in range(0,nmats-1) :
		for j in range(i+1,nmats) :
			seed_file = "%s_%s_%s_seed.ptcl" % (os.path.join(seed_path,"d"), i, j)
			if os.path.exists(seed_file):
				input_seed_file = "%s_%s_%s_seed.ptcl" % (os.path.join(input_seed_path,"d"), i, j)
				shutil.copyfile(seed_file, input_seed_file)
				g_ndoubs+=1
				doub_text += "%s %s y %s y\n" % (i, j, input_seed_file)

	isect_file_text += "%s\n" % g_ndoubs + doub_text

	f = open(fname, 'w')
	f.write(isect_file_text)
	f.close()

	g_num_isects = g_nquads+g_ntrips+g_ndoubs 

if __name__ == "__main__" :

# Check that we got the right arguments......
	if len(argv) < 2 or not(os.path.exists(argv[1])):
		print('Usage: %s [model_config] [binary_path]' % argv[0]);
		exit(1)

	MAX_SIZING_FIELD = ""
	SIZING_SCALE_VAR = ""

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

	try:
		max_procs
	except NameError:
		max_procs = 4

# Done checking arguments

	cur_dir = os.getcwd();
	os.chdir(model_output_path);

	Utils.output_path = model_output_path
	Utils.current_stage = 6

	Utils.delete_complete_log()
	Utils.delete_error_log()
	Utils.rec_running_log()
	
	start_time = time.time()

	procs = []
	dnm = "junctions"

	cwd = os.getcwd()
	if not os.path.exists(dnm) :
		os.mkdir(dnm)

	gen_isect_file( os.path.join(dnm,"m1.txt") )

	idx = 0
	lns = ["%d\n" % len(mats)]
	for i in mats :
		n = Utils.extract_root_matname(mat_names[idx])
		lns.append("%s/%s\n" % (cwd, n))
		idx = idx + 1

	f = open(os.path.join(dnm,'particle_params.txt'), 'w')
	f.writelines(lns)
	f.close()

	if MAX_SIZING_FIELD != "" :
		max_sizing_field = MAX_SIZING_FIELD
	
	#pip was not passed down to optimize-particles-system
	if SIZING_SCALE_VAR != "" :
		sizing_scale_var = SIZING_SCALE_VAR
	else:
		sizing_scale_var = 1.0

	if ROI_X != "" :
		roi_x = ROI_X
	else:
		roi_x = 0.5
	
	if ROI_Y != "" :
		roi_y = ROI_Y
	else:
		roi_y = 0.5
		
	if ROI_Z != "" :
		roi_z = ROI_Z
	else:
		roi_z = 0.5

	psys_txt = [
			"ENERGY                  radial\n",
			"NUMBER_OF_SURFACES      %d\n" % len(mats),
			"NUMBER_OF_INTERSECTIONS %d\n" % g_num_isects,
			"INTERSECTION_FILE       %s/%s/m1.txt\n" % (cwd, dnm),
			"BASE_FILE_NAME          %s/%s/particle_params.txt\n" % (cwd,dnm),
			"INIT_NUM_POINTS         5\n", #unused
			"MAX_SF                  %f\n" % max_sizing_field,
			"SIZING_SCALE            %f\n" % sizing_scale_var,
			"ROI_X           		 %f\n" % roi_x,
			"ROI_Y            		 %f\n" % roi_y,
			"ROI_Z            		 %f\n" % roi_z,
			]

	f = open(os.path.join(dnm,'psystem_input.txt'), 'w')
	f.writelines(psys_txt)
	f.close()

	# now run rendersf3d on the resulting junctions directory
	optimize_ps_cmmd = r'"%s" %s %s 0 %s' % (os.path.join(binary_path,"optimize-particle-system"), os.path.join("junctions","psystem_input.txt"), num_particle_iters, g_nquads+g_ntrips-1)
	print("**** running command: %s\n" % optimize_ps_cmmd)
	Utils.do_system(optimize_ps_cmmd,print_only,use_shell)

	# Now copy the junctions results for the quads and triples back to the
	# input_seed_path files (but need to map names!)
	# And then fork a bunch of processes -- once for each double
	# and call optimize_ps_cmd just for that interface #

	print("**** copying quad and tri junctions back into seeds folder\n")
	for i in range(0, g_nquads+g_ntrips) :
		shutil.copyfile(g_quad_tri_juncs[i], g_quad_tri_seeds[i])
	
	procs = []

	for i in range(g_nquads+g_ntrips, g_nquads+g_ntrips+g_ndoubs) :
		# We constrain the maximum number of processes that can be 
		# created at once

		optimize_ps_cmmd = r'"%s" %s %s %s %s ' % (os.path.join(binary_path,"optimize-particle-system"), os.path.join("junctions","psystem_input.txt"), num_particle_iters, i, i)
		
		print("**** parallel process %s of %s: %s\n"%(i-g_nquads-g_ntrips+1,g_ndoubs,optimize_ps_cmmd))
		
		start_proc(optimize_ps_cmmd,procs,max_procs,i-g_nquads-g_ntrips+1)

		

	wait_on_procs(procs)

	# change directory to junctions and write out files for visualization
	os.chdir("junctions")
	print("**** Change directory")

	ptcl_path = os.path.join(binary_path,"PtclToField")

	for f in (glob.glob('./particle_params.txt_*.ptcl')):
		(params_base, params_ext) = os.path.splitext(f)
		ptcl_data_cmmd = r'"%s" %s %s.pcv.fld' % (ptcl_path, f, params_base)
		Utils.do_system(ptcl_data_cmmd,print_only,use_shell)
		print("**** doing command %s" % ptcl_data_cmmd)
		ptcl_nodata_cmmd = r'"%s" -nodata %s %s.pc.fld' % (ptcl_path, f, params_base)
		Utils.do_system(ptcl_nodata_cmmd,print_only,use_shell)
		print("**** doing command %s" % ptcl_nodata_cmmd)

	if os.path.exists("particle-all.pc.fld"):
		os.remove("particle-all.pc.fld")
		
	join_fields_path = 	os.path.join(binary_path,"JoinFields")
	join_fields_nodata_cmmd = r'"%s" -indexdata particle-all.pc.fld %s' % (join_fields_path, " ".join(glob.glob('./particle_params.txt_*.pc.fld')))
	Utils.do_system(join_fields_nodata_cmmd,print_only,use_shell)

	if os.path.exists("particle-all.pcv.fld"):
		os.remove("particle-all.pcv.fld")
	join_fields_data_cmmd = r'"%s" particle-all.pcv.fld %s' % (join_fields_path, " ".join(glob.glob('./particle_params.txt_*.pcv.fld')))
	Utils.do_system(join_fields_data_cmmd,print_only,use_shell)

	os.chdir("..")

	stop_time = time.time()
	f = open("distribute-particles-runtime.txt","w")
	f.write("%lf" % (stop_time-start_time))
	f.close()

	Utils.rec_completed_log()
	os.chdir(cur_dir)
