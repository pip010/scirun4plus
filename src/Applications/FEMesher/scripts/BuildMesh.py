#!/usr/bin/env python
#
#  For more information, please see: http://software.sci.utah.edu
# 
# The MIT License
#
# Copyright (c) 2009 Scientific Computing and Imaging Institute,
# University of Utah.
#
#
# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the 
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
#

import os
import sys
from sys import argv
import getopt
import subprocess
from subprocess import Popen
import re
import math
from math import sqrt
import time

import Utils

print_only = False
use_shell = True
if sys.platform == "win32" :
	use_shell = False
scirun_processes = []

def pause_for_scirun_cmd(scirun_arg,  pause) :
	global scirun_processes

	sci_run_cmmd = "%s %s" % (os.path.join(binary_path,"scirun"), scirun_arg)
	print(sci_run_cmmd)
	try:
		p = subprocess.Popen(sci_run_cmmd,shell=True,bufsize=0,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		scirun_processes.append(p)
		if ( not pause ) :
			(stdoutdata, stderrdata) = p.communicate()
	except (OSError,):
		e = sys.exc_info()[1]
		print("Execution failed with error ", e)
		sys.exit(0)
	done = (not pause)
	while not done:
		input_str = Utils.get_input("Continue? (y/n)")
		if (input_str == "y"):
			done = True
		elif (input_str == "n"):
			print("Please close all SCIRun instances to exit...")
			# N.B.:  Volume rendering over a forwarded X11 ssh session may cause SCIRun to crash
			#        which will hang the python process while waiting for SCIRun to terminate.
			for proc in scirun_processes:
				proc.wait()
			sys.exit(0)

	#print p.pid, os.getpid()
	
	#scirun_processes.append(p)
	#time.sleep(2)
	#done = (not pause)
	#while not done:
	#    input_str = raw_input("Continue? (y/n)")
	#    if (input_str == "y"):
	#        done = True
	#    elif (input_str == "n"):
	#        print "Please close all SCIRun instances to exit..."
	 #       # N.B.:  Volume rendering over a forwarded X11 ssh session may cause SCIRun to crash
	 #       #        which will hang the python process while waiting for SCIRun to terminate.
	 #       for proc in scirun_processes:
	 #           proc.wait()
	  #      sys.exit(0)

def parse_nrrd_for_voxel_radius(fname):

	# run unu head on the input file and capture the result
	p = subprocess.Popen(["%s head %s" % (os.path.join(binary_path,"unu"), fname)],shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	Utils.rec_pid(p.pid)
	(out,err) = p.communicate() 

	spacing = ('1','1','1')

	# compose a regular expression for parsing out NRRD space directions
	tri = "\s*\(([\-\d\.]+),\s*([\-\d\.]+),\s*([\-\d\.]+)\)\s*"
	spc_re = re.compile(("space directions:%s%s%s" % (tri, tri, tri)).encode())
	res = spc_re.search(out)
	if (res != None):
		digits = res.groups()
		x_spacing = float(digits[0])+float(digits[3])+float(digits[6])
		y_spacing = float(digits[1])+float(digits[4])+float(digits[7])
		z_spacing = float(digits[2])+float(digits[5])+float(digits[8])  
		spacing = (x_spacing,y_spacing,z_spacing)
	else:
		# compose a regular expression for parsing out NRRD spacing
		pat = re.compile((r'spacings:\s*([\-\d\.]+)\s+([\-\d\.]+)\s+([\-\d\.]+)').encode())
		res = pat.search(out)
		if (res != None):
			spacing = res.groups()

	# NOTE: the next line should be removed when spacing/orientation is
	#  fixed in the sizing field and particle sytem code.  For now, it's
	#  just here as a safety precaution.  Secondary note: the above
	#  code should be enhanced so NRRD lines that start with a '#' are
	#  ignored.
	# spacing = ('1','1','1')

	radius = sqrt(pow(float(spacing[0]),2)+pow(float(spacing[1]),2)+pow(float(spacing[2]),2))

	medial_axis_point_scale = radius * 4.0/15.0
	os.environ['MEDIAL_AXIS_POINT_SCALE'] = '%f' % medial_axis_point_scale
	particles_vector_scale = radius * 17.0/15.0
	os.environ['PARTICLES_VECTOR_SCALE'] = '%f' % particles_vector_scale
	particles_point_scale = radius * 4.0/15.0
	os.environ['PARTICLES_POINT_SCALE'] = '%f' % particles_point_scale
	
#    print "fname=%s  spacing=%s %s %s\n"%(fname,spacing[0],spacing[1],spacing[2])
#    print "point scale: " + os.environ['MEDIAL_AXIS_POINT_SCALE']
#    print "vec scale: " + os.environ['PARTICLES_VECTOR_SCALE']
#    print "point scale2: " + os.environ['PARTICLES_POINT_SCALE']

#    print radius
#    sys.exit(0)

model_output_path = ""

def usage() :
	print("Usage: BuildMesh.py [-sFirstStage:LastStage,--stages=FirstStage:LastStage] [-h,--help] [-i,--interactive] [-d,--display-only] [-p, --binary-path] model_config")
	print("-s/--stages= allow a subset of the stages to be run.  The stages are:\n")
	print("1\tMake solo NRRDs from a label NRRD.")
	print("2\tExtract the material surfaces from each solo NRRD.")
	print("3\tCompute the medial axis for each material surface.")
	print("4\tCompute the sizing field (representing local feature size)\n\tfrom the medial axis.")
	print("5\tGenerate the initial sampling of the material interfaces.")
	print("6\tFrom the seeds, run a particle system to optimize the\n\tsampling of the material interfaces.")
	print("7\tGenerate a surface mesh for each of the materials, and\n\tfor all materials combined.")
	print("8\tGenerate a volume mesh for each of the materials, and for\n\tall materials combined.\n")
	print("So, for example")
	print("\tBuildMesh.py -s1:3 my_model_config.py")
	print("would run stages 1 through 3 on the model specified by my_model_config.py")
	print("-i/--interactive launches a SCIRun network at the end of each relevant")
	print("pipeline stage to help the user visualize the results.")
	print("-d/--display-only displays the output of key pipeline stages using SCIRun without running the meshing pipeline stages.")
	print("If used alone the files must already exist (the pipeline must have already been run) in order for this flag to work.")
	print("-p/--binary-path Path to executable programs used in the BioMesh3d pipeline.")
	print("If this option is not used, this path must be set in the system PATH.")
	print("-c/--client Specify the machine that will send input to and collect output (rendered images) from SCIRun")
	print("-h/--help prints this screen.")

if __name__ == '__main__' :

	args = argv
	# Check that we got the right arguments......

	try:
		opts, args = Utils.get_options(args[1:],"ihds:p:c:",["help","interactive","display-only","stages=","binary-path=","client="]);
	except (getopt.GetoptError,):
		usage()
		sys.exit(-1) 

	start_stage = 1    
	stop_stage = 8
	is_interactive = False
	display_only = False

	scripts_dir = os.path.abspath(os.path.dirname(argv[0]))
	main_bin_dir, dummy = os.path.split(scripts_dir)     
	binary_path=""
	
	# are we in an installed directory?
	binary_path=main_bin_dir

	if sys.platform == "win32":
		unu_exec = "unu.exe"
	else:
		unu_exec = "unu"
		os.setpgrp();

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
	
	print("default binary_dir = %s " % binary_path)
	thin_client = False
	scirun_thin_client_args = ""

	for o,a in opts:
		if o in ["-s","--stages"]:
			start_stage_str,stop_stage_str = a.split(":")
			start_stage = int(start_stage_str)
			stop_stage = int(stop_stage_str)
		if o in ["-i","--interactive"]:
			is_interactive = True
		if o in ["-h","--help"]:
			usage()
			sys.exit(0)
		if o in ["-d","--display-only"]:
			display_only = True
		if o in ["-p","--binary-path"]:
			binary_path  = a
		if o in ["-c","--client"]:
			thin_client = True;
			ip_address,server_port,port_range = a.split(":")
#            print ("%s %s %s\n" % (ip_address,server_port,port_range))
			scirun_thin_client_args = " -ip_address " + ip_address + " -server_port " + server_port + " -port_range " + port_range + " ";
			print ("scirun thin client args: %s" %  scirun_thin_client_args)

	pause_for_stage = (start_stage != stop_stage)
	  
	model_config = ""
	if len(args) < 1:
		print("Model configuration file is needed to run meshing pipeline.")
		getUserInput = True
		while getUserInput:
			model_config = Utils.get_input("Which model configuration file to use?  ")
			if (os.path.exists(model_config)):
				getUserInput = False
			else:
				print("\nModel configuration file "+ model_config + " does not exist.\nTo create a model_config file automatically, use\n"+os.path.join(scripts_dir,"BuildModelConfig.py")+".\n")
				input_str = Utils.get_input("Try again? (y/n)")
				if (input_str == "n"):
					print("Exiting...")
					sys.exit(0)
	else:
		model_config = args[0]
		if not(os.path.exists(model_config)):
			print("Model configuration file %s does not exist." % args[0])
			sys.exit(-1)
		
	model_config = os.path.abspath(model_config)
	
#    print(model_config)
#    str = open(model_config)
#    print(str)

	try:
		exec(open(model_config).read())
	except:
		err_message = "There are problems with the model config file: "+str(sys.exc_info()[1])
		Utils.write_error(err_message)
		print(err_message)
		sys.exit(-1)
		
	model_path, dummy = os.path.split(model_config)
	
	if not(os.path.isabs(model_input_file)) :
	  model_input_file = os.path.normpath(os.path.join(model_path, model_input_file))

	if not(os.path.isabs(model_output_path)) :
	  model_output_path = os.path.normpath(os.path.join(model_path, model_output_path))

	print("*******************************************")
	print("model_config = %s" % model_config)
	print("model_input_file = %s" % model_input_file)
	print("model_output_path = %s" % model_output_path)
	print("*******************************************")

	#python_exec = sys.executable
	#if len(python_exec) == 0:
		# try to supply a reasonable default
	python_exec = sys.executable

	if not os.path.exists(model_input_file) :
	  print("The model_input_file %s does not exist" % model_input_file)
	  sys.exit(-1)

	if not os.path.exists(model_output_path) :
		os.makedirs(model_output_path)

	# Done checking arguments

	overwrite = True
	Utils.output_path = model_output_path

	print("Running stages %s through %s" % (start_stage,stop_stage))

	parse_nrrd_for_voxel_radius(model_input_file)
	
	server_ip_address = "127.0.0.1"
	
	print("")
	print("***************")
	print("Stage 1 -- Make Solo Nrrds")
	print("***************")
	current_stage = 1

	if ((start_stage <= current_stage) and (stop_stage >= current_stage)):
		if not display_only:
			Utils.record_stage_run_file(current_stage, server_ip_address)
			Utils.current_stage = current_stage
			Utils.delete_error_log()
			Utils.delete_complete_log()
			Utils.delete_running_log()
		Utils.rec_pid(os.getpid(), overwrite) #added by rtao, "w" mode delete old pids files
		script_command = os.path.join(scripts_dir,"MakeSoloNrrd.py")
		solo_cmmd = r'"%s" -u "%s" "%s" "%s"' % (python_exec, script_command, model_config, binary_path)
		#print solo_cmmd
		if not display_only:
			Utils.do_system(solo_cmmd,print_only,use_shell)
	
	print("")
	print("***************")
	print("Stage 2 -- Extract Material Surfaces")
	print("***************")
	current_stage += 1
	if (start_stage <= current_stage and stop_stage >= current_stage):
		if not display_only:
			Utils.record_stage_run_file(current_stage, server_ip_address)
			Utils.current_stage = current_stage
			Utils.delete_error_log()
			Utils.delete_complete_log()
			Utils.delete_running_log()
		Utils.rec_pid(os.getpid(), overwrite) #added by rtao, "w" mode delete old pids files
		script_command = os.path.join(scripts_dir,"ComputeMaterialBoundary.py")
		bdry_cmmd = r'"%s" -u "%s" "%s" "%s"' % (python_exec, script_command, model_config, binary_path)
		if not display_only:
			Utils.do_system(bdry_cmmd,print_only,use_shell)

		if (is_interactive):
			os.environ['COMPUTE_MATERIAL_BOUNDARY_OUTPUT'] = '%s/isosurface-all.ts.fld' % model_output_path.replace("\\", "/")
			#print(scripts_dir)
			#pause_for_scirun_cmd(scirun_thin_client_args + "-e %sMaterialBoundaries.srn" % scripts_dir,pause_for_stage, current_stage)
			scirun_net = os.path.join(scripts_dir,"MaterialBoundaries.srn")
			pause_for_scirun_cmd(scirun_thin_client_args + "--nosplash -noeai -e \"%s\"" % scirun_net,pause_for_stage)

	print("")
	print("***************")
	print("Stage 3 -- Compute Medial Axis")
	print("***************")
	current_stage += 1
	
	if (start_stage <= current_stage and stop_stage >= current_stage):
		if not display_only:
			Utils.record_stage_run_file(current_stage, server_ip_address)
			Utils.current_stage = current_stage

		try :
			if ( constant_sizing_value == -1 ) :
				raise NameError("Invalid Sizing Value.")
	
			print(" ")
			print("Skipping -- constant sizing field requested!")
			Utils.rec_completed_log()
	
		except NameError :
	
				Utils.rec_pid(os.getpid(), overwrite) #added by rtao, "w" mode delete old pids files
			
				script_command = os.path.join(scripts_dir,"ComputeMaterialMedialAxis.py")
				ma_cmmd = r'"%s" -u "%s" "%s" "%s"' % (python_exec, script_command, model_config, binary_path)
				if not display_only:
					Utils.do_system(ma_cmmd,print_only,use_shell)
				if (is_interactive):
					os.environ['COMPUTE_MATERIAL_BOUNDARY_OUTPUT'] = '%s/isosurface-all.ts.fld' % model_output_path.replace("\\", "/")
					os.environ['COMPUTE_MATERIAL_MEDIAL_AXIS_OUTPUT'] = '%s/ma-all.pc.fld' % model_output_path.replace("\\", "/")
					scirun_net = os.path.join(scripts_dir,"MedialAxisAndBoundaries.srn")
					pause_for_scirun_cmd(scirun_thin_client_args + "--nosplash -noeai -e \"%s\"" % scirun_net, pause_for_stage)
			
	print("")
	print("***************")
	print("Stage 4 -- Compute Sizing Field")
	print("***************")
	current_stage += 1
	if (start_stage <= current_stage and stop_stage >= current_stage):
		if not display_only:
			Utils.record_stage_run_file(current_stage, server_ip_address)
			Utils.current_stage = current_stage
		Utils.rec_pid(os.getpid(), overwrite) #added by rtao, "w" mode delete old pids files
		script_command = os.path.join(scripts_dir,"ComputeSizingField.py")
		sf_cmmd = r'"%s" -u "%s" "%s" "%s"' % (python_exec, script_command, model_config, binary_path)
		if not display_only:
			Utils.do_system(sf_cmmd,print_only,use_shell)
		if (is_interactive):
			initial_mat_idx = min(1,len(mat_names)-1)
			initial_mat = Utils.extract_root_matname(mat_names[initial_mat_idx])
			os.environ['COMPUTE_SIZING_FIELD_OUTPUT'] = '%s/%s_sf.nrrd' % (model_output_path.replace("\\", "/"),initial_mat.replace("\\", "/"))
			os.environ['COMPUTE_MATERIAL_BOUNDARY_OUTPUT_SOLO'] = "%s/%s_isosurface.ts.fld" % (model_output_path.replace("\\", "/"),initial_mat.replace("\\", "/"))
			os.environ['COMPUTE_MATERIAL_MEDIAL_AXIS_OUTPUT_SOLO'] = "%s/%s_ma.pc.fld" % (model_output_path.replace("\\", "/"),initial_mat.replace("\\", "/"))
			scirun_net = os.path.join(scripts_dir,"SizingFieldMedialAxisAndBoundaries.srn")
			pause_for_scirun_cmd(scirun_thin_client_args + "--nosplash -noeai -e \"%s\"" % scirun_net, pause_for_stage)

	print("")
	print("***************")
	print("Stage 5 -- Generate Initial Material Interface Sampling")
	print("***************")
	current_stage += 1
	if (start_stage <= current_stage and stop_stage >= current_stage):
		if not display_only:
			Utils.record_stage_run_file(current_stage, server_ip_address)
			Utils.current_stage = current_stage
		Utils.rec_pid(os.getpid(), overwrite) #added by rtao, "w" mode delete old pids files

		script_command = os.path.join(scripts_dir,"GenerateSeeds.py")
		seeds_cmmd = r'"%s" -u "%s" "%s" "%s"' % (python_exec, script_command, model_config, binary_path)
		if not display_only:
			Utils.do_system(seeds_cmmd,print_only,use_shell)

	print("")
	print("***************")
	print("Stage 6 -- Run Particle System")
	print("***************")
	current_stage += 1
	if (start_stage <= current_stage and stop_stage >= current_stage):
		if not display_only:
			Utils.record_stage_run_file(current_stage, server_ip_address)
			Utils.current_stage = current_stage
		Utils.rec_pid(os.getpid(), overwrite) #added by rtao, "w" mode delete old pids files
		
		script_command = os.path.join(scripts_dir,"DistributeParticles.py")
		ps_cmmd = r'"%s" -u "%s" "%s" "%s"' % (python_exec, script_command, model_config, binary_path)
		if not display_only:
			Utils.do_system(ps_cmmd,print_only,use_shell)
		if (is_interactive):
			os.environ['DISTRIBUTE_PARTICLES_PC_OUTPUT'] = '%s/junctions/particle-all.pc.fld' % model_output_path.replace("\\", "/")
			os.environ['DISTRIBUTE_PARTICLES_PCV_OUTPUT'] = '%s/junctions/particle-all.pcv.fld' % model_output_path.replace("\\", "/")
			scirun_net = os.path.join(scripts_dir,"DistributedParticles.srn")
			pause_for_scirun_cmd(scirun_thin_client_args + "--nosplash -noeai -e \"%s\"" % scirun_net, pause_for_stage)

	print("")
	print("***************") 
	print("Stage 7 -- Generate Surface Mesh")
	print("***************")
	current_stage += 1
	if (start_stage <= current_stage and stop_stage >= current_stage):
		if not display_only:
			Utils.record_stage_run_file(current_stage, server_ip_address)
			Utils.current_stage = current_stage
		Utils.rec_pid(os.getpid(), overwrite) #added by rtao, "w" mode delete old pids files
		script_command = os.path.join(scripts_dir,"BuildMaterialInterfaceMesh.py")
		iface_cmmd = r'"%s" -u "%s" "%s" "%s"' % (python_exec, script_command, model_config, binary_path)
		if not display_only:
			Utils.do_system(iface_cmmd,print_only,use_shell)

	print("")
	print("***************")
	print("Stage 8 -- Generate Volume Mesh")
	print("***************")
	current_stage += 1
	if (start_stage <= current_stage and stop_stage >= current_stage):
		if not display_only:
			Utils.record_stage_run_file(current_stage, server_ip_address)
			Utils.current_stage = current_stage
		Utils.rec_pid(os.getpid(), overwrite) #added by rtao, "w" mode delete old pids files
		script_command = os.path.join(scripts_dir,"BuildVolumetricMesh.py")
		vol_cmmd = r'"%s" -u "%s" "%s" "%s"' % (python_exec, script_command, model_config, binary_path)
		if not display_only:
			Utils.do_system(vol_cmmd,print_only,use_shell)
		if (is_interactive):
			os.environ['BUILD_VOLUMETRIC_MESH_OUTPUT'] = '%s/junctions/particle-union.tets-labeled_transformed.fld' % model_output_path.replace("\\", "/")
			#print(os.getpid())
			scirun_net = os.path.join(scripts_dir,"VolumetricMesh.srn")
			pause_for_scirun_cmd(scirun_thin_client_args + "--nosplash -noeai -e \"%s\"" % scirun_net, pause_for_stage)

	print("")
	print("**************************")
	print("Biomesh3D completed")
	print("**************************")
	Utils.write_stage_completed(model_output_path)
