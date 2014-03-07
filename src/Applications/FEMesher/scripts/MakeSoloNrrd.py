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


if sys.version_info[0] < 3 :
	from thread import start_new_thread
else :
	from _thread import start_new_thread

from math import fabs

maxval_re = re.compile("max: ([\d\.\d]+)")
tri = "\(([\-\d\.]+),([\-\d\.]+),([\-\d\.]+)\)\s?"
spc_re = re.compile("space directions: %s%s%s" % (tri, tri, tri))

maxval = 0
spacing = None
print_only = False
use_shell = True
if sys.platform == "win32" :
	use_shell = False

stage=1

def inspect_nrrd(innrrd) :
	global maxval
	
	#print("Binary path: " + binary_path)
	cmmd = r'"%s" minmax "%s"' % (unu_path, innrrd)
	#print(cmmd)
	p = subprocess.Popen(cmmd,  shell=use_shell, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=model_output_path)
	
	Utils.rec_pid(p.pid) #added by rtao
	
	stdoutdata = p.communicate()[0]

	if (sys.version_info[0] < 3) :
		for l in string.split(stdoutdata, '\n') :
			mo = maxval_re.match(l.decode())
			if mo != None :
				maxval = int(mo.group(1))
	else :
		for l in stdoutdata.split(bytes('\n','ascii')) :
			mo = maxval_re.match(l.decode())
			if mo != None :
				maxval = int(mo.group(1))

def unorient_nrrd(innrrd) :
	nrrd_name = innrrd[:-5]
	unorient_fname = "%s_unorient.nrrd" % nrrd_name
	transform_fname = "%s_transform.tf" % nrrd_name
	unorient_cmd = r'"%s" -input "%s" -output "%s" -transform "%s"' % (os.path.join(binary_path,"UnorientNrrdAndGetTransform"), innrrd, unorient_fname, transform_fname)
	Utils.do_system(unorient_cmd,print_only,use_shell)

	print("*** Unorienting complete.")
	return unorient_fname


def pad_nrrd(innrrd,output_path) :
	outnrrd = "%s_pad.nrrd" % innrrd[:-5]
	dummy, filename  = os.path.split(outnrrd)
	outnrrd = os.path.join(output_path,filename)

	cmmd = r'"%s" pad -min -4 -4 -4 -max M+4 M+4 M+4 -b pad -v 0 -i "%s" -o "%s"' % (unu_path, innrrd, outnrrd)

	Utils.do_system(cmmd,print_only,use_shell)
	return outnrrd

def make_lut (mats,name) :
	global maxval
	
	if (not getattr(mats,'__iter__', False)) :
		mats = (mats,)

	mats = Utils.create_label_list(mats)
		
	if mats[0] == -1 : return #done by make_others_lut
	
	hit_str = "1.0\n"
	miss_str = "0.0\n"
		
	f = open("%s.lut.raw" % name, "w")
	flist = []
	
	for i in range(0, maxval + 1) :
		if i in mats :
			flist.append(hit_str)
		else :
			flist.append(miss_str)
	f.writelines(flist)
	f.close()

	cmmd = r'"%s" make -i %s.lut.raw -t float -s %d -e ascii -o %s.lut.nrrd' % (unu_path, name, maxval+1, name)
	# print cmmd
	Utils.do_system(cmmd,print_only,use_shell)

def make_others_lut (mats,name) :
	global maxval

	mats = Utils.create_label_list(mats)
	
	hit_str = "0.0\n"
	miss_str = "1.0\n"

	f = open("%s.lut.raw" % name, "w")
	flist = []

	for i in range(0, maxval + 1) :
		if i in mats:
			flist.append(hit_str)
		else:
			flist.append(miss_str)
	f.writelines(flist)
	f.close()

	cmmd = "%s make -i %s.lut.raw -t float -s %d -e ascii "\
		   "-o %s.lut.nrrd" % (unu_path, name, maxval+1, name)
	Utils.do_system(cmmd,print_only,use_shell)

def make_solo_material(idx, innrrd, name) :
	global maxval
	make_lut(idx,name)
	oname = "%s.solo.nrrd" % name
	cmmd = r'"%s" lut -m %s.lut.nrrd -min 0 -max %d -t float -i "%s" -o "%s"' % (unu_path, name, maxval, innrrd, oname)
	Utils.do_system(cmmd,print_only,use_shell)
	return oname


param_lines = []


def parallel_cmmds(cmmds) :
	procs = []
	for c in cmmds :
		if print_only :
			print(c)
		else :
			procs.append(popen2.Popen4(c))
	return procs
		
def wait_on_procs(procs) :    

	# the output from processes need to be dumped or the
	# buffer fills up and it sleeps.
	for p in procs :
		print("starting thread to dump")
		print(p)
		start_new_thread(dump_log, (p, ))
		

	#once all tight procs are done, make the vols    
	done = print_only
	count = 0
	while not done :
		time.sleep(1)
		Utils.rec_running_log()
		ndone = 0
		
		#all tight procs need to be finished before moving on.
		np = 0
		for p in procs :
			if p.poll() != None : # None means process not finished.
				ndone = ndone + 1

		if ndone == len(procs) :
			done = True
		else :
			if count % 30 == 0 :
				print("%d processes are complete of %d" % (ndone, len(procs)))
		count = count + 1


# return the cmmd to run only.
def make_tight_cmmd(innrd, name, r) :
	oname = "%s.tight.nrrd" % name
	cmmd = r'"%s" %s "%s" "%s"' % (os.path.join(binary_path,"morphsmooth"), r, innrd, oname)
	return oname, cmmd

# launch a process and return the Popen4 object to poll
def make_tight(innrd, name, r) :
	global model_output_path
	oname, cmmd = make_tight_cmmd(innrd, name, r)

	# print cmmd

	if print_only :
		print(cmmd)
	else :
		return oname, subprocess.Popen(cmmd, shell=use_shell, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=model_output_path)

def compute_tightened_labels(tightened_files) :
	cmmd = "%s %s" % (os.path.join(binary_path,"ComputeTightenedLabels"), tightened_files)
	print("Tightened command ", cmmd)
	return subprocess.Popen(cmmd, shell=use_shell, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=model_output_path)

def dump_log(p) :
	if print_only :
		return

	# Warning: The data read is buffered in memory, so do not use this method if the data size is large or unlimited.
	# See Python docs, section 18.1.2.
	(stdoutdata, stderrdata) = p.communicate()
	f = open("%d.log" % p.pid, "w")
	f.writelines(stdoutdata.decode())
	f.close()
	
if __name__ == "__main__" :

# Check that we got the right arguments......
	if len(argv) < 2 or not( os.path.exists( argv[1] ) ):
		print('Usage: %s [model_config] [binary_path]' % argv[0])
		sys.exit(1)

	
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

	unu_path = os.path.join(binary_path,"unu")

	if not(os.path.exists(model_output_path)):
		os.makedirs(model_output_path)

	Utils.output_path = model_output_path
	Utils.current_stage = 1

	curr_work_dir = os.getcwd()
	os.chdir(model_output_path)

	Utils.delete_complete_log()
	Utils.delete_error_log()
	Utils.rec_running_log()
	
	start_time = time.time()

	procs = []
	original_nrrd = model_input_file

	inspect_nrrd(original_nrrd)
	padded_nrrd = pad_nrrd(original_nrrd,model_output_path)

	nrrd = unorient_nrrd(padded_nrrd)
	os.remove(padded_nrrd)

	#create the lut for remaining materials special
	count = 0
	others_name = ''
	for i in mats :
		others_name = Utils.extract_root_matname(mat_names[count])
		count = count + 1
		if i == -1 :
			break
			

	make_others_lut (mats,"%s"% others_name)

	idx = 0
	tightened_files = ""

	for i in mats :

		n = Utils.extract_root_matname(mat_names[idx])

		r = mat_radii
		print("Working on material: %s" % n)
		idx = idx + 1
		cur = make_solo_material(i, nrrd, n)
		cur, p = make_tight(cur, n, r)

		Utils.rec_pid(p.pid)  #added by rtao, for killing the entire process tree

		procs.append(p)
		tightened_files += cur + " "
#        param_lines.append('%s\n' % (cur.rstrip(".nrrd") + "-corrected.nrrd")) 
		param_lines.append('%s\n' % cur) 

	#write param file
	nmats = len(param_lines)
	lns = []
	lns.append("%d\n" % nmats)
	for l in param_lines :
		lns.append(l)
	f = open("medial_axis_param_file.txt", 'w')
	f.writelines(lns)
	f.close()

	wait_on_procs(procs)
	procs = []
			
	p = compute_tightened_labels(tightened_files)
	Utils.rec_pid(p.pid)
	while p.poll() == None :
		Utils.rec_running_log()
		continue

	stop_time = time.time()
	f = open("make-solo-nrrd-runtime.txt","w")
	f.write("%lf" % (stop_time - start_time))
	f.close()

	Utils.rec_completed_log()
	os.chdir(curr_work_dir)

	sys.exit(0)
