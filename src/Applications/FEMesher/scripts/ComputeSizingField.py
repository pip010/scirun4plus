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

from sys import argv
import sys
import subprocess
import re
import os
import time
import glob

if sys.version_info[0] < 3 :
    from thread import start_new_thread
else :
    from _thread import start_new_thread

from math import fabs
import fileinput
import Utils

print_only = False
use_shell = True
if sys.platform == "win32" :
    use_shell = False

def assert_exists(file):
    if not(os.path.exists(file)):
        print("No such file %s. Perhaps you need to run a previous stage."%file)
        sys.exit(-1)

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

    # the output from processes need to be dumped or the
    # buffer fills up and it sleeps.
    for p in procs :
        print("starting thread to dump")
        print(p)
        start_new_thread(dump_log, (p, ))
        

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
            if p.poll() != None : # None means process not finished.
                ndone = ndone + 1

        if ndone == len(procs) :
            done = True
        else :
            if count % 30 == 0 :
                print("%d processes are complete of %d" % (ndone, len(procs)))
        count = count + 1

        

def create_sizing_field_cmmd(ma_par, ma_base, indicator) :
    try :
        cmmd = r'"%s" -c %f %s %s %d outputnrrd' % (os.path.join(binary_path,"sizingfield"), constant_sizing_value, ma_par, ma_base, indicator) # sizingfield.cxx writes all of the external files produced
    except NameError :
        cmmd  = r'"%s" %s %s %d outputnrrd' % (os.path.join(binary_path,"sizingfield"), ma_par, ma_base, indicator) # sizingfield.cxx writes all of the external files produced
    return cmmd

    
def create_sizing_field(ma_par, ma_base, indicator) :
    global model_output_path

    cmmd = create_sizing_field_cmmd(ma_par, ma_base, indicator)
    
    #print(cmmd)
    if print_only :
        print(cmmd)
    else :
        #return subprocess.Popen(cmmd, bufsize=0,shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=model_output_path)
        return subprocess.Popen(cmmd, bufsize=0,shell=use_shell, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=model_output_path)
        #return subprocess.Popen(cmmd, shell=True, stderr=subprocess.STDOUT, cwd=model_output_path)

if __name__ == "__main__" :

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

    if not(os.path.exists(model_output_path)):
        os.makedirs(model_output_path)

    Utils.output_path = model_output_path
    Utils.current_stage = 4

    curr_work_dir = os.getcwd()
    os.chdir(model_output_path)

    Utils.delete_complete_log()
    Utils.delete_error_log()
    Utils.rec_running_log()

    start_time = time.time()

    idx = 0
    procs = []

    for i in mats :
        n = Utils.extract_root_matname(mat_names[idx])
        assert_exists("medial_axis_param_file.txt")
        p = create_sizing_field("medial_axis_param_file.txt", n, idx)
        procs.append(p)
        Utils.rec_pid(p.pid)  #added by rtao, for killing the entire process tree
        idx = idx + 1
    
    wait_on_procs(procs)

    stop_time = time.time()
    f = open("compute-sizing-field-runtime.txt","w")
    f.write("%lf" % (stop_time-start_time))
    f.close()

    Utils.rec_completed_log()

    os.chdir(curr_work_dir)
