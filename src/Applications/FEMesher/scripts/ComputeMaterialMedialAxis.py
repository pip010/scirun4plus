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

if sys.version_info[0] < 3 :
    from thread import start_new_thread
else :
    from _thread import start_new_thread

from math import fabs
import Utils


maxval_re = re.compile("max: ([\d\.\d]+)")
axsizes_re = re.compile("sizes: (\d+) (\d+) (\d+)")
tri = "\(([\-\d\.]+),([\-\d\.]+),([\-\d\.]+)\)\s?"
spc_re = re.compile("space directions: %s%s%s" % (tri, tri, tri))

maxval = 0
axsizes = ()
spacing = None
print_only = False
use_shell = True
if sys.platform == "win32" :
    use_shell = False

def assert_exists(file):
    if not(os.path.exists(file)):
        print("No such file %s. Perhaps you need to run a previous stage."%file)
        exit(-1)


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
    Utils.current_stage = 3

    curr_work_dir = os.getcwd()
    os.chdir(model_output_path)

    Utils.delete_complete_log()
    Utils.delete_error_log()
    Utils.rec_running_log()
    
    field_ma_names = []

    start_time = time.time();
    idx = 0    
    for i in mats :
        n = Utils.extract_root_matname(mat_names[idx])
        p_ptcl_filename = "%s_ma.ptcl" % n
        field_ma_names.append("%s_ma.pc.fld" % n)
        iso_fname = "%s_isosurface.ts.fld" % n
        if (idx == 0):
            iso_fname = "%s_isosurface.ts.fld" % n
        assert_exists(iso_fname)
        (p_ptcl_root,p_ptcl_ext) = os.path.splitext(p_ptcl_filename)

        try:
            start_bins = initial_medial_axis_bins
        except NameError:
            start_bins = 100

        gen_ma_cmd = r'"%s" -surface %s -levels %s -start_num_bins %s -output %s.pc.fld' % (os.path.join(binary_path,"GenerateMedialAxisPoints"), iso_fname, refinement_levels, start_bins, p_ptcl_root)
        #print("**** medial axis command: %s " % gen_ma_cmd)
        Utils.do_system(gen_ma_cmd,print_only,use_shell)

        field_to_ptcl_cmmd = r'"%s" %s.pc.fld %s' % (os.path.join(binary_path,"FieldToPtcl"), p_ptcl_root, p_ptcl_filename)
        # print field_to_ptcl_cmmd
        Utils.do_system(field_to_ptcl_cmmd,print_only,use_shell)
# these next two lines may not be necessary, since .pc.fld already exists
        ptcl_to_field_cmmd = r'"%s" -nodata %s %s.pc.fld' % (os.path.join(binary_path,"PtclToField"), p_ptcl_filename, p_ptcl_root)
        Utils.do_system(ptcl_to_field_cmmd,print_only,use_shell)
        idx = idx + 1

    join_field_cmmd = r'"%s" -indexdata ma-all.pc.fld %s' % (os.path.join(binary_path,"JoinFields"), " ".join(field_ma_names))
    Utils.do_system(join_field_cmmd,print_only,use_shell)

    stop_time = time.time();
    f = open("compute-material-medial-axis-runtime.txt","w")
    f.write("%lf" % (stop_time-start_time))
    f.close()

    Utils.rec_completed_log()
    os.chdir(curr_work_dir)
