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

combined_prefix = "particle-union"

print_only = False
use_shell = True
if sys.platform == "win32" :
    use_shell = False

def assert_exists(file):
    if not(os.path.exists(file)):
        print("No such file %s. Perhaps you need to run a previous stage."%file)
        sys.exit(-1)

def do_system(cmmd) :
    if print_only :
        print(cmmd)
    else :
        ret = -1
        try:
            ret = subprocess.call(cmmd, shell=use_shell)
            if (ret < 0):
                sys.exit(ret)
        except OSError as e:
            print("Execution failed with error ", e)
            sys.exit(ret)

def make_combined_vol(d) :

    assert_exists("medial_axis_param_file.txt")

    unjittered_interface_pts = "%s/%s.fld" % (d,combined_prefix)
    assert_exists(unjittered_interface_pts)
    jittered_interface_pts = "%s/%s-jit.fld" % (d,combined_prefix)
    jit_percent = 1E-6

    cmmd = r'"%sJitterPointCloud" %s %s %s' % (binary_path,jit_percent,unjittered_interface_pts,jittered_interface_pts)
    print("Jitter Command: %s" % cmmd)
    do_system(cmmd)

#    tet_fname1 = "%s/%s.tets1.fld" % (d, combined_prefix)
    tet_fname = "%s/%s.tets.fld" % (d, combined_prefix)

    cmmd = r'"%sInterfaceWithTetGen" -cmmd_line %s -main %s -out %s' % (binary_path, tetgen_joined_vol_flags, jittered_interface_pts, tet_fname)
    do_system(cmmd)

#    ref_tet_name = "%s/%s.tets-ref.fld" % (d, combined_prefix)

#    cmmd = r'"%sInterfaceWithTetGen" -cmmd_line pYzqAa10 -main %s -out %s' % (binary_path, tet_fname, ref_tet_name)
#    do_system(cmmd)

#    cmmd = r'"%sInterfaceWithTetGen" -cmmd_line %s -main %s -out %s' % (binary_path, tetgen_joined_vol_flags2, tet_fname1, tet_fname)
#    do_system(cmmd)

    # add labels to the field
    # nrrdLabelFile tetField labeledTetField
    labeled_tet_fname = "%s/%s.tets-labeled.fld" % (d, combined_prefix)
    cmmd = r'"%sLabelTets" %s/dominant_labelmap.nrrd %s %s' % (binary_path,model_output_path,tet_fname,labeled_tet_fname)
    print("%s\n" % cmmd)
    do_system(cmmd)

    transform_with_transform(labeled_tet_fname)

def transform_with_transform(tet_file) :
    assert_exists(tet_file)
    print("%s\n" % tet_file)
    transformed_fname = "%s_transformed.fld" % tet_file[:-4]
    transf = "%s_pad_transform.tf" % model_input_file[:-5]
    dummy,filename = os.path.split(transf)
    transf = os.path.join(model_output_path,filename)
    print("Transforming %s" % tet_file)
    output_path = tet_file[:-4]
    trans_cmmd = r'"%sTransformFieldWithTransform" -input %s -transform %s -output %s' % (binary_path, tet_file, transf, transformed_fname)
    do_system(trans_cmmd)
    print("Done Transforming")

def node2ptcl(node_fname,ptcl_fname) :
    assert_exists(node_fname)
    
    node_file = open(node_fname, 'r')
    ptcl_file = open(ptcl_fname, 'w')

    line = node_file.readline()
    vals = line.split()

    n_pts = vals[0]
    
    ptcl_file.write("%s\n" % n_pts);

    for line in node_file :
        vals = line.split()
        ptcl_file.write("%s %s %s 0 0 0 0\n" % (vals[1],vals[2],vals[3]));
    node_file.close()
    ptcl_file.close()

def make_nodes_from_dir(d) :
    filename = '%s/%s.pts' % (d, combined_prefix)
    if os.path.exists(filename) :
        os.remove(filename)

    # Basically does the same thing as 'cat'
    pts_file = open('%s/%s.pts' % (d, combined_prefix), 'w')
    for f in glob.glob('%s/particle_params*.pts' % d):
        f = open(f, 'r')
        pts_file.writelines(f.readlines())
        f.close()
    pts_file.close()

    pts_fname = "%s/%s.pts" % (d,combined_prefix)
    node_fname = "%s/%s.node" % (d,combined_prefix)

    cmmd = r'"%spts2node" %s %s' % (binary_path, pts_fname, node_fname)
    do_system(cmmd)
    
    # Convert Pts format to ptcl
    ptcl_fname = "%s/%s.ptcl" % (d,combined_prefix)
    node2ptcl(node_fname,ptcl_fname)
    # Then use build-in converter to produce point cloud field
    field_fname = "%s/%s.fld" % (d,combined_prefix)
    cmmd = r'"%s/PtclToField" -nodata %s %s' % (binary_path, ptcl_fname, field_fname)
    do_system(cmmd)


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
        if binary_path[len(binary_path) - 1] != os.sep:
            binary_path += os.sep
    else:
        binary_path = ""  

# Done checking arguments

    cur_dir = os.getcwd();
    os.chdir(model_output_path);

    start_time = time.time()

    dir_name = 'junctions'
    assert_exists(dir_name)

    make_nodes_from_dir(dir_name) 
    make_combined_vol(dir_name) 

    stop_time = time.time()
    f = open("build-volumetric-mesh2-runtime.txt","w")
    f.write("%lf" % (stop_time-start_time))
    f.close()

    os.chdir(cur_dir)
