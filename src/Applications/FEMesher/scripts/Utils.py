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
# Software is furniShed to do so, subject to the following conditions:
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

import getopt
import subprocess
import sys
import time
import os

current_stage = 1
output_path = "."

def write_stage_completed(model_output_path) :
    stage_path = os.path.join(model_output_path,"mesh_state.txt")
    contents = ""
    counter = 0
    
    try:
    
        f = open(stage_path, 'r')
        for line in f :
            if counter == 0:
                contents += "completed\n"
            elif counter == 2:
                contents += str(int(time.time())) + '\n'
            else:
                contents += line
            counter += 1
    except (IOError,):
        contents += "completed\n"
        contents += "\n"
        contents += str(int(time.time())) + '\n'
        contents += "\n"
    
    txt_file = open(stage_path, 'w')
    txt_file.write(contents)
    txt_file.close()
    
def write_error(error_message):
    stage_path = os.path.join(output_path,"mesh_state.txt")
    contents = ""
    counter = 0
    
    try:
        f = open(stage_path, 'r')
        for line in f :
            if counter == 0:
                contents += "error " + error_message + "\n"
            elif counter == 2:
                contents += str(int(time.time())) + '\n'
            else:
                contents += line
            counter += 1
        f.close()
    except (IOError,):
        contents += "error " + error_message + "\n"
        contents += '\n'
        contents += str(int(time.time())) + '\n'
        contents += '\n'
        
    txt_file = open(stage_path, 'w')
    txt_file.write(contents)
    txt_file.close()

def record_stage_run_file(stage, ip_address):
    stage_path = os.path.join(output_path,"mesh_state.txt")
    counter = 0
    description = ""
    try:
        f = open(stage_path, 'r')
        for line in f :
            if counter == 3:
                description = line
            counter += 1

        f.close()
    except (IOError,):
        description = "\n"
    contents = "running " + ip_address + "\n"
    contents += str(stage) + "\n"
    contents += str(int(time.time())) + "\n"
    contents += description

    txt_file = open(stage_path, 'w')
    txt_file.write(contents)
    txt_file.close()

def get_input(str):
    if (sys.version_info[0] >= 3):
        input_str = input(str)
    else:
        input_str = raw_input(str)
    return input_str

def get_options(arglist,short_form,long_form) :
    try:
        opts, arglist = getopt.getopt(arglist,short_form,long_form)
    except (getopt.GetoptError,):
        err = sys.exc_info()[1]
        # print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        raise

    return opts, arglist

def rec_pid(pid,overwrite=False) :
    pid_log_file = 'pid_%d.txt' % current_stage
    full_pid_log_file = r'%s' % os.path.join(output_path, pid_log_file)

    #if (not os.path.exists(full_pid_log_file)) :
    #    pid_f_handle = open(full_pid_log_file, "w")
    #else :
    pid_f_handle = open(full_pid_log_file, "a")

    pid_f_handle.write("%d\n" % pid)
    pid_f_handle.close()


def rec_running_log(overwrite=False) :
    pid_log_file = 'running_stage_%d.txt' % current_stage
    full_pid_log_file = r'%s' % os.path.join(output_path, pid_log_file)

    if (overwrite) :
        pid_f_handle = open(full_pid_log_file, "w")
    else :
        pid_f_handle = open(full_pid_log_file, "a")

    pid_f_handle.write("%d\n" % ( time.time() ))
    pid_f_handle.close()


def rec_completed_log(overwrite=False) :
    completed_log_file = 'completed_stage_%d.txt' % current_stage
    full_completed_log_file = r'%s' % os.path.join(output_path, completed_log_file)

    if (overwrite) :
        completed_f_handle = open(full_completed_log_file, "w")
    else :
        completed_f_handle = open(full_completed_log_file, "a")

    completed_f_handle.write("%d\n" % ( time.time() ))
    completed_f_handle.close()

def rec_error_log(overwrite=False) :
    error_log_file = 'error_stage_%d.txt' % current_stage
    full_error_log_file = r'%s' % os.path.join(output_path, error_log_file)

    if (overwrite) :
        error_f_handle = open(full_error_log_file, "w")
    else :
        error_f_handle = open(full_error_log_file, "a")

    error_f_handle.write("%d\n" % ( time.time() ))
    error_f_handle.close()

def delete_running_log() :
    running_log_file = 'running_stage_%d.txt' % current_stage
    full_running_log_file = r'%s' % os.path.join(output_path, running_log_file)
    if ( os.path.isfile( full_running_log_file) ):
        os.remove( full_running_log_file )

def delete_error_log() :
    error_log_file = 'error_stage_%d.txt' % current_stage
    full_error_log_file = r'%s' % os.path.join(output_path, error_log_file)
    if ( os.path.isfile( full_error_log_file) ):
        os.remove( full_error_log_file )
    
def delete_complete_log() :
    complete_log_file = 'complete_stage_%d.txt' % current_stage
    full_complete_log_file = r'%s' % os.path.join(output_path,complete_log_file)
    if ( os.path.isfile( full_complete_log_file) ):
        os.remove( full_complete_log_file )

def delete_pid_log() :
    pid_log_file = 'pid_%d.txt' % current_stage
    full_pid_log_file = r'%s' % os.path.join(output_path, pid_log_file)
    if ( os.path.isfile( full_pid_log_file) ):
        os.remove( full_pid_log_file )
        
def do_system(cmmd,print_only,use_shell) :
    if print_only :
        print(cmmd)
    else :
        ret = -1
        try:
            print(cmmd)
            #p = subprocess.Popen(cmmd, shell=True)
            p = subprocess.Popen(cmmd, shell=use_shell)
            #print("pid in medial axis python", p.pid)
            rec_pid(p.pid)
            #ret = p.wait()

            while True:
                 p.poll()
                 ret = p.returncode
                 #print("return code:", ret)
                 rec_running_log()
                 if ( ret !=  None ):
                     break;
                 else:
                     time.sleep(1)  #check every 10 second
   
            #ret = subprocess.call(cmmd, shell=True)
            if (ret != 0):
                rec_error_log()
                #delete_running_log()
                write_error(r"Command '%s' failed with exit code %d" % (cmmd,ret))
                sys.exit(ret)

            if (ret == 0):
                #delete_running_log()
                rec_completed_log()
                
        except (OSError,):
            e = sys.exc_info()[1]
            print("Execution failed with error ", str(e))
            write_error(str(e))
            #delete_running_log()
            rec_error_log()
            write_error(r"Failed to run command '%s'" % cmmd)
            sys.exit(-1)

def extract_root_matname(matname) :
	# if we got a list of names, use the first name in the list
	# This handles the case where the user has combined multiple labels
	# into a single label, but wants to preserve the names, so for example:
	# ('heart',('ventricles','atria')).  We'll use only the first name in
	# naming files however.
	if (getattr(matname,'__iter__', False) and not isinstance(matname,str)) :
		matname=matname[0]
	return matname
	
def create_label_list(labels) :
	# we allow the user to group labels so that they are treated as one segment.
	# In order to determine the labels that will become part of this new label,
	# we produce a list that contains all the constituent labels.
	label_list = []
	if ((not getattr(labels,'__iter__', False))) :
		label_list.append(labels)
	else:
		for l in labels :
			label_list.extend(create_label_list(l))
	return label_list
	
