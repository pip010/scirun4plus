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

def get_options(arglist,short_form,long_form) :
    try:
        opts, arglist = getopt.getopt(arglist,short_form,long_form)
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage()
        sys.exit(-1)

    return opts, arglist


def do_system(cmmd,print_only,use_shell) :
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

