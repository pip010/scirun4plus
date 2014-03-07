#  
#  For more information, please see: http://software.sci.utah.edu
#  
#  The MIT License
#  
#  Copyright (c) 2009 Scientific Computing and Imaging Institute,
#  University of Utah.
#  
#  License for the specific language governing rights and limitations under
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
#    File   : fix_externals.py
#    Author : Martin Cole
#    Date   : Mon Mar  3 14:51:16 2008

# Use this on a new branch to fix all externals to the specific
# revision at the time the branch is made.
# 1) make a branch.
# 2) check out the branch.
# 3) run this from the top level i.e in SCIRun
# 4) run svn status to verify the changes.
# 5) svn commit your changes to the branch.

from sys import argv
import os
import re

prop_re = re.compile('([\w\.]+)\s+([\w\.:/]+)')
rev_re = re.compile('Revision:\s+(\d+)')

def find_externals(d) :
    
    for root, dirs, files in os.walk(d):
        cmmd = 'svn propget svn:externals %s' % root
        sin, std, err = os.popen3(cmmd)
        sl = std.readlines()

        if len(sl) > 0 :
            print "%s has externals property" % root
            for l in sl :
                mo = prop_re.match(l)
                if mo != None :
                    cmmd = 'svn info %s' % mo.group(2)
                    sin1, std1, err1 = os.popen3(cmmd)
                    rev = None
                    for r in std1.readlines() :
                        mo1 = rev_re.match(r)
                        if mo1 != None :
                            rev = int(mo1.group(1))
                    s = "%s -r%d %s" % (mo.group(1), rev, mo.group(2))
                    cmmd = 'svn propset svn:externals "%s" %s' % (s, root)
                    print cmmd
                    os.system(cmmd)
                    
        if '.svn' in dirs:
            dirs.remove('.svn')  # don't visit .svn directories

def restore_trunk_externals(d) :
    for root, dirs, files in os.walk(d):
        cmmd = 'svn propget svn:externals %s' % root
        sin, std, err = os.popen3(cmmd)
        sl = std.readlines()
        
        if len(sl) > 0 :
            print "%s has externals property" % root
            for l in sl :
                mo = prop_re.match(l)
                if mo != None :
                    cmmd = 'svn info %s' % mo.group(2)
 
                    s = "%s %s" % (mo.group(1), mo.group(2))
                    cmmd = 'svn propset svn:externals "%s" %s' % (s, root)
                    print cmmd
                    os.system(cmmd)

        if '.svn' in dirs:
            dirs.remove('.svn')  # don't visit .svn directories
            
def usage() :
    print ""
    print "Invalid input, Usage:"
    print "fix_externals.py <options> path"
    print "  valid options :"
    print "  -t, --trunk == remove all revisions from external references"
    exit(1)
        
if __name__ == '__main__' :
    if len(argv) > 2 :
        if argv[1] == '-t' or argv[1] == '--trunk' :
            restore_trunk_externals(argv[2])
        else :
            usage()
    else :
        find_externals(argv[1])
