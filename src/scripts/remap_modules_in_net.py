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
#    File   : remap_modules_in_net.py
#    Author : Martin Cole
#    Date   : Thu Oct 12 10:06:20 2006

import sys
import re
import os

do_port_swap = 1

# the list of module ids that match CreateFieldFromNrrdData
# gets reset on each new srn file.  This is used to swap
# first and third input port connections to this module
port_swap_input_ids = []

# the list of module ids that match SplitFieldIntoNrrdData
# gets reset on each new srn file.  This is used to swap
# first and third input port connections to this module
port_swap_output_ids = []


ws_exp = re.compile('\s')
pcm_exp = re.compile('(\w+)::(\w+)::(\w+)')
from_to_exp = re.compile('([\w:]+)\s+([\w:]+|###).*')
remap = {}

con_exp = re.compile('(.+addConnection)\s+\$(m\d+)\s+(\d+)\s+\$(m\d+)\s+(\w+)(.+)')
mod_exp = re.compile('(.+addModuleAtPosition)\s+"(\w+)"\s+"(\w+)"\s+"(\w+)"(.+)')
modid_exp = re.compile('set\s+(m\d+)\s+\[addModuleAtPosition')
def parse_net(fn) :
  global remap
  global do_port_swap
  global port_swap_input_ids
  global port_swap_output_ids

  var_tup = ()
  outfile = []
  newln = ''
  f = open(fn, 'r')
  for ln in f.readlines() :
    mo = mod_exp.match(ln)
    cmo = con_exp.match(ln)
    if mo != None :
      key = "%s::%s::%s" % (mo.group(2), mo.group(3), mo.group(4))
      add_swap_input = 0
      add_swap_output = 0

      if do_port_swap and mo.group(4) == "CreateFieldFromNrrdData" :
        add_swap_input = 1
      if do_port_swap and mo.group(4) == "SplitFieldIntoNrrdData" :
        add_swap_output = 1
        
      if remap.has_key(key) :
        nmo = pcm_exp.match(remap[key])
        newln =  '%s "%s" "%s" "%s"%s\n' % (mo.group(1), nmo.group(1),
                                            nmo.group(2), nmo.group(3),
                                            mo.group(5))
        outfile.append(newln)
        if nmo.group(3) == "CreateFieldFromNrrdData" :
          add_swap_input = 1
        if nmo.group(3) == "SplitFieldIntoNrrdData" :
          add_swap_output = 1
      else :
        outfile.append(ln)
      if add_swap_input :
        idmo = modid_exp.match(ln)
        if idmo != None :
          port_swap_input_ids.append(idmo.group(1))
      if add_swap_output :
        idmo = modid_exp.match(ln)
        if idmo != None :
          port_swap_output_ids.append(idmo.group(1))

          
    elif cmo != None :
      newln = ln
      for i in port_swap_input_ids :
        if i == cmo.group(4) :
          newtp = cmo.group(5)
          if int(cmo.group(5)) == 0 :
            newtp = "2"
          if int(cmo.group(5)) == 2 :
            newtp = "0"
          newln =  '%s $%s %s $%s %s%s\n' % (cmo.group(1), cmo.group(2),
                                             cmo.group(3), cmo.group(4),
                                             newtp, cmo.group(6))

      for o in port_swap_output_ids :
        if o == cmo.group(2) :
          newfp = cmo.group(3)
          if int(cmo.group(3)) == 0 :
            newfp = "2"
          if int(cmo.group(3)) == 2 :
            newfp = "0"
          newln =  '%s $%s %s $%s %s%s\n' % (cmo.group(1), cmo.group(2),
                                             newfp, cmo.group(4),
                                             cmo.group(5), cmo.group(6))
      outfile.append(newln)
    else :
      t = fix_old_variables(ln)
      if t[0] :
        outfile.append(t[1])
      else :
        outfile.append(ln)

  f = open(fn, 'w')
  f.writelines(outfile)


fld_exp = re.compile('(SCIRun)::(Fields)\w+::(\w+)')
        
def build_remap(fn) :
  global remap
  
  lines = []
  try:
    f = open(fn)
    lines = f.readlines()
  except:
    print "Failure reading %s" % fn
    return False
  
  lcount = 1
  for l in lines :
    mo = from_to_exp.match(l)
    if mo == None :
      print "line %d: Invalid input.\n\t(%s)" % (lcount, l[:-1])
    else :
      if mo.group(2) != '###' and mo.group(1) != mo.group(2) :
        remap[mo.group(1)] =  mo.group(2)
        mo1 = fld_exp.match(mo.group(1))
        if mo1 != None :
          bwcompat = "%s::%s::%s" % (mo1.group(1), mo1.group(2), mo1.group(3))
          remap[bwcompat] = mo.group(2)
    lcount = lcount + 1
  if len(remap) :
    return True
  return False

svnp = re.compile("\.svn")
def repair_each(arg, dirname, names) :
  mo = svnp.search(dirname)
  if mo != None :
    return # do not work in .svn dirs!

  for f in names :
    fn = dirname + os.sep + f
    if os.path.isfile(fn) and fn[-4:] == '.net' :
      parse_net(fn)

def fix_old_variables(ln) :
  var_rn_exp = re.compile('(.*)on release(.*)')
  mo = var_rn_exp.match(ln)
  if mo != None :
    newln = ln
    print "found"
    newln = "%sOn Release%s\n" % (mo.group(1), mo.group(2))
    return (True, newln)
  
  newln = ln
  return (False, newln)


if __name__ == '__main__' :
  fn = ''
  remapping = ''
  try:
    fn = sys.argv[1]
    remapping = sys.argv[2]
  except:
    print "problem with arguments"
    print "python remap_modules_in_net.py <path/to/.srn> <path/to/remapping_file>"
    sys.exit(1)

  if build_remap(remapping) :

    if os.path.isdir(fn) :
      os.path.walk(fn, repair_each, ())
      sys.exit(1)
  
    parse_net(fn)
  
  else :
    print "Bad or empty remapping table, exiting..."
