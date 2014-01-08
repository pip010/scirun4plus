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
#    File   : remap_modules_in_srn.py
#    Author : Martin Cole
#    Date   : Thu Oct 12 10:06:20 2006

from xml.dom.minidom import getDOMImplementation
from xml.dom import *
import sys
import re
import os

do_port_swap = 0

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

def get_document(fn) :
  impl = getDOMImplementation()
  doc = minidom.parse(fn)
  return doc

  
def swap_names(m) :
  global remap
  global do_port_swap
  add_swap_input = 0
  add_swap_output = 0
  
  pkg = m.attributes['package'].nodeValue
  cat = m.attributes['category'].nodeValue
  nam = m.attributes['name'].nodeValue

  if do_port_swap and nam == "CreateFieldFromNrrdData" :
    add_swap_input = 1

  if do_port_swap and nam == "SplitFieldIntoNrrdData" :
    add_swap_output = 1
    
     
  key = "%s::%s::%s" % (pkg, cat, nam)
  if remap.has_key(key) :
    mo = pcm_exp.match(remap[key])
    m.setAttribute('package', mo.group(1))
    m.setAttribute('category', mo.group(2))
    m.setAttribute('name', mo.group(3))
   
    if mo.group(3) == "CreateFieldFromNrrdData" :
      add_swap_input = 1

    if mo.group(3) == "SplitFieldIntoNrrdData" :
      add_swap_output = 1

  if add_swap_input :
    mod_id = m.attributes['id'].nodeValue
    port_swap_input_ids.append(mod_id)

  if add_swap_output :
    mod_id = m.attributes['id'].nodeValue
    port_swap_output_ids.append(mod_id)

def swap_ports(n) :
  global port_swap_input_ids
  global port_swap_output_ids

  for mod_id in port_swap_input_ids :
    if mod_id == n.attributes['to'].nodeValue :
      toport = int(n.attributes['toport'].nodeValue)
      if toport == 0:
        n.setAttribute('toport', "2")
      if toport == 2:
        n.setAttribute('toport', "0")

  for mod_id in port_swap_output_ids :
    if mod_id == n.attributes['from'].nodeValue :
      fromport = int(n.attributes['fromport'].nodeValue)
      if fromport == 0:
        n.setAttribute('fromport', "2")
      if fromport == 2:
        n.setAttribute('fromport', "0")

def repair_file(fn) :
  global port_swap_input_ids
  global port_swap_output_ids

  port_swap_input_ids = []
  port_swap_output_ids = []

  doc = get_document(fn)
  walk_modules(doc.documentElement, swap_names)
  if len(port_swap_input_ids) > 0 or len(port_swap_output_ids) > 0:
    walk_connections(doc.documentElement, swap_ports)
  fd = open(fn, "w")
  doc.writexml(fd, "", "", '')
  doc.unlink()

def recurse_modules(children, action) :
  for c in children :
    if c.nodeName == 'modules' :
      mc = c.childNodes
      for m in mc :
        if m.nodeName == 'subnet' :
          for snc in m.childNodes :
            if snc.nodeName == 'network' :
              recurse_modules(snc.childNodes, action)
        if m.nodeName == 'module' :
          action(m)


def recurse_connections(children, action) :
  global port_swap_input_ids
  global port_swap_output_ids
  for c in children :
    if c.nodeName == 'modules' :
      mc = c.childNodes
      for m in mc :
        if m.nodeName == 'subnet' :
          for snc in m.childNodes :
            if snc.nodeName == 'network' :
              # cache any module ids 
              il = []
              ol = []
              for i in port_swap_input_ids :
                il.append(i)
              for o in port_swap_output_ids :
                ol.append(o)
              # clear ids for subnet
              port_swap_input_ids = []
              port_swap_output_ids = []
              
              recurse_connections(snc.childNodes, action)
              
              #restore cached ids
              port_swap_input_ids = []
              port_swap_output_ids = []
              for i in il :
                port_swap_input_ids.append(i)
              for o in ol :
                port_swap_output_ids.append(o)
              
    if c.nodeName == 'connections' :
      conc = c.childNodes
      for n in conc :
        if n.nodeName == 'connection' :          
          action(n)

def walk_modules(doc_node, action) :
  children = doc_node.childNodes
  recurse_modules(children, action)

def walk_connections(doc_node, action) :
  children = doc_node.childNodes
  recurse_connections(children, action)

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
    if os.path.isfile(fn) and fn[-4:] == '.srn' :
      repair_file(fn)

if __name__ == '__main__' :
  fn = ''
  remapping = ''
  try:
    fn = sys.argv[1]
    remapping = sys.argv[2]
  except:
    print "WARNING:   This script will overwrite the file or files specified"
    print "           Please backup your .srn files first."
    print "Usage:"
    print "python remap_modules_in_srn.py <path/to/.srn> <path/to/remapping_file>"
    print ''
    s = '\tIf the first argument is a directory, then work recursively \non all .srn files under that directory, otherwise the argument should\nspecify a single .srn file. The remapping_file contains the mapping\nfrom old to new module info.  The file should have one remapping per\nline like so: \n\noldpackage:::oldcategory::oldmodule  newpackage::newcategory::newmodule\n\n\tSee SCIRun/src/scripts/module-remapping.txt.\n\n'
    print s
    sys.exit(1)

  if build_remap(remapping) :

    try:
      if os.path.isdir(fn) :
        os.path.walk(fn, repair_each, ())
      else :
        repair_file(fn)
    except:
      print "Error: Is %s really a .srn file?" % fn
  else :
    print "Error: Bad or empty remapping file, exiting..."
