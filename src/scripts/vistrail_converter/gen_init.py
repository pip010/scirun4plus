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
#    File   : gen_init.py
#    Author : Martin Cole
#    Date   : Tue Nov 27 14:31:45 2007


import os, sys
import re
import glob

sys.path.append(os.path.abspath(".."))
from AlgoInfo import *

#the package identifiers
vpversion = '0.9.1'
vpidentifier = 'edu.utah.sci.vistrails.scirun'
vpname = 'SCIRun'

def get_vt_core_type(t) :
    if t == 'string' :
        return "String"
    if t == 'double' :
        return "Float"
    if t == 'int' :
        return "Integer"
    if t == 'bool' :
        return "Boolean"
    
def wrap_algo(ai, cl, il, constants) :
    """
    giving an AlgoInfo instance (ai) append strings to the 2 input lists.
    cl is the list of strings that declare the class wrappings, il is the
    list that makes up the package initialize function.
    constants are datatypes, and categories. For data ports, and organizing
    wrapped modules into categories.
    """

    if (constants.count(ai.category_) == 0) :
        constants.append(ai.category_)


    #register module
    il.append("  reg.add_module(%s)\n" % ai.name_)

    #declare wrapped class
    cl.append("class %s(%s) :\n" % (ai.name_, ai.category_))
    cl.append("  def compute(self) :\n")
    cl.append("    p = sr_py.%sAlg()\n" % ai.name_)
    #handle each parameter
    for pn, typ, dv in ai.params_ :
        
        cl.append("    if self.hasInputFromPort('%s') :\n" % pn)     
        cl.append("      p.set_%s(self.getInputFromPort('%s'))\n" % (pn, pn))

        il.append("  reg.add_input_port(%s, '%s',\n" % (ai.name_, pn))
        il.append("                   ")
        t = get_vt_core_type(typ)
        il.append("(core.modules.basic_modules.%s, 'tip'), True)\n" % t)

    #handle ports
    #inputs
    call = ""
    for pn, typ in ai.inputs_ :
        if (constants.count(typ) == 0) :
            constants.append(typ)

        pnv = pn.replace(" ", "_")
        call = "%s%s, " % (call, pnv)
        if typ == "String" :
            cl.append("    %s = ''\n" % pnv)
        else :
            cl.append("    %s = 0\n" % pnv)
            
        cl.append("    if self.hasInputFromPort('%s') :\n" % pn) 
        cl.append("      %s = self.getInputFromPort('%s')\n" % (pnv, pn))
        
        il.append("  reg.add_input_port(%s, '%s',\n" % (ai.name_, pn))
        il.append("                   ")
        if (typ == "String") :
            il.append("(core.modules.basic_modules.String, 'tip'))\n")
        else:
            il.append('(%s, "%s"))\n' % (typ, typ))

    # the call to the algo
    cl.append("    results = p.execute(%s)\n" % call[:-2])

    #outputs
    i = 0
    for pn, typ in ai.outputs_ :
        if (constants.count(typ) == 0) :
            constants.append(typ)

        #pn = "%s out" % pn
        if (len(ai.outputs_) == 1):
            cl.append("    self.setResult('%s', results)\n" % pn)
        else:
            cl.append( "    self.setResult('%s', results[%s])\n" % (pn, i))
        i = i + 1
        
        il.append("  reg.add_output_port(%s, '%s',\n" % (ai.name_, pn))
        il.append("                   ")
        if (typ == "String") :
            il.append("(core.modules.basic_modules.String, 'tip'))\n")
        else:
            il.append('(%s, "%s"))\n' % (typ, typ))

    cl.append("\n")
    il.append("\n")

             

def write_init(cl, il, constants) :
    """ write_init() writes the necessary __init__.py for the package"""
    
    for new_constant in constants :
        s = "class %s(Constant):\n  def compute(self): \n"\
            "    pass\n\n" % new_constant
        cl.insert(0, s)  
        il.insert(5,"  reg.add_module(%s, abstract=True)\n" % new_constant)
    cl.insert(0, "\n")  

    # add viewercell registration
    il.append("\n\n")
    il.append("  import viewercell\n")
    il.append("  viewercell.registerSelf()\n")
    il.append("  import cm2viewcell\n")
    il.append("  cm2viewcell.registerSelf()\n")
    il.append("\n\n")
    il.append("def package_dependencies():\n")
    il.append("  return ['edu.utah.sci.vistrails.spreadsheet']\n")
    il.append("\n\n")
    il.append("def finalize():\n")
    il.append("  sr_py.terminate()\n")        
    il.append("  time.sleep(.5)\n")

    header = open("init.template").readlines()
##     for i in set(self.module_reader.imports):
##         header.append("import %s\n" % i)
    header.append("\n")

    header.append('version = "' + vpversion + '"\n')
    header.append('identifier = "' + vpidentifier + '"\n')
    header.append('name = "' + vpname + '"\n')

    outfile = open("__init__.py", "w")
    outfile.writelines(header)
    outfile.writelines(cl)
    outfile.writelines(il)
    outfile.close()
    print "__init__.py created in current directory."


def visit(all_info, dirname, names) :

    if dirname[-4:] != ".svn" :
        for n in names :
            if n[-4:] == ".xml" :
                fullname = dirname + os.sep + n
                all_info.append(parse_algo_info(fullname))


if __name__ == '__main__':
    info_path = os.path.abspath(".." + os.sep + ".." + os.sep + "Core" +
                                os.sep + "Algorithms" + os.sep +
                                "Interface" + os.sep + "Info")

    
    all_info = []
    os.path.walk(info_path, visit, all_info)

    class_lines = []
    init_lines = []
    init_lines.append("\ndef initialize(*args, **keywords):\n")
    init_lines.append("\n  env = []\n  for k in os.environ.keys() :\n")
    init_lines.append('    estr = "%s=%s" % (k, os.environ[k])\n')
    init_lines.append("    env.append(estr)\n  sr_py.init_sr_py(env)\n\n")
    init_lines.append("  reg = core.modules.module_registry\n\n")
    constants = []
    
    for ai in all_info :
        wrap_algo(ai, class_lines, init_lines, constants)

    write_init(class_lines, init_lines, constants)
