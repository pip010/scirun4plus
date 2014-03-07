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
#    File   : gen_algo_info.py
#    Author : Martin Cole
#    Date   : Thu Oct 25 10:58:16 2007

from xml.dom import minidom
import os, sys
import re
import glob

class AlgoInfo :
    """
    storage for all the information in an algorithm needed to
    produce its interface.
    
    """
    def __init__(self) :
        """
        default values for the class.
        """
        self.name_ = ''
        self.category_ = ''
        self.version_ = 0.0
        self.params_ = []
        self.inputs_ = []
        self.outputs_ = []
        
    def print_self(self) :
        """
        debugging print for the class.
        """
        print "self.name_ = %s" % self.name_
        print "self.category_ = %s" % self.category_
        print "self.version_ = %f" % self.version_
        print "self.params_:"
        for p in self.params_ :
            print "  %s, %s, %s" % p

        print "self.inputs_:"
        for i in self.inputs_ :
            print "  %s, %s" % i

        print "self.outputs_:"
        for o in self.outputs_ :
            print "  %s, %s" % o

    def sort_params(self) :
        self.params_.sort()


def parse_algo_info(fn) :
    """
    parse_algo_info(fn) :
    parses fn, which must be a path to an AlgoInfo xml file, and return
    and AlgoInfo object with all the info from the file cached within.

    """
    a = AlgoInfo()    
    doc = minidom.parse(fn)
    rn = doc.documentElement

    a.name_ = rn.getAttribute('name')
    a.category_ = rn.getAttribute('category')
    a.version_ = float(rn.getAttribute('version'))

    pnodes = rn.getElementsByTagName('parameter')
    for pn in pnodes :
        pname = pn.getAttribute('name')
        typ = ''
        dv = None
        #print pname
        tnodes = pn.getElementsByTagName('type')
        for tn in tnodes :
            for n in tn.childNodes :
                t = n.nodeName
                if t == 'int'  :
                    typ = t
                    dv = int(n.getAttribute('default'))
                elif t == 'float' or t == 'double' :
                    typ = t
                    dv = float(n.getAttribute('default'))
                elif t == 'bool'  :
                    typ = t
                    tmp = int(n.getAttribute('default'))
                    dv = "false"
                    if bool(tmp) :
                        dv = "true"
                elif t == 'string'  :
                    typ = t
                    dv = n.getAttribute('default')
                elif t == '#text' :
                    pass
                else :
                    print "ERROR: %s unhandled." % t

        a.params_.append((pname, typ, dv))

    dnodes = rn.getElementsByTagName('data')
    for dn in dnodes :
        inodes = dn.getElementsByTagName('input')
        for i in inodes :
            nm = i.getAttribute('name')
            t = i.getAttribute('type')
            if i.hasAttribute('seq') :
                t = 'sequence(%s)' % t
                nm = "%s_seq" % nm


            a.inputs_.append((nm, t))
        onodes = dn.getElementsByTagName('output')
        for o in onodes :
            nm = o.getAttribute('name')
            t = o.getAttribute('type')
            a.outputs_.append((nm, t))

    return a

if __name__ == "__main__" :
    a = parse_algo_info(sys.argv[1])
    a.print_self()
