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

# Python script for removing the first description tag from xml files.
#
# usage 'python remove_description.py foo.xml > new.xml

import sys, regex, string, regsub

rex_start = regex.compile('.*<description>.*');
rex_end = regex.compile('.*</description>.*');

def cont(file, line):
    while 1:
        n = rex_end.match(line);
        if n >= 0:
            break
        line = file.readline();
        if not line: break

def gen_newfile(filename):
    file = open(filename, 'r')
    first = 1
    while 1:
        line = file.readline()
        if not line: break

        n = rex_start.match(line)
        if n >= 0 and first:
            cont(file, line);
            first = 0
            continue

        print line,

if __name__ == '__main__':
    import sys
    gen_newfile(sys.argv[1])
        
