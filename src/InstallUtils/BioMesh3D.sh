#!/usr/bin/env bash

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

echo "\$1 = "$1
DIR=$1/Contents/Resources
VERSION=`uname -v|sed 's/\./ /g' | awk '{print $4}'`

# TIGER WE NEED TO START X11
if [[ $VERSION -lt 9 ]]; then
  export DISPLAY=:0
  if [[ -d /Applications/Utilities/X11.app ]]; then
      /usr/bin/open -a /Applications/Utilities/X11.app 
  elif [[ -d /Applications/X11.app ]]; then
      /usr/bin/open -a /Applications/X11.app 
  fi
fi

/usr/X11R6/bin/xterm -geometry 60x20 -e "env;python $DIR/bin/FEMesher/BuildMesh.py $2 || sleep 60";
echo "Finished BioMesh3D"

 



