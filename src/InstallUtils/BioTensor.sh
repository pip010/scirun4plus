#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2010 Scientific Computing and Imaging Institute,
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

# This script is just a short cut to run the BioTensor SCIRun Power App
# for OS X binary distributions

APP_BASE=`dirname "$0"`
DIR=$APP_BASE/@SCIRUN_APP_NAME@/Contents/Resources
VERSION=`uname -r | cut -d . -f 1`

# Just loading the needed packages.
export SCIRUN_LOAD_PACKAGE="SCIRun,Teem"

export SCIRUN_OBJDIR="${DIR}"
export SCIRUN_SRCDIR="${DIR}/src"
export SCIRUN_THIRDPARTY_DIR="${DIR}"
export SCIRUN_ITCL_WIDGETS="${SCIRUN_THIRDPARTY_DIR}/lib/iwidgets/scripts"
export SCIRUN_TCL_LIBRARY="${SCIRUN_THIRDPARTY_DIR}/lib/sciruntcl1.0"
export TCL_LIBRARY="${SCIRUN_THIRDPARTY_DIR}/lib/tcl"

if [ "$VERSION" -lt "9" ]; then
  # TIGER WE NEED TO START X11
  export DISPLAY=:0
  if [ -d /Applications/Utilities/X11.app ]; then
      /usr/bin/open -a /Applications/Utilities/X11.app 
  elif [ -d /Applications/X11.app ]; then
      /usr/bin/open -a /Applications/X11.app 
  fi
  echo "Starting BioTensor"
  /usr/X11R6/bin/xterm -geometry 60x20 -e "env;\"$SCIRUN_OBJDIR/bin/scirun\" $SCIRUN_TCL_LIBRARY/BioTensor.app $1 || sleep 60";
  echo "Finished BioTensor"
else
  echo "Starting BioTensor version"
  $SCIRUN_OBJDIR/bin/scirun $SCIRUN_TCL_LIBRARY/BioTensor.app $1
  echo "Finished BioTensor version"
  exit 0
fi
