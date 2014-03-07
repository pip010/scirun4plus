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

IF(WIN32)

  GET_FILENAME_COMPONENT(PROG_PATH
     "[HKEY_LOCAL_MACHINE\\Software\\Microsoft\\Windows\\CurrentVersion;ProgramFilesDir]" ABSOLUTE
  )

  GET_FILENAME_COMPONENT(PROG_PATH_x86
        "[HKEY_LOCAL_MACHINE\\Software\\Microsoft\\Windows\\CurrentVersion;ProgramFilesDir (x86)]" ABSOLUTE
  )

  FILE(GLOB INCREDIBUILD_SETUP_MATCHES "${PROG_PATH}/Xoreax" "${PROG_PATH}/Xoreax/IncrediBuild" "${PROG_PATH_x86}/Xoreax" "${PROG_PATH_x86}/Xoreax/IncrediBuild")
  FIND_PROGRAM(INCREDIBUILD_SETUP_COMMAND NAMES buildconsole.exe BuildConsole.exe PATHS ${INCREDIBUILD_SETUP_MATCHES})

  MARK_AS_ADVANCED(INCREDIBUILD_SETUP_COMMAND)

ELSE(WIN32)
  MESSAGE(SEND_ERROR "IncrediBuild is only available for Windows.")
ENDIF(WIN32)
