
#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2011 Scientific Computing and Imaging Institute,
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

##################################################
# Set sources
##################################################

SET(CORE_UTILS_SRCS
  AtomicCounter.h
  Buffer.h
  Buffer.cc
  EnumClass.h
  Exception.h
  Exception.cc
  FilesystemUtil.h
  FilesystemUtil.cc
  KeyValuePairList.h
  KeyValuePairList.cc
  Lockable.h
  Log.h
  Log.cc
  LogStreamer.h
  LogStreamer.cc
  Parallel.h
  Parallel.cc
  Runnable.h
  Runnable.cc
  ScopedCounter.h
  Singleton.h
  StackBasedVector.h
  StackVector.h
  StringUtil.h
  StringUtil.cc
  StringParser.h
  StringParser.cc
  StdVectorBuffer.h
  StdVectorBuffer.cc
  Timer.h
  Timer.cc  
)

##################################################
# Build static library
##################################################

ADD_LIBRARY(MML_Core_Utils STATIC ${CORE_UTILS_SRCS})

IF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
  SET(ENABLE_PIC ON)
ENDIF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)

IF(CMAKE_COMPILER_IS_GNUCC)
IF(ENABLE_PIC)
  SET_TARGET_PROPERTIES(MML_Core_Utils
                        PROPERTIES
                          COMPILE_FLAGS "-w -fPIC"
                       )
ELSE(ENABLE_PIC)
  SET_TARGET_PROPERTIES(MML_Core_Utils
                        PROPERTIES
                        COMPILE_FLAGS "-w")
ENDIF(ENABLE_PIC)
ENDIF(CMAKE_COMPILER_IS_GNUCC)


TARGET_LINK_LIBRARIES(MML_Core_Utils ${SCI_BOOST_LIBRARY})
