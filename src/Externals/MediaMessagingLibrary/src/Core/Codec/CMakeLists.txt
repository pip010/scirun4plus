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

SET(CORE_CODEC_SRCS
Codec.cc
Codec.h
CompressedImage3.cc
CompressedImage3.h
Decoder.cc
Decoder.h
Encoder.cc
Encoder.h
FFMPEGDecoder.cc
FFMPEGDecoder.h
Image3.cc
Image3.h
MPEGDecoder.cc
MPEGDecoder.h
MPEGEncoder.cc
MPEGEncoder.h
PNGDecoder.cc
PNGDecoder.h
PNGEncoder.cc
PNGEncoder.h
VP8Decoder.cc
VP8Decoder.h
VP8Encoder.cc
VP8Encoder.h
)

SET(CORE_CODEC_LIBRARY MML_Core_Codec CACHE INTERNAL "Codec library name." FORCE)
SET(CORE_CODEC_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE INTERNAL "MediaCodec include." FORCE)


##################################################
# Build static library
##################################################

INCLUDE_DIRECTORIES( ${CORE_CODEC_INCLUDE} ${FFMPEG_INCLUDE} ${SCI_PNG_INCLUDE} ${SCI_VPX_INCLUDE} )

ADD_LIBRARY(${CORE_CODEC_LIBRARY} STATIC ${CORE_CODEC_SRCS})

IF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
  SET(ENABLE_PIC ON)
ENDIF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)

IF(CMAKE_COMPILER_IS_GNUCC)
IF(ENABLE_PIC)
  SET_TARGET_PROPERTIES(${CORE_CODEC_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-w -fPIC"
                       )
ELSE(ENABLE_PIC)
  SET_TARGET_PROPERTIES(${CORE_CODEC_LIBRARY}
                        PROPERTIES
                        COMPILE_FLAGS "-w")
ENDIF(ENABLE_PIC)
ENDIF(CMAKE_COMPILER_IS_GNUCC)

TARGET_LINK_LIBRARIES(${CORE_CODEC_LIBRARY} MML_Core_Utils ${SCI_PNG_LIBRARY} ${SCI_BOOST_LIBRARY} ${FFMPEG_LIBRARIES} ${SCI_VPX_LIBRARY} )
