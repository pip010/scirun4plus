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

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(VERSION 2.6)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(SEG3D)

#-----------------------------------------------------------------------------
# SEG3D version number.

SET(SEG3D_VERSION_MAJOR "1" CACHE INTERNAL "Seg3D major version number." FORCE)
SET(SEG3D_VERSION_MINOR "13" CACHE INTERNAL "Seg3D minor version number." FORCE)
SET(SEG3D_VERSION_PATCH "0" CACHE INTERNAL "Seg3D patch number." FORCE)

SET(SEG3D_VERSION_STRING "${SEG3D_VERSION_MAJOR}.${SEG3D_VERSION_MINOR}.${SEG3D_VERSION_PATCH}" PARENT_SCOPE)

#-----------------------------------------------------------------------------
# Test for some required system information.
INCLUDE (CMakeBackwardCompatibilityC)

#-----------------------------------------------------------------------------
# Configure testing support.  This should be done before any
# MESSAGE(FATAL_ERROR ...) commands are invoked.
INCLUDE( CTest )
IF(BUILD_TESTING)
  ENABLE_TESTING()
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)

#----------------------------------------------------------------------------
# Output directories.
IF(NOT BUILD_FROM_TOP_LEVEL)
  SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
  SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
  #output directory for installing all include files.")
  MARK_AS_ADVANCED(CMAKE_LIBRARY_OUTPUT_DIRECTORY CMAKE_RUNTIME_OUTPUT_DIRECTORY)
ENDIF(NOT BUILD_FROM_TOP_LEVEL)

SET(SEG3D_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET (EXECUTABLE_OUTPUT_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} CACHE INTERNAL "Output directory for Seg3D executable")

# DV Find SCIRun Core stuff
SET(CMAKE_MODULE_PATH ${SEG3D_SOURCE_DIR})

IF (WITH_ITK)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/itk_temp/Source/)
ELSE(WITH_ITK)
  MESSAGE(FATAL_ERROR "Cannot build Seg3D without ITK. Please set WITH_ITK to ON.")
ENDIF(WITH_ITK)

IF(WIN32)
  IF(NOT BUILD_CORE_LIBS_SHARED)
    ADD_DEFINITIONS(-DBUILD_SCIRUN_STATIC)
  ENDIF(NOT BUILD_CORE_LIBS_SHARED)
#  This defines M_PI for windows
  ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
ENDIF(WIN32)


SET(Seg3D_SRCS 
  main.cc
  AutoviewTool.cc 
  BrushFloodFill.cc 
  BrushTool.cc 
  CLUTLevelsTool.cc
  CropTool.cc
  ArithmeticTool.cc
  FloodFillCopyTool.cc
  FloodfillTool.cc
  InhomogeneityCorrectionFilter.cc
  ITKBinaryDilateErodeImageFilterTool.cc
  ITKConfidenceConnectedImageFilterTool.cc
  ITKCurvatureAnisotropicDiffusionImageFilterTool.cc
  ITKDiscreteGaussianImageFilterTool.cc
  ITKCannyEdgeDetectionImageFilterTool.cc
  ITKGradientMagnitudeImageFilterTool.cc
  ITKHoleFillImageFilterTool.cc
  ITKNeighborhoodConnectedImageFilterTool.cc
  ITKOtsuImageFilterTool.cc
  ITKSpeedToPathGradientDescentImageFilterTool.cc
  ITKThresholdSegmentationLevelSetImageFilterTool.cc
  ITKSpeedToPathRegularStepGradientDescentImageFilterTool.cc
  ITKSpeedToPathIterateNeighborhoodImageFilterTool.cc
  KeyToolSelectorTool.cc
  LayerButton.cc
  LayerMergeTool.cc
  LineTool.cc
  NrrdVolume.cc
  Painter.cc 
  PainterSignalTargets.cc 
  PanTool.cc 
  PointerToolSelectorTool.cc
  PolylineTool.cc
  ProbeTool.cc 
  ResampleTool.cc
  SeedTool.cc
  Seg3DFrame.cc
  Seg3DVersion.cc
  Seg3DwxGuiUtils.cc
  SessionReader.cc
  SessionWriter.cc
  SliceWindow.cc
  SpeedToPathTool.cc
  StatisticsTool.cc
  ThresholdTool.cc
  WindowLevelTool.cc
  Undo.cc
  VolumeFilter.cc
  VolumeOps.cc
  VolumeSlice.cc
  WXOpenGLContext.cc
  ZoomTool.cc
  GuiCode/seg3devents.cpp
  GuiCode/brushpanel.cpp
  GuiCode/histoeqfilter.cpp
  GuiCode/intensitycorrectionfilterpanel.cpp
  GuiCode/itkcurvatureanistopricdiffusionfilter.cpp
  GuiCode/itkconfidenceconnectedfilter.cpp
  GuiCode/itkneighbourhoodconnectedfilter.cpp
  GuiCode/itkbinarydilateerodefilter.cpp
  GuiCode/itkthresholdsegmentationlevelsetfilter.cpp
  GuiCode/cropvolume.cpp
  GuiCode/arithmeticvolume.cpp
  GuiCode/cursorinformation.cpp
  GuiCode/fliptool.cpp
  GuiCode/inverttool.cpp
  GuiCode/resampletool.cpp
  GuiCode/medianfiltertool.cpp
  GuiCode/maskfilter.cpp
  GuiCode/masklogicalfilter.cpp
  GuiCode/measurementtool.cpp
  GuiCode/movescaletoolpanel.cpp
  GuiCode/optionlessfilter.cpp
  GuiCode/polylinetoolpanel.cpp
  GuiCode/thresholdtoolpanel.cpp
  GuiCode/windowleveltoolpanel.cpp
  GuiCode/itkDiscreteGaussianImageFilter.cpp
  GuiCode/itkCannyEdgeDetectionImageFilter.cpp
  GuiCode/itkspeedtopathgradientdescentfilter.cpp
  GuiCode/itkspeedtopathregularstepgradientdescentfilter.cpp
  GuiCode/itkspeedtopathiterateneighborhoodfilter.cpp
)

IF(WIN32)
  SET(Seg3D_SRCS ${Seg3D_SRCS} seg3d.rc)
ENDIF(WIN32)

SET(Seg3D_LIBS
  ${XML_LIBRARY} 
  ${M_LIBRARY} 
  ${GL_LIBRARY} 
  ${SCI_TEEM_LIBRARY}
  ${SLIVR_LIBRARY}
  ${SCI_ITK_LIBRARIES}
  Core_Algorithms_DataIO
  Core_Matlab 
  Core_Basis 
  Core_Datatypes 
  Core_Events 
  Core_Exceptions 
  Core_Geom 
  Core_Geometry 
  Core_Math 
  Core_Persistent 
  Core_Skinner 
  Core_Thread 
  Core_Util 
  Core_Volume
  ITKCommon
  ITKBasicFilters
  ITKNumerics
)

IF(X11_FOUND)
  SET(Seg3D_LIBS ${Seg3D_LIBS} ${GLU_LIBRARY})
ENDIF(X11_FOUND)

ADD_EXECUTABLE(Seg3D MACOSX_BUNDLE ${Seg3D_SRCS})
TARGET_LINK_LIBRARIES(Seg3D ${Seg3D_LIBS})


IF(wxWidgets_FOUND)
  TARGET_LINK_LIBRARIES( Seg3D ${wxWidgets_LIBRARIES})
ENDIF(wxWidgets_FOUND)

IF(APPLE AND NOT X11_FOUND)
  SET_TARGET_PROPERTIES(Seg3D PROPERTIES LINK_FLAGS "-framework Carbon")
ENDIF(APPLE AND NOT X11_FOUND)


IF(APPLE)
  SET(DATA_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Seg3D.app/Contents/MacOS/data)
ELSE(APPLE)
  IF(WIN32)
    IF(MSVC)
       # If building on Windows using Visual Studio (the recommended method), use the VS macros
       # to figure out where the executables and libraries are being generated.
      SET(DATA_PATH "\$\(TargetDir\)\\data")
   ELSE(MSVC)
      # Although Visual Studio is the recommended build tool, it may be helpful to set up
      # an alternative for other tools.
      FILE(TO_NATIVE_PATH "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/data" DATA_PATH_WINDOWS)
      SET(DATA_PATH ${DATA_PATH_WINDOWS})
    ENDIF(MSVC)
  ELSE(WIN32)
    SET(DATA_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/data)
  ENDIF(WIN32)
ENDIF(APPLE)

ADD_CUSTOM_COMMAND(TARGET Seg3D
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND}
                   ARGS -E make_directory
                   ${DATA_PATH})
                   
IF(APPLE_LEOPARD AND WITH_X11)
  SET_TARGET_PROPERTIES(Seg3D PROPERTIES LINK_FLAGS  "${GL_LINK} ${X_LINK}")
ENDIF(APPLE_LEOPARD AND WITH_X11)

ADD_CUSTOM_COMMAND(TARGET Seg3D
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND}
                   ARGS -E copy_directory 
                   ${SEG3D_SOURCE_DIR}/data/ 
                   ${DATA_PATH}/)

ADD_CUSTOM_COMMAND(TARGET Seg3D
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND}
                   ARGS -E remove_directory 
                   ${DATA_PATH}/.svn)

IF(APPLE)
  FILE(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app)
  FILE(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents)
  FILE(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/MacOS)
  FILE(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/MacOS/data)
  FILE(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/Resources)

  FILE(GLOB data_files "${SEG3D_SOURCE_DIR}/data/*.skin")
  FOREACH(data_file ${data_files})
    GET_FILENAME_COMPONENT(data_name ${data_file} NAME)
    CONFIGURE_FILE(${SEG3D_SOURCE_DIR}/data/${data_name}  ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/MacOS/data/${data_name} COPYONLY)
  ENDFOREACH(data_file ${data_files})

  FILE(GLOB data_files "${SEG3D_SOURCE_DIR}/data/*.png")
  FOREACH(data_file ${data_files})
    GET_FILENAME_COMPONENT(data_name ${data_file} NAME)
    CONFIGURE_FILE(${SEG3D_SOURCE_DIR}/data/${data_name}  ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/MacOS/data/${data_name} COPYONLY)
  ENDFOREACH(data_file ${data_files})

  FILE(GLOB data_files "${SEG3D_SOURCE_DIR}/data/*.xrc")
  FOREACH(data_file ${data_files})
    GET_FILENAME_COMPONENT(data_name ${data_file} NAME)
    CONFIGURE_FILE(${SEG3D_SOURCE_DIR}/data/${data_name}  ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/MacOS/data/${data_name} COPYONLY)
  ENDFOREACH(data_file ${data_files})

  CONFIGURE_FILE(${SEG3D_SOURCE_DIR}/OSX/Info.plist ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/Info.plist COPYONLY)
  CONFIGURE_FILE(${SEG3D_SOURCE_DIR}/OSX/Resources/Defaults.plist ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/Resources/Defaults.plist COPYONLY)
  CONFIGURE_FILE(${SEG3D_SOURCE_DIR}/OSX/Resources/Seg3D.icns ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/Resources/Seg3D.icns COPYONLY)
  CONFIGURE_FILE(${SEG3D_SOURCE_DIR}/data/splash-seg3d.png ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/Resources/splash-seg3d.png COPYONLY)
  CONFIGURE_FILE(${SEG3D_SOURCE_DIR}/data/Seg3D.xrc ${EXECUTABLE_OUTPUT_PATH}/Seg3D.app/Contents/Resources/Seg3D.xrc COPYONLY)

ENDIF(APPLE)
