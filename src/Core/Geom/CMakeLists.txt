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

# Sources of Core/Geom classes

SET(Core_Geom_SRCS 
  BBoxCache.cc		                
  ColorMapTex.cc	                
  ColorMappedNrrdTextureObj.cc          
  DirectionalLight.cc	                
  DrawInfoOpenGL.cc	                
  FreeType.cc		                
  GeomBox.cc		                
  GeomColorMap.cc	                
  GeomColorMappedNrrdTextureObj.cc
  GeomCone.cc		                
  GeomContainer.cc	                
  GeomCull.cc	                
  GeomCylinder.cc	                
  GeomDL.cc                             
  GeomDisk.cc		                
  GeomGlyph.cc
  GeomGrid.cc		                
  GeomGroup.cc		                
  GeomLine.cc		  
  GeomMaterial.cc		    		
  GeomObj.cc
  GeomPick.cc		                
  GeomPoint.cc		    		
  GeomPolyline.cc	                
  GeomQuads.cc		                
  GeomRenderMode.cc
  GeomResourceManager.cc
  GeomScene.cc		                
  GeomSphere.cc		                
  GeomSticky.cc		                
  GeomSwitch.cc		                
  GeomStippleOccluded.cc		
  GeomTexRectangle.cc	    		
  GeomText.cc           
  GeomTorus.cc		                
  GeomTransform.cc	                
  GeomTri.cc		                
  GeomTriangles.cc	                
  GeomVertexPrim.cc	                
  GeomViewerItem.cc			
  HeadLight.cc		    		
  IndexedGroup.cc	    		
  Light.cc		    		
  Lighting.cc		    		
  NrrdBitmaskOutline.cc  		
  OpenGLViewport.cc	    		
  OpenGLContext.cc
  Path.cc		    	          	    		
  PointLight.cc		    		
  SpotLight.cc		    		
  TextRenderer.cc 	    		
  FontManager.cc 	    		
  TextureObj.cc 	    		
  TexSquare.cc	    		
  View.cc		    		
  X11Lock.cc
  OSXOpenGLContext.cc
)

IF(WIN32)
  SET(Core_Geom_SRCS ${Core_Geom_SRCS} Win32OpenGLContext.cc)
ELSE(WIN32)
  SET(Core_Geom_SRCS ${Core_Geom_SRCS} X11OpenGLContext.cc)
ENDIF(WIN32)

SCIRUN_ADD_LIBRARY(Core_Geom ${Core_Geom_SRCS})


TARGET_LINK_LIBRARIES( Core_Geom 
  ${THREAD_LIBRARY} 
  ${GL_LIBRARY}
  ${SCI_ZLIB_LIBRARY}
  ${SCI_PNG_LIBRARY}
  ${SCI_TEEM_LIBRARY}
  ${SCI_FREETYPE_LIBRARY}
  ${M_LIBRARY}
  ${SCI_GLEW_LIBRARY}
  Core_Exceptions
  Core_Persistent
  Core_Geometry
  Core_Datatypes
  Core_Math
  Core_Util
  Core_Thread
)

IF(X11_FOUND)
  TARGET_LINK_LIBRARIES( Core_Geom
    ${X_LIBRARY}
  )      
ENDIF(X11_FOUND)


IF(APPLE)
  IF(NOT APPLE_LEOPARD OR NOT WITH_X11)
    TARGET_LINK_LIBRARIES( Core_Geom "-framework Carbon -framework AGL")    
  ELSE(NOT APPLE_LEOPARD OR NOT WITH_X11)
    SET_TARGET_PROPERTIES(Core_Geom PROPERTIES LINK_FLAGS  "${GL_LINK} ${X_LINK}")
  ENDIF(NOT APPLE_LEOPARD OR NOT WITH_X11)

ENDIF(APPLE)


IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Core_Geom)
ENDIF(BUILD_SHARED_LIBS)

