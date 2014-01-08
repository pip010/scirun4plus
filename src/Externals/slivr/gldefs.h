//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : gldefs.h
//    Author : Martin Cole
//    Date   : Fri Dec 29 11:26:44 2006

#if ! defined(gldefs_h)
#define gldefs_h


#include <stdio.h>
#include <GL/glew.h>

#if defined(OGL_DBG)
#  define CHECK_OPENGL_ERROR()                                          \
   {                                                                    \
    GLenum errCode;                                                     \
    int cnt = 1;                                                        \
    while ((errCode = glGetError()) != GL_NO_ERROR) {                   \
      fprintf(stderr, "OpenGL Error(%d):%s:%d '%s'\n", errCode, __FILE__,   \
              __LINE__, (const char*)gluErrorString(errCode));          \
    }                                                                   \
  }
#else
#  define CHECK_OPENGL_ERROR()     
#endif

#if defined(OGL_DBG)
#  define CHECK_FRAMEBUFFER_STATUS()                                       \
{                                                                          \
  GLenum status;                                                           \
  status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);                \
  switch(status) {                                                         \
  case GL_FRAMEBUFFER_COMPLETE_EXT:                                        \
    break;                                                                 \
  case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:                           \
    fprintf(stderr, "%s:%d - %s\n", __FILE__, __LINE__,                    \
            "GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT");                   \
    break;                                                                 \
  case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:                   \
    fprintf(stderr, "%s:%d - %s\n", __FILE__, __LINE__,                    \
            "GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT");           \
    break;                                                                 \
  case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:                           \
    fprintf(stderr, "%s:%d - %s\n", __FILE__, __LINE__,                    \
            "GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT");                   \
    break;                                                                 \
  case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:                              \
    fprintf(stderr, "%s:%d - %s\n", __FILE__, __LINE__,                    \
            "GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT");                      \
    break;                                                                 \
  case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:                          \
    fprintf(stderr, "%s:%d - %s\n", __FILE__, __LINE__,                    \
            "GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT");                  \
    break;                                                                 \
  case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:                          \
    fprintf(stderr, "%s:%d - %s\n", __FILE__, __LINE__,                    \
            "GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT");                  \
    break;                                                                 \
  case GL_FRAMEBUFFER_UNSUPPORTED_EXT:                                     \
    fprintf(stderr, "%s:%d - %s\n", __FILE__, __LINE__,                    \
            "GL_FRAMEBUFFER_UNSUPPORTED_EXT");                             \
    break;                                                                 \
  default:                                                                 \
    /* programming error; will fail on all hardware */                     \
    assert(0);                                                             \
  }                                                                        \
}                                                                          \

#else
#  define CHECK_FRAMEBUFFER_STATUS() 
#endif

#endif
