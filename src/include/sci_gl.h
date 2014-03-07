/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/

#if !defined(SCI_GL_H)
#define SCI_GL_H

// Uncomment this line to turn on verbose OPENGL ERROR checking.
//#define OGL_DBG 1
#include <slivr/gldefs.h>

#if defined(_WIN32)

#include <windows.h>

#ifndef GLAPIENTRY
#  define GLAPIENTRY APIENTRY
#endif
#undef SCISHARE
#define SCISHARE __declspec(dllimport)

#else // ! _WIN32

#undef SCISHARE
#define SCISHARE
#ifndef GLAPIENTRY
  #define GLAPIENTRY
#endif

#ifndef GLAPI
  #define GLAPI
#endif

#endif // _WIN32

#include <GL/glew.h>

#ifdef __cplusplus
  extern "C" {
#endif // __cplusplus
  extern SCISHARE int sci_glew_init();
#ifdef __cplusplus
  }
#endif // __cplusplus

#endif  // #define SCI_GL_H
