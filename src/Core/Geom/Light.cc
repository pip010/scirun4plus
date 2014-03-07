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



/*
 *  Light.cc: Base class for light sources
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   September 1994
 *
 */

#include <sci_gl.h>
#include <sci_glx.h>

#include <Core/Geom/Light.h>
#include <Core/Thread/MutexPool.h>

namespace SCIRun {

MutexPool light_lock_pool("LightObj pool", 10);

Light::Light(const std::string& name, bool on, bool transformed)
  : UsedWithLockingHandle<Mutex&>(light_lock_pool.getMutex()),
    name(name), on(on), transformed( transformed )
{
}

Light::~Light()
{
}

void
Light::opengl_reset_light(int i)
{
  float f[4];
  f[0]=0.0; f[1]=0.0; f[2]=0.0; f[3]=1.0;
  glLightfv((GLenum)(GL_LIGHT0+i), GL_AMBIENT, f);

  if ( i != 0 )
  {
    glLightfv((GLenum)(GL_LIGHT0+i), GL_DIFFUSE, f);
    glLightfv((GLenum)(GL_LIGHT0+i), GL_SPECULAR, f);
  }
  else
  {
    f[0] = 1.0; f[1]=1.0; f[2]=1.0; f[3]=1.0;
    glLightfv((GLenum)(GL_LIGHT0+i), GL_DIFFUSE, f);
    glLightfv((GLenum)(GL_LIGHT0+i), GL_SPECULAR, f);
  }
  f[0]=0.0; f[1]=0.0; f[2]=-1.0; f[3]=1.0;
  glLightfv((GLenum)(GL_LIGHT0+i), GL_POSITION, f);

  glLightfv((GLenum)(GL_LIGHT0+i), GL_SPOT_DIRECTION, f);
  f[0] = 180.0;
  glLightfv((GLenum)(GL_LIGHT0+i), GL_SPOT_CUTOFF, f);
  f[0] = 0.0;
  glLightfv((GLenum)(GL_LIGHT0+i), GL_SPOT_EXPONENT, f);
  glLightfv((GLenum)(GL_LIGHT0+i), GL_LINEAR_ATTENUATION, f);
  glLightfv((GLenum)(GL_LIGHT0+i), GL_QUADRATIC_ATTENUATION, f);
  f[0] = 1.0;
  glLightfv((GLenum)(GL_LIGHT0+i), GL_CONSTANT_ATTENUATION, f);
}


} // End namespace SCIRun


