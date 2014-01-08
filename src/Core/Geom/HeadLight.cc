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
 *  HeadLight.cc:  A Point light source
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

#include <Core/Geom/HeadLight.h>
#include <Core/Geom/View.h>

namespace SCIRun {

HeadLight::HeadLight(const std::string& name, const Color& c, bool on)
: Light(name,on,false), c(c)
{
}

HeadLight::~HeadLight()
{
}

void
HeadLight::opengl_setup(const View& /*view*/, DrawInfoOpenGL*, int& idx)
{
  if ( on )
  {
    const int i = idx++;
    float f[4];

    opengl_reset_light(i);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity(); // Never transformed: its a headlight afterall--Kurt

    f[0] = f[1] = f[3] = 0.0;
    f[2] = 1.0;
    glLightfv((GLenum)(GL_LIGHT0+i), GL_POSITION, f);
    c.get_color(f);     
    glLightfv((GLenum)(GL_LIGHT0+i), GL_DIFFUSE, f);
    glLightfv((GLenum)(GL_LIGHT0+i), GL_SPECULAR, f);

    glPopMatrix();
  }
}

} // End namespace SCIRun

