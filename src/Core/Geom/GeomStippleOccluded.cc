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
 *  GeomStippleOccluded: 
 *    Shows geometry w/ a stipple patter where it would
 *    normally be occluded because of the depth test
 *    Useful for widgets
 *
 *  Written by:
 *   McKay Davis
 *   SCI Institue
 *   University of Utah
 *   September 2005
 *
 */
 
#include <Core/Util/Debug.h>

#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomStippleOccluded.h>

namespace SCIRun {

GeomStippleOccluded::GeomStippleOccluded(GeomHandle obj, int state)
: GeomSwitch(obj, state)
{
  DEBUG_CONSTRUCTOR("GeomStippleOccluded")
}

GeomStippleOccluded::GeomStippleOccluded(const GeomStippleOccluded& copy)
: GeomSwitch(copy.child_, copy.state)
{
  DEBUG_CONSTRUCTOR("GeomStippleOccluded")
}

GeomStippleOccluded::~GeomStippleOccluded()
{
  DEBUG_DESTRUCTOR("GeomStippleOccluded")
}

GeomObj* GeomStippleOccluded::clone()
{
    return new GeomStippleOccluded(*this);
}

const GLubyte stipple_pattern[] = 
{ 
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
  0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55 };


void
GeomStippleOccluded::draw(DrawInfoOpenGL* di, Material* matl, double time)
{
  if (state && child_.get_rep()) 
  {
    child_->draw(di, matl, time);
    glDepthRange(0.05, 0.08);
    glEnable(GL_POLYGON_STIPPLE);
    glPolygonStipple(stipple_pattern);
    child_->draw(di, matl, time);
    glDepthRange(0.0, 1.0);
    glDisable(GL_POLYGON_STIPPLE);
  }
}

} // End namespace SCIRun

