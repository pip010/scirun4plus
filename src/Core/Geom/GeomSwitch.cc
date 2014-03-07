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
 *  Switch.cc:  Turn Geometry on and off
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   January 1995
 *
 */
 
#include <Core/Util/Debug.h>
#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomSwitch.h>

namespace SCIRun {

GeomSwitch::GeomSwitch(GeomHandle obj, int state)
: GeomContainer(obj), state(state)
{
  DEBUG_CONSTRUCTOR("GeomSwitch")
}

GeomSwitch::GeomSwitch(const GeomSwitch& copy)
: GeomContainer(copy), state(copy.state)
{
  DEBUG_CONSTRUCTOR("GeomSwitch")
}

GeomSwitch::~GeomSwitch()
{
  DEBUG_DESTRUCTOR("GeomSwitch")
}

GeomObj* GeomSwitch::clone()
{
  return new GeomSwitch(*this);
}

void GeomSwitch::set_state(int st)
{
  state=st;
}

int GeomSwitch::get_state()
{
  return state;
}

void GeomSwitch::get_bounds(BBox& bbox)
{
  if (state && child_.get_rep()) child_->get_bounds(bbox);
}

void
GeomSwitch::fbpick_draw(DrawInfoOpenGL* di, Material* matl, double time)
{
  if (state && child_.get_rep())
  {
    child_->fbpick_draw(di, matl, time);
  }
}

void
GeomSwitch::draw(DrawInfoOpenGL* di, Material* matl, double time)
{
  if (state && child_.get_rep())
  {
    child_->draw(di, matl, time);
  }
}

} // End namespace SCIRun

