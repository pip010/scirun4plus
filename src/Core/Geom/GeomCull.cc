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

#include <Core/Util/Debug.h>

#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomCull.h>

namespace SCIRun {

GeomCull::GeomCull(GeomHandle child, Vector *normal) :
  GeomContainer(child), normal_(0) 
{
  DEBUG_CONSTRUCTOR("GeomCull")
  if (normal) normal_ = new Vector(*normal);
}

GeomCull::GeomCull(const GeomCull &copy) :
  GeomContainer(copy), normal_(0) 
{
  DEBUG_CONSTRUCTOR("GeomCull")
  if (copy.normal_) normal_ = new Vector(*copy.normal_);
}

GeomCull::~GeomCull()
{
  DEBUG_DESTRUCTOR("GeomCull")
}

GeomObj *
GeomCull::clone() 
{
  return new GeomCull(*this);
}
  
void
GeomCull::set_normal(Vector *normal) 
{
  if (normal_) 
  {
    delete normal_;
    normal_ = 0;
  }
  
  if (normal) 
  {
    normal_ = new Vector(*normal);
  }
}

void
GeomCull::draw(DrawInfoOpenGL* di, Material* matl, double time)
{
  if (normal_)
  {
    double mat[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, mat);
    if (Dot(Vector(mat[2], mat[6], mat[10]), *normal_) < 0) return;
  }
  child_->draw(di,matl,time);
}

}
