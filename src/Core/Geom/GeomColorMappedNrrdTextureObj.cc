//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  
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

#include <Core/Util/Debug.h>
#include <Core/Geom/GeomColorMappedNrrdTextureObj.h>
#include <Core/Geom/ColorMappedNrrdTextureObj.h>

#include <Core/Geometry/BBox.h>
#include <Core/Geom/DrawInfoOpenGL.h>

#include <sci_gl.h>

using std::ostream;

namespace SCIRun {

GeomColorMappedNrrdTextureObj::GeomColorMappedNrrdTextureObj()
{
  DEBUG_CONSTRUCTOR("GeomColorMappedNrrdTextureObj")
}

GeomColorMappedNrrdTextureObj::GeomColorMappedNrrdTextureObj
(ColorMappedNrrdTextureObjHandle &cmnto)
  : GeomObj(),
    cmnto_(cmnto),
    alpha_cutoff_(0.0),
    offset_(0.0)
{  
  DEBUG_CONSTRUCTOR("GeomColorMappedNrrdTextureObj")
}

GeomColorMappedNrrdTextureObj::GeomColorMappedNrrdTextureObj
(const GeomColorMappedNrrdTextureObj &copy ) : GeomObj(copy),
                                               cmnto_(copy.cmnto_)
{
  DEBUG_CONSTRUCTOR("GeomColorMappedNrrdTextureObj")
}

GeomColorMappedNrrdTextureObj::~GeomColorMappedNrrdTextureObj()
{
  DEBUG_DESTRUCTOR("GeomColorMappedNrrdTextureObj")
  cmnto_ = 0;
}

void
GeomColorMappedNrrdTextureObj::set_alpha_cutoff(double alpha) 
{
  alpha_cutoff_ = alpha;
}

void
GeomColorMappedNrrdTextureObj::set_offset(double offset) {
  offset_ = offset;
}

GeomObj* 
GeomColorMappedNrrdTextureObj::clone() 
{
  return new GeomColorMappedNrrdTextureObj( *this );
}

void
GeomColorMappedNrrdTextureObj::draw(DrawInfoOpenGL* di, 
                                    Material* matl, double)
{
  const double old_ambient = di->ambient_scale_;
  di->ambient_scale_ = 150;

  if (!pre_draw(di, matl, 1)) {
    di->ambient_scale_ = old_ambient;
    return;
  }

  glEnable(GL_ALPHA_TEST);  
  glAlphaFunc(GL_GREATER, 0.1);
  glEnable(GL_POLYGON_OFFSET_FILL);

  // Very scientific, don't change :-o
  const Vector view(di->view_.lookat() - di->view_.eyep());
  const double d = 1.0 / view.length() - 1.0;
  glPolygonOffset(d * offset_, offset_);

  // Draw the quad with no alpha.
  cmnto_->draw_quad(true);

  glPolygonOffset(0.0, 0.0);
  glDisable(GL_POLYGON_OFFSET_FILL);
  glDisable(GL_ALPHA_TEST);  

  di->ambient_scale_ = old_ambient;

  post_draw(di);
}

void GeomColorMappedNrrdTextureObj::get_bounds( BBox& bb ) 
{
  cmnto_->get_bounds(bb);
}


} // End namespace SCIRun

