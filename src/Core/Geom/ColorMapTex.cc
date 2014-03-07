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
 *  ColorMapTex.cc: ?
 *
 *  Written by:
 *   Author: ?
 *   Department of Computer Science
 *   University of Utah
 *   Date: ?
 *
 */


#include <Core/Geom/ColorMapTex.h>
#include <Core/Geom/GeomColorMap.h>
#include <Core/Geom/GeomQuads.h>
#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/DrawInfoOpenGL.h>


namespace SCIRun {

ColorMapTex::ColorMapTex(const Point &p1, const Point &p2,
			 const Point &p3, const Point &p4,
			 ColorMapHandle cmap)
  : GeomContainer(0)
{
  GeomFastQuads *quad = new GeomFastQuads();
  const double min = cmap->getMin();
  const double max = cmap->getMax();
  quad->add (p1, min, p2, max, p3, max, p4, min);

  child_ = new GeomColorMap(quad, cmap);
}

ColorMapTex::ColorMapTex( const ColorMapTex &copy )
  : GeomContainer(copy)
{
}


ColorMapTex::~ColorMapTex()
{
}

GeomObj*
ColorMapTex::clone()
{
  return new ColorMapTex( *this );
}

void
ColorMapTex::draw(DrawInfoOpenGL *di, Material *matl, double time)
{
  if (child_.get_rep())
  {
    const bool lit = di->lighting_;
    di->lighting_ = false;
    const double diffuse = di->diffuse_scale_;
    di->diffuse_scale_ = 1.0;
    child_->draw(di, matl, time);
    di->lighting_ = lit;
    di->diffuse_scale_ = diffuse;
  }
}

} // End namespace SCIRun

