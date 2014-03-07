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
 *  BBoxCache.cc: ?
 *
 *  Written by:
 *   Author?
 *   Department of Computer Science
 *   University of Utah
 *   Date?
 *
 */

#include <Core/Geom/BBoxCache.h>

namespace SCIRun {


GeomBBoxCache::GeomBBoxCache(GeomHandle obj)
 : GeomContainer(obj), bbox_cached(0)
{

}

GeomBBoxCache::GeomBBoxCache(GeomHandle obj, const BBox &box)
  : GeomContainer(obj), bbox_cached(true)
{
  bbox.extend( box );
}

GeomObj* GeomBBoxCache::clone()
{
  return 0;
}

void GeomBBoxCache::reset_bbox()
{
  bbox_cached = false;
  GeomContainer::reset_bbox();
}

void GeomBBoxCache::get_bounds(BBox& box)
{
  if (!bbox_cached || !bbox.valid()) 
  {
    bbox.reset();
    child_->get_bounds(bbox);
    bbox_cached = true;
  }
  
  box.extend( bbox );
}

} // End namespace SCIRun
