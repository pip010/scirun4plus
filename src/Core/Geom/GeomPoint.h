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
 * Pt.h: GeomPoint objects
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   Feb 1995
 *
 */

#ifndef SCI_Geom_Point_h
#define SCI_Geom_Point_h 1

#include <Core/Geom/GeomObj.h>
#include <Core/Datatypes/Color.h>
#include <Core/Datatypes/Material.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/BBox.h>

#include <vector>

#include <Core/Geom/share.h>

namespace SCIRun {

class SCISHARE GeomPoints : public GeomObj
{
  public:
    GeomPoints();
    GeomPoints(const GeomPoints&);
    virtual ~GeomPoints();
    virtual GeomObj* clone();

    virtual void get_bounds(BBox&);

    virtual void setPointSize( double size ) { point_size_ = size; };

    virtual inline void add(const Point& p, unsigned int idx = 0) 
    {
      points_.push_back(static_cast<float>(p.x()));
      points_.push_back(static_cast<float>(p.y()));
      points_.push_back(static_cast<float>(p.z()));
      item_idx_.push_back(idx);
    }

    virtual void add(const Point& p, const MaterialHandle &c, unsigned int idx = 0);
    virtual void add(const Point& p, double index, unsigned int idx = 0);
    virtual void add(const Point& p, const MaterialHandle &c, double index, unsigned int idx = 0);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  protected:
    double point_size_;

    std::vector<float>         points_;
    std::vector<unsigned char> colors_;
    std::vector<float>         indices_;
    std::vector<unsigned int>  item_idx_;

    bool pickable;  // hack so we don't draw non-pickable pts during a pick
};


class SCISHARE GeomTranspPoints : public GeomPoints
{
  public:
    GeomTranspPoints();
    GeomTranspPoints(const GeomTranspPoints&);
    virtual ~GeomTranspPoints();
    virtual GeomObj* clone();

    virtual void Sort();
    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  protected:
    std::vector<unsigned int> xindices_;
    std::vector<unsigned int> yindices_;
    std::vector<unsigned int> zindices_;

    bool xreverse_;
    bool yreverse_;
    bool zreverse_;
};

} // End namespace SCIRun

#endif /* SCI_Geom_Point_h */
