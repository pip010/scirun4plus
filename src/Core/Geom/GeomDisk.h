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
 *  Disc.h:  Disc object
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#ifndef SCI_Geom_Disk_h
#define SCI_Geom_Disk_h 1

#include <Core/Geom/GeomObj.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>

namespace SCIRun {

class GeomDisk : public GeomObj
{
  public:
    Point cen;
    Vector n;
    double rad;
    int nu;
    int nv;

    void adjust();
    void move(const Point&, const Vector&, double, int nu=20, int nv=2);

    GeomDisk(int nu=20, int nv=2);
    GeomDisk(const Point&, const Vector&, double, int nu=20, int nv=2);
    GeomDisk(const GeomDisk&);
    virtual ~GeomDisk();

    virtual GeomObj* clone();
    virtual void get_bounds(BBox&);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  private:
    Vector v1;
    Vector v2;
    Vector zrotaxis;
    double zrotangle;
};

} // End namespace SCIRun

#endif /* SCI_Geom_Disk_h */
