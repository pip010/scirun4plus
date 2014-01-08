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
 * GeomTorus.h: Torus objects
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   January 1995
 *
 */

#ifndef SCI_Geom_Torus_h
#define SCI_Geom_Torus_h 1

#include <Core/Geom/GeomTriangles.h>

#include <Core/Datatypes/Material.h>
#include <Core/Datatypes/Color.h>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>

#include <Core/Geom/share.h>


namespace SCIRun {

class SCISHARE GeomTorus : public GeomObj
{
  public:
    GeomTorus(int nu=50, int nv=8);
    GeomTorus(const Point&, const Vector&, double, double,
              int nu=50, int nv=8);
    GeomTorus(const GeomTorus&);
    virtual ~GeomTorus();

    virtual GeomObj* clone();

    virtual void adjust();
    virtual void move(const Point&, const Vector&, double, double,
          int nu=50, int nv=8);

    virtual void get_bounds(BBox&);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  protected:
    Point center_;
    Vector axis_;
    double major_radius_;
    double minor_radius_;
    int nu_;
    int nv_;

    Vector zrotaxis_;
    double zrotangle_;
};


class SCISHARE GeomTorusArc : public GeomTorus
{
  public:
    GeomTorusArc(int nu=50, int nv=8);
    GeomTorusArc(const Point&, const Vector&, double, double, 
                 const Vector& zero, double start_angle, double arc_angle,
                 int nu=50, int nv=8);
    GeomTorusArc(const GeomTorusArc&);
    virtual ~GeomTorusArc();

    virtual GeomObj* clone();

    virtual void adjust();
    virtual void move(const Point&, const Vector&, double, double,
          const Vector& zero, double start_angle, double arc_angle,
          int nu=50, int nv=8);

    virtual void get_bounds(BBox&);
    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  protected:
    Vector zero_;
    double start_angle_;
    double arc_angle_;
    Vector yaxis_;
};

} // End namespace SCIRun

#endif /* SCI_Geom_Torus_h */
