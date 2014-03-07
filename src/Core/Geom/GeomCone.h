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
 *  Cone.h: Cone object
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#ifndef SCI_Geom_Cone_h
#define SCI_Geom_Cone_h 1

#include <Core/Geom/GeomObj.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Datatypes/Material.h>

#include <vector>

#include <Core/Geom/share.h>

namespace SCIRun {

class SCISHARE GeomCone : public GeomObj
{
  protected:
    Vector v1;
    Vector v2;
    double tilt;
    double height;
    Vector zrotaxis;
    double zrotangle;
    int nu;
    int nv;
  public:
    Point bottom;
    Point top;
    Vector axis;
    double bot_rad;
    double top_rad;

    void adjust();
    void move(const Point&, const Point&, double, double, int nu=20, int nv=1);

    GeomCone(int nu=20, int nv=1);
    GeomCone(const Point&, const Point&, double, double, int nu=20, int nv=1);
    GeomCone(const GeomCone&);
    virtual ~GeomCone();

    virtual GeomObj* clone();
    virtual void get_bounds(BBox&);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);
};

class SCISHARE GeomCappedCone : public GeomCone
{
  public:
    GeomCappedCone(int nu=20, int nv=1, int nvdisk1=1, int nvdisk2=1);
    GeomCappedCone(const Point&, const Point&, double, double, 
                   int nu=20, int nv=1, int nvdisk1=1, int nvdisk2=1);
    GeomCappedCone(const GeomCappedCone&);
    virtual ~GeomCappedCone();

    virtual GeomObj* clone();

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  private:
      int nvdisk1;
      int nvdisk2;
};


class SCISHARE GeomCones : public GeomObj
{
  protected:
    double radius_;
    int  nu_;
    std::vector<Point> points_;
    std::vector<unsigned char> colors_;
    std::vector<float> indices_;
    std::vector<double> radii_;
    
  public:
    GeomCones(int nu = 8, double radius = 1.0);
    GeomCones(const GeomCones &copy);
    virtual ~GeomCones();

    virtual GeomObj* clone();
    virtual void get_bounds(BBox&);

    bool add(const Point &p0, const Point &p1);
    bool add(const Point &p0, const Point &p1, const MaterialHandle &c);
    bool add(const Point &p0, const Point &p1, float index);

    bool add_radius(const Point &p0, const Point &p1, double r);
    bool add_radius(const Point &p0, const Point &p1,
        const MaterialHandle &c, double r);
    bool add_radius(const Point &p0, const Point &p1, float index, double r);
    void set_radius(double val) { radius_ = val; reset_bbox(); }
    void set_nu(int nu);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);
};

} // End namespace SCIRun

#endif /* SCI_Geom_Cone_h */
