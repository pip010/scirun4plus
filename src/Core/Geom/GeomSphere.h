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
 * Sphere.h: Sphere objects
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#ifndef SCI_Geom_Sphere_h
#define SCI_Geom_Sphere_h 1

#include <Core/Geom/GeomObj.h>
#include <Core/Datatypes/Material.h>
#include <Core/Geometry/Point.h>

#include <vector>

#include <Core/Geom/share.h>

namespace SCIRun {

class SCISHARE GeomSphere : public GeomObj
{
  public:
    Point cen;
    double rad;
    int nu;
    int nv;
    
    void adjust();
    void move(const Point&, double, int nu=20, int nv=10);
    void move(const Point& _cen);
      
    GeomSphere(int nu=20, int nv=10);
    GeomSphere(const Point &location, double radius, int nu=20, int nv=10);
    GeomSphere(const GeomSphere &copy);
    virtual ~GeomSphere();
    
    virtual GeomObj* clone();
    virtual void get_bounds(BBox&);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

    // This is a helper function which determins the nu and nv given an
    // approximate number of polygons desired.
    static void getnunv(const int num_polygons, int &nu, int &nv);
};


class SCISHARE GeomSuperquadric : public GeomObj
{
  public:
    GeomSuperquadric(int axis, double A, double B, int nu, int nv);
    GeomSuperquadric(const GeomSuperquadric &copy);
    virtual ~GeomSuperquadric();
    
    virtual GeomObj* clone();
    virtual void get_bounds(BBox&);
    
    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  private:
    int axis_;
    double A_, B_;
    int nu_, nv_;

    std::vector<float> points_;
    std::vector<float> normals_;
    std::vector<unsigned short> tindices_;
    std::vector<unsigned short> qindices_;

    void compute_geometry();

    GeomSuperquadric();
};


class SCISHARE GeomSpheres : public GeomObj
{
  public:
    GeomSpheres(double radius = 1.0, int nu=8, int nv=8);
    GeomSpheres(const GeomSpheres &copy);
    virtual ~GeomSpheres();
    
    virtual GeomObj* clone();
    virtual void get_bounds(BBox&);
    
    void add(const Point &center, unsigned int idx = 0);
    void add(const Point &center, const MaterialHandle &mat, unsigned int idx = 0);
    void add(const Point &center, float index, unsigned int idx = 0);

    // If radius is too small, the sphere is not added and false is returned.
    bool add_radius(const Point &cen, double radius);
    bool add_radius(const Point &cen, double radius, const MaterialHandle &mat);
    bool add_radius(const Point &cen, double radius, float index);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);
    virtual void fbpick_draw(DrawInfoOpenGL*, Material*, double time);

  private:
    std::vector<Point> centers_;
    std::vector<double> radii_;
    std::vector<unsigned char> colors_;
    std::vector<float> indices_;
    int nu_;
    int nv_;
    double global_radius_;
    std::vector<unsigned int> item_idx_;
};

} // End namespace SCIRun

#endif /* SCI_Geom_Sphere_h */
