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
 *  GeomTriangles.h: Triangle Strip object
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   June 1994
 *
 */

#ifndef CORE_GEOM_GEOMTRIANGLES_H
#define CORE_GEOM_GEOMTRIANGLES_H 1

#include <Core/Geom/GeomVertexPrim.h>

#include <vector>

#include <Core/Geom/share.h>

// SCIRun uses SLIVR's Color class
class SLIVR::Color;

namespace SCIRun {

class SCISHARE GeomTriangles : public GeomVertexPrim
{
  public:
    Array1<Vector> normals;
    GeomTriangles();
    GeomTriangles(const GeomTriangles&);
    virtual ~GeomTriangles();
    
    int size(void);
    void add(const Point&, const Point&, const Point&);
    void add(const Point&, const Vector&,
             const Point&, const Vector&,
             const Point&, const Vector&);
    void add(const Point&, const MaterialHandle&,
             const Point&, const MaterialHandle&,
             const Point&, const MaterialHandle&);
    void add(const Point&, const Color&,
             const Point&, const Color&,
             const Point&, const Color&);
    void add(const Point&, const Vector&, const MaterialHandle&,
             const Point&, const Vector&, const MaterialHandle&,
             const Point&, const Vector&, const MaterialHandle&);
    void add(GeomVertex*, GeomVertex*, GeomVertex*);
    virtual GeomObj* clone();

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  private:
    void add(const Point&);
    void add(const Point&, const Vector&);
    void add(const Point&, const MaterialHandle&);
    void add(const Point&, const Vector&, const MaterialHandle&);
    void add(GeomVertex*);
};


class SCISHARE GeomFastTriangles : public GeomObj
{
  public:
    GeomFastTriangles();
    GeomFastTriangles(const GeomFastTriangles&);
    virtual ~GeomFastTriangles();
    virtual GeomObj* clone();

    int size(void);
    void add(const Point &p0, const Point &p1, const Point &p2);
    void add(const Point &p0, const MaterialHandle &m0,
       const Point &p1, const MaterialHandle &m1,
       const Point &p2, const MaterialHandle &m2);
    void add(const Point &p0, double cindex0,
       const Point &p1, double cindex1,
       const Point &p2, double cindex2);
    void add(const Point &p0, const MaterialHandle &m0, double cindex0,
       const Point &p1, const MaterialHandle &m1, double cindex1,
       const Point &p2, const MaterialHandle &m2, double cindex2);

    void add(const Point &p0, const Vector &n0,
       const Point &p1, const Vector &n1,
       const Point &p2, const Vector &n2);
    void add(const Point &p0, const Vector &n0, const MaterialHandle &m0,
       const Point &p1, const Vector &n1, const MaterialHandle &m1,
       const Point &p2, const Vector &n2, const MaterialHandle &m2);
    void add(const Point &p0, const Vector &n0, double cindex0,
       const Point &p1, const Vector &n1, double cindex1,
       const Point &p2, const Vector &n2, double cindex2);
    void add(const Point &p0, const Vector &n0,
       const MaterialHandle &m0, double cindex0,
       const Point &p1, const Vector &n1,
       const MaterialHandle &m1, double cindex1,
       const Point &p2, const Vector &n2,
       const MaterialHandle &m2, double cindex2);

    void add( const std::vector< std::vector< std::pair<Point, Vector> > >& tristrips );
    void add( const std::vector< std::vector< std::pair<Point, Vector> > >& tristrips,
        const MaterialHandle &mat );
    void add( const std::vector< std::vector< std::pair<Point, Vector> > >& tristrips,
        double cindex );
    void add( const std::vector< std::vector< std::pair<Point, Vector> > >& tristrips,
        const MaterialHandle &mat, double cindex );

    virtual void get_bounds(BBox& bb);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);
    virtual void drawVertexData(DrawInfoOpenGL*);
    virtual void fbpick_draw(DrawInfoOpenGL*, Material*, double time);

  protected:
    std::vector<float>         points_;
    std::vector<unsigned char> colors_;
    std::vector<float>         indices_;
    std::vector<float>         normals_;
    std::vector<float>         face_normals_;
    MaterialHandle             material_;
};


class SCISHARE GeomTranspTriangles : public GeomFastTriangles
{
  public:
    GeomTranspTriangles();
    GeomTranspTriangles(const GeomTranspTriangles&);
    virtual ~GeomTranspTriangles();
    virtual GeomObj* clone();

    void Sort();

    virtual void draw(DrawInfoOpenGL*, Material*, double time);
    virtual void fbpick_draw(DrawInfoOpenGL*, Material*, double time);

  protected:
    std::vector<unsigned int> xlist_;
    std::vector<unsigned int> ylist_;
    std::vector<unsigned int> zlist_;

    bool xreverse_;
    bool yreverse_;
    bool zreverse_;
};


class SCISHARE GeomFastTrianglesTwoSided : public GeomObj
{
  public:
    GeomFastTrianglesTwoSided();
    GeomFastTrianglesTwoSided(const GeomFastTrianglesTwoSided&);
    virtual ~GeomFastTrianglesTwoSided();
    virtual GeomObj* clone();

    int size(void);
    void add(const Point &p0, const Point &p1, const Point &p2);
    void add(const Point &p0, const MaterialHandle &m0, const MaterialHandle &k0,
       const Point &p1, const MaterialHandle &m1, const MaterialHandle &k1,
       const Point &p2, const MaterialHandle &m2, const MaterialHandle &k2);
    void add(const Point &p0, double cindex0, double dindex0,
       const Point &p1, double cindex1, double dindex1,
       const Point &p2, double cindex2, double dindex2);
    void add(const Point &p0,
       const MaterialHandle &m0, double cindex0,
       const MaterialHandle &k0, double dindex0,
       const Point &p1,
       const MaterialHandle &m1, double cindex1,
       const MaterialHandle &k1, double dindex1,
       const Point &p2,
       const MaterialHandle &m2, double cindex2,
       const MaterialHandle &k2, double dindex2);

    void add(const Point &p0, const Vector &n0,
       const Point &p1, const Vector &n1,
       const Point &p2, const Vector &n2);
    void add(const Point &p0, const Vector &n0, const MaterialHandle &m0, const MaterialHandle &k0,
       const Point &p1, const Vector &n1, const MaterialHandle &m1, const MaterialHandle &k1,
       const Point &p2, const Vector &n2, const MaterialHandle &m2, const MaterialHandle &k2);     
    void add(const Point &p0, const Vector &n0, double i0, double j0,
       const Point &p1, const Vector &n1, double i1, double j1,
       const Point &p2, const Vector &n2, double i2, double j2);
    void add(const Point &p0, const Vector &n0,
       const MaterialHandle &m0, double cindex0,
       const MaterialHandle &k0, double dindex0,
       const Point &p1, const Vector &n1,
       const MaterialHandle &m1, double cindex1,
       const MaterialHandle &k1, double dindex1,
       const Point &p2, const Vector &n2,
       const MaterialHandle &m2, double cindex2,
       const MaterialHandle &k2, double dindex2);

    virtual void get_bounds(BBox& bb);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);
    virtual void drawVertexData(DrawInfoOpenGL*);

  protected:
    std::vector<float> points_;
    std::vector<unsigned char> colors_;
    std::vector<unsigned char> colors2_;  
    std::vector<float> indices_;
    std::vector<float> indices2_;
    std::vector<float> normals_;
    std::vector<float> face_normals_;
    MaterialHandle material_;
    MaterialHandle material2_;
};


class SCISHARE GeomTranspTrianglesTwoSided : public GeomFastTrianglesTwoSided
{
  public:
    GeomTranspTrianglesTwoSided();
    GeomTranspTrianglesTwoSided(const GeomTranspTrianglesTwoSided&);
    virtual ~GeomTranspTrianglesTwoSided();
    virtual GeomObj* clone();

    void Sort();

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  protected:
    std::vector<unsigned int> xlist_;
    std::vector<unsigned int> ylist_;
    std::vector<unsigned int> zlist_;

    bool xreverse_;
    bool yreverse_;
    bool zreverse_;
};


class SCISHARE GeomTrianglesP: public GeomObj
{
  public:
    GeomTrianglesP();
    virtual ~GeomTrianglesP();

    void SetColor(double rr, double gg, double bb) 
    {
      has_color_ = 1;
      r_ = rr; g_ = gg; b_ = bb;
    };

    int size(void);

    int add(const Point&, const Point&, const Point&);
      
    // below is a virtual function - makes some things easier...
    virtual int vadd(const Point&, const Point&, const Point&);

    void reserve_clear(int);   // reserves storage... and clears

    virtual GeomObj* clone();

    virtual void get_bounds(BBox&);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

    virtual void get_triangles( Array1<float> &v);
    
  protected:
    Array1<float> points_;
    Array1<float> normals_;

    int has_color_;
    double r_, g_, b_;  // actual color values...
};


class SCISHARE GeomTrianglesPT1d : public GeomTrianglesP
{
  public:
    GeomTrianglesPT1d();
    virtual ~GeomTrianglesPT1d();

    int add(const Point&, const Point&, const Point&,
      const float&, const float&, const float&);
    
    unsigned char* cmap;

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  protected:
    Array1<float> scalars_;
};


class SCISHARE GeomTranspTrianglesP : public GeomTrianglesP
{
  public:
    GeomTranspTrianglesP();
    GeomTranspTrianglesP(double);
    virtual ~GeomTranspTrianglesP();

    // below computes the x/y/z center points as well...
    // should only use below...
    virtual int vadd(const Point&, const Point&, const Point&);

    // function below sorts the polygons...
    void Sort();

    // function below merges in another "list" - also clears
    //void MergeStuff(GeomTranspTrianglesP*);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  protected:
    double alpha_;

    std::vector<std::pair<float, unsigned int> > xlist_;
    std::vector<std::pair<float, unsigned int> > ylist_;
    std::vector<std::pair<float, unsigned int> > zlist_;

    bool sorted_;
};


class SCISHARE GeomTrianglesPC: public GeomTrianglesP
{
  public:
    GeomTrianglesPC();
    virtual ~GeomTrianglesPC();

    int add(const Point&, const Color&,
            const Point&, const Color&,
            const Point&, const Color&);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  private:
    Array1<float> colors_;
};


class SCISHARE GeomTrianglesVP: public GeomObj
{
  public:
    GeomTrianglesVP();
    virtual ~GeomTrianglesVP();

    int size(void);

    int add(const Point&, const Vector&, 
            const Point&, const Vector&, 
            const Point&, const Vector&);
      
    void reserve_clear(int);   // reserves storage... and clears

    virtual GeomObj* clone();

    virtual void get_bounds(BBox&);

    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  protected:
    Array1<float> points_;
    Array1<float> normals_;
};


class SCISHARE GeomTrianglesVPC: public GeomTrianglesVP
{
  public:
    GeomTrianglesVPC();
    virtual ~GeomTrianglesVPC();

    int add(const Point&, const Vector&, const Color&,
            const Point&, const Vector&, const Color&,
            const Point&, const Vector&, const Color&);
    
    virtual void draw(DrawInfoOpenGL*, Material*, double time);

  private:
    Array1<float> colors_;
};

} // End namespace SCIRun

#endif /* SCI_Geom_Triangles_h */
