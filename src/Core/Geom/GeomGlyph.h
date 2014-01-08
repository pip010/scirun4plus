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
 * GeomGlyph.h: Glyph objects
 *
 *  Written by:
 *   Allen R. Sanderson
 *   Scientific Computing & Imaging Inst
 *   University of Utah
 *   March 2007
 *
 */

#ifndef SCI_Geom_Glyph_h
#define SCI_Geom_Glyph_h 1

#include <Core/Datatypes/Color.h>
#include <Core/Datatypes/Material.h>
#include <Core/Geom/GeomTriangles.h>
#include <Core/Geom/GeomQuads.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Transform.h>

#include <Core/Math/MiscMath.h>
#include <Core/Math/TrigTable.h>

#include <iostream>
#include <vector>

#include <Core/Geom/share.h>

namespace SCIRun {

class SCISHARE GeomGlyphBase
{
public:
  typedef std::vector< std::pair< Point, Vector > > QuadStrip;

  GeomGlyphBase() 
  {
    tables_.resize(21);
    for (size_t j=0;j<21;j++)
    {
      tables_[j].build_table(j,0.0,2*M_PI);
    }
  }
  
  virtual ~GeomGlyphBase() {}

  virtual GeomObj *getObj() { return 0; };

  // Arrow
  template <class T>
  void add_arrow(const Point&, const T& t,
		 double radius, double length,
		 int nu=20, int nv=2);
  
  template <class T>
  void add_arrow(const Point&, const T& t,
		 double radius, double length,
		 const MaterialHandle &mat,
		 int nu=20, int nv=2);
  
  template <class T>
  void add_arrow(const Point&, const T& t,
		 double radius, double length,
		 const double cindex,
		 int nu=20, int nv=2);

  template <class T>
  void add_arrow(const Point&, const T& t,
		 double radius, double length,
		 const MaterialHandle &mat,
		 const double cindex,
		 int nu=20, int nv=2);

  // Box
  template <class T>
  void add_box(const Point& center, const T& t,
	       double x_side, double y_side, double z_side,
	       const bool RGB = false);

  template <class T>
  void add_box(const Point& center, const T& t,
	       double x_side, double y_side, double z_side,
	       const MaterialHandle &mat);

  template <class T>
  void add_box(const Point& center, const T& t,
	       double x_side, double y_side, double z_side,
	       const double cindex);

  template <class T>
  void add_box(const Point& center, const T& t,
	       double x_side, double y_side, double z_side,
	       const MaterialHandle &mat, const double cindex);

  // Capped Cylinder
  template <class T>
  void add_capped_cylinder(const Point&, const T& t,
			   double radius1, double radius2, double length,
			   int nu=20, int nv=2);
  
  template <class T>
  void add_capped_cylinder(const Point&, const T& t,
			   double radius1, double radius2, double length,
			   const MaterialHandle &mat,
			   int nu=20, int nv=2);
  
  template <class T>
  void add_capped_cylinder(const Point&, const T& t,
			   double radius1, double radius2, double length,
			   const double cindex,
			   int nu=20, int nv=2);

  template <class T>
  void add_capped_cylinder(const Point&, const T& t,
			   double radius1, double radius2, double length,
			   const MaterialHandle &mat, const double cindex,
			   int nu=20, int nv=2);

  // Cylinder
  template <class T>
  void add_cylinder(const Point&, const T& t,
		    double radius1, double radius2, double length,
		    int nu=20, int nv=2);
  
  template <class T>
  void add_cylinder(const Point&, const T& t,
		    double radius1, double radius2, double length,
		    const MaterialHandle &mat,
		    int nu=20, int nv=2);

  template <class T>
  void add_cylinder(const Point&, const T& t,
		    double radius1, double radius2, double length,
		    const MaterialHandle &mat0, const MaterialHandle &mat1,
		    int nu=20, int nv=25);
  
  template <class T>
  void add_cylinder(const Point&, const T& t,
		    double radius1, double radius2, double length,
		    const double cindex,
		    int nu=20, int nv=2);

  template <class T>
  void add_cylinder(const Point&, const T& t,
		    double radius1, double radius2, double length,
		    const double cindex0, const double cindex1,
		    int nu=20, int nv=25);

  template <class T>
  void add_cylinder(const Point&, const T& t,
		    double radius1, double radius2, double length,
		    const MaterialHandle &mat, const double cindex,
		    int nu=20, int nv=2);

  template <class T>
  void add_cylinder(const Point&, const T& t,
		    double radius1, double radius2, double length,
		    const MaterialHandle &mat0, const double cindex0,
		    const MaterialHandle &mat1, const double cindex1,
		    int nu=20, int nv=25);

  // Ellipsoid
  template <class T, class S>
  void add_ellipsoid(const Point&, const T& t, S scales,
		     int nu=20, int nv=20, int half=0);
  
  template <class T, class S>
  void add_ellipsoid(const Point&, const T& t, S scales,
		     const MaterialHandle &mat,
		     int nu=20, int nv=20, int half=0);

  template <class T, class S>
  void add_ellipsoid(const Point&, const T& t, S scales,
		     const double cindex,
		     int nu=20, int nv=20, int half=0);

  template <class T, class S>
  void add_ellipsoid(const Point&, const T& t, S scales,
		     const MaterialHandle &mat,
		     const double cindex,
		     int nu=20, int nv=20, int half=0);

  // Helix
  template <class T>
  void add_helix(const Point&, const T& t,
		 double major_radius1, double major_radius2,
		 double minor_radius, double length,
		 unsigned int pitch,
		 int nu=20, int nv=8);
  
  template <class T>
  void add_helix(const Point&, const T& t,
		 double major_radius1, double major_radius2,
		 double minor_radius, double length,
		 unsigned int pitch,
		 const MaterialHandle &mat,
		 int nu=20, int nv=8);

  template <class T>
  void add_helix(const Point&, const T& t,
		 double major_radius1, double major_radius2,
		 double minor_radius, double length,
		 unsigned int pitch,
		 const MaterialHandle &mat0, const MaterialHandle &mat1,
		 int nu=20, int nv=8);
  
  template <class T>
  void add_helix(const Point&, const T& t,
		 double major_radius1, double major_radius2,
		 double minor_radius, double length,
		 unsigned int pitch,
		 const double cindex,
		 int nu=20, int nv=8);

  template <class T>
  void add_helix(const Point&, const T& t,
		 double major_radius1, double major_radius2,
		 double minor_radius, double length,
		 unsigned int pitch,
		 const double cindex0, const double cindex1,
		 int nu=20, int nv=8);

  template <class T>
  void add_helix(const Point&, const T& t,
		 double major_radius1, double major_radius2,
		 double minor_radius, double length,
		 unsigned int pitch,
		 const MaterialHandle &mat, const double cindex,
		 int nu=20, int nv=8);

  template <class T>
  void add_helix(const Point&, const T& t,
		 double major_radius1, double major_radius2,
		 double minor_radius, double length,
		 unsigned int pitch,
		 const MaterialHandle &mat0, const double cindex0,
		 const MaterialHandle &mat1, const double cindex1,
		 int nu=20, int nv=8);

  // Sphere
  void add_sphere(const Point&, double radius,
		 int nu=20, int nv=20, int half=0);
  
  void add_sphere(const Point&, double radius,
		 const MaterialHandle &mat,
		 int nu=20, int nv=20, int half=0);

  void add_sphere(const Point&, double radius,
		 const double cindex,
		 int nu=20, int nv=20, int half=0);

  void add_sphere(const Point&, double radius,
		 const MaterialHandle &mat,
		 const double cindex,
		 int nu=20, int nv=20, int half=0);

  // Superquadric
  template <class T>
  void add_superquadric(const Point&, const T& t,
			Vector scales, double A, double B, int axis, 
			int nu=20, int nv=8);
  
  template <class T>
  void add_superquadric(const Point&, const T& t,
			Vector scales, double A, double B, int axis, 
			const MaterialHandle &mat,
			int nu=20, int nv=8);

  template <class T>
  void add_superquadric(const Point&, const T& t,
			Vector scales, double A, double B, int axis, 
			const double cindex,
			int nu=20, int nv=8);
  
  template <class T>
  void add_superquadric(const Point&, const T& t,
			Vector scales, double A, double B, int axis, 
			const MaterialHandle &mat,
			const double cindex,
			int nu=20, int nv=8);

  // Torus
  template <class T>
  void add_torus(const Point&, const T& t,
		 double major_radius, double minor_radius,
		 int nu=20, int nv=8);
  
  template <class T>
  void add_torus(const Point&, const T& t,
		 double major_radius, double minor_radius,
		 const MaterialHandle &mat,
		 int nu=20, int nv=8);

  template <class T>
  void add_torus(const Point&, const T& t,
		 double major_radius, double minor_radius,
		 const double cindex,
		 int nu=20, int nv=8);

  template <class T>
  void add_torus(const Point&, const T& t,
		 double major_radius, double minor_radius,
		 const MaterialHandle &mat,
		 const double cindex,
		 int nu=20, int nv=8);

public:
  static void gen_transforms(const Point& center,
			     const Vector& normal,
			     Transform& trans,
			     Transform& rotate);


  static void gen_transforms(const Point& center,
			     const Transform& transform,
			     Transform& trans,
			     Transform& rotate);


protected:

  inline double spow(double e, double x)
  {
    // This for round off of very small numbers.
    if( Abs( e ) < 1.0e-6)
      e = 0.0;

    if (e < 0.0)
    {
      return (double)(Pow(Abs(e), x) * -1.0);
    }
    else
    {
      return (double)(Pow(e, x));
    }
  }

  virtual void add_strip(const std::vector<QuadStrip>& ) {};
  
  
  virtual void add_strip(const std::vector<QuadStrip>&, const MaterialHandle& ) {};


  virtual void add_strip(const std::vector<QuadStrip>&, double) {};

  virtual void add_strip(const std::vector<QuadStrip>&,
                         const MaterialHandle &, double) {};


  template <class T>
  void gen_box(const Point& center, const T& t, 
	       double x_side, double y_side, double z_side,
	       std::vector< QuadStrip >& quadstrips);

  template <class T>
  void gen_capped_cylinder(const Point& center, const T& t,
			   double radius1, double radius2, double length,
			   int nu, int nv,
			   std::vector< QuadStrip >& quadstrips);

  template <class T>
  void gen_cylinder(const Point& center, const T& t,
		    double radius1, double radius2, double length,
		    int nu, int nv,
		    std::vector< QuadStrip >& quadstrips);

  template <class T>
  void gen_disk(const Point& center, const T& t,
		double radius, double thickness,
		int nu, int nv, bool reverse,
		std::vector< QuadStrip >& quadstrips);

  template <class T, class S>
  void gen_ellipsoid(const Point& center, const T& t, S scales,
		     int nu, int nv, int half,
		     std::vector< QuadStrip >& quadstrips );

  template <class T>
  void gen_helix(const Point& center, const T& t,
		 double major_radius1, double major_radius2,
		 double minor_radius, double length,
		 unsigned int pitch,
		 int nu, int nv,
		 std::vector< QuadStrip >& quadstrips );

  template <class T>
  void gen_superquadric(const Point& center, const T& t,
			Vector scales, double A, double B, int axis,
			int nu, int nv,
			std::vector< QuadStrip >& quadstrips );
  
  template <class T>
  void gen_torus(const Point& center, const T& t,
		 double major_radius, double minor_radius,
		 int nu, int nv,
		 std::vector< QuadStrip >& quadstrips);
  
  protected:
    // So we only compute it once instead for each glyph:
    std::vector<SinCosTable> tables_;
};


class SCISHARE GeomGlyph : public GeomGlyphBase, public GeomFastQuads
{
  public:
    GeomGlyph();
    GeomGlyph(const GeomGlyph&);
    virtual ~GeomGlyph();

    virtual GeomObj* clone();

    virtual GeomObj *getObj() { return (GeomObj*) this; };

  protected:
    virtual void add_strip( const std::vector< QuadStrip >&
          quadstrips ) 
    {
      add( quadstrips );
    }
    
    virtual void add_strip( const std::vector< QuadStrip >&
          quadstrips, const MaterialHandle &mat ) 
    {
      add( quadstrips, mat );
    }

    virtual void add_strip( const std::vector< QuadStrip >&
          quadstrips, double cindex ) 
    {
      add( quadstrips, cindex );
    }

    virtual void add_strip( const std::vector< QuadStrip >&
          quadstrips,
          const MaterialHandle &mat, double cindex ) 
    {
      add( quadstrips, mat, cindex );
    }
};


class SCISHARE GeomTranspGlyph : public GeomGlyphBase, public GeomTranspQuads
{
  public:
    GeomTranspGlyph();
    GeomTranspGlyph(const GeomTranspGlyph&);
    virtual ~GeomTranspGlyph();

    virtual GeomObj* clone();

    virtual GeomObj *getObj() { return (GeomObj*) this; };

  protected:
    virtual void add_strip( const std::vector< QuadStrip >&
          quadstrips ) 
    {
      add( quadstrips );
    }
    
    
    virtual void add_strip( const std::vector< QuadStrip >&
          quadstrips, const MaterialHandle &mat ) 
    {
      add( quadstrips, mat );
    }
    
    virtual void add_strip( const std::vector< QuadStrip >&
          quadstrips, double cindex ) 
    {
      add( quadstrips, cindex );
    }

    virtual void add_strip( const std::vector< QuadStrip >&
          quadstrips,
          const MaterialHandle &mat, double cindex ) 
    {
      add( quadstrips, mat, cindex );
    }
};


template <class T>
void GeomGlyphBase::gen_box(const Point& center, const T& t, 
			    double x_side, double y_side, double z_side,
			    std::vector< QuadStrip >& quadstrips )
{
  double half_x_side = x_side * 0.5;
  double half_y_side = y_side * 0.5;
  double half_z_side = z_side * 0.5;

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  //Draw the Box
  Point p1 = trans * Point(-half_x_side,  half_y_side,  half_z_side);
  Point p2 = trans * Point(-half_x_side,  half_y_side, -half_z_side);
  Point p3 = trans * Point( half_x_side,  half_y_side,  half_z_side);
  Point p4 = trans * Point( half_x_side,  half_y_side, -half_z_side);

  Point p5 = trans * Point(-half_x_side, -half_y_side,  half_z_side);
  Point p6 = trans * Point(-half_x_side, -half_y_side, -half_z_side);
  Point p7 = trans * Point( half_x_side, -half_y_side,  half_z_side);
  Point p8 = trans * Point( half_x_side, -half_y_side, -half_z_side);

  Vector v1 = rotate * Vector(half_x_side, 0,            0);
  Vector v2 = rotate * Vector(0,           half_y_side,  0);
  Vector v3 = rotate * Vector(0,           0,            half_z_side);

  Vector v4 = rotate * Vector(-half_x_side, 0,           0);
  Vector v5 = rotate * Vector(0,           -half_y_side, 0);
  Vector v6 = rotate * Vector(0,           0,            -half_z_side);

  QuadStrip quadstrip1;
  QuadStrip quadstrip2;
  QuadStrip quadstrip3;
  QuadStrip quadstrip4;
  QuadStrip quadstrip5;
  QuadStrip quadstrip6;

  // +X
  quadstrip1.push_back( std::make_pair(p7, v1) );
  quadstrip1.push_back( std::make_pair(p8, v1) );
  quadstrip1.push_back( std::make_pair(p3, v1) );
  quadstrip1.push_back( std::make_pair(p4, v1) );

  // +Y
  quadstrip2.push_back( std::make_pair(p3, v2) );
  quadstrip2.push_back( std::make_pair(p4, v2) );
  quadstrip2.push_back( std::make_pair(p1, v2) );
  quadstrip2.push_back( std::make_pair(p2, v2) );

  // +Z
  quadstrip3.push_back( std::make_pair(p5, v3) );
  quadstrip3.push_back( std::make_pair(p7, v3) );
  quadstrip3.push_back( std::make_pair(p1, v3) );
  quadstrip3.push_back( std::make_pair(p3, v3) );

  // -X
  quadstrip4.push_back( std::make_pair(p1, v4) );
  quadstrip4.push_back( std::make_pair(p2, v4) );
  quadstrip4.push_back( std::make_pair(p5, v4) );
  quadstrip4.push_back( std::make_pair(p6, v4) );

  // -Y
  quadstrip5.push_back( std::make_pair(p5, v5) );
  quadstrip5.push_back( std::make_pair(p6, v5) );
  quadstrip5.push_back( std::make_pair(p7, v5) );
  quadstrip5.push_back( std::make_pair(p8, v5) );

  // -Z
  quadstrip6.push_back( std::make_pair(p2, v6) );
  quadstrip6.push_back( std::make_pair(p4, v6) );
  quadstrip6.push_back( std::make_pair(p6, v6) );
  quadstrip6.push_back( std::make_pair(p8, v6) );
  
  quadstrips.push_back( quadstrip1 );
  quadstrips.push_back( quadstrip2 );
  quadstrips.push_back( quadstrip3 );
  quadstrips.push_back( quadstrip4 );
  quadstrips.push_back( quadstrip5 );
  quadstrips.push_back( quadstrip6 );
}



template <class T>
void GeomGlyphBase::gen_capped_cylinder(const Point& center, const T& t,
					double radius1, double radius2, double length,
					int nu, int nv,
					std::vector< QuadStrip >& quadstrips )
{
//  SinCosTable tab1(nu,  0, 2*Pi);

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  // Draw the cylinder
  gen_cylinder(center, t, radius1, radius2, length, nu, nv, quadstrips );

  // Draw the first cap
  if( radius1 > 0 )
    gen_disk(center, t, 0.0, radius1, nu, nv, true, quadstrips);

  // Draw the second cap
  Vector offset = rotate * Vector(0,0,1);
  offset.safe_normalize();
  offset *= length;

  if( radius2 > 0 )
    gen_disk(center+offset, t, 0.0, radius2, nu, nv, false, quadstrips);
}


template <class T>
void GeomGlyphBase::gen_cylinder(const Point& center, const T& t,
				 double radius1, double radius2, double length,
				 int nu, int nv,
				 std::vector< QuadStrip >& quadstrips )
{
  nu++; //Bring nu to expected value for shape.

  if ( nu > 20) nu = 20;
  SinCosTable& tab1 = tables_[nu];

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  // Draw the cylinder
  double dz = length / (float) nv;
  double dr = (radius2-radius1) / (float) nv;

  for (int v=0; v<nv; v++)
  {
    double z1 = dz * (float) v; 
    double z2 = z1 + dz;

    double r1 = radius1 + dr * (float) v; 
    double r2 = r1 + dr;

    QuadStrip quadstrip;

    for (int u=0; u<nu; u++)
    {
      double nx = tab1.sin(u);
      double ny = tab1.cos(u);

      double x1 = r1 * nx;
      double y1 = r1 * ny;

      double x2 = r2 * nx;
      double y2 = r2 * ny;

      double nx1 = length * nx;
      double ny1 = length * ny;
      
      Point p1 = trans * Point(x1, y1, z1);
      Point p2 = trans * Point(x2, y2, z2);

      Vector v1 = rotate * Vector(nx1, ny1, -dr);
      v1.safe_normalize();

      quadstrip.push_back( std::make_pair(p1, v1) );
      quadstrip.push_back( std::make_pair(p2, v1) );
    }

    quadstrips.push_back( quadstrip );
  }
}


template <class T>
void GeomGlyphBase::gen_disk(const Point& center, const T& t,
			     double radius1, double radius2,
			     int nu, int nv, bool reverse,
			     std::vector< QuadStrip >& quadstrips)
{
  nu++; //Bring nu to expected value for shape.

  if ( nu > 20) nu = 20;
  SinCosTable& tab1 = tables_[nu];

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

 // Draw the disk
  Vector n1 = trans * Vector(0, 0, 1) * (reverse ? -1.0 : 1.0);

  double z1 = 0;

  double dr = (radius2-radius1) / (float) nv;

  for (int v=0; v<nv; v++)
  {
    double r1 = radius1 + dr * (float) (v+1); 
    double r2 = radius1 + dr * (float) v;

    QuadStrip quadstrip1;

    for (int u=0; u<nu; u++)
    {
      double nx = tab1.sin(u);
      double ny = tab1.cos(u);

      double x1 = r1 * nx;
      double y1 = r1 * ny;

      double x2 = r2 * nx;
      double y2 = r2 * ny;

      Point p1 = trans * Point(x1, y1, z1);
      Point p2 = trans * Point(x2, y2, z1);

      if( reverse )
      {
        quadstrip1.push_back( std::make_pair(p2, n1) );
        quadstrip1.push_back( std::make_pair(p1, n1) );
      }
      else
      {
        quadstrip1.push_back( std::make_pair(p1, n1) );
        quadstrip1.push_back( std::make_pair(p2, n1) );
      }
    }

    quadstrips.push_back( quadstrip1 );
  }
}


template <class T, class S>
void GeomGlyphBase::gen_ellipsoid(const Point& center, const T& t,
				  S scales,
				  int nu, int nv, int half,
				  std::vector< QuadStrip >& quadstrips )
{
  nu++; //Bring nu to expected value for shape.

  double start=0, stop=M_PI;

  // Half ellipsoid criteria.
  if( half == -1) start = M_PI/2.0;
  if( half ==  1) stop  = M_PI/2.0;
  if( half != 0 ) nv /= 2;

  // Should only happen when doing half ellipsoids.
  if( nv < 2 ) nv = 2;

  SinCosTable tab1(nu, 0, 2*M_PI);
  SinCosTable tab2(nv, start, stop);
  
  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  trans.post_scale ( Vector(1.0,1.0,1.0) * scales );
  rotate.post_scale( Vector(1.0,1.0,1.0) / scales );

  // Draw the ellipsoid
  for (int v=0; v<nv-1; v++)
  {
    double nr1 = tab2.sin(v+1);
    double nr2 = tab2.sin(v  );
    
    double nz1 = tab2.cos(v+1);
    double nz2 = tab2.cos(v  );
    
    QuadStrip quadstrip;

    for (int u=0; u<nu; u++)
    {
      double nx = tab1.sin(u);
      double ny = tab1.cos(u);
	    
      double x1 = nr1 * nx;
      double y1 = nr1 * ny;
      double z1 = nz1;

      double x2 = nr2 * nx;
      double y2 = nr2 * ny;
      double z2 = nz2;
	    
      Point p1 = trans * Point(x1, y1, z1);
      Point p2 = trans * Point(x2, y2, z2);

      Vector v1 = rotate * Vector( x1, y1, z1 );
      Vector v2 = rotate * Vector( x2, y2, z2 );

      v1.safe_normalize();
      v2.safe_normalize();

      quadstrip.push_back( std::make_pair(p1, v1) );
      quadstrip.push_back( std::make_pair(p2, v2) );
    }
	
    quadstrips.push_back( quadstrip );
  }
}


template <class T>
void GeomGlyphBase::gen_helix(const Point& center, const T& t,
			      double major_radius1, double major_radius2,
			      double minor_radius, double length,
			      unsigned int pitch,
			      int nu, int nv,
			      std::vector< QuadStrip >& quadstrips )
{
  nu++; //Bring nu to expected value for shape.

  SinCosTable tab1(nu, 0, 2*M_PI);
  SinCosTable tab2(nv, 0, 2*M_PI, minor_radius);

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  // Draw the helix

  // The total number of sections is based on the pitch times the
  // number of sections in each rotation which is nu-1.
  double nsections = (double) (pitch * (nu-1));

  double dz = length / nsections;
  double dr = (major_radius2-major_radius1) / nsections;

  for (int v=0; v < nv-1; v++)
  {
    Vector base1 = Vector(tab2.sin(v  ), 0, tab2.cos(v  ) );
    Vector base2 = Vector(tab2.sin(v+1), 0, tab2.cos(v+1) );

    QuadStrip quadstrip;
      
    for (unsigned int p=0; p<pitch; p++)
    {
      unsigned int full = (p == pitch-1) ? 0 : 1;
      
      for (unsigned int u=0; u<nu-full; u++)
      {
        double section = p * (nu-1) + u;

        double zz = dz * section;
        double rad = major_radius1 + dr * section;

        // Base point for this section
        Point base( rad, 0, zz );

        // Cord length of the section - this is used to get the
        // angle between the section and the pitch.
        double dl = (2.0*rad + dr) * sin( 2.0*M_PI/(nu-1.0) );

        // Rotation about x to get the correct pitch.
        Transform pitch_rotate;
        pitch_rotate.post_rotate(atan2(dz,dl), Vector(1,0,0));

        // Rotation about z to place the normal into the correct xy
        // location.
        Transform z_rotate;
        z_rotate.load_identity();

        z_rotate.set_mat_val( 0, 0,  tab1.sin(u));
        z_rotate.set_mat_val( 0, 1,  tab1.cos(u));
        z_rotate.set_mat_val( 1, 0,  tab1.cos(u));
        z_rotate.set_mat_val( 1, 1, -tab1.sin(u));

        // Tilt the two base normals to the correct pitch. Then rotate
        // into the correct xy location.
        Vector n1 = z_rotate * (pitch_rotate * base1);
        Vector n2 = z_rotate * (pitch_rotate * base2);

        // Rotate the centerline into the correct xy location.
        Point centerline = z_rotate * base;

        // Offest from the centerline then the final world
        // rotation/translation.
        Point p1 = trans * (centerline + n1);
        Point p2 = trans * (centerline + n2);

        // Final world rotation.
        Vector v1 = rotate * -n1;
        Vector v2 = rotate * -n2;
        
        v1.safe_normalize();
        v2.safe_normalize();
        
        quadstrip.push_back( std::make_pair(p1, v1) );
        quadstrip.push_back( std::make_pair(p2, v2) );
      }
    }

    quadstrips.push_back( quadstrip );
  }

}


template <class T>
void GeomGlyphBase::gen_superquadric(const Point& center, const T& t,
				     Vector scales, double A, double B,
				     int axis,			     
				     int nu, int nv,
				     std::vector< QuadStrip >& quadstrips )
{
  nu++; //Bring nu to expected value for shape.

  SinCosTable tab1(nu, 0, 2*M_PI);
  SinCosTable tab2(nv, 0,   M_PI);

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  trans.post_scale( Vector(1.0,1.0,1.0) * scales);
  rotate.post_scale(Vector(1.0,1.0,1.0) / scales);

  double nr[2];
  double nz[2];

  Point point;
  Vector normal;

  for (int v=0; v < nv-1; v++)
  {
    nr[0] = tab2.sin(v+1);
    nr[1] = tab2.sin(v);

    nz[0] = tab2.cos(v+1);
    nz[1] = tab2.cos(v);

    QuadStrip quadstrip;

    for (int u=0; u<nu; u++)
    {
      double nx = tab1.sin(u);
      double ny = tab1.cos(u);

      for( unsigned int i=0; i<2; i++ )
      {
        const double x = spow(nr[i], B) * spow(nx, A);
        const double y = spow(nr[i], B) * spow(ny, A);
        const double z = spow(nz[i], B);
        
        const float nnx = spow(nr[i], 2.0-B) * spow(nx, 2.0-A);
        const float nny = spow(nr[i], 2.0-B) * spow(ny, 2.0-A);
        const float nnz = spow(nz[i], 2.0-B);
        
        Point  point  = trans  * Point( x, y, z );
        Vector normal = rotate * Vector( nnx, nny, nnz );

        normal.safe_normalize();

        quadstrip.push_back( std::make_pair(point, normal) );
      }
    }

    quadstrips.push_back( quadstrip );
  }
}

template <class T>
void GeomGlyphBase::gen_torus(const Point& center, const T& t,
			      double major_radius, double minor_radius,
			      int nu, int nv,
			      std::vector< QuadStrip >& quadstrips )
{
  nu++; //Bring nu to expected value for shape.

  SinCosTable tab1(nu, 0, 2*M_PI);
  SinCosTable tab2(nv, 0, 2*M_PI, minor_radius);

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  // Draw the torus
  for (int v=0; v<nv-1; v++)
  {
    double z1 = tab2.cos(v+1);
    double z2 = tab2.cos(v);

    double nr1 = tab2.sin(v+1);
    double nr2 = tab2.sin(v);

    double r1 = major_radius + nr1;
    double r2 = major_radius + nr2;

    QuadStrip quadstrip;

    for (int u=0; u<nu; u++)
    {
      double nx = tab1.sin(u);
      double ny = tab1.cos(u);

      double x1 = r1 * nx;
      double y1 = r1 * ny;

      double x2 = r2 * nx;
      double y2 = r2 * ny;

      Point p1 = trans * Point(x1, y1, z1);
      Point p2 = trans * Point(x2, y2, z2);

      Vector v1 = rotate * Vector(nr1*nx, nr1*ny, z1);
      Vector v2 = rotate * Vector(nr2*nx, nr2*ny, z2);

      v1.safe_normalize();
      v2.safe_normalize();

      quadstrip.push_back( std::make_pair(p1, v1) );
      quadstrip.push_back( std::make_pair(p2, v2) );
    }

    quadstrips.push_back( quadstrip );
  }
}


// Arrow
template <class T>
void GeomGlyphBase::add_arrow(const Point& center, const T& t,
			      double radius, double length,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  double ratio = 2.0;

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  Vector offset = rotate * Vector(0,0,1);
  offset.safe_normalize();
  offset *= length * ratio;

  gen_cylinder(center, t,
	       radius/10.0, radius/10.0, length*ratio, nu, nv, quadstrips);

  gen_cylinder(center+offset, t,
	       radius, 0.0, length, nu, nv, quadstrips);

  add_strip( quadstrips );
}

template <class T>
void GeomGlyphBase::add_arrow(const Point& center, const T& t,
			      double radius, double length,
			      const MaterialHandle &mat,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  double ratio = 2.0;

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  Vector offset = rotate * Vector(0,0,1);
  offset.safe_normalize();
  offset *= length * ratio;

  gen_cylinder(center, t,
	       radius/10.0, radius/10.0, length*ratio, nu, nv, quadstrips);

  gen_cylinder(center+offset, t,
	       radius, 0.0, length, nu, nv, quadstrips);

  add_strip( quadstrips, mat );
}

template <class T>
void GeomGlyphBase::add_arrow(const Point& center, const T& t,
			      double radius, double length,
			      const double cindex,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  double ratio = 2.0;

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  Vector offset = rotate * Vector(0,0,1);
  offset.safe_normalize();
  offset *= length * ratio;

  gen_cylinder(center, t,
	       radius/10.0, radius/10.0, length*ratio, nu, nv, quadstrips);

  gen_cylinder(center+offset, t,
	       radius, 0.0, length, nu, nv, quadstrips);

  add_strip( quadstrips, cindex );
}

template <class T>
void GeomGlyphBase::add_arrow(const Point& center, const T& t,
			      double radius, double length,
			      const MaterialHandle &mat,
			      const double cindex,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  double ratio = 2.0;

  Transform trans;
  Transform rotate;
  gen_transforms( center, t, trans, rotate );

  Vector offset = rotate * Vector(0,0,1);
  offset.safe_normalize();
  offset *= length * ratio;

  gen_cylinder(center, t,
	       radius/10.0, radius/10.0, length*ratio, nu, nv, quadstrips);

  gen_cylinder(center+offset, t,
	       radius, 0.0, length, nu, nv, quadstrips);

  add_strip( quadstrips, mat, cindex );
}

// Box
template <class T>
void GeomGlyphBase::add_box(const Point& center, const T& t,
			    double x_side, double y_side, double z_side,
			    const bool RGB)
{
  std::vector< QuadStrip > quadstrips;

  gen_box(center, t, x_side, y_side, z_side, quadstrips);

  // Instead of the default color make the boxes RGB.
  if( RGB )
  {
    MaterialHandle mat[3];

    mat[0] = new Material();
    mat[0]->transparency = 1.0;
    mat[0]->diffuse = Color( 1.0, 0.0, 0.0 );
    
    mat[1] = new Material();
    mat[1]->transparency = 1.0;
    mat[1]->diffuse = Color( 0.0, 1.0, 0.0 );
    
    mat[2] = new Material();
    mat[2]->transparency = 1.0;
    mat[2]->diffuse = Color( 0.0, 0.0, 1.0 );
    
    for( unsigned int i=0; i<quadstrips.size()/2; i++ ) {

      std::vector< QuadStrip > quadstrips_rgb;

      quadstrips_rgb.push_back( quadstrips[i] );
      quadstrips_rgb.push_back( quadstrips[i+3] );

      add_strip( quadstrips_rgb, mat[i] );
    }
  } else
    add_strip( quadstrips );
}


template <class T>
void GeomGlyphBase::add_box(const Point& center, const T& t,
			    double x_side, double y_side, double z_side,
			    const MaterialHandle &mat)
{
  std::vector< QuadStrip > quadstrips;

  gen_box(center, t, x_side, y_side, z_side, quadstrips);

  add_strip( quadstrips, mat );
}

template <class T>
void GeomGlyphBase::add_box(const Point& center, const T& t,
			    double x_side, double y_side, double z_side,
			    const double cindex)
{
  std::vector< QuadStrip > quadstrips;

  gen_box(center, t, x_side, y_side, z_side, quadstrips);

  add_strip( quadstrips, cindex );
}

template <class T>
void GeomGlyphBase::add_box(const Point& center, const T& t,
			    double x_side, double y_side, double z_side,
			    const MaterialHandle &mat, const double cindex)
{
  std::vector< QuadStrip > quadstrips;

  gen_box(center, t, x_side, y_side, z_side, quadstrips);

  add_strip( quadstrips, mat, cindex );
}



// Capped Cylinder
template <class T>
void GeomGlyphBase::add_capped_cylinder(const Point& center, const T& t,
					double radius1, double radius2,
					double length,
					int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_capped_cylinder(center, t,
		      radius1, radius2, length, nu, nv, quadstrips);

  add_strip( quadstrips );
}


template <class T>
void GeomGlyphBase::add_capped_cylinder(const Point& center, const T& t,
					double radius1, double radius2,
					double length,
					const MaterialHandle &mat,
					int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_capped_cylinder(center, t,
		      radius1, radius2, length, nu, nv, quadstrips);

  add_strip( quadstrips, mat );
}


template <class T>
void GeomGlyphBase::add_capped_cylinder(const Point& center, const T& t,
					double radius1, double radius2,
					double length,
					const double cindex,
					int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_capped_cylinder(center, t,
		      radius1, radius2, length, nu, nv, quadstrips);

  add_strip( quadstrips, cindex );
}

template <class T>
void GeomGlyphBase::add_capped_cylinder(const Point& center, const T& t,
					double radius1, double radius2,
					double length,
					const MaterialHandle &mat,
					const double cindex,
					int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_capped_cylinder(center, t,
		      radius1, radius2, length, nu, nv, quadstrips);

  add_strip( quadstrips, mat, cindex );
}

// Cylinder
template <class T>
void GeomGlyphBase::add_cylinder(const Point& center, const T& t,
				 double radius1, double radius2,
				 double length,
				 int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_cylinder(center, t,
	       radius1, radius2, length, nu, nv, quadstrips);

  add_strip( quadstrips );
}


template <class T>
void GeomGlyphBase::add_cylinder(const Point& center, const T& t,
				 double radius1, double radius2,
				 double length,
				 const MaterialHandle &mat,
				 int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_cylinder(center, t,
	       radius1, radius2, length, nu, nv, quadstrips);

  add_strip( quadstrips, mat );
}

template <class T>
void GeomGlyphBase::add_cylinder(const Point& center, const T& t,
				 double radius1, double radius2,
				 double length,
				 const MaterialHandle &mat0,
				 const MaterialHandle &mat1,
				 int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_cylinder(center, t,
	       radius1, radius2, length, nu, nv, quadstrips);

  MaterialHandle mat = new Material();

  double dt = (mat1->transparency - mat0->transparency) / (double) (nv-1);
  Color  dc = (mat1->diffuse      - mat0->diffuse     ) / (double) (nv-1);

  for (int v=0; v<nv; v++)
  {
    mat->transparency = mat0->transparency + dt * (double) v;
    mat->diffuse      = mat0->diffuse      + dc * (double) v;

    std::vector< QuadStrip > quadstrips_d;

    quadstrips_d.push_back( quadstrips[v] );

    add_strip( quadstrips_d, mat );
  }
}


template <class T>
void GeomGlyphBase::add_cylinder(const Point& center, const T& t,
				 double radius1, double radius2,
				 double length,
				 const double cindex,
				 int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_cylinder(center, t,
	       radius1, radius2, length, nu, nv, quadstrips);

  add_strip( quadstrips, cindex );
}

template <class T>
void GeomGlyphBase::add_cylinder(const Point& center, const T& t,
				 double radius1, double radius2,
				 double length,
				 const MaterialHandle &mat,
				 const double cindex,
				 int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_cylinder(center, t,
	       radius1, radius2, length, nu, nv, quadstrips);

  add_strip( quadstrips, mat, cindex );
}

template <class T>
void GeomGlyphBase::add_cylinder(const Point& center, const T& t,
				 double radius1, double radius2,
				 double length,
				 const double cindex0, const double cindex1,
				 int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_cylinder(center, t,
	       radius1, radius2, length, nu, nv, quadstrips);

  double cindex;

  double di = (cindex1 - cindex0) / (double) (nv-1);

  for (int v=0; v<nv; v++)
  {
    cindex = cindex0 + di * (double) v;

    std::vector< QuadStrip > quadstrips_d;

    quadstrips_d.push_back( quadstrips[v] );

    add_strip( quadstrips_d, cindex );
  }
}

template <class T>
void GeomGlyphBase::add_cylinder(const Point& center, const T& t,
				 double radius1, double radius2,
				 double length,
				 const MaterialHandle &mat0,
				 const double cindex0,
				 const MaterialHandle &mat1,
				 const double cindex1,
				 int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_cylinder(center, t,
	       radius1, radius2, length, nu, nv, quadstrips);

  MaterialHandle mat = new Material();

  double dt = (mat1->transparency - mat0->transparency) / (double) (nv-1);
  Color  dc = (mat1->diffuse      - mat0->diffuse     ) / (double) (nv-1);

  double cindex;

  double di = (cindex1 - cindex0) / (double) nv;

  for (int v=0; v<nv; v++) {
    cindex = cindex0 + di * (double) v;

    mat->transparency = mat0->transparency + dt * (double) v;
    mat->diffuse      = mat0->diffuse      + dc * (double) v;

    std::vector< QuadStrip > quadstrips_d;

    quadstrips_d.push_back( quadstrips[v] );

    add_strip( quadstrips_d, mat, cindex );
  }
}


// Ellipsoid
template <class T, class S>
void GeomGlyphBase::add_ellipsoid(const Point& center, const T& t,
				  S scales,
				  int nu, int nv, int half)
{
  std::vector< QuadStrip > quadstrips;

  gen_ellipsoid(center, t, scales, nu, nv, half, quadstrips);

  add_strip( quadstrips );
}


template <class T, class S>
void GeomGlyphBase::add_ellipsoid(const Point& center, const T& t,
				  S scales,
				  const MaterialHandle &mat,
				  int nu, int nv, int half)
{
  std::vector< QuadStrip > quadstrips;

  gen_ellipsoid(center, t, scales, nu, nv, half, quadstrips);

  add_strip( quadstrips, mat );
}


template <class T, class S>
void GeomGlyphBase::add_ellipsoid(const Point& center, const T& t,
				  S scales,
				  const double cindex,
				  int nu, int nv, int half)
{
  std::vector< QuadStrip > quadstrips;
  
  gen_ellipsoid(center, t, scales, nu, nv, half, quadstrips);

  add_strip( quadstrips, cindex );
}

template <class T, class S>
void GeomGlyphBase::add_ellipsoid(const Point& center, const T& t,
				  S scales,
				  const MaterialHandle &mat,
				  const double cindex,
				  int nu, int nv, int half)
{
  std::vector< QuadStrip > quadstrips;

  gen_ellipsoid(center, t, scales, nu, nv, half, quadstrips);

  add_strip( quadstrips, mat, cindex );
}


// Helix
template <class T>
void GeomGlyphBase::add_helix(const Point& center, const T& t,
			      double major_radius1, double major_radius2,
			      double minor_radius, double length,
			      unsigned int pitch,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_helix(center, t, major_radius1, major_radius2,
	    minor_radius, length, pitch, nu, nv, quadstrips);

  add_strip( quadstrips );
}


template <class T>
void GeomGlyphBase::add_helix(const Point& center, const T& t,
			      double major_radius1, double major_radius2,
			      double minor_radius, double length,
			      unsigned int pitch,
			      const MaterialHandle &mat,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_helix(center, t, major_radius1, major_radius2,
	    minor_radius, length, pitch, nu, nv, quadstrips);

  add_strip( quadstrips, mat );
}

template <class T>
void GeomGlyphBase::add_helix(const Point& center, const T& t,
			      double major_radius1, double major_radius2,
			      double minor_radius, double length,
			      unsigned int pitch,
			      const MaterialHandle &mat0,
			      const MaterialHandle &mat1,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_helix(center, t, major_radius1, major_radius2,
	    minor_radius, length, pitch, nu, nv, quadstrips);

  MaterialHandle mat = new Material();

  double dt = (mat1->transparency - mat0->transparency) / (double) (nv-1);
  Color  dc = (mat1->diffuse      - mat0->diffuse     ) / (double) (nv-1);

  for (unsigned int v=0; v<nv; v++) {
    mat->transparency = mat0->transparency + dt * (double) v;
    mat->diffuse      = mat0->diffuse      + dc * (double) v;

    std::vector< QuadStrip > quadstrips_d;

    quadstrips_d.push_back( quadstrips[v] );

    add_strip( quadstrips_d, mat );
  }
}


template <class T>
void GeomGlyphBase::add_helix(const Point& center, const T& t,
			      double major_radius1, double major_radius2,
			      double minor_radius, double length,
			      unsigned int pitch,
			      const double cindex,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_helix(center, t, major_radius1, major_radius2,
	    minor_radius, length, pitch, nu, nv, quadstrips);

  add_strip( quadstrips, cindex );
}

template <class T>
void GeomGlyphBase::add_helix(const Point& center, const T& t,
			      double major_radius1, double major_radius2,
			      double minor_radius, double length,
			      unsigned int pitch,
			      const MaterialHandle &mat,
			      const double cindex,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_helix(center, t, major_radius1, major_radius2,
	    minor_radius, length, pitch, nu, nv, quadstrips);

  add_strip( quadstrips, mat, cindex );
}

template <class T>
void GeomGlyphBase::add_helix(const Point& center, const T& t,
			      double major_radius1, double major_radius2,
			      double minor_radius, double length,
			      unsigned int pitch,
			      const double cindex0, const double cindex1,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_helix(center, t, major_radius1, major_radius2,
	    minor_radius, length, pitch, nu, nv, quadstrips);

  double cindex;

  double di = (cindex1 - cindex0) / (double) (nv-1);

  for (unsigned int v=0; v<nv; v++) {
    cindex = cindex0 + di * (double) v;

    std::vector< QuadStrip > quadstrips_d;

    quadstrips_d.push_back( quadstrips[v] );

    add_strip( quadstrips_d, cindex );
  }
}

template <class T>
void GeomGlyphBase::add_helix(const Point& center, const T& t,
			      double major_radius1, double major_radius2,
			      double minor_radius, double length,
			      unsigned int pitch,
			      const MaterialHandle &mat0,
			      const double cindex0,
			      const MaterialHandle &mat1,
			      const double cindex1,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_helix(center, t, major_radius1, major_radius2,
	    minor_radius, length, pitch, nu, nv, quadstrips);

  MaterialHandle mat = new Material();

  double dt = (mat1->transparency - mat0->transparency) / (double) (nv-1);
  Color  dc = (mat1->diffuse      - mat0->diffuse     ) / (double) (nv-1);

  double cindex;

  double di = (cindex1 - cindex0) / (double) nv;

  for (unsigned int v=0; v<nv; v++) 
  {
    cindex = cindex0 + di * (double) v;

    mat->transparency = mat0->transparency + dt * (double) v;
    mat->diffuse      = mat0->diffuse      + dc * (double) v;

    std::vector< QuadStrip > quadstrips_d;

    quadstrips_d.push_back( quadstrips[v] );

    add_strip( quadstrips_d, mat, cindex );
  }
}


// Superquadric
template <class T>
void GeomGlyphBase::add_superquadric(const Point& center, const T& t,
				     Vector scales, double A, double B, 
				     int axis,
				     int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_superquadric(center, t, scales, A, B, axis, nu, nv, quadstrips);

  add_strip( quadstrips );
}

template <class T>
void GeomGlyphBase::add_superquadric(const Point& center, const T& t,
				     Vector scales, double A, double B, 
				     int axis,
				     const MaterialHandle &mat,
				     int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_superquadric(center, t, scales, A, B, axis, nu, nv, quadstrips);

  add_strip( quadstrips, mat );
}

template <class T>
void GeomGlyphBase::add_superquadric(const Point& center, const T& t,
				     Vector scales, double A, double B, 
				     int axis,
				     const double cindex,
				     int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_superquadric(center, t, scales, A, B, axis, nu, nv, quadstrips);

  add_strip( quadstrips, cindex );
}

template <class T>
void GeomGlyphBase::add_superquadric(const Point& center, const T& t,
				     Vector scales, double A, double B, 
				     int axis,
				     const MaterialHandle &mat,
				     const double cindex,
				     int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_superquadric(center, t, scales, A, B, axis, nu, nv, quadstrips);

  add_strip( quadstrips, mat, cindex );
}

// Torus
template <class T>
void GeomGlyphBase::add_torus(const Point& center, const T& t,
			      double major_radius, double minor_radius,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_torus(center, t, major_radius, minor_radius, nu, nv, quadstrips);

  add_strip( quadstrips );
}

template <class T>
void GeomGlyphBase::add_torus(const Point& center, const T& t,
			      double major_radius, double minor_radius,
			      const MaterialHandle &mat,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_torus(center, t, major_radius, minor_radius, nu, nv, quadstrips);

  add_strip( quadstrips, mat );
}

template <class T>
void GeomGlyphBase::add_torus(const Point& center, const T& t,
			      double major_radius, double minor_radius,
			      const double cindex,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_torus(center, t, major_radius, minor_radius, nu, nv, quadstrips);

  add_strip( quadstrips, cindex );
}

template <class T>
void GeomGlyphBase::add_torus(const Point& center, const T& t,
			      double major_radius, double minor_radius,
			      const MaterialHandle &mat,
			      const double cindex,
			      int nu, int nv)
{
  std::vector< QuadStrip > quadstrips;

  gen_torus(center, t, major_radius, minor_radius, nu, nv, quadstrips);

  add_strip( quadstrips, mat, cindex );
}


} // End namespace SCIRun

#endif /* SCI_Geom_Glyph_h */
