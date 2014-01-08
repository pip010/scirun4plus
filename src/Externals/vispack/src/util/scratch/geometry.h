// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: geometry.h,v 1.1.1.1 2003/02/12 16:51:54 whitaker Exp $



#ifndef	MATH_GEOMETRY_H
#define	MATH_GEOMETRY_H

#include <math.h>
#include "util/mathutil.h"
//#include "matrix/matrix.h"
class VISVISMatrix;

#ifndef ERROR
#define ERROR(a) printf("ERROR: %s \n", a);
#endif

// ---------------------------------------------------------------------------

enum VISCoordIndex			// used to specify rotations
{
    X = 0,
    Y = 1,
    Z = 2
};

typedef float GeomReal;

class Point2;
class Point3;
class Point4;
class VISMatrix4;
class Transform;
class Quaternion;


class Point2
{
  public:
    Point2()	{ x(0); y(0); }
    Point2(GeomReal xx, GeomReal yy) { x(xx); y(yy); }
    Point2(const GeomReal *p) 	{ x(p[0]);  y(p[1]); }
    Point2(const Point2& p)	{ operator=(p); }
    
    // provide access to indivual components.;
    GeomReal	x()	const	{ return _d[0]; }
    GeomReal	y()	const	{ return _d[1]; }
    
    void	x(GeomReal x)	{ _d[0] = x; }
    void	y(GeomReal y)	{ _d[1] = y; }
    
    // Allow indexing on point objects
    // Range checking is performed when debugging is enabled
    GeomReal&	at(const unsigned int i) {
#ifdef DEBUG
	if (i > 1)
	    ERROR("Range check failure: Point2: operator[]");
#endif
	return _d[i];
    }
    const GeomReal&	operator[](const unsigned i) const {
#ifdef DEBUG
	if (i > 1)
	    ERROR("Range check failure: Point2: operator[]");
#endif
	return _d[i];
    }
    
    boolean operator==(const Point2& p) const {
	return (p.x() == x() && p.y() == y());
    }
    boolean operator!=(const Point2& p) const	{ return !operator==(p); }

    Point2& operator*=(GeomReal s) { _d[0] *= s;  _d[1] *= s; return *this; }
    Point2& operator/=(GeomReal s) { _d[0] /= s;  _d[1] /= s; return *this; }
    Point2& operator+=(const Point2& p)	{ _d[0] += p.x(); _d[1] += p.y();
					  return *this; }
    Point2& operator-=(const Point2& p)	{ _d[0] -= p.x(); _d[1] -= p.y();
					  return *this; }
    // misc
    void normalize();
    Point2 normal() const;
    Point2 dominant() const;

    Point2& operator=(const Point2& p) { x(p.x()); y(p.y()); return *this; }
    Point2& operator=(const GeomReal* p){ x(p[0]); y(p[1]); return *this; }

    void print(ostream& ostr = cout) const;
  protected:
    GeomReal	_d[2];
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Inline functions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline GeomReal length2(const Point2& p)
{
    return ( sqr(p.x()) + sqr(p.y()) );
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline GeomReal length(const Point2& p)
{
    return sqrt(length2(p));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point2 operator-(const Point2& p)
{
    return Point2(-p.x(), -p.y());
}

// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

inline Point2 operator+(const Point2& a, const Point2& b)
{
    Point2 r(a); r += b; return r;
}

// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

inline Point2 operator-(const Point2& a, const Point2& b)
{
    Point2 r(a); r -= b; return r;
}

// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

inline Point2 operator*(GeomReal s, const Point2& a)
{
    Point2 r(a); r *= s; return r;
}

// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

inline Point2 operator*(const Point2& a, GeomReal s)
{
    Point2 r(a); r *= s; return r;
}

// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

inline Point2 operator/(const Point2& a, GeomReal s)
{
    Point2 r(a); r /= s; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline GeomReal operator*(const Point2& a, const Point2& b)
{
    GeomReal f;
    
    f  = a.x() * b.x();
    f += a.y() * b.y();
    
    return f;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point2	average(const Point2& p1,const Point2& p2)
{
    return Point2((p1.x() + p2.x())/2.0,
		  (p1.y() + p2.y())/2.0);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point2	average(const Point2& p1,const Point2& p2,const Point2& p3)
{
    return Point2((p1.x() + p2.x() + p3.x())/3.0,
		  (p1.y() + p2.y() + p3.y())/3.0);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point2	average(const Point2& p1,const Point2& p2,
			const Point2& p3,const Point2& p4)
{
    return Point2((p1.x() + p2.x() + p3.x() + p4.x())/4.0,
		  (p1.y() + p2.y() + p3.y() + p4.y())/4.0);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline ostream& operator<<(ostream& ostr, const Point2& p)
{
    p.print(ostr);
    return ostr;
}


// geometric coordinate in 3D space;
class Point3		
{
  protected:
    GeomReal	_d[3];
  public:
    Point3()	{ x(0); y(0); z(0); }
    Point3(GeomReal xx, GeomReal yy, GeomReal zz) { x(xx);  y(yy);  z(zz); }
    Point3(const GeomReal *p) 	{ x(p[0]);  y(p[1]);  z(p[2]); }
    Point3(const Point3& p)	{ operator=(p); }
    Point3(const Point4& p)	{ operator=(p); }
    
    // provide access to indivual components.;
    GeomReal	x()	const	{ return _d[0]; }
    GeomReal	y()	const	{ return _d[1]; }
    GeomReal	z()	const	{ return _d[2]; }
    
    void	x(GeomReal x)	{ _d[0] = x; }
    void	y(GeomReal y)	{ _d[1] = y; }
    void	z(GeomReal z)	{ _d[2] = z; }
    
    // Allow indexing on point objects
    // Range checking is performed when debugging is enabled
    GeomReal&	at(const unsigned int i) {
#ifdef DEBUG
	if (i > 2)
	    ERROR("Range check failure: Point3: operator[]");
#endif
	return _d[i];
    }
    const GeomReal&	operator[](const unsigned i) const {
#ifdef DEBUG
	if (i > 2)
	    ERROR("Range check failure: Point3: operator[]");
#endif
	return _d[i];
    }
    
    boolean operator==(const Point3& p) const {
	return (p.x() == x() && p.y() == y() && p.z() == z());
    }
    boolean operator!=(const Point3& p) const	{ return !operator==(p); }

    Point3& operator*=(GeomReal s) { _d[0] *= s;  _d[1] *= s;  _d[2] *= s; return *this; }
    Point3& operator/=(GeomReal s) { _d[0] /= s;  _d[1] /= s;  _d[2] /= s; return *this; }
    Point3& operator+=(const Point3& p)	{ _d[0] += p.x(); _d[1] += p.y(); _d[2] += p.z();
					  return *this; }
    Point3& operator-=(const Point3& p)	{ _d[0] -= p.x(); _d[1] -= p.y(); _d[2] -= p.z();
					  return *this; }
    // Point self-multiplications
    Point3& operator*=(const VISMatrix4&);
    Point3& operator*=(const Quaternion&);
    
    // rotations
    void rotate_only(const Transform&);

    // scalings
    void scale(GeomReal);
    void scale(GeomReal,GeomReal,GeomReal);

    // misc
    void normalize();
    Point3 normal() const;
    Point3 dominant() const;

    // general transform(s)
    Point3& transform(const Transform* t);
    Point3& transform(const Transform& t);

    Point3& operator=(const Point3& p) { x(p.x()); y(p.y()); z(p.z()); return *this; }
    Point3& operator=(const Point4&);
    Point3& operator=(const GeomReal* p){ x(p[0]); y(p[1]); z(p[2]); return *this; }

    void print(ostream& ostr = cout) const;
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Point3 operator declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Transform operator|(const Point3&,const Transform&);
Transform operator&(const Point3&,const Transform&);
Point3 normal(const Point3&);
Point3 dominant(Point3&);
Point3 operator*(const VISMatrix4&,const Point3&);
Point3 operator*(const Point3&,const VISMatrix4&);
Point3 operator*(const Point3&,const Quaternion&);
Point3 operator*(const Quaternion&,const Point3&);
Point3 rotate_only(const Point3&,const Transform&);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Inline functions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline GeomReal length2(const Point3& p)
{
    return ( sqr(p.x()) + sqr(p.y()) + sqr(p.z()) );
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline GeomReal length(const Point3& p)
{
    return sqrt(length2(p));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point3 operator-(const Point3& p)
{
    return Point3(-p.x(), -p.y(), -p.z());
}

// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

inline Point3 operator+(const Point3& a, const Point3& b)
{
    Point3 r(a); r += b; return r;
}

// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

inline Point3 operator-(const Point3& a, const Point3& b)
{
    Point3 r(a); r -= b; return r;
}

// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

inline Point3 operator*(GeomReal s, const Point3& a)
{
    Point3 r(a); r *= s; return r;
}

// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

inline Point3 operator*(const Point3& a, GeomReal s)
{
    Point3 r(a); r *= s; return r;
}

// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

inline Point3 operator/(const Point3& a, GeomReal s)
{
    Point3 r(a); r /= s; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline GeomReal operator*(const Point3& a, const Point3& b)
{
    GeomReal f;
    
    f  = a.x() * b.x();
    f += a.y() * b.y();
    f += a.z() * b.z();
    
    return f;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point3 operator^(const Point3& a, const Point3& b)
{
    return Point3(
	   (a.y() * b.z()) - (a.z() * b.y()),
	   (a.z() * b.x()) - (a.x() * b.z()),
	   (a.x() * b.y()) - (a.y() * b.x()));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point3	average(const Point3& p1,const Point3& p2)
{
    return Point3((p1.x() + p2.x())/2.0,
		  (p1.y() + p2.y())/2.0,
		  (p1.z() + p2.z())/2.0);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point3	average(const Point3& p1,const Point3& p2,const Point3& p3)
{
    return Point3((p1.x() + p2.x() + p3.x())/3.0,
		  (p1.y() + p2.y() + p3.y())/3.0,
		  (p1.z() + p2.z() + p3.z())/3.0);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point3	average(const Point3& p1,const Point3& p2,
			const Point3& p3,const Point3& p4)
{
    return Point3((p1.x() + p2.x() + p3.x() + p4.x())/4.0,
		  (p1.y() + p2.y() + p3.y() + p4.y())/4.0,
		  (p1.z() + p2.z() + p3.z() + p4.z())/4.0);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point3	scale(const Point3& p,GeomReal f)
{
    Point3 r(p);
    r.scale(f);
    return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point3	scale(const Point3& p,GeomReal fx,GeomReal fy,GeomReal fz)
{
    Point3 r(p);
    r.scale(fx,fy,fz);
    return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline ostream& operator<<(ostream& ostr, const Point3& p)
{
    p.print(ostr);
    return ostr;
}


// Homogenous point in 3D space of the form x, y, z, w

class Point4 	
{
  protected:
    GeomReal _d[4];
    
  public:
    Point4 () { x(0); y(0); z(0); w(1); }
    Point4 (GeomReal xx, GeomReal yy, GeomReal zz, GeomReal ww = 1.0) {
	x(xx); y(yy); z(zz); w(ww); }
    Point4 (const GeomReal *p) { x(p[0]); y(p[1]); z(p[2]); w(p[3]); }
    Point4(const Point4& p) { operator=(p); }
    Point4(const Point3& p) { operator=(p); }
    
    GeomReal x()	const	{ return _d[0];}
    GeomReal y()	const	{ return _d[1];}
    GeomReal z()	const	{ return _d[2];}
    GeomReal w()	const	{ return _d[3];}
    
    void x(GeomReal x)	{ _d[0] = x;}
    void y(GeomReal y)	{ _d[1] = y;}
    void z(GeomReal z)	{ _d[2] = z;}
    void w(GeomReal w)	{ _d[3] = w;}
    
    boolean operator==(const Point3&) const;
    boolean operator!=(const Point3& p) const	{ return !operator==(p); }
    boolean operator==(const Point4&) const;
    boolean operator!=(const Point4& p) const	{ return !operator==(p); }
    
    // Scales and adds;
    Point4& operator*=(GeomReal s) { _d[0] *= s;  _d[1] *= s;  _d[2] *= s;
				  return *this; }
    Point4& operator+=(const Point3&);
    Point4& operator-=(const Point3&);
    Point4& operator*=(const VISMatrix4&);
    
    GeomReal&	at(const unsigned int i) {
#ifdef DEBUG
	if (i > 3)
	    ERROR("Range check failure: Point4: at()");
#endif
	return _d[i];
    }
    const GeomReal&	operator[](const unsigned int i) const {
#ifdef DEBUG
	if (i > 3)
	    ERROR("Range check failure: Point4: operator[]");
#endif
	return _d[i];
    }

    // other useful operations
    void homogenize();
    
    // rotations
    void rotate_only(const Transform&);

    // general transformations
    Point4& transform(const Transform* t);
    Point4& transform(const Transform& t);

    Point4& operator=(const Point4& p) { x(p.x()); y(p.y()); z(p.z()); w(p.w());
					 return *this; }
    Point4& operator=(const Point3& p) { x(p.x()); y(p.y()); z(p.z()); w(1.0);
					 return *this; }
    Point4& operator=(const GeomReal* p){ x(p[0]); y(p[1]); z(p[2]); w(p[3]);
					 return *this; }

    void print(ostream& ostr = cout) const;
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Point4 operator declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Point4 operator*(const VISMatrix4&,const Point4&);
Point4 rotate_only(const Point4&,const Transform&);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Inline functions

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point4 operator-(const Point4& p)
{
    return Point4 (-p.x(), -p.y(), -p.z(), p.w());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point4 operator+(const Point4& a, const Point3& b)
{
    Point4 r(a); r += b; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point4 operator+(const Point3& a, const Point4& b)
{
    Point4 r(b); r += a; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point4 operator-(const Point4& a, const Point3& b)
{
    Point4 r(a); r -= b; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point4 operator-(const Point3& a, const Point4& b)
{
    Point4 r; r = -(b-a); return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point4 operator*(GeomReal s, const Point4& a)
{
    Point4 r(a); r *= s; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point4 operator*(const Point4& a, GeomReal s)
{
    Point4 r(a); r *= s; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point4 operator*(const Point4& p, const VISMatrix4& a)
{	
    // (row vector) * VISMatrix4
    Point4 r(p); r *= a; return r;
}


class VISVector3 : public Point3
{
  public:
    VISVector3() : Point3() {}
    VISVector3(GeomReal x, GeomReal y, GeomReal z) : Point3(x,y,z) {}
    VISVector3(const GeomReal *p) : Point3(p) {}
    VISVector3(const Point3& p) : Point3(p) {}
    VISVector3(const VISVector3& d) { operator=(d); }
    VISVector3(const Point4& d) { operator=(d); }

    // Point transform operators have been defined in geometry.C
    VISVector3& operator*=(const VISMatrix4&);

    VISVector3& transform(const Transform* t);
    VISVector3& transform(const Transform& t);
    VISVector3  cross(const VISVector3& v) const;
    GeomReal dot(const VISVector3& v) const;
    
    VISVector3& operator=(const VISVector3& d) {
	x(d.x()); y(d.y()); z(d.z()); return *this; }
    VISVector3& operator=(const Point4& d) {
	x(d.x()/d.w()); y(d.y()/d.w()); z(d.z()/d.w()); return *this; }
    VISVector3& operator=(const Point3& p) {
	x(p.x()); y(p.y()); z(p.z()); return *this; }
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// VISVector3 operator declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

VISVector3 operator*(const VISMatrix4&,const VISVector3&);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Inline functions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISVector3 operator*(const VISVector3& d,const VISMatrix4& a)
{
    // (row vector) * VISMatrix4
    VISVector3 nd(d); nd *= a; return nd;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline ostream& operator<<(ostream& ostr, const VISVector3& d)
{
    d.print(ostr);
    return ostr;
}


class Direct3 : public Point3
{
  public:
    Direct3() : Point3() {}
    Direct3(GeomReal x, GeomReal y, GeomReal z) : Point3(x,y,z) { normalize(); }
    Direct3(const GeomReal *p) : Point3(p) { normalize(); }
    Direct3(const Point3& p) : Point3(p) { normalize(); }
    Direct3(const Direct3& d) { operator=(d); }

    // Point transform operators have been defined in geometry.C
    Direct3& operator*=(const VISMatrix4&);

    Direct3& transform(const Transform* t);
    Direct3& transform(const Transform& t);

    Direct3& operator=(const Direct3& d) { x(d.x()); y(d.y()); z(d.z()); return *this; }
    Direct3& operator=(const Point3& p) { x(p.x()); y(p.y()); z(p.z());
					 normalize(); return *this; }

    void print(ostream& ostr = cout) const { Point3::print(ostr); }
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Direct3 operator declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Direct3 operator*(const VISMatrix4&,const Direct3&);
/*
Direct3 operator*(const Direct3&,const VISMatrix4&);
*/

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Inline functions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Direct3 operator*(const Direct3& d,const VISMatrix4& a)
{
    // (row vector) * VISMatrix4
    Direct3 nd(d); nd *= a; return nd;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline ostream& operator<<(ostream& ostr, const Direct3& d)
{
    d.print(ostr);
    return ostr;
}


class VISBBox
{
  protected:
    Point3	_pmin,_pmax;
  public:
    VISBBox() { _pmin = Point3(0,0,0); _pmax = Point3(0,0,0); }
    VISBBox(GeomReal xmin,GeomReal xmax,
	    GeomReal ymin,GeomReal ymax,
	    GeomReal zmin,GeomReal zmax) {
	_pmin = Point3(xmin,ymin,zmin);
	_pmax = Point3(xmax,ymax,zmax);
    }
    VISBBox(const VISBBox& b) { operator=(b); }
    VISBBox(const Point3& pmin,const Point3& pmax);

    boolean operator==(const VISBBox&) const;
    boolean operator!=(const VISBBox& b) const {
	return !operator==(b);
    }
    
    GeomReal	xmin()	const	{ return _pmin.x(); }
    GeomReal	xmax()	const	{ return _pmax.x(); }
    GeomReal	ymin()	const	{ return _pmin.y(); }
    GeomReal	ymax()	const	{ return _pmax.y(); }
    GeomReal	zmin()	const	{ return _pmin.z(); }
    GeomReal	zmax()	const	{ return _pmax.z(); }
    const Point3& pmin() const { return _pmin; }
    const Point3& pmax() const { return _pmax; }

    Point3 center() const;

    VISBBox& scale(GeomReal f) { operator*=(f); return *this; }

    VISBBox& merge(const VISBBox&);
    VISBBox& merge(const Point3&);

    VISBBox& transform(const Transform* t);
    VISBBox& transform(const Transform& t);

    VISBBox& operator+=(const VISBBox& bb) { merge(bb); return *this; }
    VISBBox& operator*=(const VISMatrix4&);
    VISBBox& operator*=(GeomReal);

    VISBBox& operator=(const VISBBox&);

    void print(ostream& ostr = cout) const;
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// VISBBox operator declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/*
VISBBox merge(const VISBBox&,const VISBBox&);
VISBBox operator+(const VISBBox& b1,const VISBBox& b2);
VISBBox operator*(const VISBBox& bbox,const Transform& t);
*/

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Inline functions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISBBox merge(const VISBBox& b1,const VISBBox& b2)
{
    VISBBox	result(b1);
    result.merge(b2);
    return result;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISBBox operator+(const VISBBox& b1,const VISBBox& b2)
{
    return merge(b1,b2);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISBBox operator*(const VISBBox& bbox,const VISMatrix4& m)
{
    VISBBox result(bbox);
    result *= m;
    return result;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISBBox operator*(const VISBBox& bbox,GeomReal f)
{
    VISBBox result(bbox);
    result *= f;
    return result;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISBBox operator*(GeomReal f,const VISBBox& bbox)
{
    VISBBox result(bbox);
    result *= f;
    return result;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline ostream& operator<<(ostream& ostr, const VISBBox& bb)
{
    bb.print(ostr);
    return ostr;
}


class Angle3
{
  protected:
    GeomReal	_a[3];
  public:
    Angle3()	{ x(0); y(0); z(0); }
    Angle3(GeomReal xx, GeomReal yy, GeomReal zz) { x(xx);  y(yy);  z(zz); }
    Angle3(const GeomReal *a) 	{ x(a[0]);  y(a[1]);  z(a[2]); }
    Angle3(const Angle3& a)	{ operator=(a); }
    
    // provide access to indivual components.;
    GeomReal	x()	const	{ return _a[0]; }
    GeomReal	y()	const	{ return _a[1]; }
    GeomReal	z()	const	{ return _a[2]; }
    
    void	x(GeomReal x)	{ _a[0] = x; }
    void	y(GeomReal y)	{ _a[1] = y; }
    void	z(GeomReal z)	{ _a[2] = z; }

    boolean operator==(const Angle3& a) const {
	return (a.x() == x() && a.y() == y() && a.z() == z());
    }
    boolean operator!=(const Angle3& a) const {
	return !operator==(a);
    }
    
    // Allow indexing on point objects
    // Range checking is performed when debugging is enabled
    GeomReal&	at(const unsigned int i) {
#ifdef DEBUG
	if (i > 2)
	    ERROR("Range check failure: Angle3: operator[]");
#endif
	return _a[i];
    }
    const GeomReal&	operator[](const unsigned i) const {
#ifdef DEBUG
	if (i > 2)
	    ERROR("Range check failure: Angle3: operator[]");
#endif
	return _a[i];
    }

    Angle3& operator=(const Angle3& a) { x(a.x()); y(a.y()); z(a.z()); return *this; }
    Angle3& operator=(const GeomReal* a){ x(a[0]);  y(a[1]);  z(a[2]); return *this; }
    
    void print(ostream& ostr = cout) const;
};


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/*    Remove this from patchobj.C when you check in geometry.h - DB  */
inline GeomReal operator*(const Point4& p1,const Point4& p2)
{
    GeomReal	f;

    f  = p1.x() * p2.x();
    f += p1.y() * p2.y();
    f += p1.z() * p2.z();
    f += p1.w() * p2.w();

    return f;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline ostream& operator<<(ostream& ostr, const Point4& p)
{
    p.print(ostr);
    return ostr;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline GeomReal length2(const Point4& p)
{
    return ( sqr(p.x()) + sqr(p.y()) + sqr(p.z()) * sqr(p.w()) );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline GeomReal length(const Point4& p)
{
    return sqrt( length2(p) );
}


typedef	GeomReal	M4Rep[4][4];	// Internal rep of a VISMatrix4

class VISMatrix4		// Abstract class representing a transform matrix;
{
  protected:
    M4Rep _m;
  public:
    VISMatrix4() { }
    VISMatrix4(GeomReal a1, GeomReal a2, GeomReal a3, GeomReal a4,
	    GeomReal b1, GeomReal b2, GeomReal b3, GeomReal b4,
	    GeomReal c1, GeomReal c2, GeomReal c3, GeomReal c4,
	    GeomReal d1, GeomReal d2, GeomReal d3, GeomReal d4);
    VISMatrix4(const M4Rep m)    { operator=(m); }
    VISMatrix4(const VISMatrix4& m) { operator=(m); }
    VISMatrix4(const Quaternion& q) { operator=(q); }
    
    VISMatrix4& zero();
    VISMatrix4& ident();

    boolean	operator==(const VISMatrix4&) const;
    boolean	operator!=(const VISMatrix4& m) const {
	return !operator==(m);
    }

    GeomReal& at(const unsigned i,const unsigned j) {
#ifdef DEBUG
	if (i> 3 || j > 3) {
	    ERROR("Range check failure: VISMatrix4: at()");
	}
#endif
	changed();
	return _m[i][j];
    }
    const GeomReal&	itemAt(const unsigned int i,const unsigned j) const {
#ifdef DEBUG
	if (i > 3 || j > 3)
	    ERROR("Range check failure: VISMatrix4: itemAt()");
#endif
	return _m[i][j];
    }
    const GeomReal* operator[](const unsigned i) const {
#ifdef DEBUG
	if (i> 3) {
	    ERROR("Range check failure: VISMatrix4: operator[]");
	}
#endif
	return (GeomReal*)(_m[i]);
    }
    
    // adds; note scales are not defined for VISMatrix4 (this would be
    // a no-op in Homogeneous Coordinates);
    void operator+=(const VISMatrix4&);
    void operator-=(const VISMatrix4&);
    void operator*=(const VISMatrix4&);
    void operator*=(GeomReal);
    
    friend VISMatrix4 inverse(const VISMatrix4&);


    void mulRow(unsigned,GeomReal);
    void swapRows(unsigned,unsigned);
    void subRows(unsigned,unsigned,GeomReal);

    // This is how we can extract values from the matrix
    void assignTo(M4Rep m) const {
	for (unsigned i=0; i<4; ++i)
	    for (unsigned j=0; j<4; ++j)
		m[i][j] = _m[i][j];
    }

    VISMatrix4& operator=(const M4Rep m) {
	memcpy(_m,m,sizeof(GeomReal)*16); return *this; }
    VISMatrix4& operator=(const VISMatrix4& m) {
	memcpy(_m,m._m,sizeof(GeomReal)*16); return *this; }
    VISMatrix4& operator=(const Quaternion& q);

    // set values close to zero to zero
    void print(ostream& ostr = cout) const;

    virtual void changed() {}
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// VISMatrix operator declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

VISMatrix4 transpose(const VISMatrix4&);
/*
VISMatrix4 operator+(const VISMatrix4&,const VISMatrix4&);
VISMatrix4 operator-(const VISMatrix4&,const VISMatrix4&);
VISMatrix4 operator-(const VISMatrix4&);
*/
VISMatrix4 operator*(const VISMatrix4&,const VISMatrix4&);
// VISMatrix4 operator*(const VISMatrix4&,GeomReal);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Inline functions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISMatrix4& VISMatrix4::ident()
{
    this->zero();
    at(0,0) = at(1,1) = at(2,2) = at(3,3) = 1.0;
    return *this;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline void VISMatrix4::operator*=(const VISMatrix4& m)
{
    operator=((*this) * m);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISMatrix4 operator+(const VISMatrix4& a, const VISMatrix4& b)
{
    VISMatrix4 r(a); r += b; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISMatrix4 operator-(const VISMatrix4& a, const VISMatrix4& b)
{
    VISMatrix4 r(a); r -= b; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISMatrix4 operator-(const VISMatrix4& m)
{
    VISMatrix4 r; r.zero(); r -= m; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISMatrix4 operator*(const VISMatrix4& a, GeomReal b)
{
    VISMatrix4 r(a); r *= b; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline ostream& operator<<(ostream& ostr, const VISMatrix4& m)
{
    m.print(ostr);
    return ostr;
}


class Transform : public VISMatrix4
{
  public:
    typedef int FlagField;
    enum { NoneValid =			0x0};
    enum { All =			0x7FFFFFFF};
    enum { Valid =			0x1 };
    enum { Identity =			0x2 };
    enum { Mirrored =			0x4 };

    Transform() { }
    Transform(GeomReal a1, GeomReal a2, GeomReal a3, GeomReal a4,
		 GeomReal b1, GeomReal b2, GeomReal b3, GeomReal b4,
		 GeomReal c1, GeomReal c2, GeomReal c3, GeomReal c4,
		 GeomReal d1, GeomReal d2, GeomReal d3, GeomReal d4) :
	VISMatrix4(a1,a2,a3,a4,
		b1,b2,b3,b4,
		c1,c2,c3,c4,
		d1,d2,d3,d4) { }
    Transform(const M4Rep m)			: VISMatrix4(m) { }
    Transform(const VISMatrix4& m)		: VISMatrix4(m) { }
    Transform(const Transform& t) 	: VISMatrix4(t) { }
    Transform(const Quaternion& q) 		: VISMatrix4(q) { }
    
// this is part of our transition out of this shit RTW 3-13-99    
    Transform(const VISVISMatrix &T);

    Transform& operator|=(const Point3&);

    boolean	operator==(const Transform&) const;
    boolean	operator!=(const Transform& t) const {
	return !operator==(t);
    }

    // Scale operators
    friend Transform scaling(GeomReal s);
    friend Transform scaling(GeomReal x, GeomReal y, GeomReal z);
    Transform& operator&=(const Point3&);
    Transform& operator&=(GeomReal s);

    boolean decompose(Transform& scale,
		      Transform& shear,
		      Transform& rot,
		      Transform& trans,
		      Transform& persp) const;
    boolean decompose(GeomReal& sx,GeomReal& sy,GeomReal& sz,
		      GeomReal& sxy,GeomReal& sxz,GeomReal& syz,
		      GeomReal& rx,GeomReal& ry,GeomReal& rz,
		      GeomReal& tx,GeomReal& ty,GeomReal& tz,
		      GeomReal& px,GeomReal& py,GeomReal& pz,GeomReal& pw) const;
    boolean decompose_rotation(GeomReal&,GeomReal&,GeomReal&) const;
    boolean decompose_translation(GeomReal&,GeomReal&,GeomReal&) const;

    GeomReal det3() const;
    Transform& rotate_only();
    Transform& translate_only();
    
    Transform& operator=(const Transform& m) { 
	VISMatrix4::operator=(m); return *this; }
    Transform& operator=(const VISMatrix4& m) { 
	VISMatrix4::operator=(m); return *this; }
    Transform& operator=(const Quaternion& q) {
	VISMatrix4::operator=(q); return *this; }

    virtual void print(ostream& ostr = cout) const { VISMatrix4::print(ostr); }
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Transform operator declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
// Translate operators;
//Transform operator|(const Transform&, const Point3&);
    
// Scale operators;
//Transform operator&(const Transform&, const Point3&);
//Transform operator&(const Transform&, GeomReal);

// rotations
Transform rotation(const VISCoordIndex, GeomReal);
Transform rotation(const Point3&, GeomReal);
Transform rotation(const Direct3&, const VISCoordIndex = X);
Transform perspective(GeomReal);
Transform rotate_only(const Transform&);
Transform translate_only(const Transform&);
GeomReal det3(const Transform&);

// Translate operator
Transform operator|(const Point3& p, const Transform& t);
Transform translation(const Point3&);
Transform translation (GeomReal x, GeomReal y, GeomReal z);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Inline functions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Transform operator|(const Transform& t, const Point3& p)
{
    Transform r(t); r |= p; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Transform scaling(GeomReal s)
{
    return scaling(s, s, s);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Transform operator&(const Transform& t, const Point3& p)
{
    Transform r(t); r &= p; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Transform& Transform::operator&=(GeomReal s)
{
    (*this) &= Point3(s,s,s); return *this;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Transform operator& (const Transform& t, GeomReal s)
{
    Transform r(t); r &= s; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Transform operator&(GeomReal s, const Transform& t)
{
    Point3 n(s,s,s); Transform r = n & t; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline ostream& operator<<(ostream& ostr, const Transform& t)
{
    t.print(ostr); return ostr;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point3& Point3::transform(const Transform& t)
{
    (*this) *= t; return *this;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Direct3& Direct3::transform(const Transform& t)
{
    (*this) *= t; return *this;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline VISVector3& VISVector3::transform(const Transform& t)
{
    (*this) *= t; return *this;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Point4& Point4::transform(const Transform& t)
{
    (*this) *= t; return *this;
}


class Quaternion
{
  protected:
    GeomReal _d[4];
  public:
    Quaternion() {_d[0] = _d[1] = _d[2] = 0; _d[3] = 1;}
    Quaternion(GeomReal x, GeomReal y, GeomReal z, GeomReal w) {
	_d[0] = x; _d[1] = y; _d[2] = z; _d[3] = w;}
    Quaternion(const GeomReal *q) {
	_d[0] = q[0]; _d[1] = q[1];
	_d[2] = q[2]; _d[3] = q[3];
    }
    Quaternion(const Point3& p) { operator=(p); }	// direction, magnitude
    Quaternion(const Quaternion& q) { operator=(q); }
    Quaternion(const Angle3& a) { operator=(a); }
    Quaternion(const Transform& m) { operator=(m); }

    inline boolean operator==(const Quaternion&) const;
    boolean operator!=(const Quaternion& q) const {
	return !operator==(q);
    }
    
    GeomReal x() const {return _d[0];}
    GeomReal y() const {return _d[1];}
    GeomReal z() const {return _d[2];}
    GeomReal w() const {return _d[3];}
    
    void x(GeomReal x) {_d[0] = x;}
    void y(GeomReal y) {_d[1] = y;}
    void z(GeomReal z) {_d[2] = z;}
    void w(GeomReal w) {_d[3] = w;}

    GeomReal&	at(const unsigned int i) {
#ifdef DEBUG
	if (i > 3)
	    VISHandler::panic ("Range check failure",
			       "Quaternion",
			       "at()");
#endif
	return _d[i];
    }
    const GeomReal&	itemAt(const unsigned int i) const {
#ifdef DEBUG
	if (i > 3)
	    VISHandler::panic("Range check failure",
			      "Quaternion",
			      "itemAt()");
#endif
	return _d[i];
    }
    const GeomReal& operator[] (const unsigned int i) const {
#ifdef DEBUG
	if (i > 3)
	    VISHandler::panic ("Range check failure",
			       "Quaternion",
			       "operator[]");
#endif
	return _d[i];
    }

    Quaternion& operator+=(const Quaternion&);
    Quaternion& operator-=(const Quaternion&);
    Quaternion& operator*=(const Quaternion&);
    Quaternion& operator/=(const Quaternion&);

    Quaternion& operator=(const Quaternion& q) { x(q.x()); y(q.y());
					         z(q.z()); w(q.w());
					         return *this; }
    Quaternion& operator=(const Angle3&);
    Quaternion& operator=(const Transform&);
    Quaternion& operator=(const Point3&);
    Quaternion normal() const;		// Returns the normal of this
    void normalize();

    void print(ostream& ostr = cout) const;
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Quaternion operator declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//Quaternion operator+(const Quaternion&,const Quaternion&);
//Quaternion operator-(const Quaternion&,const Quaternion&);
//Quaternion operator*(const Quaternion&,const Quaternion&);
Quaternion inverse(const Quaternion&);
//Quaternion operator/(const Quaternion&,const Quaternion&);
Quaternion quaternion(const Point3&,const GeomReal&);
//GeomReal length(const Quaternion&);
//GeomReal length2(const Quaternion&);
Transform rotation(const Quaternion&);
//Quaternion normal(const Quaternion&);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Inline functions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline boolean	Quaternion::operator==(const Quaternion& q) const {
    return (x() == q.x() && y() == q.y() && z() == q.z() && w() == q.w());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline GeomReal length2(const Quaternion& q)
{
    return ( sqr(q.x()) + sqr(q.y()) + sqr(q.z()) + sqr(q.w()) );
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline GeomReal length(const Quaternion& q)
{
    return sqrt(length2(q));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Quaternion operator+(const Quaternion& a, const Quaternion& b)
{
    Quaternion r(a); r += b; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Quaternion operator-(const Quaternion& a, const Quaternion& b)
{
    Quaternion r(a); r -= b; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Quaternion operator*(const Quaternion& a, const Quaternion& b)
{
    Quaternion r(a); r *= b; return r;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Quaternion& Quaternion::operator/=(const Quaternion& q)
{
    *this *= inverse(q);
    return *this;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Quaternion operator/(const Quaternion& a, const Quaternion& b)
{
    return (a * inverse(b));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline Quaternion normal(const Quaternion& q)
{
    return q.normal();
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline ostream& operator<<(ostream& ostr, const Quaternion& q)
{
    q.print(ostr);
    return(ostr);
}





// ---------------------------------------------------------------------------

#endif
