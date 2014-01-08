
// this is to help transition out of this shit RTW 3-13-99
#include "matrix/matrix.h"
#include "util/geometry.h"
#include "util/mathutil.h"



// ---------------------------------------------------------------------------
// Point2
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

void Point2::print(ostream& ostr) const
{
    // { X, Y };
    ostr << "[ " << form("%6g", x())
	 << ", " << form("%6g", y())
	 << " ]";
}

// ---------------------------------------------------------------------------

Point2 normal(const Point2& p)
{
    return p.normal();
}

// ---------------------------------------------------------------------------

Point2 Point2::normal() const
{
    Point2 r = *this;
    GeomReal len = sqrt(r.x()*r.x() + r.y()*r.y());

    if (len == 0.0)
    {
	r.x(1);
	r.y(0);
    }
    else
    {
	len = 1.0/len;
	r.x( r.x() * len );
	r.y( r.y() * len );
    }
    return r;
}

// ---------------------------------------------------------------------------

void Point2::normalize()
{
    GeomReal len = sqrt(sqr(x()) + sqr(y()));

    if (len == 0.0)	
    {
	x(1.0);
	y(0.0);
    }
    else
    {
	len = 1.0/len;
	x( x() * len );
	y( y() * len );
    }
}

// ---------------------------------------------------------------------------

Point2 Point2::dominant() const
{
    GeomReal xa,ya;
    xa = fabs(x());
    ya = fabs(y());

    if (ya > xa)
	return Point2(0,y());
    else
	return Point2(x(),0);
}

// ---------------------------------------------------------------------------

Point2 dominant(Point2& p)
{
    return p.dominant();
}

// ---------------------------------------------------------------------------
// Point3
// ---------------------------------------------------------------------------

Point3&	Point3::operator=(const Point4& p)
{
    if (p.w() != 0.0) {
	x( p.x()/p.w() );
	y( p.y()/p.w() );
	z( p.z()/p.w() );
    }
    else {	// What else can you do?
	x( 0.0 );
	y( 0.0 );
	z( 0.0 );
    }

    return *this;
}

// ---------------------------------------------------------------------------


void Point3::print(ostream& ostr) const
{
    // { X, Y, Z };
    ostr << "[ " << form("%6g", x())
	 << ", " << form("%6g", y())
	 << ", " << form("%6g", z())
	 << " ]";
}

// ---------------------------------------------------------------------------

Point3 normal(const Point3& p)
{
    return p.normal();
}

// ---------------------------------------------------------------------------

Point3 Point3::normal() const
{
    Point3 r = *this;
    GeomReal len = sqrt(r.x()*r.x() + r.y()*r.y() + r.z()*r.z());

    if (len == 0.0)
    {
	r.x(1);
	r.y(0);
	r.z(0);
    }
    else
    {
	len = 1.0/len;
	r.x( r.x() * len );
	r.y( r.y() * len );
	r.z( r.z() * len );
    }
    return r;
}

// ---------------------------------------------------------------------------

Point3 Point3::dominant() const
{
    GeomReal xa,ya,za;
    xa = fabs(x());
    ya = fabs(y());
    za = fabs(z());

    if (za > xa && za > ya)
	return Point3(0,0,z());
    else if (ya > xa && ya > za)
	return Point3(0,y(),0);
    else
	return Point3(x(),0,0);
}

// ---------------------------------------------------------------------------

Point3 dominant(Point3& p)
{
    return p.dominant();
}

// ---------------------------------------------------------------------------
// Point4
// ---------------------------------------------------------------------------

void Point4::homogenize()
{
    if (w() != 0.0) {
	x( x()/w() );
	y( y()/w() );
	z( z()/w() );
	w( 1.0 );
    }
}

// ---------------------------------------------------------------------------

Point4& Point4::operator+=(const Point3& p)
{
    _d[0] += p.x() * w();
    _d[1] += p.y() * w();
    _d[2] += p.z() * w();
    return *this;
}

// ---------------------------------------------------------------------------

Point4& Point4::operator-=(const Point3& p)
{
    _d[0] -= p.x() * w();
    _d[1] -= p.y() * w();
    _d[2] -= p.z() * w();
    return *this;
}

// ---------------------------------------------------------------------------

void Point4::print(ostream& ostr) const
{
    // { X, Y, Z, W };
    ostr<< "{ " << form("%6g", x())
	<< ", " << form("%6g", y())
	<< ", " << form("%6g", z())
	<< ", " << form("%6g", w())
	<< " }";
}

// ---------------------------------------------------------------------------
// Quaternion
// ---------------------------------------------------------------------------
// Notes: There are certain inconsistencies here when compared to the
// implementations found in some of the books: Graphics Gems and the Shoemake
// paper.  The problem is when a quaternion is applied to a point as a rotation,
// then the system handedness comes in questions.  We apply points on the left-
// hand-side to convert them, which is opposite from other systems, which do
// it on the right.  Thus a quaternion multiplication, which represents the
// concatenation of two transforms, needs to be performed backwards so that
// the "order" is preserved for a left-handed system.
//
// ---------------------------------------------------------------------------
// This builds a quaternion from a vector p, whose magnitude is the
// angle in radians and whose direction is axis of rotation.
// Normalization should be done here, becuase p may not be normalized.

Quaternion& Quaternion::operator=(const Point3& p)
{
    GeomReal	angle = length(p);
    GeomReal	c = cos(angle/2.0);
    Point3	q = sin(angle/2.0) * p.normal();

    x( q.x() );
    y( q.y() );
    z( q.z() );
    w( c );
    return *this;
}

// ---------------------------------------------------------------------------
// Calculates a quaternion from three Euler angles
// Based on Shoemake's SIGGRAPH '85 paper. - DB

Quaternion& Quaternion::operator=(const Angle3& a)
{
    GeomReal sinax, cosax, sinay, cosay, sinaz, cosaz;
    
    sinax = sin(a.x()/2);
    cosax = cos(a.x()/2);
    sinay = sin(a.y()/2);
    cosay = cos(a.y()/2);
    sinaz = sin(a.z()/2);
    cosaz = cos(a.z()/2);

    _d[0] = sinax * cosay * cosaz - cosax * sinay * sinaz;
    _d[1] = cosax * sinay * cosaz + sinax * cosay * sinaz;
    _d[2] = cosax * cosay * sinaz - sinax * sinay * cosaz;
    _d[3] = cosax * cosay * cosaz + sinax * sinay * sinaz;

    return *this;
}

// ---------------------------------------------------------------------------

Quaternion& Quaternion::operator=(const Transform& m)
{
    GeomReal sx,sy,sz;
    GeomReal sxy,sxz,syz;
    GeomReal rx,ry,rz;
    GeomReal tx,ty,tz;
    GeomReal px,py,pz,pw;

    if (m.decompose(sx,sy,sz,sxy,sxz,syz,rx,ry,rz,tx,ty,tz,px,py,pz,pw)) {
	operator=(Angle3(rx,ry,rz));
    }

    return *this;
}

// ---------------------------------------------------------------------------

Quaternion& Quaternion::operator+=(const Quaternion& q)
{
    _d[0] += q.x();
    _d[1] += q.y();
    _d[2] += q.z();
    _d[3] += q.w();
    return *this;
}

// ---------------------------------------------------------------------------

Quaternion& Quaternion::operator-=(const Quaternion& q)
{
    _d[0] -= q.x();
    _d[1] -= q.y();
    _d[2] -= q.z();
    _d[3] -= q.w();
    return *this;
}

// ---------------------------------------------------------------------------
// Note: this operation here involves an implicit order-swap from that
// in the book because we pre-multiply points to do rotations, while
// the Graphics Gems and apparently the others post-multiply.  So
// we need to implicitly swap the two values.

Quaternion& Quaternion::operator*=(const Quaternion& q)
{
    GeomReal xx,yy,zz,ww;

    xx =   x()*q.w() - y()*q.z() + z()*q.y() + w()*q.x();
    yy =   x()*q.z() + y()*q.w() - z()*q.x() + w()*q.y();
    zz = - x()*q.y() + y()*q.x() + z()*q.w() + w()*q.z();
    ww = - x()*q.x() - y()*q.y() - z()*q.z() + w()*q.w();
    /*
    xx =   w()*q.x() - z()*q.y() + y()*q.z() + x()*q.w();
    yy =   z()*q.x() + w()*q.y() - x()*q.z() + y()*q.w();
    zz = - y()*q.x() + x()*q.y() + w()*q.z() + z()*q.w();
    ww = - x()*q.x() - y()*q.y() - z()*q.z() + w()*q.w();
    */
    x( xx );  y( yy );  z( zz );  w( ww );
    return *this;
}

// ---------------------------------------------------------------------------

void Quaternion::normalize()
{
    GeomReal len = length(*this);
    if (len <= 0.0) {
	x(0);
	y(0);
	z(0);
	w(1);
    }
    else {
	len = 1.0/len;
	x( x() * len );
	y( y() * len );
	z( z() * len );
	w( w() * len );
    }
}

// ---------------------------------------------------------------------------

Quaternion inverse(const Quaternion& q)
{
    GeomReal f = 1.0 / length2(q);

    Quaternion r(-q.x()*f, -q.y()*f, -q.z()*f, q.w()*f);
    return r;
}

// ---------------------------------------------------------------------------
// I made this a member function, but the old routine is still accessible

Quaternion Quaternion::normal() const
{
    Quaternion r(*this);
    r.normalize();
    return r;
}

// ---------------------------------------------------------------------------

Quaternion quaternion(const Point3& v, const GeomReal& theta)
{
    Quaternion q;
    GeomReal len = length(v);

    if (len == 0.0)
    {
	q.x(0);
	q.y(0);
	q.z(0);
	q.w(1);
    }
    else
    {
	GeomReal sl = sin(theta/2.)/len;
	q.x(v.x()*sl);
	q.y(v.y()*sl);
	q.z(v.z()*sl);
	q.w(cos(theta/2.));
    }
    return q;
}

// ---------------------------------------------------------------------------

void Quaternion::print(ostream& ostr) const
{
    // { X, Y, Z, W };
    ostr << "{ " << form("%6g", x())
	 << ", " << form("%6g", y())
	 << ", " << form("%6g", z())
	 << ", " << form("%6g", w())
	 << " }";
}

// ---------------------------------------------------------------------------
// Point3
// ---------------------------------------------------------------------------

// Does a direction transform (row vector) * Matrix4

void Point3::normalize()
{
    GeomReal len = sqrt(sqr(x()) + sqr(y()) + sqr(z()));

    if (len == 0.0)	
    {
	x(1.0);
	y(0.0);
	z(0.0);
    }
    else
    {
	len = 1.0/len;
	x( x() * len );
	y( y() * len );
	z( z() * len );
    }
}

// ---------------------------------------------------------------------------
// These are in here to keep C++ from trying other conversions
// that screw up because there are too many related classes (Point4, Point3,
// Direct3, etc).  If I don't explicitly give a routine C++ might do something
// wrong or at least complain about it.		-erose

Point3& Point3::operator*=(const Matrix4& m)
{
    Point4 np = *this;
    *this = (np * m);
    return *this;
}

// ---------------------------------------------------------------------------
// This is a different order than the book because we have a pre-multiply
// system instead of a post-multiply system

Point3& Point3::operator*=(const Quaternion& q)
{
    Quaternion qp(x(),y(),z(),0);
    Quaternion r = inverse(q) * qp * q;
    operator=(Point3(r.x(),r.y(),r.z()));
    return *this;
}

// ---------------------------------------------------------------------------

Point3 operator*(const Point3& p,const Matrix4& m)
{
    // (row vector) * Matrix4
    Point4 np = p; np *= m;		// Convert into Point4, *
    Point3 np3 = np; return np3;	// Convert back
}

// ---------------------------------------------------------------------------

Point3 operator*(const Matrix4& m,const Point3& p)
{
    Point4 p4(p);
    Point4 np = m * p4;
    Point3 np3 = np; return np3;
}

// ---------------------------------------------------------------------------

Point3 operator*(const Point3& p,const Quaternion& q)
{
    Point3 result(p);
    result *= q;
    return result;
}

// ---------------------------------------------------------------------------
// This is the correct may to apply rotations to points -- apply the point
// first -- this is a different order than the book because we have a pre-multiply
// system instead of a post-multiply system


Point3 operator*(const Quaternion& q,const Point3& p)
{
    Quaternion qp(p.x(),p.y(),p.z(),0);
    Quaternion r = q * qp * inverse(q);
    return Point3(r.x(),r.y(),r.z());
}

// ---------------------------------------------------------------------------

void Point3::scale(GeomReal f)
{
    operator*=(f);
}

// ---------------------------------------------------------------------------

void Point3::scale(GeomReal fx,GeomReal fy,GeomReal fz)
{
    x( x() * fx );
    y( y() * fy );
    z( z() * fz );
}

// ---------------------------------------------------------------------------
// This applies a transform to a point the correct way...

Point3& Point3::transform(const Transform* t)
{
#ifdef DEBUG
    if (t == NULL) {
	ERROR("NULL transform: geometry: Point3::transform(Transform*)");
	return *this;
    }
    else
#endif
	return this->transform(*t);
}

// ---------------------------------------------------------------------------
// Vector3
// ---------------------------------------------------------------------------

Vector3& Vector3::operator*=(const Matrix4& m)
{
    // (row vector) * Matrix4
    Vector3 d(*this);

    for (unsigned i=0; i<3; i++) {
	at(i) = (m[0][i] * d[0] +
		 m[1][i] * d[1] +
		 m[2][i] * d[2]);
    }

    return *this;
}

// ---------------------------------------------------------------------------

Vector3 operator*(const Matrix4& m,const Vector3& d)
{
    // Matrix4 * (column vector);
    Vector3 result;

    for (unsigned i=0; i<4; i++) {
	result.at(i) = (m[i][0] * d[0] +
			m[i][1] * d[1] +
			m[i][2] * d[2]);
    }
    return result;
}

// standard cross product

Vector3 Vector3::cross(const Vector3& v) const
{
    Vector3 result;

    result.at(0) = (_d[1] * v[2] -
		    _d[2] * v[1]);

    result.at(1) = (_d[2] * v[0] -
		    _d[0] * v[2]);

    result.at(2) = (_d[1] * v[2] -
		    _d[2] * v[1]);
    
    return result;
}

GeomReal Vector3::dot(const Vector3& v) const
{
    float total = 0;
    for (int i = 0; i < 3; i++)
	{
	    total += v[i]*_d[i];
	}
    
    return total;
}

// ---------------------------------------------------------------------------

Vector3& Vector3::transform(const Transform* t)
{
#ifdef DEBUG
    if (t == NULL) {
	ERROR("NULL transform geometry Vector3::transform(const Transform*)");
	return *this;
    }
    else
#endif
	return this->transform(*t);
}

// ---------------------------------------------------------------------------
// Direct3
// ---------------------------------------------------------------------------

Direct3& Direct3::operator*=(const Matrix4& m)
{
    // (row vector) * Matrix4
    Direct3 d(*this);

    for (unsigned i=0; i<3; i++) {
	at(i) = (m[0][i] * d[0] +
		 m[1][i] * d[1] +
		 m[2][i] * d[2]);
    }

    normalize();
    return *this;
}

// ---------------------------------------------------------------------------

Direct3 operator*(const Matrix4& m,const Direct3& d)
{
    // Matrix4 * (column vector);
    Direct3 result;

    for (unsigned i=0; i<4; i++) {
	result.at(i) = (m[i][0] * d[0] +
			m[i][1] * d[1] +
			m[i][2] * d[2]);
    }
    result.normalize();
    return result;
}

// ---------------------------------------------------------------------------

Direct3& Direct3::transform(const Transform* t)
{
#ifdef DEBUG
    if (t == NULL) {
	ERROR("NULL transform: geometry: Direct3::transform(const Transform*)");
	return *this;
    }
    else
#endif
	return this->transform(*t);
}

// ---------------------------------------------------------------------------
// GfxBBox
// ---------------------------------------------------------------------------

IrisBBox::IrisBBox(const Point3& pmin,const Point3& pmax)
{
    _pmin = pmin;
    _pmax = pmin;
    merge(pmax);
}

// ---------------------------------------------------------------------------

boolean	IrisBBox::operator==(const IrisBBox& b) const
{
    return (_pmin == b.pmin() && _pmax == b.pmax());
}

// ---------------------------------------------------------------------------

IrisBBox& IrisBBox::merge(const IrisBBox& b)
{
    if (b.xmin() < xmin()) _pmin.x(b.xmin());
    if (b.xmax() > xmax()) _pmax.x(b.xmax());
    if (b.ymin() < ymin()) _pmin.y(b.ymin());
    if (b.ymax() > ymax()) _pmax.y(b.ymax());
    if (b.zmin() < zmin()) _pmin.z(b.zmin());
    if (b.zmax() > zmax()) _pmax.z(b.zmax());

    return *this;
}

// ---------------------------------------------------------------------------

IrisBBox& IrisBBox::merge(const Point3& p)
{
    if (p.x() < xmin()) _pmin.x(p.x());
    if (p.x() > xmax()) _pmax.x(p.x());
    if (p.y() < ymin()) _pmin.y(p.y());
    if (p.y() > ymax()) _pmax.y(p.y());
    if (p.z() < zmin()) _pmin.z(p.z());
    if (p.z() > zmax()) _pmax.z(p.z());

    return *this;
}

// ---------------------------------------------------------------------------

IrisBBox& IrisBBox::operator=(const IrisBBox& b)
{
    _pmin = b.pmin();
    _pmax = b.pmax();

    return *this;
}

// ---------------------------------------------------------------------------

IrisBBox& IrisBBox::operator*=(const Matrix4& m)
{
    unsigned	i;
    Point3	corners[8];

    corners[0] = Point3(xmin(),ymin(),zmin());
    corners[1] = Point3(xmax(),ymin(),zmin());
    corners[2] = Point3(xmin(),ymax(),zmin());
    corners[3] = Point3(xmax(),ymax(),zmin());
    corners[4] = Point3(xmin(),ymin(),zmax());
    corners[5] = Point3(xmax(),ymin(),zmax());
    corners[6] = Point3(xmin(),ymax(),zmax());
    corners[7] = Point3(xmax(),ymax(),zmax());

    for (i=0; i<8; ++i) {
	corners[i] *= m;
    }

    // Construct the new bbox based on this one
    IrisBBox	result(corners[0],corners[0]);

    for (i=1; i<8; ++i) {
	result.merge(corners[i]);
    }

    operator=(result);
    return *this;
}

// ---------------------------------------------------------------------------

Point3 IrisBBox::center() const
{
    return (_pmin + _pmax)/2.0;
}

// ---------------------------------------------------------------------------

IrisBBox& IrisBBox::operator*=(GeomReal f)
{
//    GeomReal dx,dy,dz;
//    dx = (xmax() - xmin()) * f;
//    dy = (ymax() - ymin()) * f;
//    dz = (zmax() - zmin()) * f;
    
    _pmin.x( _pmin.x() - f);
    _pmax.x( _pmax.x() + f);
    _pmin.y( _pmin.y() - f);
    _pmax.y( _pmax.y() + f);
    _pmin.z( _pmin.z() - f);
    _pmax.z( _pmax.z() + f);
    
    return *this;
}

// ---------------------------------------------------------------------------

IrisBBox& IrisBBox::transform(const Transform& t)
{
    unsigned	i;
    Point3	corners[8];

    corners[0] = Point3(xmin(),ymin(),zmin());
    corners[1] = Point3(xmax(),ymin(),zmin());
    corners[2] = Point3(xmin(),ymax(),zmin());
    corners[3] = Point3(xmax(),ymax(),zmin());
    corners[4] = Point3(xmin(),ymin(),zmax());
    corners[5] = Point3(xmax(),ymin(),zmax());
    corners[6] = Point3(xmin(),ymax(),zmax());
    corners[7] = Point3(xmax(),ymax(),zmax());

    for (i=0; i<8; ++i) {
	corners[i].transform(t);
    }

    // Construct the new bbox based on this one
    IrisBBox	result(corners[0],corners[0]);

    for (i=1; i<8; ++i) {
	result.merge(corners[i]);
    }

    operator=(result);
    return *this;
}

// ---------------------------------------------------------------------------

IrisBBox& IrisBBox::transform(const Transform* t)
{
#ifdef DEBUG
    if (t == NULL) {
	ERROR("NULL transform: geometry: IrisBBox::transform(Transform*)");
	return *this;
    }
    else
#endif
	return this->transform(*t);
}

// ---------------------------------------------------------------------------

void IrisBBox::print(ostream& ostr) const
{
    ostr << "["
	 << xmin() << " " << xmax() << " "
	 << ymin() << " " << ymax() << " "
	 << zmin() << " " << zmax() << "]";
}


// ---------------------------------------------------------------------------
// Angle3
// ---------------------------------------------------------------------------


void Angle3::print(ostream& ostr) const
{
    ostr << "[" << form("%6g", x())
	 << ", " << form("%6g", y())
	 << ", " << form("%6g", z())
	 << "]";
}

// ---------------------------------------------------------------------------
// Point4
// ---------------------------------------------------------------------------

boolean	Point4::operator== (const Point4& p) const
{
    if (p.w() == 0 && w() == 0)
	// Both points are at infinity and conisidered equivalent
	return TRUE;
    else if ((p.w() == 0 && w() != 0) || (p.w() != 0 && w() == 0))
	// One of the points is at infinity but the other isn't
	return FALSE;

    // Compare the normalised values of each point
    return (Point3(*this) == Point3(p));    
}

// ---------------------------------------------------------------------------

boolean	Point4::operator== (const Point3& p) const
{
    if (w() == 0)
	// This is at infinity which a point3 cannot be so...
	return FALSE;

    // Compare the normalised values of each point
    return (Point3(*this) == p);
}

// ---------------------------------------------------------------------------

Point4& Point4::operator*=(const Matrix4& m)
{
    // (row vector) * Matrix4
    Point4 p(*this);

    for (unsigned i=0; i<4; i++) {
	this->at(i) = (m[0][i] * p[0] +
		       m[1][i] * p[1] +
		       m[2][i] * p[2] +
		       m[3][i] * p[3]);
    }

    return *this;
}

// ---------------------------------------------------------------------------

Point4 operator*(const Matrix4& m, const Point4& p)
{
    // Matrix4 * (column vector);
    Point4 result;

    for (unsigned i=0; i<4; i++) {
	result.at(i) = (m[i][0] * p[0] +
			m[i][1] * p[1] +
			m[i][2] * p[2] +
			m[i][3] * p[3]);
    }
    return result;
}

// ---------------------------------------------------------------------------
// This applies a transform to a point the correct way...

Point4& Point4::transform(const Transform* t)
{
#ifdef DEBUG
    if (t == NULL) {
	ERROR("NULL transform: geometry: Point4::transform(const Transform*)");
	return *this;
    }
    else
#endif
	return this->transform(*t);
}

// ---------------------------------------------------------------------------
// Matrix4
// ---------------------------------------------------------------------------

Matrix4::Matrix4 (GeomReal a1, GeomReal a2, GeomReal a3, GeomReal a4,
		  GeomReal b1, GeomReal b2, GeomReal b3, GeomReal b4,
		  GeomReal c1, GeomReal c2, GeomReal c3, GeomReal c4,
		  GeomReal d1, GeomReal d2, GeomReal d3, GeomReal d4)
{
    _m[0][0] = a1; _m[0][1] = a2; _m[0][2] = a3; _m[0][3] = a4; 
    _m[1][0] = b1; _m[1][1] = b2; _m[1][2] = b3; _m[1][3] = b4; 
    _m[2][0] = c1; _m[2][1] = c2; _m[2][2] = c3; _m[2][3] = c4; 
    _m[3][0] = d1; _m[3][1] = d2; _m[3][2] = d3; _m[3][3] = d4; 
}

// ---------------------------------------------------------------------------

boolean Matrix4::operator==(const Matrix4& m) const
{
    for (int i = 0; i < 4; i++)
    {
	if (m.itemAt(0,i) != itemAt(0,i) ||
	    m.itemAt(1,i) != itemAt(1,i) ||
	    m.itemAt(2,i) != itemAt(2,i) ||
	    m.itemAt(3,i) != itemAt(3,i))
	    return FALSE;
    }
    return TRUE;
    
}

// ---------------------------------------------------------------------------

Matrix4	transpose(const Matrix4& a)
{
    Matrix4 r;

    for (int i = 0; i < 4; i++)
    {
	r.at(0,i) = a[i][0];
	r.at(1,i) = a[i][1];
	r.at(2,i) = a[i][2];
	r.at(3,i) = a[i][3];
    }
    return r;
}

// ---------------------------------------------------------------------------

// Note - merely clearing the memory with some sort of "set" command is
// 	  actually not portable - since other GeomReal-point representations
//	  might not have 0.0 == (all bits zero).  Its true for IEEE, though.

Matrix4& Matrix4::zero()
{
    for (unsigned i=0; i<4; ++i) {
	for (unsigned j=0; j<4; ++j) {
	    at(i,j) = 0.0;
	}
    }
    return *this;
}

// ---------------------------------------------------------------------------

void Matrix4::operator+=(const Matrix4& a)
{
    for (int i = 0; i < 4; i++)
    {
	at(i,0) += a[i][0];
	at(i,1) += a[i][1];
	at(i,2) += a[i][2];
	at(i,3) += a[i][3];
    }
}

// ---------------------------------------------------------------------------

void Matrix4::operator-=(const Matrix4& a)
{
    for (int i = 0; i < 4; i++)
    {
	at(i,0) -= a[i][0];
	at(i,1) -= a[i][1];
	at(i,2) -= a[i][2];
	at(i,3) -= a[i][3];
    }
}

// ---------------------------------------------------------------------------

Matrix4 operator*(const Matrix4& a, const Matrix4& b)
{
    Matrix4	result;

    for (register unsigned i = 0; i < 4; i++)
 	for (register unsigned j = 0; j < 4; j++)
	    result.at(i,j) =
		(a[i][0] * b[0][j]) +
		(a[i][1] * b[1][j]) +
		(a[i][2] * b[2][j]) +
		(a[i][3] * b[3][j]);

    return result;
}

// ---------------------------------------------------------------------------

void Matrix4::operator*=(GeomReal s)
{
    for (register unsigned i = 0; i < 4; i++)
 	for (register unsigned j = 0; j < 4; j++)
	    at(i,j) *= s;
}

// ---------------------------------------------------------------------------

Matrix4& Matrix4::operator=(const Quaternion& q)
{
    operator=(rotation(q));
    return *this;
}

// ---------------------------------------------------------------------------
// Various routines used for matrix inverstion
// ---------------------------------------------------------------------------

inline void swap(GeomReal& x1,GeomReal& x2)
{
    GeomReal temp = x1;
    x1 = x2;
    x2 = temp;
}

// ---------------------------------------------------------------------------

void Matrix4::mulRow(unsigned row,GeomReal fac)
{
    for (unsigned j=0; j<4; ++j) {
	_m[row][j] *= fac;
    }
}

// ---------------------------------------------------------------------------

void Matrix4::swapRows(unsigned r1,unsigned r2)
{
    for (unsigned j=0; j<4; ++j) {
	::swap(_m[r1][j],_m[r2][j]);
    }
}

// ---------------------------------------------------------------------------

void Matrix4::subRows(unsigned r1,unsigned r2,GeomReal fac)
{
    for (unsigned j=0; j<4; ++j) {
	_m[r2][j] -= _m[r1][j] * fac;
    }
}

// ---------------------------------------------------------------------------

Matrix4	inverse(const Matrix4& a)
{
    float det,f,cmx;
    int r,c,wr;
    Matrix4 mat(a),imt;

    imt.ident();
    det = 1.0;

    for (c=0; c<4; c++) {
	f = 0.0;
	cmx = fabs(mat[c][c]);
	wr = c;
	for (r=c+1; r<4; r++) {
	    if(fabs(mat[r][c]) > cmx) {
		cmx = fabs(mat[r][c]);
		wr = r;
	    }
	}
	if (mat[wr][c] == 0.0) {
	    // Singular matrix - sorry
	    ERROR("Singular matrix inverted!: geom :inverse()");
	    imt.ident();
	    return(imt);
	}
	if (wr != c) {
	    mat.swapRows(wr,c);
	    imt.swapRows(wr,c);
	    det = (-det);
	}
	f = 1.0/mat[c][c];
	det /= f;
	mat.mulRow(c,f);
	imt.mulRow(c,f);
	for(r=c+1; r<4; r++) {
	    f = mat[r][c];
	    mat.subRows(c,r,f);
	    imt.subRows(c,r,f);
	}
    }
    for (c=3; c>=0; c--) {
	for (r=(c-1); r>=0; r--) {
	    f = mat[r][c];
	    mat.subRows(c,r,f);
	    imt.subRows(c,r,f);
	}
    }

    return imt;
}

// ----------------------------------------------------------------

void Matrix4::print(ostream& ostr) const
{
    ostr << "\n";
    for (int i = 0; i < 4; i++)
	ostr << "\t{ " << form("%6g", _m[i][0])
	     << ", " << form("%6g", _m[i][1])
	     << ", " << form("%6g", _m[i][2])
	     << ", " << form("%6g", _m[i][3])
	     << " }\n";
}

// ---------------------------------------------------------------------------
// Transform
// ---------------------------------------------------------------------------

boolean Transform::operator==(const Transform& m) const
{
    for (int i = 0; i < 4; i++)
    {
	if (m.itemAt(0,i) != itemAt(0,i) ||
	    m.itemAt(1,i) != itemAt(1,i) ||
	    m.itemAt(2,i) != itemAt(2,i) ||
	    m.itemAt(3,i) != itemAt(3,i))
	    return FALSE;
    }
    return TRUE;
    
}


Transform::Transform(const IrisMatrix& m)
{
    int i, j;
    if ((m.width() == 4)&&(m.height() == 4))
	{
	    for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		    at(j, i) = m.itemAt(j, i);
	}
    else
	cout << "matrix to transform get bad matrix" << m;
}    



// ---------------------------------------------------------------------------
// determinant of 3x3 rotation/scale matrix

GeomReal Transform::det3() const
{
    return (itemAt(0,0) * itemAt(1,1) * itemAt(2,2) +
	    itemAt(0,1) * itemAt(1,2) * itemAt(2,0) +
	    itemAt(0,2) * itemAt(1,0) * itemAt(2,1) -
	    itemAt(0,2) * itemAt(1,1) * itemAt(2,0) -
	    itemAt(0,0) * itemAt(1,2) * itemAt(2,1) -
	    itemAt(0,1) * itemAt(1,0) * itemAt(2,2));
}

// ---------------------------------------------------------------------------

GeomReal det3(const Transform& t)
{
    return t.det3();
}

// ---------------------------------------------------------------------------

Transform& Transform::rotate_only()
{
    at(3,0) = 0.0;
    at(3,1) = 0.0;
    at(3,2) = 0.0;
    at(0,3) = 0.0;
    at(1,3) = 0.0;
    at(2,3) = 0.0;
    at(3,3) = 1.0;
    return *this;
}

// ---------------------------------------------------------------------------

Transform& Transform::translate_only()
{
    for (int i=0; i<3; ++i)
	for (int j=0; j<4; ++j)
	    at(i,j) = 0.0;
    at(0,0) = 1.0;
    at(1,1) = 1.0;
    at(2,2) = 1.0;
    at(3,3) = 1.0;
    return *this;
}

// ---------------------------------------------------------------------------

Transform rotate_only(const Transform& t)
{
    Transform result = t;
    result.rotate_only();
    return result;
}

// ---------------------------------------------------------------------------

Transform translate_only(const Transform& t)
{
    Transform result = t;
    result.translate_only();
    return result;
}

// ---------------------------------------------------------------------------
// Translate operators

Transform translation(GeomReal x, GeomReal y, GeomReal z)
{
    Transform r;

    r.ident();

    r.at(3,0) = x;
    r.at(3,1) = y;
    r.at(3,2) = z;

    return r;
}

// ---------------------------------------------------------------------------

Transform translation(const Point3& p)
{
    Transform r = translation(p.x(), p.y(), p.z());
    return r;
}

// ---------------------------------------------------------------------------

Transform& Transform::operator|=(const Point3& p)
{
    for(int i = 0; i < 4; i++)
    {
	GeomReal  r = itemAt(i,3);

	at(i,0) += p.x() * r;
	at(i,1) += p.y() * r;
	at(i,2) += p.z() * r;
    }
    return *this;
}

// ---------------------------------------------------------------------------
// (Translate(P3) + Transform) doesn't commute

Transform operator|(const Point3& p, const Transform& t)
{
    Transform r(t);

    for(int i = 0; i < 4; i++) {
	r.at(3,i) =
	    (p.x() * t[0][i]) +
	    (p.y() * t[1][i]) +
	    (p.z() * t[2][i]) +
	    t[3][i];
    }

    return r;
}

// ---------------------------------------------------------------------------
// Scale operators

Transform scaling(GeomReal x, GeomReal y, GeomReal z)
{
    Transform r;
    r.zero();
    r.at(0,0) = x;
    r.at(1,1) = y;
    r.at(2,2) = z;
    r.at(3,3) = 1.0;
    return r;
}

// ---------------------------------------------------------------------------

Transform& Transform::operator&=(const Point3& p)
{
    for (int i = 0; i < 4; i++) {
	at(i,0) *= p.x();
	at(i,1) *= p.y();
	at(i,2) *= p.z();
    }
    return *this;
}

// ---------------------------------------------------------------------------

Transform operator&(const Point3& p, const Transform& t)
{
    // (Scale & Transform) doesn't commute

    Transform	r(t);

    for (unsigned i=0; i<3; i++)
    {
	r.at(i,0) = t[i][0] * p[i];
	r.at(i,1) = t[i][1] * p[i];
	r.at(i,2) = t[i][2] * p[i];
	r.at(i,3) = t[i][3] * p[i];
    }
    return r;
}

// ---------------------------------------------------------------------------
// This functions is not yet completed!

boolean Transform::decompose(Transform& scale,
				Transform& shear,
				Transform& rot,
				Transform& trans,
				Transform& persp) const
{
    GeomReal sx,sy,sz;
    GeomReal sxy,sxz,syz;
    GeomReal rx,ry,rz;
    GeomReal tx,ty,tz;
    GeomReal px,py,pz,pw;
    
    scale.ident();
    shear.ident();
    rot.ident();
    trans.ident();
    persp.ident();

    if (!decompose(sx,sy,sz,sxy,sxz,syz,rx,ry,rz,tx,ty,tz,px,py,pz,pw))
	return FALSE;

    scale.at(0,0) = sx;
    scale.at(1,1) = sy;
    scale.at(2,2) = sz;

    shear.at(1,0) = sxy;
    shear.at(2,0) = sxz;
    shear.at(2,1) = syz;

    return TRUE;
}

// ---------------------------------------------------------------------------
// The code and idea come from _Graphic Gems_, all rights reserved, etc.

boolean Transform::decompose(GeomReal& sx,GeomReal& sy,GeomReal& sz,
				GeomReal& sxy,GeomReal& sxz,GeomReal& syz,
				GeomReal& rx,GeomReal& ry,GeomReal& rz,
				GeomReal& tx,GeomReal& ty,GeomReal& tz,
				GeomReal& px,GeomReal& py,GeomReal& pz,GeomReal& pw) const
{
    unsigned	i;

    // Clear values to start
    sx = 0; sy = 0; sz = 0;
    sxy = 0; sxz = 0; syz = 0;
    rx = 0; ry = 0; rz = 0;
    tx = 0; ty = 0; tz = 0;
    px = 0; py = 0; pz = 0; pw = 0;

    // First check the degeneracy conditions:
    // The product of the [4,4] element and upper
    // 3x3 determinant must not be zero
    if (det3() * itemAt(3,3) == 0.0)
	return FALSE;

    // Next we need to find the perspective transformations for the
    // matrix.  This ivolves solving the equation
    //
    //  [ M(0,3)     [ M(0,0)  M(0,1)  M(0,2)  0   [ px
    //    M(1,3)   =   M(1,0)  M(1,1)  M(1,2)  0     py
    //    M(2,3)       M(2,0)  M(2,1)  M(2,2)  0     pz
    //	  M(3,3) ]     M(3,0)  M(3,1)  M(3,2)  1 ]   pw ]
    //
    Transform	m(*this);
    Point4		v(m[0][3],m[1][3],m[2][3],m[3][3]);
    
    m.at(0,3) = m.at(1,3) = m.at(2,3) = 0;
    m.at(3,3) = 1;
    
    Transform	minv = inverse(m);
    Point4		p = minv * v;

    // Copy out the perspective transform factors
    px = p.x();
    py = p.y();
    pz = p.z();
    pw = p.w();

    // Clear out the perspective areas
    for (i=0; i<3; ++i) {
	m.at(i,3) = 0;
    }
    m.at(3,3) = 1;	// Set the last one back to 1

    // Copy out the transform factors
    m.decompose_translation(tx,ty,tz);
    for (i=0; i<3; ++i) {
	m.at(3,i) = 0.0;
    }

    // What is now left is only an upper 3x3 matrix (m), with the
    // extra ones zero
    // Next step - isolate the shear/scale factors
    Point3	row0(m[0]);
    Point3	row1(m[1]);
    Point3	row2(m[2]);

    sx = length(row0);		// known
    m.mulRow(0,1.0/sx);		// normalize row 0

    sxy = row0 * row1;
    m.subRows(0,1,sxy);
    row1 = m[1];		// reset row1 (its changed)
    sy = length(row1);
    m.mulRow(1,1.0/sy);		// normalize row1
    sxy/= sy;			// final value
    
    sxz = row0 * row2;
    m.subRows(0,2,sxz);
    row2 = m[2];		// reset row2

    row1 = m[1];
    syz = row1 * row2;
    m.subRows(1,2,syz);
    row2 = m[2];		// reset row2

    sz = length(row2);
    m.mulRow(2,1.0/sz);		// normalize row2
    sxy /= sz;
    syz /= sz;

    // If the determinant of the resulting matrix is -1, we have
    // to negate it and the three scaling factors
    GeomReal	det = m.det3();
					// rint is the only function availablt in Solaris
    if (int(rint(det)) == -1) {		// Kludge, == -1 float fails (eps problem)
	m = -m;
	sx = -sx;
	sy = -sy;
	sz = -sz;
    }

    m.decompose_rotation(rx,ry,rz);

    // Now for some interpretation
    // Ideally: remove rotates by 180 and push into scales
    
    return TRUE;
}

// ---------------------------------------------------------------------------
// Gets angle from the sin and cos of it
//
// Due to round-off errors in the calling method decompose_rotation()
// the arguments sin and cos maybe infinitesimal smaller than -1.0 or
// +1.0, in which case acos() and asin() return NaN or even just 0.0!
// So we need to check this.                        [Dieter 01-09-95]
//

static GeomReal find_angle(GeomReal s,GeomReal c)
{
    if (c >=0){
	if( c > 1.0 ) c = 1.0;
	if (s >= 0)		// 0 .. pi/2
	    return acos(c);
	else{			// -pi/2 .. 0
	    if( s < -1.0 ) s = -1.0;
	    return asin(s);
	}
    }
    else{
	if( c < -1.0 ) c = -1.0;
	if (s >= 0)		// pi/2 .. pi
	    return acos(c);
	else			// pi .. 3pi/2
	    return M_PI + acos(-c);
    }
}

// ---------------------------------------------------------------------------
// You must give this a rotation matrix!
// Decomposes a rotation matrix into component angles about x,y,z
// I am not sure how to treat numerical instabilities that occur
// when taking arc-sines and arc-cosines.
// This has been improved somewhat and works for the "normal" cases
// but there are still some potential problems lurking in here.
//	--Chris (with technical assistance from Ross)

boolean Transform::decompose_rotation(GeomReal& x_angle,
					 GeomReal& y_angle,
					 GeomReal& z_angle) const
{
    GeomReal	cosbeta;
    
    x_angle = 0;
    y_angle = 0;
    z_angle = 0;

    cosbeta = sqrt(1 - (itemAt(0,2)*itemAt(0,2)));
    y_angle = asin( -itemAt(0,2) );
    
    if (cosbeta == 0.0) {
	x_angle = find_angle(itemAt(1,0),itemAt(1,1));
	z_angle = 0;
    }
    else {
	x_angle = find_angle( itemAt(1,2)/cosbeta, itemAt(2,2)/cosbeta );
	z_angle = find_angle( itemAt(0,1)/cosbeta, itemAt(0,0)/cosbeta );
    }

    return TRUE;
}

// ---------------------------------------------------------------------------
// Just return the translation component of a Transform -DB

boolean Transform::decompose_translation(GeomReal& tx,
				            GeomReal& ty,
				            GeomReal& tz) const
{
    // Copy out the translation components
    tx = itemAt(3,0);
    ty = itemAt(3,1);
    tz = itemAt(3,2);
    
    return TRUE;
}

// ---------------------------------------------------------------------------
// Rotations

void Point3::rotate_only(const Transform& t)
{
    // (row vector) * Matrix4
    GeomReal xx,yy,zz;

    xx = x()*t[0][0] + y()*t[1][0] + z()*t[2][0];
    yy = x()*t[0][1] + y()*t[1][1] + z()*t[2][1];
    zz = x()*t[0][2] + y()*t[1][2] + z()*t[2][2];

    x(xx); y(yy); z(zz);
}

// ---------------------------------------------------------------------------

void Point4::rotate_only(const Transform& t)
{
    // (row vector) * Matrix4
    GeomReal xx,yy,zz;

    xx = x()*t[0][0] + y()*t[1][0] + z()*t[2][0];
    yy = x()*t[0][1] + y()*t[1][1] + z()*t[2][1];
    zz =  x()*t[0][2] + y()*t[1][2] + z()*t[2][2];

    x(xx); y(yy); z(zz);
}

// ---------------------------------------------------------------------------

Point3 rotate_only(const Point3& p,const Transform& t)
{
    Point3 r(p);
    r.rotate_only(t);
    return r;
}

// ---------------------------------------------------------------------------

Point4 rotate_only(const Point4& p,const Transform& t)
{
    Point4 r(p);
    r.rotate_only(t);
    return r;
}

// ---------------------------------------------------------------------------

Transform rotation(const IrisCoordIndex i, GeomReal theta)
{
    Transform r;
    GeomReal s = sin(theta), c = cos(theta);
    int j = ((int)i + 1) % 3, k = ((int)i + 2) % 3;

    r.ident();
    r.at(j,j) = c;
    r.at(k,k) = c;
    r.at(j,k) = s;
    r.at(k,j) = -s;
    return r;
}

// ---------------------------------------------------------------------------

Transform rotation(const Point3& p, GeomReal theta)
{
    // rotate theta radians about the direction vector p;
    // returns the transform Rv*Rz*inv(Ra) where Rv aligns the z axis
    // with p and Rz rotates by theta about the z axis.
    // Rv = Ra * Rb where Ra is a rotation about x and Rb a rotation about y
    // ca, sa, cb and -sb are the sin and cos values for Ra and Rb
    //
    Transform v, rz, r;
    Point3 n = normal(p);
    GeomReal y, z, ca, cb, sa, sb;

    y = n.y();
    z = n.z();
    cb = sqrt(y*y + z*z);
    if (cb == 0.0) {
	r = rotation(X, theta*n.x());
    }
    else {
	ca = z/cb;
	sa = y/cb;
	sb = n.x();
	//r.zero();	// Not needed since r is assigned to
	v.ident();
	v.at(0,0) = cb;
	v.at(0,2) = sb;
	v.at(1,0) = -sa*sb;
	v.at(1,1) = ca;
	v.at(1,2) = sa*cb;
	v.at(2,0) = -ca*sb;
	v.at(2,1) = -sa;
	v.at(2,2) = ca*cb;
	v.at(3,3) = 1.0;
	rz = rotation(Z,theta);
	r = v * rz;
	r *= transpose(v); // for rotations, transpose = inverse;
    }
    return r;
}

// ---------------------------------------------------------------------------
// This constructs a rotation that takes the arbitrary vector v and produces
// a rotation matrix from it that rotates the coordinate system to align
// on its axis.  Three columns are computed, c[0], which is just v, and c[1..2],
// which are computed.

Transform rotation(const Direct3& v,const IrisCoordIndex axis)
{
    Transform r;
    GeomReal z;
    Point3 c[3];		// c1 is v, or will be
    int i,j;
    int	iv,i1,i2;
    GeomReal csign = 1;

    r.ident();

    z = sqrt(1 - sqr(v.z()));

    switch (axis) {
      case X:
	iv = 0;
	i1 = 1;
	i2 = 2;
	//csign = -1;
	break;
      case Y:
	iv = 1;
	i1 = 0;
	i2 = 2;
	csign = -1;
	break;
      case Z:
	iv = 2;
	i1 = 0;
	i2 = 1;
	//csign = -1;
	break;
    }
    c[iv] = v;
    if (z == 0) {
	// c[iv] == [0,0,1] since z == 0
	c[i1] = csign * Point3(1,0,0);
	c[i2] = Point3(0,1,0);
    }
    else {
	c[i1] = csign * Point3(-v.y()/z,
			       v.x()/z,
			       0);
	c[i2] = Point3((v.x() * v.z())/z,
		       -(v.y() * v.z())/z,
		       z);
    }

    // Insert the vectors into the rotation matrix
    for (i=0; i<3; ++i) {		// row
	for (j=0; j<3; ++j) {		// col
	    r.at(i,j) = c[j][i];
	}
    }

    return r;
}

// ---------------------------------------------------------------------------

Transform rotation(const Quaternion& q)
{
    Transform r;
    Quaternion n = normal(q);
    GeomReal x = n.x();
    GeomReal y = n.y();
    GeomReal z = n.z();
    GeomReal w = n.w();
    GeomReal xx = 2*x*x;
    GeomReal yy = 2*y*y;
    GeomReal zz = 2*z*z;
    //GeomReal ww = 2*w*w;
    GeomReal xy = 2*x*y;
    GeomReal xz = 2*x*z;
    GeomReal xw = 2*x*w;
    GeomReal yz = 2*y*z;
    GeomReal yw = 2*y*w;
    GeomReal zw = 2*z*w;

    r.ident();
    r.at(0,0) = 1.0 - yy - zz;
    r.at(0,1) = xy + zw;
    r.at(0,2) = xz - yw;
    r.at(1,0) = xy - zw;
    r.at(1,1) = 1.0 - xx - zz;
    r.at(1,2) = yz + xw;
    r.at(2,0) = xz + yw;
    r.at(2,1) = yz - xw;
    r.at(2,2) = 1.0 - xx - yy;
    return r;
}

// ---------------------------------------------------------------------------

Transform perspective(GeomReal fratio)
{
    Transform r;

    r.ident();
    r.at(2,3) = fratio;

    return r;
}

// ---------------------------------------------------------------------------
