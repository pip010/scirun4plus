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
 *  Cylinder.h: Cylinder Object
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#include <Core/Util/Debug.h>

#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomCylinder.h>
#include <Core/Geom/GeomTri.h>

#include <Core/Math/MiscMath.h>

#include <iostream>

namespace SCIRun {

GeomCylinder::GeomCylinder(int nu, int nv)
: GeomObj(), bottom(0,0,0), top(0,0,1), rad(1), nu(nu), nv(nv)
{
  DEBUG_CONSTRUCTOR("GeomCylinder")
  adjust();
}

GeomCylinder::GeomCylinder(const Point& bottom, const Point& top,
			   double rad, int nu, int nv)
: GeomObj(), bottom(bottom), top(top), rad(rad), nu(nu), nv(nv)
{
  DEBUG_CONSTRUCTOR("GeomCylinder")
  adjust();
}

void 
GeomCylinder::move(const Point& _bottom, const Point& _top,
			double _rad, int _nu, int _nv)
{
  bottom=_bottom;
  top=_top;
  rad=_rad;
  nu=_nu;
  nv=_nv;
  adjust();
}

GeomCylinder::GeomCylinder(const GeomCylinder& copy)
: GeomObj(copy), v1(copy.v1), v2(copy.v2), bottom(copy.bottom), top(copy.top),
  rad(copy.rad), nu(copy.nu), nv(copy.nv)
{
  DEBUG_CONSTRUCTOR("GeomCylinder")
  adjust();
}

GeomCylinder::~GeomCylinder()
{
  DEBUG_DESTRUCTOR("GeomCylinder")
}

void 
GeomCylinder::adjust()
{
  axis=top-bottom;
  height=axis.length();
	axis.find_orthogonal(v1, v2);

  v1*=rad;
  v2*=rad;

  Vector z(0,0,1);
  if(Abs(axis.y())+Abs(axis.x()) < 1.e-5)
  {
    // Only in x-z plane...
    zrotaxis=Vector(0,-1,0);
  } 
  else 
  {
    zrotaxis=Cross(axis, z);
    zrotaxis.safe_normalize();
  }
  double cangle=Dot(z, axis)/height;
  zrotangle=-acos(cangle);
}

GeomObj* 
GeomCylinder::clone()
{
  return new GeomCylinder(*this);
}

void 
GeomCylinder::get_bounds(BBox& bb)
{
  bb.extend_disk(bottom, axis, rad);
  bb.extend_disk(top, axis, rad);
}

void
GeomCylinder::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (height < 1.e-6 || rad < 1.e-6)return;
  if (!pre_draw(di, matl, 1)) return;
  glPushMatrix();
  glTranslated(bottom.x(), bottom.y(), bottom.z());
  glRotated(RtoD(zrotangle), zrotaxis.x(), zrotaxis.y(), zrotaxis.z());
  di->polycount_+=2*(nu-1)*(nv-1);
  gluCylinder(di->qobj_, rad, rad, height, nu, nv);
  glPopMatrix();
  post_draw(di);
}

void
GeomCappedCylinder::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (height < 1.e-6 || rad < 1.e-6)return;
  if (!pre_draw(di, matl, 1)) return;
  glPushMatrix();
  glTranslated(bottom.x(), bottom.y(), bottom.z());
  glRotated(RtoD(zrotangle), zrotaxis.x(), zrotaxis.y(), zrotaxis.z());
  di->polycount_+=2*(nu-1)*(nv-1);
  gluCylinder(di->qobj_, rad, rad, height, nu, nv);
  // Bottom endcap
  di->polycount_+=2*(nu-1)*(nvdisk-1);
  gluDisk(di->qobj_, 0, rad, nu, nvdisk);
  // Top endcap
  glTranslated(0, 0, height);
  di->polycount_+=2*(nu-1)*(nvdisk-1);
  gluDisk(di->qobj_, 0, rad, nu, nvdisk);
  glPopMatrix();
  post_draw(di);
}

GeomCylinders::GeomCylinders(int nu, double r)
  : radius_(r),
    nu_(nu)
{
  DEBUG_CONSTRUCTOR("GeomCylinders")
}

GeomCylinders::GeomCylinders(const GeomCylinders& copy)
  : GeomObj(copy),
    radius_(copy.radius_),
    nu_(copy.nu_),
    points_(copy.points_),
    colors_(copy.colors_)
{
  DEBUG_CONSTRUCTOR("GeomCylinders")
}

GeomCylinders::~GeomCylinders()
{
  DEBUG_DESTRUCTOR("GeomCylinders")
}

GeomObj* 
GeomCylinders::clone()
{
  return new GeomCylinders(*this);
}

void
GeomCylinders::get_bounds(BBox& bb)
{
  for (unsigned int i = 0; i < points_.size(); i+=2)
  {
    Vector axis(points_[i] - points_[i+1]);
    bb.extend_disk(points_[i], axis, radius_);
    bb.extend_disk(points_[i], axis, radius_);
  }
}

void
GeomCylinders::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  di->polycount_+=points_.size() * nu_ * 2;

  const bool texturing =
    di->using_cmtexture_ && indices_.size() == points_.size();
  if (texturing)
  {
    glColor3d(di->diffuse_scale_, di->diffuse_scale_, di->diffuse_scale_);

    glEnable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_1D, di->cmtexture_);
  }

  const bool coloring = colors_.size() == points_.size() * 4;

  float tabx[41];
  float taby[41];
  for (int j=0; j<nu_; j++)
  {
    tabx[j] = sin(2.0 * M_PI * j / nu_);
    taby[j] = cos(2.0 * M_PI * j / nu_);
  }
  tabx[nu_] = tabx[0];
  taby[nu_] = taby[0];

  for (unsigned int i=0; i < points_.size(); i+=2)
  {
    Vector v0(points_[i+1] - points_[i+0]);

    Vector v1, v2;
    v0.find_orthogonal(v1, v2);
    if (v0.length2() < 1e-5) //numeric_limits<float>::epsilon())
      continue;
    v1 *= radius_;
    v2 *= radius_;

    float matrix[16];
    matrix[0] = v1.x();
    matrix[1] = v1.y();
    matrix[2] = v1.z();
    matrix[3] = 0.0;
    matrix[4] = v2.x();
    matrix[5] = v2.y();
    matrix[6] = v2.z();
    matrix[7] = 0.0;
    matrix[8] = v0.x();
    matrix[9] = v0.y();
    matrix[10] = v0.z();
    matrix[11] = 0.0;
    matrix[12] = points_[i].x();
    matrix[13] = points_[i].y();
    matrix[14] = points_[i].z();
    matrix[15] = 1.0;

    glPushMatrix();
    glMultMatrixf(matrix);

    glBegin(GL_QUAD_STRIP);
    for (int k=0; k<nu_+1; k++)
    {
      glNormal3f(tabx[k], taby[k], 0.0);

      if (coloring) { glColor3ubv(&(colors_[i*4])); }
      if (texturing) { glTexCoord1f(indices_[i]); }
      glVertex3f(tabx[k], taby[k], 0.0);

      if (coloring) { glColor3ubv(&(colors_[(i+1)*4])); }
      if (texturing) { glTexCoord1f(indices_[i+1]); }
      glVertex3f(tabx[k], taby[k], 1.0);
    }
    glEnd();

    glPopMatrix();
  }

  glDisable(GL_TEXTURE_1D);

  post_draw(di);
}

static unsigned char
COLOR_FTOB(double v)
{
  const int inter = (int)(v * 255 + 0.5);
  if (inter > 255) return 255;
  if (inter < 0) return 0;
  return (unsigned char)inter;
}

bool
GeomCylinders::add(const Point& p1, const Point& p2)
{
  if ((p1 - p2).length2() > 1.0e-12)
  {
    points_.push_back(p1);
    points_.push_back(p2);
    return true;
  }
  return false;
}

bool
GeomCylinders::add(const Point& p1, const MaterialHandle &c1,
		   const Point& p2, const MaterialHandle &c2)
{
  if ((p1 - p2).length2() > 1.0e-12)
  {
    points_.push_back(p1);
    points_.push_back(p2);
    
    const unsigned char r0 = COLOR_FTOB(c1->diffuse.r());
    const unsigned char g0 = COLOR_FTOB(c1->diffuse.g());
    const unsigned char b0 = COLOR_FTOB(c1->diffuse.b());
    const unsigned char a0 = COLOR_FTOB(c1->transparency);

    colors_.push_back(r0);
    colors_.push_back(g0);
    colors_.push_back(b0);
    colors_.push_back(a0);

    const unsigned char r1 = COLOR_FTOB(c2->diffuse.r());
    const unsigned char g1 = COLOR_FTOB(c2->diffuse.g());
    const unsigned char b1 = COLOR_FTOB(c2->diffuse.b());
    const unsigned char a1 = COLOR_FTOB(c2->transparency);

    colors_.push_back(r1);
    colors_.push_back(g1);
    colors_.push_back(b1);
    colors_.push_back(a1);
    return true;
  }
  return false;
}

bool
GeomCylinders::add(const Point& p1, float index1,
		   const Point& p2, float index2)
{
  if ((p1 - p2).length2() > 1.0e-12)
  {
    points_.push_back(p1);
    points_.push_back(p2);
    indices_.push_back(index1);
    indices_.push_back(index2);
    return true;
  }
  return false;
}

void
GeomCylinders::set_nu_nv(int nu, int /*nv*/)
{
  if (nu < 3) { nu_ = 3; }
  if (nu > 20) { nu_ = 20; }
  else { nu_ = nu; }
}

GeomCappedCylinder::GeomCappedCylinder(int nu, int nv, int nvdisk)
: GeomCylinder(nu, nv), nvdisk(nvdisk)
{
  DEBUG_CONSTRUCTOR("GeomCappedCylinder")
}

GeomCappedCylinder::GeomCappedCylinder(const Point& bottom, const Point& top,
				       double rad, int nu, int nv, int nvdisk)
: GeomCylinder(bottom, top, rad, nu, nv), nvdisk(nvdisk)
{
  DEBUG_CONSTRUCTOR("GeomCappedCylinder")
}

GeomCappedCylinder::GeomCappedCylinder(const GeomCappedCylinder& copy)
: GeomCylinder(copy), nvdisk(copy.nvdisk)
{
  DEBUG_CONSTRUCTOR("GeomCappedCylinder")
}

GeomCappedCylinder::~GeomCappedCylinder()
{
  DEBUG_DESTRUCTOR("GeomCappedCylinder")
}

GeomObj* 
GeomCappedCylinder::clone()
{
  return new GeomCappedCylinder(*this);
}

// Multiple Capped Geometry.

GeomCappedCylinders::GeomCappedCylinders(int nu, double r)
  : GeomCylinders(nu, r)
{
  DEBUG_CONSTRUCTOR("GeomCappedCylinders")
}

GeomCappedCylinders::GeomCappedCylinders(const GeomCappedCylinders& copy)
  : GeomCylinders(copy)
{
  DEBUG_CONSTRUCTOR("GeomCappedCylinders")
}

GeomCappedCylinders::~GeomCappedCylinders()
{
  DEBUG_DESTRUCTOR("GeomCappedCylinders")
}

void
GeomCappedCylinders::add_radius(const Point& p1, const Point& p2, double r)
{
  if (add(p1, p2)) { radii_.push_back(r); }
}

void
GeomCappedCylinders::add_radius(const Point& p1, const MaterialHandle &c1,
				const Point& p2, const MaterialHandle &c2,
				double r)
{
  if (add(p1, c1, p2, c2)) { radii_.push_back(r); }
}

void
GeomCappedCylinders::add_radius(const Point& p1, float index1,
				const Point& p2, float index2, double r)
{
  if (add(p1, index1, p2, index2)) { radii_.push_back(r); }
}

GeomObj* GeomCappedCylinders::clone()
{
  return new GeomCappedCylinders(*this);
}

void
GeomCappedCylinders::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  di->polycount_+=points_.size() * nu_ * 2;

  const bool texturing =
    di->using_cmtexture_ && indices_.size() == points_.size();
  if (texturing)
  {
    glColor3d(di->diffuse_scale_, di->diffuse_scale_, di->diffuse_scale_);

    glEnable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_1D, di->cmtexture_);
  }

  const bool coloring = colors_.size() == points_.size() * 4;
  const bool use_local_radii = radii_.size() == points_.size()/2;

  float tabx[41];
  float taby[41];
  for (int j=0; j<nu_; j++)
  {
    tabx[j] = sin(2.0 * M_PI * j / nu_);
    taby[j] = cos(2.0 * M_PI * j / nu_);
  }
  tabx[nu_] = tabx[0];
  taby[nu_] = taby[0];

  for (unsigned int i=0; i < points_.size(); i+=2)
  {
    int k;
    Vector v0(points_[i+1] - points_[i+0]);
    Vector v1, v2;
    v0.find_orthogonal(v1, v2);
    if (use_local_radii)
    {
      v1 *= radii_[i/2];
      v2 *= radii_[i/2];
    }
    else
    {
      v1 *= radius_;
      v2 *= radius_;
    }

    float matrix[16];
    matrix[0] = v1.x();
    matrix[1] = v1.y();
    matrix[2] = v1.z();
    matrix[3] = 0.0;
    matrix[4] = v2.x();
    matrix[5] = v2.y();
    matrix[6] = v2.z();
    matrix[7] = 0.0;
    matrix[8] = v0.x();
    matrix[9] = v0.y();
    matrix[10] = v0.z();
    matrix[11] = 0.0;
    matrix[12] = points_[i].x();
    matrix[13] = points_[i].y();
    matrix[14] = points_[i].z();
    matrix[15] = 1.0;

    glPushMatrix();
    glMultMatrixf(matrix);

    glBegin(GL_QUAD_STRIP);
    for (k=0; k<nu_+1; k++)
    {
      glNormal3f(tabx[k], taby[k], 0.0);

      if (coloring) { glColor3ubv(&(colors_[i*4])); }
      if (texturing) { glTexCoord1f(indices_[i]); }
      glVertex3f(tabx[k], taby[k], 0.0);

      if (coloring) { glColor3ubv(&(colors_[(i+1)*4])); }
      if (texturing) { glTexCoord1f(indices_[i+1]); }
      glVertex3f(tabx[k], taby[k], 1.0);
    }
    glEnd();

    // Bottom cap
    if (coloring) { glColor3ubv(&(colors_[i*4])); }
    if (texturing) { glTexCoord1f(indices_[i]); }
    glNormal3f(0.0, 0.0, -1.0);
    glBegin(GL_TRIANGLE_FAN);
    glVertex3f(0.0, 0.0, 0.0);
    for (k = 0; k < nu_+1; k++)
    {
      glVertex3f(tabx[k], taby[k], 0.0);
    }
    glEnd();

    // Top cap
    if (coloring) { glColor3ubv(&(colors_[(i+1)*4])); }
    if (texturing) { glTexCoord1f(indices_[i+1]); }
    glNormal3f(0.0, 0.0, 1.0);
    glBegin(GL_TRIANGLE_FAN);
    glVertex3f(0.0, 0.0, 1.0);
    for (k = nu_; k >= 0; k--)
    {
      glVertex3f(tabx[k], taby[k], 1.0);
    }
    glEnd();

    glPopMatrix();
  }

  glDisable(GL_TEXTURE_1D);

  post_draw(di);
}

} // End namespace SCIRun

