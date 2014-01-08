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

#include <Core/Util/Debug.h>

#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomCone.h>
#include <Core/Geom/GeomTri.h>

#include <Core/Math/MiscMath.h>
#include <iostream>

namespace SCIRun {

GeomCone::GeomCone(int nu, int nv)
  : GeomObj(),
    nu(nu), nv(nv),
    bottom(0,0,0), top(0,0,1),
    bot_rad(1), top_rad(0)
{
  DEBUG_CONSTRUCTOR("GeomCone")
  adjust();
}

GeomCone::GeomCone(const Point& bottom, const Point& top,
		   double bot_rad, double top_rad, int nu, int nv)
: GeomObj(),
  nu(nu), nv(nv),
  bottom(bottom), top(top),
  bot_rad(bot_rad), top_rad(top_rad)
{
  DEBUG_CONSTRUCTOR("GeomCone")
  adjust();
}

void GeomCone::move(const Point& _bottom, const Point& _top,
		    double _bot_rad, double _top_rad, int _nu, int _nv)
{
  bottom=_bottom;
  top=_top;
  bot_rad=_bot_rad;
  top_rad=_top_rad;
  nu=_nu;
  nv=_nv;
  adjust();
}

GeomCone::GeomCone(const GeomCone& copy)
  : GeomObj(),
    v1(copy.v1), v2(copy.v2),
    nu(copy.nu), nv(copy.nv),
    bottom(copy.bottom), top(copy.top),
    bot_rad(copy.bot_rad), top_rad(copy.top_rad)
{
  DEBUG_CONSTRUCTOR("GeomCone")
  adjust();
}

GeomCone::~GeomCone()
{
  DEBUG_DESTRUCTOR("GeomCone")
}

GeomObj* GeomCone::clone()
{
  return new GeomCone(*this);
}

void 
GeomCone::adjust()
{
  axis=top-bottom;
  height=axis.length();
  if(height < 1.e-6)
  {
    std::cerr << "Degenerate Cone!\n";
  } 
  else 
  {
    axis.find_orthogonal(v1, v2);
  }
  tilt=(bot_rad-top_rad)/axis.length2();
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

void 
GeomCone::get_bounds(BBox& bb)
{
  bb.extend_disk(bottom, axis, bot_rad);
  bb.extend_disk(top, axis, top_rad);
}

void
GeomCone::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (height < 1.e-6 || (bot_rad < 1.e-6 && top_rad < 1.e-6)) return;
  if (!pre_draw(di, matl, 1)) return;
  glPushMatrix();
  glTranslated(bottom.x(), bottom.y(), bottom.z());
  glRotated(RtoD(zrotangle), zrotaxis.x(), zrotaxis.y(), zrotaxis.z());
  di->polycount_+=2*(nu-1)*(nv-1);
  gluCylinder(di->qobj_, bot_rad, top_rad, height, nu, nv);
  glPopMatrix();
  post_draw(di);
}

// Capped Geometry

GeomCappedCone::GeomCappedCone(int nu, int nv, int nvdisk1, int nvdisk2)
: GeomCone(nu, nv), nvdisk1(nvdisk1), nvdisk2(nvdisk2)
{
  DEBUG_CONSTRUCTOR("GeomCappedCone")
}

GeomCappedCone::GeomCappedCone(const Point& bottom, const Point& top,
			       double bot_rad, double top_rad, int nu, int nv,
			       int nvdisk1, int nvdisk2)
: GeomCone(bottom, top, bot_rad, top_rad, nu, nv), nvdisk1(nvdisk1),
  nvdisk2(nvdisk2)
{
  DEBUG_CONSTRUCTOR("GeomCappedCone")
}

GeomCappedCone::GeomCappedCone(const GeomCappedCone& copy)
: GeomCone(copy), nvdisk1(copy.nvdisk1), nvdisk2(copy.nvdisk2)
{
  DEBUG_CONSTRUCTOR("GeomCappedCone")
}

GeomCappedCone::~GeomCappedCone()
{
  DEBUG_DESTRUCTOR("GeomCappedCone")
}

GeomObj* 
GeomCappedCone::clone()
{
  return new GeomCappedCone(*this);
}

void
GeomCappedCone::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (height < 1.e-6 || (bot_rad < 1.e-6 && top_rad < 1.e-6))return;
  if (!pre_draw(di, matl, 1)) return;
  glPushMatrix();
  glTranslated(bottom.x(), bottom.y(), bottom.z());
  glRotated(RtoD(zrotangle), zrotaxis.x(), zrotaxis.y(), zrotaxis.z());
  di->polycount_+=2*(nu-1)*(nv-1);
  gluCylinder(di->qobj_, bot_rad, top_rad, height, nu, nv);
  if (bot_rad > 1.e-6)
  {
    // Bottom endcap
    di->polycount_+=2*(nu-1)*(nvdisk1-1);
    gluDisk(di->qobj_, 0, bot_rad, nu, nvdisk1);
  }
  if (top_rad > 1.e-6)
  {
    // Top endcap
    glTranslated(0, 0, height);
    di->polycount_+=2*(nu-1)*(nvdisk2-1);
    gluDisk(di->qobj_, 0, top_rad, nu, nvdisk2);
  }
  glPopMatrix();
  post_draw(di);
}

// GeomCones, accelerated for many objects.

GeomCones::GeomCones(int nu, double r)
  : radius_(r),
    nu_(nu)
{
  DEBUG_CONSTRUCTOR("GeomCones")
}

GeomCones::GeomCones(const GeomCones& copy)
  : GeomObj(copy),
    radius_(copy.radius_),
    nu_(copy.nu_),
    points_(copy.points_),
    colors_(copy.colors_)
{
  DEBUG_CONSTRUCTOR("GeomCones")
}

GeomCones::~GeomCones()
{
  DEBUG_DESTRUCTOR("GeomCones")
}

GeomObj* GeomCones::clone()
{
  return new GeomCones(*this);
}

void
GeomCones::get_bounds(BBox& bb)
{
  for (unsigned int i = 0; i < points_.size(); i+=2)
  {
    Vector axis(points_[i] - points_[i+1]);
    bb.extend_disk(points_[i], axis, radius_);
    bb.extend_disk(points_[i+1], axis, 0);
  }
}

void
GeomCones::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  di->polycount_ += points_.size() * nu_ / 2;

  const bool texturing =
    di->using_cmtexture_ && indices_.size() == points_.size() / 2;
  if (texturing)
  {
    glColor3d(di->diffuse_scale_, di->diffuse_scale_, di->diffuse_scale_);

    glEnable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_1D, di->cmtexture_);
  }

  const bool coloring = colors_.size() == points_.size() * 2;
  const bool use_local_radii = radii_.size() == points_.size()/2;

  const float nz0 = 1.0/6.0;
  const float nzm = 1.0/sqrt(1.0 + 1.0 + nz0*nz0);
  const float nz = nz0 * nzm;
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

    if (coloring) { glColor3ubv(&(colors_[i*2])); }
    if (texturing) { glTexCoord1f(indices_[i/2]); }

    glBegin(GL_QUAD_STRIP);
    for (int k = 0; k <= nu_; k++)
    {
      glNormal3f(tabx[k]*nzm, taby[k]*nzm, nz);
      glVertex3f(tabx[k], taby[k], 0.0);
      glVertex3f(0.0, 0.0, 1.0);
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
GeomCones::add(const Point& p1, const Point& p2)
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
GeomCones::add(const Point& p1, const Point &p2, const MaterialHandle &c)
{
  if ((p1 - p2).length2() > 1.0e-12)
  {
    points_.push_back(p1);
    points_.push_back(p2);
    
    const unsigned char r0 = COLOR_FTOB(c->diffuse.r());
    const unsigned char g0 = COLOR_FTOB(c->diffuse.g());
    const unsigned char b0 = COLOR_FTOB(c->diffuse.b());
    const unsigned char a0 = COLOR_FTOB(c->transparency);

    colors_.push_back(r0);
    colors_.push_back(g0);
    colors_.push_back(b0);
    colors_.push_back(a0);

    return true;
  }
  return false;
}

bool
GeomCones::add(const Point& p1, const Point& p2, float index)
{
  if ((p1 - p2).length2() > 1.0e-12)
  {
    points_.push_back(p1);
    points_.push_back(p2);
    indices_.push_back(index);
    return true;
  }
  return false;
}

bool
GeomCones::add_radius(const Point& p1, const Point& p2, double r)
{
  if (add(p1, p2))
  {
    radii_.push_back(r);
    return true;
  }
  return false;
}

bool
GeomCones::add_radius(const Point& p1, const Point &p2,
		      const MaterialHandle &c, double r)
{
  if (add(p1, p2, c))
  {
    radii_.push_back(r);
    return true;
  }
  return false;
}

bool
GeomCones::add_radius(const Point& p1, const Point& p2, float index, double r)
{
  if (add(p1, p2, index))
  {
    radii_.push_back(r);
    return true;
  }
  return false;
}

void
GeomCones::set_nu(int nu)
{
  if (nu < 3) { nu_ = 3; }
  if (nu > 40) { nu_ = 40; }
  else { nu_ = nu; }
}

} // End namespace SCIRun

