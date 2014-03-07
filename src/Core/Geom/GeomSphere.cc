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
 * GeomSphere.cc: Sphere objects
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
#include <Core/Geom/GeomSphere.h>
#include <Core/Geom/GeomTri.h>

#include <Core/Math/MiscMath.h>

namespace SCIRun {

GeomSphere::GeomSphere(int nu, int nv)
: GeomObj(), cen(0,0,0), rad(1), nu(nu), nv(nv)
{
  DEBUG_CONSTRUCTOR("GeomSphere")
  adjust();
}

GeomSphere::GeomSphere(const Point& cen, double rad, int nu, int nv)
: GeomObj(), cen(cen), rad(rad), nu(nu), nv(nv)
{
  DEBUG_CONSTRUCTOR("GeomSphere")
  adjust();
}

void GeomSphere::move(const Point& _cen, double _rad, int _nu, int _nv)
{
  cen=_cen;
  rad=_rad;
  nu=_nu;
  nv=_nv;
  adjust();
}

void GeomSphere::move(const Point& _cen) 
{
  cen = _cen;
  adjust();
}

GeomSphere::GeomSphere(const GeomSphere& copy)
: GeomObj(copy), cen(copy.cen), rad(copy.rad), nu(copy.nu), nv(copy.nv)
{
  DEBUG_CONSTRUCTOR("GeomSphere")
  adjust();
}

GeomSphere::~GeomSphere()
{
  DEBUG_DESTRUCTOR("GeomSphere")
}

void GeomSphere::adjust()
{
}

GeomObj* GeomSphere::clone()
{
  return new GeomSphere(*this);
}

void GeomSphere::get_bounds(BBox& bb)
{
  bb.extend(cen, rad);
}

void
GeomSphere::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (rad < 1.e-6)return;
  if (!pre_draw(di, matl, 1)) return;
  glPushMatrix();

  glTranslated(cen.x(), cen.y(), cen.z());
  di->polycount_+=2*(nu-1)*(nv-1);

  gluSphere(di->qobj_, rad, nu, nv);

  glPopMatrix();
  post_draw(di);
}

void GeomSphere::getnunv(int num_polygons, int &nu, int &nv) 
{
#define MIN_POLYS 8
#define MAX_POLYS 400
#define MIN_NU 4
#define MAX_NU 20
#define MIN_NV 2
#define MAX_NV 20
  // calculate the spheres nu,nv based on the number of polygons
  float t = (num_polygons - MIN_POLYS)/float(MAX_POLYS - MIN_POLYS);
  nu = int(MIN_NU + t*(MAX_NU - MIN_NU)); 
  nv = int(MIN_NV + t*(MAX_NV - MIN_NV));
}

static inline float
spow(double e, double x)
{
  // This for round off of very small numbers as raising them to a
  // power though gives correct results the input is really erroneous.
  if( fabs( e ) < 1.0e-6)
    e = 0.0;

  if (e < 0.0)
  {
    return (float)(pow(fabs(e), x) * -1.0);
  }
  else
  {
    return (float)(pow(e, x));
  }
}

GeomSuperquadric::GeomSuperquadric()
{
  DEBUG_CONSTRUCTOR("GeomSuperquadric")
}

GeomSuperquadric::GeomSuperquadric(int axis, double A, double B,
				   int nu, int nv)
  : axis_(axis),
    A_(A),
    B_(B),
    nu_(nu),
    nv_(nv)
{
  DEBUG_CONSTRUCTOR("GeomSuperquadric")
  nu_ = Max(2, nu_);
  nv_ = Max(3, nv_);
  // TODO:  ASSERT nu_ * nv_ not a short overflow.

  compute_geometry();
}

void
GeomSuperquadric::compute_geometry()
{
  int ti, pi;

  // Skip the first which all are 1,0,0 and the last which all are
  // 1,0,0. The last is implicitly skipped by dividing by nv_ rather
  // than nv-1. They are built "maunally" as triangle fans.
  for (pi = 1; pi < nv_; pi++)
  {
    const double p = (pi / (double) (nv_)) * M_PI;

    // The last is implicitly skipped by dividing by nu_ rather than
    // nu_-1. It is built by adding the first in at the beginning and end.
    for (ti = 0; ti < nu_; ti++)
    {
      const double t = (ti / (double)nu_) * (2.0 * M_PI);
      const double x = spow(cos(p), B_);
      const double y = spow(cos(t), A_) * spow(sin(p), B_);
      const double z = spow(sin(t), A_) * spow(sin(p), B_);

      const double xxb = pow(x*x, 1/B_);
      const double yya = pow(y*y, 1/A_);
      const double zza = pow(z*z, 1/A_);
      const double R = pow(yya + zza, (A_/B_)-1);
      const float nx = 2*xxb/(B_*x + 1.0e-6);
      const float ny = 2*R*yya/(B_*y + 1.0e-6);
      const float nz = 2*R*zza/(B_*z + 1.0e-6);
      const float nn = 1.0 / sqrt(nx * nx + ny * ny + nz * nz);
      
      switch(axis_) 
      {
      case 0:
        points_.push_back(x);
        points_.push_back(y);
        points_.push_back(z);
        normals_.push_back(nx * nn);
        normals_.push_back(ny * nn);
        normals_.push_back(nz * nn);
        break;
      case 1:
        points_.push_back(z);
        points_.push_back(x);
        points_.push_back(y);
        normals_.push_back(nz * nn);
        normals_.push_back(nx * nn);
        normals_.push_back(ny * nn);
        break;
      case 2:
        points_.push_back(y);
        points_.push_back(z);
        points_.push_back(x);
        normals_.push_back(ny * nn);
        normals_.push_back(nz * nn);
        normals_.push_back(nx * nn);
        break;
      }
    }
  }

  // Add north pole.
  switch(axis_) 
  {
  case 0:
    points_.push_back(1.0);
    points_.push_back(0.0);
    points_.push_back(0.0);
    normals_.push_back(1.0);
    normals_.push_back(0.0);
    normals_.push_back(0.0);
    break;
  case 1:
    points_.push_back(0.0);
    points_.push_back(1.0);
    points_.push_back(0.0);
    normals_.push_back(0.0);
    normals_.push_back(1.0);
    normals_.push_back(0.0);
    break;
  case 2: default:
    points_.push_back(0.0);
    points_.push_back(0.0);
    points_.push_back(1.0);
    normals_.push_back(0.0);
    normals_.push_back(0.0);
    normals_.push_back(1.0);
    break;
  }

  // Add south pole.
  switch(axis_) 
  {
  case 0:
    points_.push_back(-1.0);
    points_.push_back(0.0);
    points_.push_back(0.0);
    normals_.push_back(-1.0);
    normals_.push_back(0.0);
    normals_.push_back(0.0);
    break;
  case 1:
    points_.push_back(0.0);
    points_.push_back(-1.0);
    points_.push_back(0.0);
    normals_.push_back(0.0);
    normals_.push_back(-1.0);
    normals_.push_back(0.0);
    break;
  case 2: default:
    points_.push_back(0.0);
    points_.push_back(0.0);
    points_.push_back(-1.0);
    normals_.push_back(0.0);
    normals_.push_back(0.0);
    normals_.push_back(-1.0);
    break;
  }

  // North pole cap
  tindices_.push_back(points_.size()/3 - 2);
  for (ti = 0; ti < nu_; ti++)
  {
    tindices_.push_back(ti);
  }
  tindices_.push_back(0); // Wrap back around
    
  // South pole cap.
  tindices_.push_back(points_.size()/3 - 1);
  for (ti = nu_; ti >= 0; ti--)
  {
    tindices_.push_back((nv_-2) * nu_ + (ti%nu_));
  }

  // Equatorial region
  for (pi = 0; pi < nv_-2; pi++)
  {
    for (ti=0; ti < nu_; ti++)
    {
      qindices_.push_back(pi * nu_ + ti);
      qindices_.push_back((pi+1) * nu_ + ti);
    }
    qindices_.push_back(pi * nu_);     // Wrap back around
    qindices_.push_back((pi+1) * nu_); // Wrap back around
  }
}

GeomSuperquadric::GeomSuperquadric(const GeomSuperquadric& copy)
  : GeomObj(copy)
{  
  DEBUG_CONSTRUCTOR("GeomSuperquadric")
}

GeomSuperquadric::~GeomSuperquadric()
{
  DEBUG_DESTRUCTOR("GeomSuperquadric")
}

GeomObj* GeomSuperquadric::clone()
{
  return new GeomSuperquadric(*this);
}

void GeomSuperquadric::get_bounds(BBox& bb)
{
  bb.extend(Point(1.0, 1.0, 1.0));
  bb.extend(Point(-1.0, -1.0, -1.0));
}

void
GeomSuperquadric::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  di->polycount_ += nu_ * nv_;

  if (points_.size())
  {
    glVertexPointer(3, GL_FLOAT, 0, &(points_.front()));
    glEnableClientState(GL_VERTEX_ARRAY);
  }
  else
  {
    glDisableClientState(GL_VERTEX_ARRAY);
  }

  if (normals_.size())
  {
    glNormalPointer(GL_FLOAT, 0, &(normals_.front()));
    glEnableClientState(GL_NORMAL_ARRAY);
  }
  else
  {
    glDisableClientState(GL_NORMAL_ARRAY);
  }

  glDisableClientState(GL_COLOR_ARRAY);

  if (tindices_.size()) {
    glDrawElements(GL_TRIANGLE_FAN, nu_ + 2, GL_UNSIGNED_SHORT,
                   &(tindices_[0]));

    glDrawElements(GL_TRIANGLE_FAN, nu_ + 2, GL_UNSIGNED_SHORT,
                   &(tindices_[nu_+2]));
  }

  if (qindices_.size()) {
    for (int pi = 0; pi < nv_-2; pi++)
    {
      glDrawElements(GL_QUAD_STRIP, (nu_+1)*2, GL_UNSIGNED_SHORT,
                     &(qindices_[pi * (nu_+1) * 2]));
    }
  }
  post_draw(di);
}

GeomSpheres::GeomSpheres(double radius, int nu, int nv)
  : GeomObj(),
    nu_(nu),
    nv_(nv),
    global_radius_(radius)
{
  DEBUG_CONSTRUCTOR("GeomSpheres")
}

GeomSpheres::GeomSpheres(const GeomSpheres& copy)
  : GeomObj(copy),
    centers_(copy.centers_),
    radii_(copy.radii_),
    colors_(copy.colors_),
    indices_(copy.indices_),
    nu_(copy.nu_),
    nv_(copy.nv_),
    global_radius_(copy.global_radius_)
{
  DEBUG_CONSTRUCTOR("GeomSpheres")
}

GeomSpheres::~GeomSpheres()
{
  DEBUG_DESTRUCTOR("GeomSpheres")
}

GeomObj *
GeomSpheres::clone()
{
  return new GeomSpheres(*this);
}

void
GeomSpheres::get_bounds(BBox& bb)
{
  const bool ugr = !(radii_.size() == centers_.size());
  for (unsigned int i=0; i < centers_.size(); i++)
  {
    bb.extend(centers_[i], ugr?global_radius_:radii_[i]);
  }
}

static unsigned char
COLOR_FTOB(double v)
{
  const int inter = (int)(v * 255 + 0.5);
  if (inter > 255) return 255;
  if (inter < 0) return 0;
  return (unsigned char)inter;
}


void
GeomSpheres::add(const Point &center, unsigned int idx)
{
  centers_.push_back(center);
  item_idx_.push_back(idx);
}


void
GeomSpheres::add(const Point &center, const MaterialHandle &mat,
		 unsigned int idx)
{
  add(center);
  const unsigned char r0 = COLOR_FTOB(mat->diffuse.r());
  const unsigned char g0 = COLOR_FTOB(mat->diffuse.g());
  const unsigned char b0 = COLOR_FTOB(mat->diffuse.b());
  const unsigned char a0 = COLOR_FTOB(mat->transparency);
  colors_.push_back(r0);
  colors_.push_back(g0);
  colors_.push_back(b0);
  colors_.push_back(a0);
  item_idx_.push_back(idx);
}


void
GeomSpheres::add(const Point &center, float index, unsigned int idx)
{
  add(center);
  indices_.push_back(index);
  item_idx_.push_back(idx);
}


bool
GeomSpheres::add_radius(const Point &c, double r)
{
  if (r < 1.0e-6) { return false; }
  centers_.push_back(c);
  radii_.push_back(r);
  return true;
}

bool
GeomSpheres::add_radius(const Point &c, double r, const MaterialHandle &mat)
{
  if (r < 1.0e-6) { return false; }
  add_radius(c, r);
  const unsigned char r0 = COLOR_FTOB(mat->diffuse.r());
  const unsigned char g0 = COLOR_FTOB(mat->diffuse.g());
  const unsigned char b0 = COLOR_FTOB(mat->diffuse.b());
  const unsigned char a0 = COLOR_FTOB(mat->transparency);
  colors_.push_back(r0);
  colors_.push_back(g0);
  colors_.push_back(b0);
  colors_.push_back(a0);
  return true;
}

bool
GeomSpheres::add_radius(const Point &c, double r, float index)
{
  if (r < 1.0e-6) { return false; }
  add_radius(c, r);
  indices_.push_back(index);
  return true;
}


void
GeomSpheres::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  const bool ulr = radii_.size() == centers_.size();
  if (!ulr && global_radius_ < 1.0e-6) { return; }

  if (!pre_draw(di, matl, 1)) return;

  di->polycount_ += 2 * (nu_-1) * (nv_-1) * centers_.size();

  const bool using_texture =
    di->using_cmtexture_ && indices_.size() == centers_.size();
  if (using_texture)
  {
    glColor3d(di->diffuse_scale_, di->diffuse_scale_, di->diffuse_scale_);

    glEnable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_1D, di->cmtexture_);
  }

  const bool using_color = centers_.size() == colors_.size() / 4;

  glMatrixMode(GL_MODELVIEW);
  for (unsigned int i=0; i < centers_.size(); i++)
  {
    if (using_texture) { glTexCoord1f(indices_[i]); }
    if (using_color) { glColor3ubv(&(colors_[i*4])); }

    glPushMatrix();

    glTranslated(centers_[i].x(), centers_[i].y(), centers_[i].z());
    gluSphere(di->qobj_, ulr?radii_[i]:global_radius_, nu_, nv_);
        
    glPopMatrix();
  }

  glDisable(GL_TEXTURE_1D);

  post_draw(di);
}


void
GeomSpheres::fbpick_draw(DrawInfoOpenGL* di, Material* matl, double)
{
  const bool ulr = radii_.size() == centers_.size();
  if (!ulr && global_radius_ < 1.0e-6) { return; }

  if (!pre_draw(di, matl, 1)) return;

  di->polycount_ += 2 * (nu_-1) * (nv_-1) * centers_.size();

  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glMatrixMode(GL_MODELVIEW);

  unsigned char cols[4];
  for (unsigned int i = 0; i < centers_.size(); i++)
  {
    unsigned char r, g, b, a;
    // we use 0 as a non index, so increment all indecies by one,
    // here while we decrement after reading..
    
    unsigned int idx = item_idx_[i] + 1;
    //just encode the index in the rgb channels
    // this type of encoding can only handle indeces 2^24
    if ((unsigned int)(idx & 0xff000000) != 0) 
    { 
      std::cerr << "Error: encoded index is > 2^24: " << idx << std::endl;
      std::cerr << "       Picking may fail." << std::endl;
    }
    const unsigned int rmask = 0x00ff0000;
    const unsigned int gmask = 0x0000ff00;
    const unsigned int bmask = 0x000000ff;
    r = (idx & rmask) >> 16;
    g = (idx & gmask) >> 8;
    b = (idx & bmask);
    a = 255;
    
    cols[0] = r;
    cols[1] = g;
    cols[2] = b;
    cols[3] = a;

    glColor3ubv(&(cols[0]));

    glPushMatrix();

    glTranslated(centers_[i].x(), centers_[i].y(), centers_[i].z());
    gluSphere(di->qobj_, ulr?radii_[i]:global_radius_, nu_, nv_);

    glPopMatrix();
  }

  post_draw(di);
}

} // End namespace SCIRun

