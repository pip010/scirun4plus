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
 * Torus.cc: Torus objects
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
#include <Core/Geom/GeomTorus.h>
#include <Core/Geom/GeomTri.h>
#include <Core/Geometry/Transform.h>

#include <Core/Math/MiscMath.h>
#include <Core/Math/TrigTable.h>

namespace SCIRun {

GeomTorus::GeomTorus(int nu_, int nv_)
: GeomObj(), center_(0,0,0), axis_(0,0,1),
  major_radius_(1), minor_radius_(.1),
  nu_(nu_), nv_(nv_)
{
  DEBUG_CONSTRUCTOR("GeomTorus")
  adjust();
}

GeomTorus::GeomTorus(const Point& center, const Vector& axis,
		     double major_radius, double minor_radius,
		     int nu, int nv)
: GeomObj(), center_(center), axis_(axis),
  major_radius_(major_radius), minor_radius_(minor_radius),
  nu_(nu), nv_(nv)
{
  DEBUG_CONSTRUCTOR("GeomTorus")
  adjust();
}

GeomTorus::GeomTorus(const GeomTorus& copy)
: GeomObj(copy), center_(copy.center_), axis_(copy.axis_),
  major_radius_(copy.major_radius_), minor_radius_(copy.minor_radius_),
  nu_(copy.nu_), nv_(copy.nv_)
{
  DEBUG_CONSTRUCTOR("GeomTorus")
  adjust();
}

GeomTorus::~GeomTorus()
{
  DEBUG_DESTRUCTOR("GeomTorus")
}

GeomObj* 
GeomTorus::clone()
{
  return new GeomTorus(*this);
}

void 
GeomTorus::move(const Point& center, const Vector& axis,
		     double major_radius, double minor_radius,
		     int nu, int nv)
{
  center_ = center;
  axis_ = axis;
  major_radius_ = major_radius;
  minor_radius_ = minor_radius;
  nu_ = nu;
  nv_ = nv;
  adjust();
}

void 
GeomTorus::adjust()
{
  axis_.normalize();

  Vector z(0,0,1);
  if(Abs(axis_.y())+Abs(axis_.x()) < 1.e-5) 
  {
    // Only in x-z plane...
    zrotaxis_=Vector(0,-1,0);
  } 
  else 
  {
    zrotaxis_=Cross(axis_, z);
    zrotaxis_.normalize();
  }
  double cangle = Dot(z, axis_);
  zrotangle_ = -acos(cangle);
}

void 
GeomTorus::get_bounds(BBox& bb)
{
  bb.extend_disk(center_-axis_*minor_radius_, axis_, major_radius_+minor_radius_);
  bb.extend_disk(center_+axis_*minor_radius_, axis_, major_radius_+minor_radius_);
}

void
GeomTorus::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;
  glPushMatrix();
  glTranslated(center_.x(), center_.y(), center_.z());
  glRotated(RtoD(zrotangle_), zrotaxis_.x(), zrotaxis_.y(), zrotaxis_.z());
  di->polycount_+=2*(nu_-1)*(nv_-1);

  // Draw the torus
  SinCosTable tab1(nu_, 0, 2*M_PI);
  SinCosTable tab2 (nv_, 0, 2*M_PI, minor_radius_);
  SinCosTable tab2n(nv_, 0, 2*M_PI, minor_radius_);
  int u,v;
  switch(di->get_drawtype())
  {
  case DrawInfoOpenGL::WireFrame:
    for (u=0; u<nu_; u++)
    {
      double rx=tab1.sin(u);
      double ry=tab1.cos(u);
      glBegin(GL_LINE_LOOP);
      for (v=1; v<nv_; v++)
      {
        double z=tab2.cos(v);
        double rad=major_radius_+tab2.sin(v);
        double x=rx*rad;
        double y=ry*rad;
        glVertex3d(x, y, z);
      }
      glEnd();
    }
    for (v=1; v<nv_; v++)
    {
      double z=tab2.cos(v);
      double rr=tab2.sin(v);
      glBegin(GL_LINE_LOOP);
      for (u=1; u<nu_; u++)
      {
        double rad=major_radius_+rr;
        double x=tab1.sin(u)*rad;
        double y=tab1.cos(u)*rad;
        glVertex3d(x, y, z);
      }
      glEnd();
    }
    break;
  case DrawInfoOpenGL::Flat:
    for (v=0; v<nv_-1; v++)
    {
      double z1=tab2.cos(v);
      double rr1=tab2.sin(v);
      double z2=tab2.cos(v+1);
      double rr2=tab2.sin(v+1);
      double nr=-tab2n.sin(v);
      double nz=-tab2n.cos(v);
      glBegin(GL_TRIANGLE_STRIP);
      for (u=0; u<nu_; u++)
      {
        double r1=major_radius_+rr1;
        double r2=major_radius_+rr2;
        double xx=tab1.sin(u);
        double yy=tab1.cos(u);
        double x1=xx*r1;
        double y1=yy*r1;
        double x2=xx*r2;
        double y2=yy*r2;
        glNormal3d(nr*xx, nr*yy, nz);
        glVertex3d(x1, y1, z1);
        glVertex3d(x2, y2, z2);
      }
      glEnd();
    }
    break;
  case DrawInfoOpenGL::Gouraud:
    for (v=0; v<nv_-1; v++)
    {
      double z1=tab2.cos(v);
      double rr1=tab2.sin(v);
      double z2=tab2.cos(v+1);
      double rr2=tab2.sin(v+1);
      double nr1=-tab2n.sin(v);
      double nr2=-tab2n.sin(v+1);
      double nz1=-tab2n.cos(v);
      double nz2=-tab2n.cos(v+1);
      glBegin(GL_TRIANGLE_STRIP);
      for (u=0; u<nu_; u++)
      {
        double r1=major_radius_+rr1;
        double r2=major_radius_+rr2;
        double xx=tab1.sin(u);
        double yy=tab1.cos(u);
        double x1=xx*r1;
        double y1=yy*r1;
        double x2=xx*r2;
        double y2=yy*r2;
        glNormal3d(nr1*xx, nr1*yy, nz1);
        glVertex3d(x1, y1, z1);
        glNormal3d(nr2*xx, nr2*yy, nz2);
        glVertex3d(x2, y2, z2);
      }
      glEnd();
    }
    break;      
  }
  glPopMatrix();
  post_draw(di);
}

// GeomTorusArc

GeomTorusArc::GeomTorusArc(int nu_, int nv_)
  : GeomTorus(nu_, nv_), zero_(0,1,0), arc_angle_(M_PI)
{
}

GeomTorusArc::GeomTorusArc(const Point& center, const Vector& axis,
			   double major_radius, double minor_radius,
			   const Vector& zero,
			   double start_angle, double arc_angle,
			   int nu, int nv)
  : GeomTorus(center, axis, major_radius, minor_radius, nu, nv),
    zero_(zero), start_angle_(start_angle), arc_angle_(arc_angle)
{
}

GeomTorusArc::GeomTorusArc(const GeomTorusArc& copy)
  : GeomTorus(copy),
    zero_(copy.zero_),
    start_angle_(copy.start_angle_),
    arc_angle_(copy.arc_angle_)
{
}

GeomTorusArc::~GeomTorusArc()
{
}

GeomObj* 
GeomTorusArc::clone()
{
  return new GeomTorusArc(*this);
}

void 
GeomTorusArc::move(const Point& center, const Vector& axis,
			double major_radius, double minor_radius,
			const Vector& zero,
			double start_angle, double arc_angle,
			int nu, int nv)
{
  center_ = center;
  axis_ = axis;
  major_radius_ = major_radius;
  minor_radius_ = minor_radius;
  nu_ = nu;
  nv_ = nv;
  zero_ = zero;
  start_angle_ = start_angle;
  arc_angle_ = arc_angle;
  adjust();
}


void 
GeomTorusArc::adjust()
{
  axis_.normalize();
  zero_.normalize();
  yaxis_ = Cross(axis_, zero_);
}

void 
GeomTorusArc::get_bounds(BBox& bb)
{
  bb.extend_disk(center_-axis_*minor_radius_, axis_, major_radius_+minor_radius_);
  bb.extend_disk(center_+axis_*minor_radius_, axis_, major_radius_+minor_radius_);
}

void
GeomTorusArc::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;
  glPushMatrix();
  glTranslated(center_.x(), center_.y(), center_.z());
  double matrix[16];
  matrix[0]=zero_.x(); matrix[1]=zero_.y(); matrix[2]=zero_.z(); matrix[3]=0;
  matrix[4]=yaxis_.x();matrix[5]=yaxis_.y();matrix[6]=yaxis_.z();matrix[7]=0;
  matrix[8]=axis_.x(); matrix[9]=axis_.y(); matrix[10]=axis_.z();matrix[11]=0;
  matrix[12]=0;        matrix[13]=0;        matrix[14]=0;        matrix[15]=1;
  glMultMatrixd(matrix);
  di->polycount_+=2*(nu_-1)*(nv_-1);

  // Draw the torus
  double a1=start_angle_;
  double a2=start_angle_-arc_angle_;
  if (a1 > a2)
  {
    double tmp=a1;
    a1=a2;
    a2=tmp;
  }
  SinCosTable tab1 (nu_, a1, a2);
  SinCosTable tab2 (nv_, 0, 2*M_PI, minor_radius_);
  SinCosTable tab2n(nv_, 0, 2*M_PI, minor_radius_);
  int u,v;
  switch(di->get_drawtype())
  {
  case DrawInfoOpenGL::WireFrame:
    {
      double srx=tab1.sin(0);
      double sry=tab1.cos(0);
      glBegin(GL_LINE_LOOP);
      for (v=1; v<nv_; v++)
      {
        double sz=tab2.cos(v);
        double srad=major_radius_+tab2.sin(v);
        double sx=srx*srad;
        double sy=sry*srad;
        glVertex3d(sx, sy, sz);
        glVertex3d(srx*major_radius_, sry*major_radius_, 0);
      }
      glEnd();

      srx=tab1.sin(nu_-1);
      sry=tab1.cos(nu_-1);
      glBegin(GL_LINE_LOOP);
      for (v=1; v<nv_; v++)
      {
        double sz=tab2.cos(v);
        double srad=major_radius_+tab2.sin(v);
        double sx=srx*srad;
        double sy=sry*srad;
        glVertex3d(sx, sy, sz);
        glVertex3d(srx*major_radius_, sry*major_radius_, 0);
      }
      glEnd();
        
      for (u=0; u<nu_; u++)
      {
        double rx=tab1.sin(u);
        double ry=tab1.cos(u);
        glBegin(GL_LINE_LOOP);
        for (v=1; v<nv_; v++)
        {
          double z=tab2.cos(v);
          double rad=major_radius_+tab2.sin(v);
          double x=rx*rad;
          double y=ry*rad;
          glVertex3d(x, y, z);
        }
        glEnd();
      }
      for (v=1; v<nv_; v++)
      {
        double z=tab2.cos(v);
        double rr=tab2.sin(v);
        glBegin(GL_LINE_LOOP);
        for (u=1; u<nu_; u++)
        {
          double rad=major_radius_+rr;
          double x=tab1.sin(u)*rad;
          double y=tab1.cos(u)*rad;
          glVertex3d(x, y, z);
        }
        glEnd();
      }
    }
    break;
  case DrawInfoOpenGL::Flat:
    for (v=0; v<nv_-1; v++)
    {
      double z1=tab2.cos(v);
      double rr1=tab2.sin(v);
      double z2=tab2.cos(v+1);
      double rr2=tab2.sin(v+1);
      glBegin(GL_TRIANGLE_STRIP);
      for (u=0; u<nu_; u++)
      {
        double r1=major_radius_+rr1;
        double r2=major_radius_+rr2;
        double xx=tab1.sin(u);
        double yy=tab1.cos(u);
        double x1=xx*r1;
        double y1=yy*r1;
        double x2=xx*r2;
        double y2=yy*r2;
        glVertex3d(x1, y1, z1);
        glVertex3d(x2, y2, z2);
      }
      glEnd();
    }
    break;
  case DrawInfoOpenGL::Gouraud:
    for (v=0; v<nv_-1; v++)
    {
      double z1=tab2.cos(v);
      double rr1=tab2.sin(v);
      double z2=tab2.cos(v+1);
      double rr2=tab2.sin(v+1);
      double nr=-tab2n.sin(v);
      double nz=-tab2n.cos(v);
      glBegin(GL_TRIANGLE_STRIP);
      for (u=0; u<nu_; u++)
      {
        double r1=major_radius_+rr1;
        double r2=major_radius_+rr2;
        double xx=tab1.sin(u);
        double yy=tab1.cos(u);
        double x1=xx*r1;
        double y1=yy*r1;
        double x2=xx*r2;
        double y2=yy*r2;
        glNormal3d(nr*xx, nr*yy, nz);
        glVertex3d(x1, y1, z1);
        glVertex3d(x2, y2, z2);
      }
      glEnd();
    }
    break;
  }
  glPopMatrix();
  post_draw(di);
}

} // End namespace SCIRun


