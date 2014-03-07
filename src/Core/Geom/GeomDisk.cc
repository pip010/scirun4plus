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
 *  GeomDisk.h:  Disk object
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#include <Core/Util/Debug.h>
#include <Core/Geom/GeomDisk.h>
#include <Core/Geom/GeomTri.h>
#include <Core/Geom/DrawInfoOpenGL.h>

#include <Core/Math/MiscMath.h>

namespace SCIRun {

GeomDisk::GeomDisk(int nu, int nv)
: GeomObj(), n(0,0,1), rad(1), nu(nu), nv(nv)
{
  DEBUG_CONSTRUCTOR("GeomDisk")
  adjust();
}

GeomDisk::GeomDisk(const Point& cen, const Vector& n,
		   double rad, int nu, int nv)
: GeomObj(), cen(cen), n(n), rad(rad), nu(nu), nv(nv)
{
  DEBUG_CONSTRUCTOR("GeomDisk")
  adjust();
}

void 
GeomDisk::move(const Point& _cen, const Vector& _n,
		    double _rad, int _nu, int _nv)
{
  cen=_cen;
  n=_n;
  rad=_rad;
  nu=_nu;
  nv=_nv;
  adjust();
}

GeomDisk::GeomDisk(const GeomDisk& copy)
  : GeomObj(), cen(copy.cen), n(copy.n),
    rad(copy.rad), nu(copy.nu), nv(copy.nv),
    v1(copy.v1), v2(copy.v2)
{
  DEBUG_CONSTRUCTOR("GeomDisk")
  adjust();
}

GeomDisk::~GeomDisk()
{
  DEBUG_DESTRUCTOR("GeomDisk")
}

GeomObj* 
GeomDisk::clone()
{
  return new GeomDisk(*this);
}

void 
GeomDisk::adjust()
{
  n.find_orthogonal(v1, v2);
  n.normalize();
  Vector z(0,0,1);
  if(Abs(n.y()) < 1.e-5)
  {
    // Only in x-z plane...
    zrotaxis=Vector(0,-1,0);
  } 
  else 
  {
    zrotaxis=Cross(n, z);
    zrotaxis.normalize();
  }
  double cangle=Dot(z, n);
  zrotangle=-acos(cangle);
}

void 
GeomDisk::get_bounds(BBox& bb)
{
  bb.extend_disk(cen, n, rad);
}

void
GeomDisk::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (rad < 1.e-6) return;
  if (!pre_draw(di, matl, 1)) return;
  glPushMatrix();
  glTranslated(cen.x(), cen.y(), cen.z());
  glRotated(RtoD(zrotangle), zrotaxis.x(), zrotaxis.y(), zrotaxis.z());
  di->polycount_+=2*(nu-1)*(nv-1);
  gluDisk(di->qobj_, 0, rad, nu, nv);
  glPopMatrix();
  post_draw(di);
}

} // End namespace SCIRun


