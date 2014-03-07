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
 *  GeomQuads.cc: Fast Quads object
 *
 *  Written by:
 *   Michael Callahan
 *   Department of Computer Science
 *   University of Utah
 *   May 2003
 *
 */

#include <Core/Util/Debug.h>

#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomQuads.h>
#include <Core/Geom/GeomTri.h>

#include <functional>
#include <algorithm>

namespace SCIRun {

static bool
pair_less(const std::pair<double, unsigned int> &a,
	  const std::pair<double, unsigned int> &b)
{
  return a.first < b.first;
}

GeomFastQuads::GeomFastQuads()
  : material_(0)
{
  DEBUG_CONSTRUCTOR("GeomFastQuads")
}

GeomFastQuads::GeomFastQuads(const GeomFastQuads& copy)
  : GeomObj(copy), 
    points_(copy.points_),
    colors_(copy.colors_),
    normals_(copy.normals_),
    material_(0)
{
  DEBUG_CONSTRUCTOR("GeomFastQuads")
}

GeomFastQuads::~GeomFastQuads()
{
  DEBUG_DESTRUCTOR("GeomFastQuads")
}

GeomObj*
GeomFastQuads::clone()
{
  return new GeomFastQuads(*this);
}

int
GeomFastQuads::size()
{
  return points_.size() / 12;
}

void
GeomFastQuads::get_bounds(BBox& bb)
{
  for(size_t i=0;i<points_.size();i+=3)
  {
    bb.extend(Point(points_[i+0], points_[i+1], points_[i+2]));
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
GeomFastQuads::add(const Point &p0, const Point &p1,
		   const Point &p2, const Point &p3)
{
  // Assume planar, use first three to compute normal.
  Vector n(Cross(p1-p0, p2-p0));
  if (n.length2() < (p1-p0).length2()*1e-5) n = Cross(p2-p1,p3-p1);
  add(p0, n, p1, n, p2, n, p3, n);
}

void
GeomFastQuads::add(const Point &p0, const MaterialHandle &m0,
		   const Point &p1, const MaterialHandle &m1,
		   const Point &p2, const MaterialHandle &m2,
		   const Point &p3, const MaterialHandle &m3)
{
  // Assume planar, use first three to compute normal.
  Vector n(Cross(p1-p0, p2-p0));
  if (n.length2() < (p1-p0).length2()*1e-5) n = Cross(p2-p1,p3-p1);
  add(p0, n, m0, p1, n, m1, p2, n, m2, p3, n, m3);
}

void
GeomFastQuads::add(const Point &p0, double i0,
		   const Point &p1, double i1,
		   const Point &p2, double i2,
		   const Point &p3, double i3)
{
  // Assume planar, use first three to compute normal.
  Vector n(Cross(p1-p0, p2-p0));
  if (n.length2() < (p1-p0).length2()*1e-5) n = Cross(p2-p1,p3-p1);
  add(p0, n, i0, p1, n, i1, p2, n, i2, p3, n, i3);
}

void
GeomFastQuads::add(const Point &p0, const MaterialHandle &m0, double i0,
		   const Point &p1, const MaterialHandle &m1, double i1,
		   const Point &p2, const MaterialHandle &m2, double i2,
		   const Point &p3, const MaterialHandle &m3, double i3)
{
  // Assume planar, use first three to compute normal.
  Vector n(Cross(p1-p0, p2-p0));
  if (n.length2() < (p1-p0).length2()*1e-5) n = Cross(p2-p1,p3-p1);
  add(p0, n, m0, i0, p1, n, m1, i1, p2, n, m2, i2, p3, n, m3, i3);
}

void
GeomFastQuads::add(const Point &p0, const Vector &n0,
		   const Point &p1, const Vector &n1,
		   const Point &p2, const Vector &n2,
		   const Point &p3, const Vector &n3)
{
  points_.push_back(p0.x());
  points_.push_back(p0.y());
  points_.push_back(p0.z());

  points_.push_back(p1.x());
  points_.push_back(p1.y());
  points_.push_back(p1.z());

  points_.push_back(p2.x());
  points_.push_back(p2.y());
  points_.push_back(p2.z());

  points_.push_back(p3.x());
  points_.push_back(p3.y());
  points_.push_back(p3.z());


  normals_.push_back(n0.x());
  normals_.push_back(n0.y());
  normals_.push_back(n0.z());

  normals_.push_back(n1.x());
  normals_.push_back(n1.y());
  normals_.push_back(n1.z());

  normals_.push_back(n2.x());
  normals_.push_back(n2.y());
  normals_.push_back(n2.z());

  normals_.push_back(n3.x());
  normals_.push_back(n3.y());
  normals_.push_back(n3.z());
}

void
GeomFastQuads::add(const Point &p0, const Vector &n0,
		   const MaterialHandle &m0,
		   const Point &p1, const Vector &n1,
		   const MaterialHandle &m1,
		   const Point &p2, const Vector &n2,
		   const MaterialHandle &m2,
		   const Point &p3, const Vector &n3,
		   const MaterialHandle &m3)
{
  add(p0, n0, p1, n1, p2, n2, p3, n3);

  colors_.push_back(COLOR_FTOB(m0->diffuse.r()));
  colors_.push_back(COLOR_FTOB(m0->diffuse.g()));
  colors_.push_back(COLOR_FTOB(m0->diffuse.b()));
  colors_.push_back(COLOR_FTOB(m0->transparency));

  colors_.push_back(COLOR_FTOB(m1->diffuse.r()));
  colors_.push_back(COLOR_FTOB(m1->diffuse.g()));
  colors_.push_back(COLOR_FTOB(m1->diffuse.b()));
  colors_.push_back(COLOR_FTOB(m1->transparency));

  colors_.push_back(COLOR_FTOB(m2->diffuse.r()));
  colors_.push_back(COLOR_FTOB(m2->diffuse.g()));
  colors_.push_back(COLOR_FTOB(m2->diffuse.b()));
  colors_.push_back(COLOR_FTOB(m2->transparency));

  colors_.push_back(COLOR_FTOB(m3->diffuse.r()));
  colors_.push_back(COLOR_FTOB(m3->diffuse.g()));
  colors_.push_back(COLOR_FTOB(m3->diffuse.b()));
  colors_.push_back(COLOR_FTOB(m3->transparency));

  material_ = m0;
}

void
GeomFastQuads::add(const Point &p0, const Vector &n0, double i0,
		   const Point &p1, const Vector &n1, double i1,
		   const Point &p2, const Vector &n2, double i2,
		   const Point &p3, const Vector &n3, double i3)
{
  add(p0, n0, p1, n1, p2, n2, p3, n3);

  indices_.push_back(i0);
  indices_.push_back(i1);
  indices_.push_back(i2);
  indices_.push_back(i3);
}

void
GeomFastQuads::add(const Point &p0, const Vector &n0,
		   const MaterialHandle &mat0, double i0,
		   const Point &p1, const Vector &n1,
		   const MaterialHandle &mat1, double i1,
		   const Point &p2, const Vector &n2,
		   const MaterialHandle &mat2, double i2,
		   const Point &p3, const Vector &n3,
		   const MaterialHandle &mat3, double i3)
{
  add(p0, n0, mat0, p1, n1, mat1, p2, n2, mat2, p3, n3, mat3);

  indices_.push_back(i0);
  indices_.push_back(i1);
  indices_.push_back(i2);
  indices_.push_back(i3);
}

void
GeomFastQuads::drawVertexData(DrawInfoOpenGL* di)
{
  glShadeModel(GL_FLAT);
  glDisable(GL_NORMALIZE);

  if (points_.size()) 
  {
    glVertexPointer(3, GL_FLOAT, 0, &(points_.front()));
    glEnableClientState(GL_VERTEX_ARRAY);
  }
  else 
  {
    glDisableClientState(GL_VERTEX_ARRAY);
  }

  if (colors_.size())
  {
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(colors_.front()));
    glEnableClientState(GL_COLOR_ARRAY);
  }
  else
  {
    glDisableClientState(GL_COLOR_ARRAY);
  }

  if (di->using_cmtexture_ && indices_.size() == points_.size() / 3 &&
      indices_.size())
  {
    glTexCoordPointer(1, GL_FLOAT, 0, &(indices_[0]));
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glColor3d(di->diffuse_scale_, di->diffuse_scale_, di->diffuse_scale_);

    glEnable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_1D, di->cmtexture_);
  }
  else
  {
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  }

  if (di->currently_lit_ && normals_.size())
  {
    glEnable(GL_NORMALIZE);

    glNormalPointer(GL_FLOAT, 0, &(normals_.front()));

    glEnableClientState(GL_NORMAL_ARRAY);
    if (di->get_drawtype() != DrawInfoOpenGL::Flat)
    {
      glShadeModel(GL_SMOOTH);
    }
  }
  else
  {
    glDisableClientState(GL_NORMAL_ARRAY);
  }

  if (material_.get_rep()) { di->set_material(material_.get_rep()); }
}

void
GeomFastQuads::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  drawVertexData( di );
  glDrawArrays(GL_QUADS, 0, points_.size()/3);

  glDisableClientState(GL_NORMAL_ARRAY);
  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
  glDisable(GL_TEXTURE_1D);

  post_draw(di);
}

GeomTranspQuads::GeomTranspQuads()
  : xreverse_(false),
    yreverse_(false),
    zreverse_(false)
{
  DEBUG_CONSTRUCTOR("GeomTranspQuads")
}

GeomTranspQuads::GeomTranspQuads(const GeomTranspQuads& copy)
  : GeomFastQuads(copy),
    xlist_(copy.xlist_),
    ylist_(copy.ylist_),
    zlist_(copy.zlist_),
    xreverse_(copy.xreverse_),
    yreverse_(copy.yreverse_),
    zreverse_(copy.zreverse_)
{
  DEBUG_CONSTRUCTOR("GeomTranspQuads")
}

GeomTranspQuads::~GeomTranspQuads()
{
  DEBUG_DESTRUCTOR("GeomTranspQuads")
}

GeomObj* 
GeomTranspQuads::clone()
{
  return new GeomTranspQuads(*this);
}

void
GeomTranspQuads::SortPolys()
{
  const unsigned int vsize = points_.size() / 12;
  if (xlist_.size() == vsize*4) return;

  xreverse_ = false;
  yreverse_ = false;
  zreverse_ = false;

  std::vector<std::pair<float, unsigned int> > tmp(vsize);
  unsigned int i;

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*12+0] + points_[i*12+3] +
      points_[i*12+6] + points_[i*12+9];
    tmp[i].second = i*4;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  xlist_.resize(vsize*4);
  for (i=0; i < vsize; i++)
  {
    xlist_[i*4+0] = tmp[i].second + 0;
    xlist_[i*4+1] = tmp[i].second + 1;
    xlist_[i*4+2] = tmp[i].second + 2;
    xlist_[i*4+3] = tmp[i].second + 3;
  }

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*12+1] + points_[i*12+4] +
      points_[i*12+7] + points_[i*12+10];
    tmp[i].second = i*4;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  ylist_.resize(vsize*4);
  for (i=0; i < vsize; i++)
  {
    ylist_[i*4+0] = tmp[i].second + 0;
    ylist_[i*4+1] = tmp[i].second + 1;
    ylist_[i*4+2] = tmp[i].second + 2;
    ylist_[i*4+3] = tmp[i].second + 3;
  }

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*12+2] + points_[i*12+5] +
      points_[i*12+8] + points_[i*12+11];
    tmp[i].second = i*4;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  zlist_.resize(vsize*4);
  for (i=0; i < vsize; i++)
  {
    zlist_[i*4+0] = tmp[i].second + 0;
    zlist_[i*4+1] = tmp[i].second + 1;
    zlist_[i*4+2] = tmp[i].second + 2;
    zlist_[i*4+3] = tmp[i].second + 3;
  }
}

void
GeomTranspQuads::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  drawVertexData( di );

  SortPolys();
  get_view(di);

  std::vector<unsigned int> &clist =
    (di->axis_==0)?xlist_:((di->axis_==1)?ylist_:zlist_);

  bool &reverse =
    (di->axis_==0)?xreverse_:((di->axis_==1)?yreverse_:zreverse_);

  if (di->dir_ ==  1 &&  reverse ||
      di->dir_ == -1 && !reverse)
  {
    std::reverse(clist.begin(), clist.end());
    reverse = !reverse;
  }

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glFrontFace(reverse?GL_CW:GL_CCW);

  if (clist.size())
    glDrawElements(GL_QUADS, clist.size(), GL_UNSIGNED_INT, &(clist.front()));

  glFrontFace(GL_CCW);

  glDisableClientState(GL_NORMAL_ARRAY);
  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_BLEND);

  post_draw(di);
}

GeomFastQuadsTwoSided::GeomFastQuadsTwoSided()
  : material_(0), material2_(0)
{
  DEBUG_CONSTRUCTOR("GeomFastQuadsTwoSided")
}

GeomFastQuadsTwoSided::GeomFastQuadsTwoSided(const GeomFastQuadsTwoSided& copy)
  : GeomObj(copy),
    points_(copy.points_),
    colors_(copy.colors_),
    colors2_(copy.colors2_),
    normals_(copy.normals_),
    material_(copy.material_),
    material2_(copy.material2_)
{
  DEBUG_CONSTRUCTOR("GeomFastQuadsTwoSided")
}

GeomFastQuadsTwoSided::~GeomFastQuadsTwoSided()
{
  DEBUG_DESTRUCTOR("GeomFastQuadsTwoSided")
}

GeomObj*
GeomFastQuadsTwoSided::clone()
{
  return new GeomFastQuadsTwoSided(*this);
}

int
GeomFastQuadsTwoSided::size()
{
  return points_.size() / 12;
}

void
GeomFastQuadsTwoSided::get_bounds(BBox& bb)
{
  for(unsigned int i=0;i<points_.size();i+=3)
  {
    bb.extend(Point(points_[i+0], points_[i+1], points_[i+2]));
  }
}

void
GeomFastQuadsTwoSided::add(const Point &p0, const Point &p1,
			   const Point &p2, const Point &p3)
{
  // Assume planar, use first three to compute normal.
  Vector n(Cross(p1-p0, p2-p0));
  add(p0, n, p1, n, p2, n, p3, n);
}

void
GeomFastQuadsTwoSided::add(const Point &p0, const MaterialHandle &m0, const MaterialHandle &k0,
			   const Point &p1, const MaterialHandle &m1, const MaterialHandle &k1,
			   const Point &p2, const MaterialHandle &m2, const MaterialHandle &k2,
			   const Point &p3, const MaterialHandle &m3, const MaterialHandle &k3)
{
  // Assume planar, use first three to compute normal.
  Vector n(Cross(p1-p0, p2-p0));
  add(p0, n, m0, k0, p1, n, m1, k1, p2, n, m2, k2, p3, n, m3, k3);
}

void
GeomFastQuadsTwoSided::add(const Point &p0, double i0, double j0,
			   const Point &p1, double i1, double j1,
			   const Point &p2, double i2, double j2,
			   const Point &p3, double i3, double j3)
{
  // Assume planar, use first three to compute normal.
  Vector n(Cross(p1-p0, p2-p0));
  add(p0, n, i0, j0, p1, n, i1, j1, p2, n, i2, j2, p3, n, i3, j3);
}

void
GeomFastQuadsTwoSided::add(const Point &p0,
			   const MaterialHandle &m0, double i0,
			   const MaterialHandle &k0, double j0,
			   const Point &p1,
			   const MaterialHandle &m1, double i1,
			   const MaterialHandle &k1, double j1,
			   const Point &p2,
			   const MaterialHandle &m2, double i2,
			   const MaterialHandle &k2, double j2,
			   const Point &p3,
			   const MaterialHandle &m3, double i3,
			   const MaterialHandle &k3, double j3 )
{
  // Assume planar, use first three to compute normal.
  Vector n(Cross(p1-p0, p2-p0));
  add(p0, n, m0, i0, k0, j0, p1, n, m1, i1, k1, j1,
      p2, n, m2, i2, k2, j2, p3, n, m3, i3, k3, j3);
}

void
GeomFastQuadsTwoSided::add(const Point &p0, const Vector &n0,
                           const Point &p1, const Vector &n1,
                           const Point &p2, const Vector &n2,
                           const Point &p3, const Vector &n3)
{
  points_.push_back(p0.x());
  points_.push_back(p0.y());
  points_.push_back(p0.z());

  points_.push_back(p1.x());
  points_.push_back(p1.y());
  points_.push_back(p1.z());

  points_.push_back(p2.x());
  points_.push_back(p2.y());
  points_.push_back(p2.z());

  points_.push_back(p3.x());
  points_.push_back(p3.y());
  points_.push_back(p3.z());


  normals_.push_back(n0.x());
  normals_.push_back(n0.y());
  normals_.push_back(n0.z());

  normals_.push_back(n1.x());
  normals_.push_back(n1.y());
  normals_.push_back(n1.z());

  normals_.push_back(n2.x());
  normals_.push_back(n2.y());
  normals_.push_back(n2.z());

  normals_.push_back(n3.x());
  normals_.push_back(n3.y());
  normals_.push_back(n3.z());
}

void
GeomFastQuadsTwoSided::add(const Point &p0, const Vector &n0,
			   const MaterialHandle &m0, const MaterialHandle &k0,
			   const Point &p1, const Vector &n1,
			   const MaterialHandle &m1, const MaterialHandle &k1,
			   const Point &p2, const Vector &n2,
			   const MaterialHandle &m2, const MaterialHandle &k2,
			   const Point &p3, const Vector &n3,
			   const MaterialHandle &m3, const MaterialHandle &k3)
{
  add(p0, n0, p1, n1, p2, n2, p3, n3);

  colors_.push_back(COLOR_FTOB(m0->diffuse.r()));
  colors_.push_back(COLOR_FTOB(m0->diffuse.g()));
  colors_.push_back(COLOR_FTOB(m0->diffuse.b()));
  colors_.push_back(COLOR_FTOB(m0->transparency));

  colors_.push_back(COLOR_FTOB(m1->diffuse.r()));
  colors_.push_back(COLOR_FTOB(m1->diffuse.g()));
  colors_.push_back(COLOR_FTOB(m1->diffuse.b()));
  colors_.push_back(COLOR_FTOB(m1->transparency));

  colors_.push_back(COLOR_FTOB(m2->diffuse.r()));
  colors_.push_back(COLOR_FTOB(m2->diffuse.g()));
  colors_.push_back(COLOR_FTOB(m2->diffuse.b()));
  colors_.push_back(COLOR_FTOB(m2->transparency));

  colors_.push_back(COLOR_FTOB(m3->diffuse.r()));
  colors_.push_back(COLOR_FTOB(m3->diffuse.g()));
  colors_.push_back(COLOR_FTOB(m3->diffuse.b()));
  colors_.push_back(COLOR_FTOB(m3->transparency));

  colors2_.push_back(COLOR_FTOB(k0->diffuse.r()));
  colors2_.push_back(COLOR_FTOB(k0->diffuse.g()));
  colors2_.push_back(COLOR_FTOB(k0->diffuse.b()));
  colors2_.push_back(COLOR_FTOB(k0->transparency));

  colors2_.push_back(COLOR_FTOB(k1->diffuse.r()));
  colors2_.push_back(COLOR_FTOB(k1->diffuse.g()));
  colors2_.push_back(COLOR_FTOB(k1->diffuse.b()));
  colors2_.push_back(COLOR_FTOB(k1->transparency));

  colors2_.push_back(COLOR_FTOB(k2->diffuse.r()));
  colors2_.push_back(COLOR_FTOB(k2->diffuse.g()));
  colors2_.push_back(COLOR_FTOB(k2->diffuse.b()));
  colors2_.push_back(COLOR_FTOB(k2->transparency));

  colors2_.push_back(COLOR_FTOB(k3->diffuse.r()));
  colors2_.push_back(COLOR_FTOB(k3->diffuse.g()));
  colors2_.push_back(COLOR_FTOB(k3->diffuse.b()));
  colors2_.push_back(COLOR_FTOB(k3->transparency));

  material_  = m0;
  material2_ = k0;
}

void
GeomFastQuadsTwoSided::add(const Point &p0, const Vector &n0, double i0, double j0,
			   const Point &p1, const Vector &n1, double i1, double j1,
			   const Point &p2, const Vector &n2, double i2, double j2,
			   const Point &p3, const Vector &n3, double i3, double j3)
{
  add(p0, n0, p1, n1, p2, n2, p3, n3);

  indices_.push_back(i0);
  indices_.push_back(i1);
  indices_.push_back(i2);
  indices_.push_back(i3);

  indices2_.push_back(j0);
  indices2_.push_back(j1);
  indices2_.push_back(j2);
  indices2_.push_back(j3);
}

void
GeomFastQuadsTwoSided::add(const Point &p0, const Vector &n0,
			   const MaterialHandle &m0, double i0,
			   const MaterialHandle &k0, double j0,
			   const Point &p1, const Vector &n1,
			   const MaterialHandle &m1, double i1,
			   const MaterialHandle &k1, double j1,
			   const Point &p2, const Vector &n2,
			   const MaterialHandle &m2, double i2,
			   const MaterialHandle &k2, double j2,
			   const Point &p3, const Vector &n3,
			   const MaterialHandle &m3, double i3,
			   const MaterialHandle &k3, double j3)
{
  add(p0, n0, m0, k0, p1, n1, m1, k1, p2, n2, m2, k2, p3, n3, m3, k3);

  indices_.push_back(i0);
  indices_.push_back(i1);
  indices_.push_back(i2);
  indices_.push_back(i3);

  indices2_.push_back(j0);
  indices2_.push_back(j1);
  indices2_.push_back(j2);
  indices2_.push_back(j3);
}

void
GeomFastQuadsTwoSided::drawVertexData(DrawInfoOpenGL* di)
{
  di->polycount_ += size();

  glShadeModel(GL_FLAT);
  glDisable(GL_NORMALIZE);

  if (points_.size()) 
  {
    glVertexPointer(3, GL_FLOAT, 0, &(points_.front()));
    glEnableClientState(GL_VERTEX_ARRAY);
  }
  else 
  {
    glDisableClientState(GL_VERTEX_ARRAY);
  }

  if (colors_.size())
  {
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(colors_.front()));
    glEnableClientState(GL_COLOR_ARRAY);
  }
  else
  {
    glDisableClientState(GL_COLOR_ARRAY);
  }

  if (di->using_cmtexture_ && indices_.size() == points_.size() / 3 &&
      indices_.size())
  {
    glTexCoordPointer(1, GL_FLOAT, 0, &(indices_[0]));
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glColor3d(di->diffuse_scale_, di->diffuse_scale_, di->diffuse_scale_);

    glEnable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_1D, di->cmtexture_);
  }
  else
  {
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  }

  if (di->currently_lit_ && normals_.size())
  {
    glEnable(GL_NORMALIZE);

    glNormalPointer(GL_FLOAT, 0, &(normals_.front()));

    glEnableClientState(GL_NORMAL_ARRAY);
    if (di->get_drawtype() != DrawInfoOpenGL::Flat)
    {
      glShadeModel(GL_SMOOTH);
    }
  }
  else
  {
    glDisableClientState(GL_NORMAL_ARRAY);
  }

  if (material_.get_rep()) { di->set_material(material_.get_rep()); }
}
  
void
GeomFastQuadsTwoSided::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  drawVertexData( di );

  bool cullenable = false;
  if (glIsEnabled(GL_CULL_FACE)) cullenable = true;

  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glDrawArrays(GL_QUADS, 0, points_.size()/3);

  glCullFace(GL_FRONT);
 
  if (material2_.get_rep()) { di->set_material(material2_.get_rep()); }

  if (colors2_.size())
  {
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(colors2_.front()));
    glEnableClientState(GL_COLOR_ARRAY);
  }
  else
  {
    glDisableClientState(GL_COLOR_ARRAY);
  }
  
  if (di->using_cmtexture_ && indices_.size() == points_.size() / 3 &&
      indices2_.size())
  {
    glTexCoordPointer(1, GL_FLOAT, 0, &(indices2_[0]));
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glColor3d(di->diffuse_scale_, di->diffuse_scale_, di->diffuse_scale_);

    glEnable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_1D, di->cmtexture_);
  }
  else
  {
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  }
  
  glDrawArrays(GL_QUADS, 0, points_.size()/3);
    
  if(!(cullenable)) glDisable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  
  glDisableClientState(GL_NORMAL_ARRAY);
  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
  glDisable(GL_TEXTURE_1D);

  post_draw(di);
}

GeomTranspQuadsTwoSided::GeomTranspQuadsTwoSided()
  : xreverse_(false),
    yreverse_(false),
    zreverse_(false)
{
}

GeomTranspQuadsTwoSided::GeomTranspQuadsTwoSided(const GeomTranspQuadsTwoSided& copy)
  : GeomFastQuadsTwoSided(copy),
    xlist_(copy.xlist_),
    ylist_(copy.ylist_),
    zlist_(copy.zlist_),
    xreverse_(copy.xreverse_),
    yreverse_(copy.yreverse_),
    zreverse_(copy.zreverse_)
{
}

GeomTranspQuadsTwoSided::~GeomTranspQuadsTwoSided()
{
}


GeomObj* GeomTranspQuadsTwoSided::clone()
{
  return new GeomTranspQuadsTwoSided(*this);
}

void
GeomTranspQuadsTwoSided::SortPolys()
{
  const unsigned int vsize = points_.size() / 12;
  if (xlist_.size() == vsize*4) return;

  xreverse_ = false;
  yreverse_ = false;
  zreverse_ = false;

  std::vector<std::pair<float, unsigned int> > tmp(vsize);
  unsigned int i;

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*12+0] + points_[i*12+3] +
      points_[i*12+6] + points_[i*12+9];
    tmp[i].second = i*4;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  xlist_.resize(vsize*4);
  for (i=0; i < vsize; i++)
  {
    xlist_[i*4+0] = tmp[i].second + 0;
    xlist_[i*4+1] = tmp[i].second + 1;
    xlist_[i*4+2] = tmp[i].second + 2;
    xlist_[i*4+3] = tmp[i].second + 3;
  }

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*12+1] + points_[i*12+4] +
      points_[i*12+7] + points_[i*12+10];
    tmp[i].second = i*4;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  ylist_.resize(vsize*4);
  for (i=0; i < vsize; i++)
  {
    ylist_[i*4+0] = tmp[i].second + 0;
    ylist_[i*4+1] = tmp[i].second + 1;
    ylist_[i*4+2] = tmp[i].second + 2;
    ylist_[i*4+3] = tmp[i].second + 3;
  }

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*12+2] + points_[i*12+5] +
      points_[i*12+8] + points_[i*12+11];
    tmp[i].second = i*4;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  zlist_.resize(vsize*4);
  for (i=0; i < vsize; i++)
  {
    zlist_[i*4+0] = tmp[i].second + 0;
    zlist_[i*4+1] = tmp[i].second + 1;
    zlist_[i*4+2] = tmp[i].second + 2;
    zlist_[i*4+3] = tmp[i].second + 3;
  }
}

void
GeomTranspQuadsTwoSided::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  drawVertexData( di );

  SortPolys();
  get_view(di);

  std::vector<unsigned int> &clist =
    (di->axis_==0)?xlist_:((di->axis_==1)?ylist_:zlist_);

  bool &reverse =
    (di->axis_==0)?xreverse_:((di->axis_==1)?yreverse_:zreverse_);

  if (di->dir_ ==  1 &&  reverse ||
      di->dir_ == -1 && !reverse)
  {
    std::reverse(clist.begin(), clist.end());
    reverse = !reverse;
  }

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glFrontFace(reverse?GL_CW:GL_CCW);

  bool cullenable = false;
  if (glIsEnabled(GL_CULL_FACE)) cullenable = true;

  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  
  if (clist.size())
    glDrawElements(GL_QUADS, clist.size(), GL_UNSIGNED_INT, &(clist.front()));

  glCullFace(GL_FRONT);
 
  if (material2_.get_rep()) { di->set_material(material2_.get_rep()); }

  if (colors2_.size())
  {
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(colors2_.front()));
    glEnableClientState(GL_COLOR_ARRAY);
  }
  else
  {
    glDisableClientState(GL_COLOR_ARRAY);
  }
  
  if (di->using_cmtexture_ && indices_.size() == points_.size() / 3 &&
      indices2_.size())
  {
    glTexCoordPointer(1, GL_FLOAT, 0, &(indices2_[0]));
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glColor3d(di->diffuse_scale_, di->diffuse_scale_, di->diffuse_scale_);

    glEnable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_1D, di->cmtexture_);
  }
  else
  {
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  }
  
  if (clist.size())
    glDrawElements(GL_QUADS, clist.size(), GL_UNSIGNED_INT, &(clist.front()));

  if(!(cullenable)) glDisable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  
  glFrontFace(GL_CCW);
    
  glDisableClientState(GL_NORMAL_ARRAY);
  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_BLEND);

  post_draw(di);
}

} // End namespace SCIRun

