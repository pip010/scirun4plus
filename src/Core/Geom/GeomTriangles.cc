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
 *  GeomTriangles.cc: Triangle Strip object
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   June 1995
 *
 */

#include <Core/Util/Debug.h>

#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomTriangles.h>
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

GeomTriangles::GeomTriangles()
{
  DEBUG_CONSTRUCTOR("GeomTriangles")
}

GeomTriangles::GeomTriangles(const GeomTriangles& copy)
: GeomVertexPrim(copy)
{
  DEBUG_CONSTRUCTOR("GeomTriangles")
}

GeomTriangles::~GeomTriangles() 
{
  DEBUG_DESTRUCTOR("GeomTriangles")
}

void
GeomTriangles::add(const Point& p1, const Point& p2, const Point& p3)
{
  Vector n(Cross(p2-p1, p3-p1));
  normals.add(n);
  GeomVertexPrim::add(p1);
  GeomVertexPrim::add(p2);
  GeomVertexPrim::add(p3);
}

int
GeomTriangles::size(void)
{
  return verts.size();
}

void
GeomTriangles::add(const Point& p1, const Vector& v1,
                   const Point& p2, const Vector& v2,
                   const Point& p3, const Vector& v3)
{
  Vector n(Cross(p2-p1, p3-p1));
  normals.add(n);
  GeomVertexPrim::add(p1, v1);
  GeomVertexPrim::add(p2, v2);
  GeomVertexPrim::add(p3, v3);
}

void
GeomTriangles::add(const Point& p1, const MaterialHandle& m1,
                   const Point& p2, const MaterialHandle& m2,
                   const Point& p3, const MaterialHandle& m3)
{
  Vector n(Cross(p2-p1, p3-p1));
  normals.add(n);
  GeomVertexPrim::add(p1, m1);
  GeomVertexPrim::add(p2, m2);
  GeomVertexPrim::add(p3, m3);
}

void
GeomTriangles::add(const Point& p1, const Color& c1,
                   const Point& p2, const Color& c2,
                   const Point& p3, const Color& c3)
{
  Vector n(Cross(p2-p1, p3-p1));
  normals.add(n);
  GeomVertexPrim::add(p1, c1);
  GeomVertexPrim::add(p2, c2);
  GeomVertexPrim::add(p3, c3);
}

void
GeomTriangles::add(const Point& p1, const Vector& v1, 
                   const MaterialHandle& m1, const Point& p2, 
                   const Vector& v2, const MaterialHandle& m2,
                   const Point& p3, const Vector& v3, 
                   const MaterialHandle& m3)
{
  Vector n(Cross(p2-p1, p3-p1));
  normals.add(n);
  GeomVertexPrim::add(p1, v1, m1);
  GeomVertexPrim::add(p2, v2, m2);
  GeomVertexPrim::add(p3, v3, m3);
}

void
GeomTriangles::add(GeomVertex* v1, GeomVertex* v2, GeomVertex* v3)
{
  Vector n(Cross(v3->p - v1->p, v2->p - v1->p));
  normals.add(n);
  GeomVertexPrim::add(v1);
  GeomVertexPrim::add(v2);
  GeomVertexPrim::add(v3);
}

GeomObj*
GeomTriangles::clone()
{
  return new GeomTriangles(*this);
}

void
GeomTriangles::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;
  if (verts.size() <= 2)
    return;
  di->polycount_+=verts.size()/3;
  if (di->currently_lit_)
  {
    glEnable(GL_NORMALIZE);

    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
      {
        for (int i=0;i<verts.size();i+=3)
        {
          glBegin(GL_LINE_LOOP);
          glNormal3d(normals[i/3].x(), normals[i/3].y(),
                     normals[i/3].z());
          verts[i]->emit_all(di);
          verts[i+1]->emit_all(di);
          verts[i+2]->emit_all(di);
          glEnd();
        }
      }
      break;
    case DrawInfoOpenGL::Flat:
      {
        glBegin(GL_TRIANGLES);
        for (int i=0;i<verts.size();i+=3)
        {
          glNormal3d(normals[i/3].x(), normals[i/3].y(),
                     normals[i/3].z());
          verts[i]->emit_point(di);
          verts[i+1]->emit_point(di);
          verts[i+2]->emit_all(di);
        }
        glEnd();
      }
      break;
    case DrawInfoOpenGL::Gouraud:
      {
        glBegin(GL_TRIANGLES);
        for (int i=0;i<verts.size();i+=3)
        {
          glNormal3d(normals[i/3].x(), normals[i/3].y(),
                     normals[i/3].z());
          verts[i]->emit_all(di);
          verts[i+1]->emit_all(di);
          verts[i+2]->emit_all(di);
        }
        glEnd();
      }
      break;
    }
  }
  else
  {
    glDisable(GL_NORMALIZE);

    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
      {
        for (int i=0;i<verts.size();i+=3)
        {
          glBegin(GL_LINE_LOOP);
          verts[i]->emit_all(di);
          verts[i+1]->emit_all(di);
          verts[i+2]->emit_all(di);
          glEnd();
        }
      }
      break;
    case DrawInfoOpenGL::Flat:
      {
        glBegin(GL_TRIANGLES);
        for (int i=0;i<verts.size();i+=3)
        {
          verts[i]->emit_point(di);
          verts[i+1]->emit_point(di);
          verts[i+2]->emit_all(di);
        }
        glEnd();
      }
      break;
    case DrawInfoOpenGL::Gouraud:
      {
        glBegin(GL_TRIANGLES);
        for (int i=0;i<verts.size();i+=3)
        {
          verts[i]->emit_all(di);
          verts[i+1]->emit_all(di);
          verts[i+2]->emit_all(di);
        }
        glEnd();
      }
      break;
    }
    glEnable(GL_NORMALIZE);
  }
  post_draw(di);
}

GeomFastTriangles::GeomFastTriangles() :
  material_(0)
{
  DEBUG_CONSTRUCTOR("GeomFastTriangles")
}

GeomFastTriangles::GeomFastTriangles(const GeomFastTriangles& copy) :
  GeomObj(copy),
  points_(copy.points_),
  colors_(copy.colors_),
  normals_(copy.normals_),
  face_normals_(copy.face_normals_),
  material_(0)
{
  DEBUG_CONSTRUCTOR("GeomFastTriangles")
}

GeomFastTriangles::~GeomFastTriangles()
{
  DEBUG_DESTRUCTOR("GeomFastTriangles")
}

GeomObj*
GeomFastTriangles::clone()
{
  return new GeomFastTriangles(*this);
}

int
GeomFastTriangles::size()
{
  return points_.size() / 9;
}

void
GeomFastTriangles::get_bounds(BBox& bb)
{
  for(unsigned int i=0;i<points_.size();i+=3)
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
GeomFastTriangles::add(const Point &p0,
		       const Point &p1,
		       const Point &p2)
{
  Vector n(Cross(p1-p0, p2-p0));

  points_.push_back(p0.x());
  points_.push_back(p0.y());
  points_.push_back(p0.z());

  points_.push_back(p1.x());
  points_.push_back(p1.y());
  points_.push_back(p1.z());

  points_.push_back(p2.x());
  points_.push_back(p2.y());
  points_.push_back(p2.z());

  face_normals_.push_back(n.x());
  face_normals_.push_back(n.y());
  face_normals_.push_back(n.z());

  face_normals_.push_back(n.x());
  face_normals_.push_back(n.y());
  face_normals_.push_back(n.z());

  face_normals_.push_back(n.x());
  face_normals_.push_back(n.y());
  face_normals_.push_back(n.z());
}


void
GeomFastTriangles::add(const Point &p0, const MaterialHandle &m0,
		       const Point &p1, const MaterialHandle &m1,
		       const Point &p2, const MaterialHandle &m2)
{
  add(p0, p1, p2);

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

  material_ = m0;
}


void
GeomFastTriangles::add(const Point &p0, double i0,
		       const Point &p1, double i1,
		       const Point &p2, double i2)
{
  add(p0, p1, p2);

  indices_.push_back(i0);
  indices_.push_back(i1);
  indices_.push_back(i2);
}


void
GeomFastTriangles::add(const Point &p0,
		       const MaterialHandle &mat0, double i0,
		       const Point &p1,
		       const MaterialHandle &mat1, double i1,
		       const Point &p2,
		       const MaterialHandle &mat2, double i2)
{
  add(p0, mat0, p1, mat1, p2, mat2);

  indices_.push_back(i0);
  indices_.push_back(i1);
  indices_.push_back(i2);
}


void
GeomFastTriangles::add(const Point &p0, const Vector &n0,
		       const Point &p1, const Vector &n1,
		       const Point &p2, const Vector &n2)
{
  add(p0, p1, p2);

  normals_.push_back(n0.x());
  normals_.push_back(n0.y());
  normals_.push_back(n0.z());

  normals_.push_back(n1.x());
  normals_.push_back(n1.y());
  normals_.push_back(n1.z());

  normals_.push_back(n2.x());
  normals_.push_back(n2.y());
  normals_.push_back(n2.z());
}

void
GeomFastTriangles::add(const Point &p0, const Vector &n0,
		       const MaterialHandle &m0,
		       const Point &p1, const Vector &n1,
		       const MaterialHandle &m1,
		       const Point &p2, const Vector &n2,
		       const MaterialHandle &m2)
{
  add(p0, m0, p1, m1, p2, m2);

  normals_.push_back(n0.x());
  normals_.push_back(n0.y());
  normals_.push_back(n0.z());

  normals_.push_back(n1.x());
  normals_.push_back(n1.y());
  normals_.push_back(n1.z());

  normals_.push_back(n2.x());
  normals_.push_back(n2.y());
  normals_.push_back(n2.z());
}


void
GeomFastTriangles::add(const Point &p0, const Vector &n0, double i0,
		       const Point &p1, const Vector &n1, double i1,
		       const Point &p2, const Vector &n2, double i2)
{
  add(p0, i0, p1, i1, p2, i2);

  normals_.push_back(n0.x());
  normals_.push_back(n0.y());
  normals_.push_back(n0.z());

  normals_.push_back(n1.x());
  normals_.push_back(n1.y());
  normals_.push_back(n1.z());

  normals_.push_back(n2.x());
  normals_.push_back(n2.y());
  normals_.push_back(n2.z());
}


void
GeomFastTriangles::add(const Point &p0, const Vector &n0,
		       const MaterialHandle &mat0, double i0,
		       const Point &p1, const Vector &n1,
		       const MaterialHandle &mat1, double i1,
		       const Point &p2, const Vector &n2,
		       const MaterialHandle &mat2, double i2)
{
  add(p0, n0, mat0, p1, n1, mat1, p2, n2, mat2);

  indices_.push_back(i0);
  indices_.push_back(i1);
  indices_.push_back(i2);
}


void
GeomFastTriangles::add( const std::vector< std::vector< std::pair<Point, Vector> > >&
			tristrips )
{
  for(size_t i=0; i<tristrips.size(); i++ ) 
  {
    const std::vector< std::pair<Point, Vector> >& tristrip = tristrips[i];

    for(size_t j=0; j<tristrip.size()-2; j+=2 ) 
    {
      add(tristrip[j  ].first, tristrip[j  ].second,
        tristrip[j+1].first, tristrip[j+1].second,
        tristrip[j+2].first, tristrip[j+2].second);
      
      add( tristrip[j+3].first, tristrip[j+3].second,
        tristrip[j+2].first, tristrip[j+2].second,
        tristrip[j+1].first, tristrip[j+1].second );
    }
  }
}


void
GeomFastTriangles::add( const std::vector< std::vector< std::pair<Point, Vector> > >&
			tristrips,
			const MaterialHandle &mat )
{

  for(size_t i=0; i<tristrips.size(); i++ ) 
  {
    const std::vector< std::pair<Point, Vector> >& tristrip = tristrips[i];

    for(size_t j=0; j<tristrip.size()-2; j+=2 ) 
    {
      add( tristrip[j  ].first, tristrip[j  ].second, mat,
        tristrip[j+1].first, tristrip[j+1].second, mat,
        tristrip[j+2].first, tristrip[j+2].second, mat );
      
      add( tristrip[j+3].first, tristrip[j+3].second, mat,
        tristrip[j+2].first, tristrip[j+2].second, mat,
        tristrip[j+1].first, tristrip[j+1].second, mat );
    }
  }
}


void
GeomFastTriangles::add( const std::vector< std::vector< std::pair<Point, Vector> > >&
			tristrips,
			const MaterialHandle &mat,
			double cindex )
{
  for(size_t i=0; i<tristrips.size(); i++ ) 
  {
    const std::vector< std::pair<Point, Vector> >& tristrip = tristrips[i];

    for(size_t j=0; j<tristrip.size()-2; j+=2 ) 
    {
      add(tristrip[j  ].first, tristrip[j  ].second, mat, cindex,
        tristrip[j+1].first, tristrip[j+1].second, mat, cindex,
        tristrip[j+2].first, tristrip[j+2].second, mat, cindex);
      
      add( tristrip[j+3].first, tristrip[j+3].second, mat, cindex,
        tristrip[j+2].first, tristrip[j+2].second, mat, cindex,
        tristrip[j+1].first, tristrip[j+1].second, mat, cindex );
    }
  }
}


void
GeomFastTriangles::add( const std::vector< std::vector< std::pair<Point, Vector> > >&
			tristrips,
			double cindex )
{
  for(size_t i=0; i<tristrips.size(); i++ ) 
  {
    const std::vector<std::pair<Point, Vector> >& tristrip = tristrips[i];

    for(size_t j=0; j<tristrip.size()-2; j+=2 ) 
    {
      add(tristrip[j  ].first, tristrip[j  ].second, cindex,
        tristrip[j+1].first, tristrip[j+1].second, cindex,
        tristrip[j+2].first, tristrip[j+2].second, cindex);
      
      add( tristrip[j+3].first, tristrip[j+3].second, cindex,
        tristrip[j+2].first, tristrip[j+2].second, cindex,
        tristrip[j+1].first, tristrip[j+1].second, cindex );
    }
  }
}

void
GeomFastTriangles::drawVertexData(DrawInfoOpenGL* di)
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

  if (di->currently_lit_ && (normals_.size() || face_normals_.size()))
  {
    glEnable(GL_NORMALIZE);

    if (di->get_drawtype() == DrawInfoOpenGL::Flat ||
        normals_.size() < face_normals_.size())
    {
      glNormalPointer(GL_FLOAT, 0, &(face_normals_.front()));
    }
    else
    {
      glNormalPointer(GL_FLOAT, 0, &(normals_.front()));
    }

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
GeomFastTriangles::fbpick_draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_TEXTURE_3D);
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  glShadeModel(GL_FLAT);
  
  unsigned int nverts = points_.size();  // number of vertices.
  std::vector<unsigned char> cols(nverts * 4); // 4 bytes per vertex.
  for (unsigned int v = 0; v < nverts * 4; v+=4) 
  {
    unsigned int face_idx = v / 12;
    unsigned char r, g, b, a;
    
    unsigned int idx = face_idx + 1;
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
    
    cols[v] = r;
    cols[v+1] = g;
    cols[v+2] = b;
    cols[v+3] = a;
  }
  glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(cols.front()));
  glEnableClientState(GL_COLOR_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, &(points_.front()));
  glEnableClientState(GL_VERTEX_ARRAY);
  glDrawArrays(GL_TRIANGLES, 0, points_.size()/3);
}

void
GeomFastTriangles::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  drawVertexData( di );
  glDrawArrays(GL_TRIANGLES, 0, points_.size()/3);

  glDisableClientState(GL_NORMAL_ARRAY);
  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
  glDisable(GL_TEXTURE_1D);

  post_draw(di);
}

GeomFastTrianglesTwoSided::GeomFastTrianglesTwoSided() :
 material_(0), material2_(0)
{
  DEBUG_CONSTRUCTOR("GeomFastTrianglesTwoSided")
}

GeomFastTrianglesTwoSided::GeomFastTrianglesTwoSided(const GeomFastTrianglesTwoSided& copy)
  : GeomObj(copy),
    points_(copy.points_),
    colors_(copy.colors_),
    colors2_(copy.colors2_),
    indices_(copy.indices_),
    indices2_(copy.indices2_),
    normals_(copy.normals_),
    face_normals_(copy.face_normals_)
{
  DEBUG_CONSTRUCTOR("GeomFastTrianglesTwoSided")
}

GeomFastTrianglesTwoSided::~GeomFastTrianglesTwoSided()
{
  DEBUG_DESTRUCTOR("GeomFastTrianglesTwoSided")
}

GeomObj*
GeomFastTrianglesTwoSided::clone()
{
  return new GeomFastTrianglesTwoSided(*this);
}

int
GeomFastTrianglesTwoSided::size()
{
  return points_.size() / 9;
}

void
GeomFastTrianglesTwoSided::get_bounds(BBox& bb)
{
  for(unsigned int i=0;i<points_.size();i+=3)
  {
    bb.extend(Point(points_[i+0], points_[i+1], points_[i+2]));
  }
}

void
GeomFastTrianglesTwoSided::add(const Point &p0,
			       const Point &p1,
			       const Point &p2)
{
  Vector n(Cross(p1-p0, p2-p0));

  points_.push_back(p0.x());
  points_.push_back(p0.y());
  points_.push_back(p0.z());

  points_.push_back(p1.x());
  points_.push_back(p1.y());
  points_.push_back(p1.z());

  points_.push_back(p2.x());
  points_.push_back(p2.y());
  points_.push_back(p2.z());

  face_normals_.push_back(n.x());
  face_normals_.push_back(n.y());
  face_normals_.push_back(n.z());

  face_normals_.push_back(n.x());
  face_normals_.push_back(n.y());
  face_normals_.push_back(n.z());

  face_normals_.push_back(n.x());
  face_normals_.push_back(n.y());
  face_normals_.push_back(n.z());
}


void
GeomFastTrianglesTwoSided::add(const Point &p0, const Vector &n0,
			       const Point &p1, const Vector &n1,
			       const Point &p2, const Vector &n2)
{
  add(p0, p1, p2);

  normals_.push_back(n0.x());
  normals_.push_back(n0.y());
  normals_.push_back(n0.z());

  normals_.push_back(n1.x());
  normals_.push_back(n1.y());
  normals_.push_back(n1.z());

  normals_.push_back(n2.x());
  normals_.push_back(n2.y());
  normals_.push_back(n2.z());
}


void
GeomFastTrianglesTwoSided::add(const Point &p0, double i0, double j0,
			       const Point &p1, double i1, double j1,
			       const Point &p2, double i2, double j2)
{
  add(p0, p1, p2);

  indices_.push_back(i0);
  indices_.push_back(i1);
  indices_.push_back(i2);

  indices2_.push_back(j0);
  indices2_.push_back(j1);
  indices2_.push_back(j2);
}

void
GeomFastTrianglesTwoSided::add(const Point &p0, const MaterialHandle &m0, const MaterialHandle &k0,
			       const Point &p1, const MaterialHandle &m1, const MaterialHandle &k1,
			       const Point &p2, const MaterialHandle &m2, const MaterialHandle &k2)
{
  add(p0, p1, p2);

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

  material_ = m0;
  material2_ = k0;
}


void
GeomFastTrianglesTwoSided::add(const Point &p0, const Vector &n0,
			       const MaterialHandle &m0, const MaterialHandle &k0,
			       const Point &p1, const Vector &n1,
			       const MaterialHandle &m1, const MaterialHandle &k1,
			       const Point &p2, const Vector &n2,
			       const MaterialHandle &m2, const MaterialHandle &k2)
{
  add(p0, m0, k0, p1, m1, k1, p2, m2, k2);

  normals_.push_back(n0.x());
  normals_.push_back(n0.y());
  normals_.push_back(n0.z());

  normals_.push_back(n1.x());
  normals_.push_back(n1.y());
  normals_.push_back(n1.z());

  normals_.push_back(n2.x());
  normals_.push_back(n2.y());
  normals_.push_back(n2.z());
}


void
GeomFastTrianglesTwoSided::add(const Point &p0, const Vector &n0, double i0, double j0,
			       const Point &p1, const Vector &n1, double i1, double j1,
			       const Point &p2, const Vector &n2, double i2, double j2)
{
  add(p0, i0, j0, p1, i1, j1, p2, i2, j2);

  normals_.push_back(n0.x());
  normals_.push_back(n0.y());
  normals_.push_back(n0.z());

  normals_.push_back(n1.x());
  normals_.push_back(n1.y());
  normals_.push_back(n1.z());

  normals_.push_back(n2.x());
  normals_.push_back(n2.y());
  normals_.push_back(n2.z());
}


void
GeomFastTrianglesTwoSided::drawVertexData(DrawInfoOpenGL* di)
{
  di->polycount_ += size();

  glShadeModel(GL_FLAT);
  glDisable(GL_NORMALIZE);

  if (points_.size()) {
    glVertexPointer(3, GL_FLOAT, 0, &(points_.front()));
    glEnableClientState(GL_VERTEX_ARRAY);
  }
  else {
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

  if (di->currently_lit_ && (normals_.size() || face_normals_.size()))
  {
    glEnable(GL_NORMALIZE);

    if (di->get_drawtype() == DrawInfoOpenGL::Flat ||
        normals_.size() < face_normals_.size())
    {
      glNormalPointer(GL_FLOAT, 0, &(face_normals_.front()));
    }
    else
    {
      glNormalPointer(GL_FLOAT, 0, &(normals_.front()));
    }

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
GeomFastTrianglesTwoSided::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  drawVertexData( di );

  bool cullenable = false;
  if (glIsEnabled(GL_CULL_FACE)) cullenable = true;

  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glDrawArrays(GL_TRIANGLES, 0, points_.size()/3);
  
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

  if (di->using_cmtexture_ && indices_.size() == points_.size() / 3  &&
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

  glDrawArrays(GL_TRIANGLES, 0, points_.size()/3);

  if (!(cullenable)) glDisable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  glDisableClientState(GL_NORMAL_ARRAY);
  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
  glDisable(GL_TEXTURE_1D);

  post_draw(di);
}

GeomTranspTriangles::GeomTranspTriangles()
  : xreverse_(false),
    yreverse_(false),
    zreverse_(false)
{
  DEBUG_CONSTRUCTOR("GeomTranspTriangles")
}

GeomTranspTriangles::GeomTranspTriangles(const GeomTranspTriangles& copy)
  : GeomFastTriangles(copy),
    xlist_(copy.xlist_),
    ylist_(copy.ylist_),
    zlist_(copy.zlist_),
    xreverse_(copy.xreverse_),
    yreverse_(copy.yreverse_),
    zreverse_(copy.zreverse_)
{
  DEBUG_CONSTRUCTOR("GeomTranspTriangles")
}

GeomTranspTriangles::~GeomTranspTriangles()
{
  DEBUG_DESTRUCTOR("GeomTranspTriangles")
}


GeomObj* GeomTranspTriangles::clone()
{
    return new GeomTranspTriangles(*this);
}

void
GeomTranspTriangles::Sort()
{
  const unsigned int vsize = points_.size() / 9;
  if (xlist_.size() == vsize*3) return;

  xreverse_ = false;
  yreverse_ = false;
  zreverse_ = false;

  std::vector<std::pair<float, unsigned int> > tmp(vsize);
  unsigned int i;

  for (i=0; i<vsize; i++)
  {
    tmp[i].first = points_[i*9+0] + points_[i*9+3] + points_[i*9+6];
    tmp[i].second = i*3;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  xlist_.resize(vsize*3);
  for (i=0; i<vsize; i++)
  {
    xlist_[i*3+0] = tmp[i].second + 0;
    xlist_[i*3+1] = tmp[i].second + 1;
    xlist_[i*3+2] = tmp[i].second + 2;
  }

  for (i=0; i<vsize;i++)
  {
    tmp[i].first = points_[i*9+1] + points_[i*9+4] + points_[i*9+7];
    tmp[i].second = i*3;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  ylist_.resize(vsize*3);
  for (i=0; i<vsize; i++)
  {
    ylist_[i*3+0] = tmp[i].second + 0;
    ylist_[i*3+1] = tmp[i].second + 1;
    ylist_[i*3+2] = tmp[i].second + 2;
  }

  for (i=0; i<vsize; i++)
  {
    tmp[i].first = points_[i*9+2] + points_[i*9+5] + points_[i*9+8];
    tmp[i].second = i*3;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  zlist_.resize(vsize*3);
  for (i=0; i<vsize; i++)
  {
    zlist_[i*3+0] = tmp[i].second + 0;
    zlist_[i*3+1] = tmp[i].second + 1;
    zlist_[i*3+2] = tmp[i].second + 2;
  }
}


GeomTranspTrianglesTwoSided::GeomTranspTrianglesTwoSided()
  : xreverse_(false),
    yreverse_(false),
    zreverse_(false)
{
  DEBUG_CONSTRUCTOR("GeomTranspTrianglesTwoSided")
}


GeomTranspTrianglesTwoSided::GeomTranspTrianglesTwoSided(const GeomTranspTrianglesTwoSided& copy)
  : GeomFastTrianglesTwoSided(copy),
    xlist_(copy.xlist_),
    ylist_(copy.ylist_),
    zlist_(copy.zlist_),
    xreverse_(copy.xreverse_),
    yreverse_(copy.yreverse_),
    zreverse_(copy.zreverse_)
{
  DEBUG_CONSTRUCTOR("GeomTranspTrianglesTwoSided")
}

GeomTranspTrianglesTwoSided::~GeomTranspTrianglesTwoSided()
{
  DEBUG_DESTRUCTOR("GeomTranspTrianglesTwoSided")
}


GeomObj* GeomTranspTrianglesTwoSided::clone()
{
  return new GeomTranspTrianglesTwoSided(*this);
}

void
GeomTranspTrianglesTwoSided::Sort()
{
  const unsigned int vsize = points_.size() / 9;
  if (xlist_.size() == vsize*3) return;

  xreverse_ = false;
  yreverse_ = false;
  zreverse_ = false;

  std::vector<std::pair<float, unsigned int> > tmp(vsize);
  unsigned int i;

  for (i=0; i<vsize; i++)
  {
    tmp[i].first = points_[i*9+0] + points_[i*9+3] + points_[i*9+6];
    tmp[i].second = i*3;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  xlist_.resize(vsize*3);
  for (i=0; i<vsize; i++)
  {
    xlist_[i*3+0] = tmp[i].second + 0;
    xlist_[i*3+1] = tmp[i].second + 1;
    xlist_[i*3+2] = tmp[i].second + 2;
  }

  for (i=0; i<vsize;i++)
  {
    tmp[i].first = points_[i*9+1] + points_[i*9+4] + points_[i*9+7];
    tmp[i].second = i*3;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  ylist_.resize(vsize*3);
  for (i=0; i<vsize; i++)
  {
    ylist_[i*3+0] = tmp[i].second + 0;
    ylist_[i*3+1] = tmp[i].second + 1;
    ylist_[i*3+2] = tmp[i].second + 2;
  }

  for (i=0; i<vsize; i++)
  {
    tmp[i].first = points_[i*9+2] + points_[i*9+5] + points_[i*9+8];
    tmp[i].second = i*3;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  zlist_.resize(vsize*3);
  for (i=0; i<vsize; i++)
  {
    zlist_[i*3+0] = tmp[i].second + 0;
    zlist_[i*3+1] = tmp[i].second + 1;
    zlist_[i*3+2] = tmp[i].second + 2;
  }
}

void
GeomTranspTriangles::fbpick_draw(DrawInfoOpenGL* di, Material* matl, double)
{

  if (!pre_draw(di, matl, 1)) return;
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_TEXTURE_3D);
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  glShadeModel(GL_FLAT);
  unsigned int nverts = points_.size();  // number of vertices.
  std::vector<unsigned char> cols(nverts * 4); // 4 bytes per vertex.
  for (unsigned int v = 0; v < nverts * 4; v+=4) 
  {
    unsigned int face_idx = v / 12;
    unsigned char r, g, b, a;

    unsigned int idx = face_idx + 1;
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
    a = 255;    cols[v] = r;
    cols[v+1] = g;
    cols[v+2] = b;
    cols[v+3] = a;
  }
  glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(cols.front()));
  glEnableClientState(GL_COLOR_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, &(points_.front()));
  glEnableClientState(GL_VERTEX_ARRAY);
  glDrawArrays(GL_TRIANGLES, 0, points_.size()/3);
}

void
GeomTranspTriangles::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  drawVertexData( di );

  Sort();
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
    glDrawElements(GL_TRIANGLES, clist.size(), GL_UNSIGNED_INT, &(clist[0]));

  glFrontFace(GL_CCW);

  glDisableClientState(GL_NORMAL_ARRAY);
  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_BLEND);

  post_draw(di);
}

void
GeomTranspTrianglesTwoSided::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 1)) return;

  drawVertexData( di );

  Sort();
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
    glDrawElements(GL_TRIANGLES, clist.size(), GL_UNSIGNED_INT, &(clist[0]));

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

  if (di->using_cmtexture_ && indices_.size() == points_.size() / 3  &&
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
    glDrawElements(GL_TRIANGLES, clist.size(), GL_UNSIGNED_INT, &(clist[0]));

  if (!(cullenable)) glDisable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  glFrontFace(GL_CCW);

  glDisableClientState(GL_NORMAL_ARRAY);
  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_BLEND);

  post_draw(di);
}

GeomTrianglesPT1d::GeomTrianglesPT1d()
:GeomTrianglesP(),cmap(0)
{
  DEBUG_CONSTRUCTOR("GeomTrianglesPT1d")
}

GeomTrianglesPT1d::~GeomTrianglesPT1d()
{
  DEBUG_DESTRUCTOR("GeomTrianglesPT1d")
}

int 
GeomTrianglesPT1d::add(const Point& p1,const Point& p2,const Point& p3,
			   const float& f1,const float& f2,const float& f3)
{
  if (GeomTrianglesP::add(p1,p2,p3)) 
  {
    scalars_.add(f1);
    scalars_.add(f2);
    scalars_.add(f3);
    return 1;
  }
  return 0;
}


void
GeomTrianglesPT1d::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  //  return;
  if (!pre_draw(di,matl,1)) return;
  di->polycount_ += size();

  if (cmap)
  { // use 1D texturing.
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glEnable(GL_TEXTURE_1D);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glTexImage1D(GL_TEXTURE_1D,0,4,
                 256,0,GL_RGBA,GL_UNSIGNED_BYTE,
                 cmap);
    glColor4f(1,1,1,1);
    glMatrixMode(GL_TEXTURE);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glAlphaFunc(GL_GREATER,0.0); // exactly 0 means draw nothing.
    glEnable(GL_ALPHA_TEST);
  }
  else
  {
    std::cerr << "No color map!\n";
    return; // don't draw if no color map.
  }

  if (di->currently_lit_)
  {
    float *pts=&points_[0];
    float *nrmls=&normals_[0];
    float *sclrs=&scalars_[0];

    int niter=size();
    glBegin(GL_TRIANGLES);
    while (niter--)
    {
      glNormal3fv(nrmls);
      nrmls+=3;
      glTexCoord1fv(sclrs);
      sclrs++;
      glVertex3fv(pts);
      pts += 3;
      glVertex3fv(pts);
      pts+=3;
      glVertex3fv(pts);
      pts+=3;
    }
    glEnd();
  }
  else
  { // no normals.
    float *pts=&points_[0];
    float *sclrs=&scalars_[0];

    int niter=size();
    glBegin(GL_TRIANGLES);
    while (niter--)
    {
      glTexCoord1fv(sclrs);
      sclrs++;
      glVertex3fv(pts);
      pts += 3;
      glVertex3fv(pts);
      pts+=3;
      glVertex3fv(pts);
      pts+=3;
    }
    glEnd();

  }

  if (cmap)
  {
    glDisable(GL_ALPHA_TEST);
    glDisable(GL_TEXTURE_1D);
  }
  post_draw(di);
}

GeomTranspTrianglesP::GeomTranspTrianglesP()
  : alpha_(0.2), sorted_( false )
{
  DEBUG_CONSTRUCTOR("GeomTranspTrianglesP")
}

GeomTranspTrianglesP::GeomTranspTrianglesP(double aval)
  : alpha_(aval), sorted_( false )
{
  DEBUG_CONSTRUCTOR("GeomTranspTrianglesP")
}

GeomTranspTrianglesP::~GeomTranspTrianglesP()
{
  DEBUG_DESTRUCTOR("GeomTranspTrianglesP")
}

int
GeomTranspTrianglesP::vadd(const Point& p1, 
			   const Point& p2,
			   const Point& p3)
{
  if (add(p1,p2,p3)) 
  {
    const unsigned int index = xlist_.size();
    const Vector center = (p1.vector()+p2.vector()+p3.vector())*(1.0/3.0);
    xlist_.push_back(std::pair<float, unsigned int>(center.x(), index));
    ylist_.push_back(std::pair<float, unsigned int>(center.y(), index));
    zlist_.push_back(std::pair<float, unsigned int>(center.z(), index));
    return 1;
  } 
  return 0;
}


void
GeomTranspTrianglesP::Sort()
{
  if (!sorted_) 
  {
    std::sort(xlist_.begin(), xlist_.end(), pair_less);
    std::sort(ylist_.begin(), ylist_.end(), pair_less);
    std::sort(zlist_.begin(), zlist_.end(), pair_less);

    sorted_ = true;
  }
}


// grows points, normals and centers...
#if 0
void GeomTranspTrianglesP::MergeStuff(GeomTranspTrianglesP* other)
{
  points_.resize ( points_.size() + other->points_.size());
  normals_.resize(normals_.size() + other->normals_.size());

  xc.resize(xc.size() + other->xc.size());
  yc.resize(yc.size() + other->yc.size());
  zc.resize(zc.size() + other->zc.size());
  
  int start = points.size()-other->points.size();
  int i;
  for(i=0;i<other->points.size();i++) {
    points[i+start] = other->points[i];
  }
  start = normals.size() - other->normals.size();
  for(i=0;i<other->normals.size();i++) {
    normals[i+start] = other->normals[i];
  }

  start = xc.size() - other->xc.size();
  for(i=0;i<other->xc.size();i++) {
    xc[start +i] = other->xc[i];
    yc[start +i] = other->yc[i];
    zc[start +i] = other->zc[i];
  }
  other->points_.resize(0);
  other->normals_.resize(0);
  other->xc.resize(0);
  other->yc.resize(0);
  other->zc.resize(0);
}
#endif


void
GeomTranspTrianglesP::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!size())
  {
    return;
  }

  if (!pre_draw(di,matl,1)) return; // yes, this is lit.

  di->polycount_ += size();

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

  if (!has_color_)
  {
    glColor4f(0,1.0,0.0,alpha_);
  }
  else
  {
    glColor4f(r_,g_,b_,alpha_);
  }

  glDepthMask(GL_FALSE); // no zbuffering for now.
  glEnable(GL_NORMALIZE);

  Sort();
  get_view(di);

  const std::vector<std::pair<float, unsigned int> > &clist =
    (di->axis_==0)?xlist_:((di->axis_==1)?ylist_:zlist_);

  const int sort_dir = (di->dir_<0)?-1:1;
  const unsigned int sort_start = (sort_dir>0)?0:(clist.size()-1);
  unsigned int i;
  unsigned int ndone = 0;
  glBegin(GL_TRIANGLES);
  if (di->currently_lit_)
  {
    for (i = sort_start ; ndone < clist.size(); ndone++, i += sort_dir)
    {
      const unsigned int nindex = clist[i].second * 3;
      const unsigned int pindex = nindex * 3;
      glNormal3fv(&normals_[nindex]);
      glVertex3fv(&points_[pindex+0]);
      glVertex3fv(&points_[pindex+3]);
      glVertex3fv(&points_[pindex+6]);
    }
  }
  else
  {
    for (i = sort_start; ndone < clist.size(); ndone++, i += sort_dir)
    {
      const int pindex = clist[i].second * 9;

      glVertex3fv(&points_[pindex+0]);
      glVertex3fv(&points_[pindex+3]);
      glVertex3fv(&points_[pindex+6]);
    }
  }
  glEnd();

  glDepthMask(GL_TRUE); // turn zbuff back on.
  glDisable(GL_BLEND);
  post_draw(di);
}

GeomTrianglesP::GeomTrianglesP()
:has_color_(0)
{
  DEBUG_CONSTRUCTOR("GeomTrianglesP")
    // don't really need to do anythin...
}

GeomTrianglesP::~GeomTrianglesP()
{
  DEBUG_DESTRUCTOR("GeomTrianglesP")
}

void 
GeomTrianglesP::get_triangles( Array1<float> &v)
{
  int end = v.size();
  v.grow (points_.size());
  for (int i=0; i<points_.size(); i++)
    v[end+i] = points_[i];
}

int 
GeomTrianglesP::size(void)
{
  return points_.size()/9;
}

void 
GeomTrianglesP::reserve_clear(int n)
{
  points_.setsize(n*9);
  normals_.setsize(n*3);

  points_.remove_all();
  normals_.remove_all();
}

int 
GeomTrianglesP::add(const Point& p1, const Point& p2, const Point& p3)
{
  Vector n(Cross(p2-p1, p3-p1));

  int idx=normals_.size();
  normals_.grow(3);
  normals_[idx+0]=n.x();
  normals_[idx+1]=n.y();
  normals_[idx+2]=n.z();

  idx=points_.size();
  points_.grow(9);
  points_[idx+0]=p1.x();
  points_[idx+1]=p1.y();
  points_[idx+2]=p1.z();
  points_[idx+3]=p2.x();
  points_[idx+4]=p2.y();
  points_[idx+5]=p2.z();
  points_[idx+6]=p3.x();
  points_[idx+7]=p3.y();
  points_[idx+8]=p3.z();
  return 1;
}

// below is just a virtual function...
int 
GeomTrianglesP::vadd(const Point& p1, const Point& p2, const Point& p3)
{
  Vector n(Cross(p2-p1, p3-p1));

  int idx=normals_.size();
  normals_.grow(3);
  normals_[idx]=n.x();
  normals_[idx+1]=n.y();
  normals_[idx+2]=n.z();

  idx=points_.size();
  points_.grow(9);
  points_[idx]=p1.x();
  points_[idx+1]=p1.y();
  points_[idx+2]=p1.z();
  points_[idx+3]=p2.x();
  points_[idx+4]=p2.y();
  points_[idx+5]=p2.z();
  points_[idx+6]=p3.x();
  points_[idx+7]=p3.y();
  points_[idx+8]=p3.z();
  return 1;
}

GeomObj* GeomTrianglesP::clone()
{
  return new GeomTrianglesP(*this);
}

void GeomTrianglesP::get_bounds(BBox& box)
{
    for(int i=0;i<points_.size();i+=3)
box.extend(Point(points_[i],points_[i+1],points_[i+2]));
}

void
GeomTrianglesP::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (points_.size() == 0)
    return;

  // DAVE: Hack for 3d texture mapping
  if (!pre_draw(di,matl,1)) return;

  di->polycount_ += size();

  if (di->currently_lit_)
  {
    glEnable(GL_NORMALIZE);

    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
    case DrawInfoOpenGL::Flat:
    case DrawInfoOpenGL::Gouraud:
      { 
        float *pts   = &points_[0];
        float *nrmls = &normals_[0];
        int niter = size();
        glBegin(GL_TRIANGLES);
        while (niter--)
        {
          glNormal3fv(nrmls);
          nrmls+=3;
          glVertex3fv(pts);
          pts += 3;
          glVertex3fv(pts);
          pts+=3;
          glVertex3fv(pts);
          pts+=3;
        }
        glEnd();
      }
                
      break;
    }
  }
  else
  { // lights are off, don't emit the normals
    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
    case DrawInfoOpenGL::Flat:
    case DrawInfoOpenGL::Gouraud:
      { 
        float *pts = &points_[0];
        int niter = size();
        glBegin(GL_TRIANGLES);
        while (niter--)
        {
          glVertex3fv(pts);
          pts += 3;
          glVertex3fv(pts);
          pts+=3;
          glVertex3fv(pts);
          pts+=3;
        }
        glEnd();
      }
                
      break;
    }
  }
  post_draw(di);
}

GeomTrianglesPC::GeomTrianglesPC()
{
  DEBUG_CONSTRUCTOR("GeomTrianglesPC")
    // don't really need to do anythin...
}

GeomTrianglesPC::~GeomTrianglesPC()
{
  DEBUG_DESTRUCTOR("GeomTrianglesPC")
}

int 
GeomTrianglesPC::add(const Point& p1, const Color& c1,
			const Point& p2, const Color& c2,
			const Point& p3, const Color& c3)
{
  if (GeomTrianglesP::add(p1,p2,p3)) 
  {
    colors_.add(c1.r());
    colors_.add(c1.g());
    colors_.add(c1.b());

    colors_.add(c2.r());
    colors_.add(c2.g());
    colors_.add(c2.b());

    colors_.add(c3.r());
    colors_.add(c3.g());
    colors_.add(c3.b());
    return 1;
  }

  return 0;
}


void
GeomTrianglesPC::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (points_.size() == 0)
    return;
  if (!pre_draw(di,matl,1)) return;

  di->polycount_ += size();

  if (di->currently_lit_)
  {
    glEnable(GL_NORMALIZE);

    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
    case DrawInfoOpenGL::Flat:
    case DrawInfoOpenGL::Gouraud:
      { 
        float *pts = &points_[0];
        float *nrmls = &normals_[0];
        float *clrs = &colors_[0];
        int niter = size();
        glBegin(GL_TRIANGLES);
        while (niter--)
        {
          glNormal3fv(nrmls);
          nrmls+=3;

          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts += 3;

          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts+=3;

          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts+=3;
        }
        glEnd();
      }
      break;
    }
    glEnable(GL_NORMALIZE);
  }
  else
  { // lights are off, don't emit the normals
    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
    case DrawInfoOpenGL::Flat:
    case DrawInfoOpenGL::Gouraud:
      { 
        float *pts = &points_[0];
        float *clrs = &colors_[0];
        int niter = size();
        glBegin(GL_TRIANGLES);
        while (niter--)
        {
          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts += 3;

          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts+=3;

          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts+=3;
        }
        glEnd();
      }
      break;
    }
  }
  post_draw(di);
}

GeomTrianglesVP::GeomTrianglesVP()
{
  DEBUG_CONSTRUCTOR("GeomTrianglesVP")
    // don't really need to do anythin...
}

GeomTrianglesVP::~GeomTrianglesVP()
{
  DEBUG_DESTRUCTOR("GeomTrianglesVP")
}

int 
GeomTrianglesVP::size(void)
{
  return points_.size()/9;
}

void
GeomTrianglesVP::reserve_clear(int n)
{
  int np = points_.size()/9;
  int delta = n - np;

  points_.remove_all();
  normals_.remove_all();

  if (delta > 0) 
  {
    points_.grow(delta);
    normals_.grow(delta);
  }
	
}

int 
GeomTrianglesVP::add(const Point& p1, const Vector &v1,
			 const Point& p2, const Vector &v2,	
			 const Point& p3, const Vector &v3)
{
  int idx=normals_.size();
  normals_.grow(9);
  normals_[idx]=v1.x();
  normals_[idx+1]=v1.y();
  normals_[idx+2]=v1.z();
  normals_[idx+3]=v2.x();
  normals_[idx+4]=v2.y();
  normals_[idx+5]=v2.z();
  normals_[idx+6]=v3.x();
  normals_[idx+7]=v3.y();
  normals_[idx+8]=v3.z();

  idx=points_.size();
  points_.grow(9);
  points_[idx]=p1.x();
  points_[idx+1]=p1.y();
  points_[idx+2]=p1.z();
  points_[idx+3]=p2.x();
  points_[idx+4]=p2.y();
  points_[idx+5]=p2.z();
  points_[idx+6]=p3.x();
  points_[idx+7]=p3.y();
  points_[idx+8]=p3.z();
  return 1;
}

GeomObj* 
GeomTrianglesVP::clone()
{
  return new GeomTrianglesVP(*this);
}

void 
GeomTrianglesVP::get_bounds(BBox& box)
{
  for(int i=0;i<points_.size();i+=3)
	{
    box.extend(Point(points_[i],points_[i+1],points_[i+2]));
  }
}


void
GeomTrianglesVP::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (points_.size() == 0)
    return;
  if (!pre_draw(di,matl,1)) return;

  di->polycount_ += size();

  if (di->currently_lit_)
  {
    glEnable(GL_NORMALIZE);

    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
    case DrawInfoOpenGL::Flat:
    case DrawInfoOpenGL::Gouraud:
      { 
        float *pts = &points_[0];
        float *nrmls = &normals_[0];
        int niter = size();
        glBegin(GL_TRIANGLES);
        while (niter--)
        {
          glNormal3fv(nrmls);
          nrmls+=3;
          glVertex3fv(pts);
          pts += 3;
          glNormal3fv(nrmls);
          nrmls+=3;
          glVertex3fv(pts);
          pts+=3;
          glNormal3fv(nrmls);
          nrmls+=3;
          glVertex3fv(pts);
          pts+=3;
        }
        glEnd();
      }
                
      break;
    }
  }
  else { // lights are off, don't emit the normals
    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
    case DrawInfoOpenGL::Flat:
    case DrawInfoOpenGL::Gouraud:
      { 
        float *pts = &points_[0];
        int niter = size();
        glBegin(GL_TRIANGLES);
        while (niter--)
        {
          glVertex3fv(pts);
          pts += 3;
          glVertex3fv(pts);
          pts+=3;
          glVertex3fv(pts);
          pts+=3;
        }
        glEnd();
      }
      break;
    }
  }
  post_draw(di);
}

GeomTrianglesVPC::GeomTrianglesVPC()
{
  DEBUG_CONSTRUCTOR("GeomTrianglesVPC")

    // don't really need to do anythin...
}

GeomTrianglesVPC::~GeomTrianglesVPC()
{
  DEBUG_DESTRUCTOR("GeomTrianglesVPC")
}

int 
GeomTrianglesVPC::add(const Point& p1, const Vector &v1, const Color& c1,
			  const Point& p2, const Vector &v2, const Color& c2,
			  const Point& p3, const Vector &v3, const Color& c3)
{
  if (GeomTrianglesVP::add(p1,v1,p2,v2,p3,v3)) 
  {
    colors_.add(c1.r());
    colors_.add(c1.g());
    colors_.add(c1.b());

    colors_.add(c2.r());
    colors_.add(c2.g());
    colors_.add(c2.b());

    colors_.add(c3.r());
    colors_.add(c3.g());
    colors_.add(c3.b());
    return 1;
  }

  return 0;
}

void
GeomTrianglesVPC::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (points_.size() == 0)
    return;
  if (!pre_draw(di,matl,1)) return;

  di->polycount_ += size();

  if (di->currently_lit_)
  {
    glEnable(GL_NORMALIZE);

    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
    case DrawInfoOpenGL::Flat:
    case DrawInfoOpenGL::Gouraud:
      { 
        float *pts = &points_[0];
        float *nrmls = &normals_[0];
        float *clrs = &colors_[0];
        int niter = size();
        glBegin(GL_TRIANGLES);
        while (niter--)
        {
          glNormal3fv(nrmls);
          nrmls+=3;
          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts += 3;

          glNormal3fv(nrmls);
          nrmls+=3;
          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts+=3;

          glNormal3fv(nrmls);
          nrmls+=3;
          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts+=3;
        }
        glEnd();
      }
                
      break;
    }
  }
  else { // lights are off, don't emit the normals
    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
    case DrawInfoOpenGL::Flat:
    case DrawInfoOpenGL::Gouraud:
      { 
        float *pts = &points_[0];
        float *clrs = &colors_[0];
        int niter = size();
        glBegin(GL_TRIANGLES);
        while (niter--)
        {
          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts += 3;

          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts+=3;

          glColor3fv(clrs);
          clrs+=3;
          glVertex3fv(pts);
          pts+=3;
        }
        glEnd();
      }
      break;
    }
  }
  post_draw(di);
}

} // End namespace SCIRun

