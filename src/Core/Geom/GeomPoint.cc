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
 * GeomPoint.cc: Points objects
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#include <Core/Util/Debug.h>
#include <Core/Geom/GeomPoint.h>
#include <Core/Geom/DrawInfoOpenGL.h>

#include <algorithm>

namespace SCIRun {

GeomPoints::GeomPoints(const GeomPoints &copy)
  : GeomObj(copy),
    point_size_(copy.point_size_), 
    points_(copy.points_), 
    pickable(copy.pickable)
{
  DEBUG_CONSTRUCTOR("GeomPoints")
}

GeomPoints::GeomPoints()
  : point_size_(1.0), pickable(false)
{
  DEBUG_CONSTRUCTOR("GeomPoints")
}

GeomPoints::~GeomPoints()
{
  DEBUG_DESTRUCTOR("GeomPoints")
}

GeomObj* 
GeomPoints::clone()
{
  return new GeomPoints(*this);
}

void
GeomPoints::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 0)) { return; }

  di->polycount_+=points_.size()/3;

  if( point_size_ > 0.0 )
  {
    glPointSize(point_size_);
  }

  if (di->pickmode_)
  {
    if (pickable)
    {
      glLoadName(0);
      float* p=&points_[0];
      for (unsigned int i=0; i<points_.size(); i+=3)
      {
        glLoadName(i/3);
        glBegin(GL_POINTS);
        glVertex3fv(p);
        glEnd();
        p+=3;
      }
    }
  }
  else
  {
    if (points_.size()) 
    {
      glVertexPointer(3, GL_FLOAT, 0, &(points_[0]));
      glEnableClientState(GL_VERTEX_ARRAY);
    }
    else 
    {
      glDisableClientState(GL_VERTEX_ARRAY);
    }
    if (colors_.size())
    {
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(colors_[0]));
      glEnableClientState(GL_COLOR_ARRAY);
    }
    else
    {
      glDisableClientState(GL_COLOR_ARRAY);
    }

    if (di->using_cmtexture_ && indices_.size() == points_.size() / 3 && indices_.size())
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

    glDrawArrays(GL_POINTS, 0, points_.size()/3);
  }

  glDisable(GL_TEXTURE_1D);

  // HACK set point size back to default
  // our scenegraph needs more graceful control of such state.
  glPointSize(di->point_size_);

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

void
GeomPoints::add(const Point &p, const MaterialHandle &m, unsigned int idx)
{
  add(p);
  
  colors_.push_back(COLOR_FTOB(m->diffuse.r()));
  colors_.push_back(COLOR_FTOB(m->diffuse.g()));
  colors_.push_back(COLOR_FTOB(m->diffuse.b()));
  colors_.push_back(COLOR_FTOB(m->transparency));
  item_idx_.push_back(idx);
}

void
GeomPoints::add(const Point& p, double index, unsigned int idx)
{
  add(p);
  indices_.push_back(index);
  item_idx_.push_back(idx);
}

void
GeomPoints::add(const Point& p, const MaterialHandle &m, double index,
		unsigned int idx)
{
  add(p, m, idx);
  indices_.push_back(index);
}

void
GeomPoints::get_bounds(BBox& bb)
{
  for (unsigned int i=0; i<points_.size(); i+=3)
  {
    bb.extend(Point(points_[i], points_[i+1], points_[i+2]));
  }
}

GeomTranspPoints::GeomTranspPoints()
  : GeomPoints(),
    xreverse_(false),
    yreverse_(false),
    zreverse_(false)
{
  DEBUG_CONSTRUCTOR("GeomTransPoints")
}


GeomTranspPoints::GeomTranspPoints(const GeomTranspPoints &copy)
  : GeomPoints(copy),
    xindices_(copy.xindices_),
    yindices_(copy.yindices_),
    zindices_(copy.zindices_),
    xreverse_(copy.xreverse_),
    yreverse_(copy.yreverse_),
    zreverse_(copy.zreverse_)
{
  DEBUG_CONSTRUCTOR("GeomTransPoints")
}

GeomTranspPoints::~GeomTranspPoints()
{
  DEBUG_DESTRUCTOR("GeomTransPoints")
}

GeomObj*
GeomTranspPoints::clone()
{
  return new GeomTranspPoints(*this);
}

void
GeomTranspPoints::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 0)) { return; }

  di->polycount_+=points_.size()/3;

  if( point_size_ > 0.0 )
  {
    glPointSize(point_size_);
  }

  Sort();
  get_view(di);

  std::vector<unsigned int> &clist =
    (di->axis_==0)?xindices_:((di->axis_==1)?yindices_:zindices_);

  bool &reverse =
    (di->axis_==0)?xreverse_:((di->axis_==1)?yreverse_:zreverse_);

  if (points_.size()) 
  {
    glVertexPointer(3, GL_FLOAT, 0, &(points_[0]));
    glEnableClientState(GL_VERTEX_ARRAY);
  }
  else 
  {
    glDisableClientState(GL_VERTEX_ARRAY);
  }

  if (colors_.size())
  {
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(colors_[0]));
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

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  if (di->dir_ ==  1 &&  reverse ||
      di->dir_ == -1 && !reverse)
  {
    std::reverse(clist.begin(), clist.end());
    reverse = !reverse;
  }

  if (clist.size())
    glDrawElements(GL_POINTS, clist.size(), GL_UNSIGNED_INT, &(clist[0]));

  glDisable(GL_BLEND);
  glDisable(GL_TEXTURE_1D);

  // HACK set point size back to default
  // our scenegraph needs more graceful control of such state.
  glPointSize(di->point_size_);

  post_draw(di);
}

static bool
pair_less(const std::pair<float, unsigned int> &a,
	  const std::pair<float, unsigned int> &b)
{
  return a.first < b.first;
}

void
GeomTranspPoints::Sort()
{
  const unsigned int vsize = points_.size() / 3;
  if (xindices_.size() == vsize) return;
    
  xreverse_ = false;
  yreverse_ = false;
  zreverse_ = false;
    
  std::vector<std::pair<float, unsigned int> > tmp(vsize);
  unsigned int i;
    
  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*3+0];
    tmp[i].second = i;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);
  
  xindices_.resize(vsize);
  for (i=0; i < vsize; i++)
  {
    xindices_[i] = tmp[i].second;
  }

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*3+1];
    tmp[i].second = i;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  yindices_.resize(vsize);
  for (i=0; i < vsize; i++)
  {
    yindices_[i] = tmp[i].second;
  }

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*3+2];
    tmp[i].second = i;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  zindices_.resize(vsize);
  for (i=0; i < vsize; i++)
  {
    zindices_[i] = tmp[i].second;
  }
}

} // End namespace SCIRun


