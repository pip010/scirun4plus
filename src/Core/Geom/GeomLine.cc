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
 *  GeomLine.cc:  Line object
 *
 *  Written by:
 *   Steven G. Parker & David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#ifdef _WIN32
#pragma warning(disable:4291) // quiet the visual C++ compiler
#endif

#include <Core/Util/Debug.h>
#include <Core/Geom/GeomLine.h>
#include <Core/Geom/DrawInfoOpenGL.h>

#include <algorithm>
#include <stdlib.h>

namespace SCIRun {

GeomLine::GeomLine(const Point& p1, const Point& p2) : 
  GeomObj(), 
  p1(p1), 
  p2(p2),
  line_width_(0.0)
{
  DEBUG_CONSTRUCTOR("GeomLine")
}

GeomLine::GeomLine(const GeomLine& copy) : 
  GeomObj(), 
  p1(copy.p1), 
  p2(copy.p2),
  line_width_(0.0)
{
  DEBUG_CONSTRUCTOR("GeomLine")
}

GeomLine::~GeomLine()
{
  DEBUG_DESTRUCTOR("GeomLine")
}

GeomObj* GeomLine::clone()
{    
  return new GeomLine(*this);
}

void GeomLine::get_bounds(BBox& bb)
{
  bb.extend(p1);
  bb.extend(p2);
}

void
GeomLine::setLineWidth(float val) 
{
  line_width_ = val;
}

void
GeomLine::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 0)) return;
  di->polycount_++;
  // Set line width. Set it

  if( line_width_ > 0.0 )
  {
    glLineWidth(line_width_);
  }

  glBegin(GL_LINE_STRIP);
  glVertex3d(p1.x(), p1.y(), p1.z());
  glVertex3d(p2.x(), p2.y(), p2.z());
  glEnd();

  // HACK set line width back to default
  // our scenegraph needs more graceful control of such state.
  glLineWidth(di->line_width_);
  post_draw(di);
}

GeomLines::GeomLines()
  : line_width_(0.0)
{
  DEBUG_CONSTRUCTOR("GeomLines")
}

GeomLines::GeomLines(const GeomLines& copy)
  : GeomObj(copy),
    line_width_(copy.line_width_),
    points_(copy.points_),
    colors_(copy.colors_),
    indices_(copy.indices_)
{
  DEBUG_CONSTRUCTOR("GeomLines")
}

GeomLines::~GeomLines()
{
  DEBUG_DESTRUCTOR("GeomLines")
}

GeomObj* 
GeomLines::clone()
{
  return new GeomLines(*this);
}

void 
GeomLines::get_bounds(BBox& bb)
{
  for(size_t i=0;i<points_.size();i+=3)
  {
    bb.extend(Point(points_[i+0], points_[i+1], points_[i+2]));
  }
}

void
GeomLines::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 0)) return;

  di->polycount_+= points_.size()/6;

  if( line_width_ > 0.0 )
    glLineWidth(line_width_);

  if (points_.size()) {
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

  glDrawArrays(GL_LINES, 0, points_.size()/3);

  glLineWidth(di->line_width_);

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

void
GeomLines::add(const Point& p1, const Point& p2)
{
  points_.push_back(p1.x());
  points_.push_back(p1.y());
  points_.push_back(p1.z());
  points_.push_back(p2.x());
  points_.push_back(p2.y());
  points_.push_back(p2.z());
}

void
GeomLines::add(const Point& p1, const MaterialHandle &c1,
	       const Point& p2, const MaterialHandle &c2)
{
  add(p1, p2);

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
}

void
GeomLines::add(const Point& p0, double cindex0,
	       const Point& p1, double cindex1)
{
  add(p0, p1);
  
  indices_.push_back(cindex0);
  indices_.push_back(cindex1);
}

void
GeomLines::add(const Point& p0, const MaterialHandle &c0, double cindex0,
	       const Point& p1, const MaterialHandle &c1, double cindex1)
{
  add(p0, c0, p1, c1);
  
  indices_.push_back(cindex0);
  indices_.push_back(cindex1);
}

GeomTranspLines::GeomTranspLines()
  : GeomLines(),
    xreverse_(false),
    yreverse_(false),
    zreverse_(false)
{
  DEBUG_CONSTRUCTOR("GeomTranspLines")
}

GeomTranspLines::GeomTranspLines(const GeomTranspLines& copy)
  : GeomLines(copy),
    xindices_(copy.xindices_),
    yindices_(copy.yindices_),
    zindices_(copy.zindices_),
    xreverse_(copy.xreverse_),
    yreverse_(copy.yreverse_),
    zreverse_(copy.zreverse_)
{
  DEBUG_CONSTRUCTOR("GeomTranspLines")
}

GeomTranspLines::~GeomTranspLines()
{
  DEBUG_DESTRUCTOR("GeomTranspLines")
}

GeomObj* GeomTranspLines::clone()
{
  return new GeomTranspLines(*this);
}

void
GeomTranspLines::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 0)) return;

  Sort();
  get_view(di);

  std::vector<unsigned int> &clist =
    (di->axis_==0)?xindices_:((di->axis_==1)?yindices_:zindices_);

  bool &reverse =
    (di->axis_==0)?xreverse_:((di->axis_==1)?yreverse_:zreverse_);

  di->polycount_+=points_.size()/6;

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  if( line_width_ > 0.0 )
    glLineWidth(line_width_);

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

  if (di->dir_ ==  1 &&  reverse ||
      di->dir_ == -1 && !reverse)
  {
    std::reverse(clist.begin(), clist.end());
    reverse = !reverse;
  }

  if (clist.size())
    glDrawElements(GL_LINES, clist.size(), GL_UNSIGNED_INT, &(clist.front()));

  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);

  glLineWidth(di->line_width_);

  glDisable(GL_BLEND);
  glDisable(GL_TEXTURE_1D);

  post_draw(di);
}

static bool
pair_less(const std::pair<float, unsigned int> &a,
	  const std::pair<float, unsigned int> &b)
{
  return a.first < b.first;
}
 
void
GeomTranspLines::Sort()
{
  const size_t vsize = points_.size() / 6;
  if (xindices_.size() == vsize*2) return;
  
  xreverse_ = false;
  yreverse_ = false;
  zreverse_ = false;

  std::vector<std::pair<float, unsigned int> > tmp(vsize);
  unsigned int i;

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*6+0] + points_[i*6+3];
    tmp[i].second = i*6;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  xindices_.resize(vsize*2);
  for (i=0; i < vsize; i++)
  {
    xindices_[i*2+0] = tmp[i].second / 3;
    xindices_[i*2+1] = tmp[i].second / 3 + 1;
  }

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*6+1] + points_[i*6+4];
    tmp[i].second = i*6;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  yindices_.resize(vsize*2);
  for (i=0; i < vsize; i++)
  {
    yindices_[i*2+0] = tmp[i].second / 3;
    yindices_[i*2+1] = tmp[i].second / 3 + 1;
  }

  for (i = 0; i < vsize;i++)
  {
    tmp[i].first = points_[i*6+2] + points_[i*6+5];
    tmp[i].second = i*6;
  }
  std::sort(tmp.begin(), tmp.end(), pair_less);

  zindices_.resize(vsize*2);
  for (i=0; i < vsize; i++)
  {
    zindices_[i*2+0] = tmp[i].second / 3;
    zindices_[i*2+1] = tmp[i].second / 3 + 1;
  }

}

TexGeomLines::TexGeomLines()
  : mutex_("TexGeomLines mutex"),
    alpha_(1.0),
    tmapid_(0),
    tex_per_seg_(1)
{
  DEBUG_CONSTRUCTOR("TexGeomLines")
}

TexGeomLines::TexGeomLines(const TexGeomLines& copy)
  : GeomObj(copy),
    mutex_("TexGeomLines mutex"), 
    points_(copy.points_)
{
  DEBUG_CONSTRUCTOR("TexGeomLines")
}

TexGeomLines::~TexGeomLines()
{
  DEBUG_DESTRUCTOR("TexGeomLines")
}

GeomObj* TexGeomLines::clone()
{
  return new TexGeomLines(*this);
}

void 
TexGeomLines::get_bounds(BBox& bb)
{
  for(int i=0; i<points_.size(); i++)
    bb.extend(points_[i]);
}

void
TexGeomLines::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 0)) return;  // lighting is turned off here.
  di->polycount_ += points_.size()/2;

  static Vector view;  // shared by all of them - should be in ViewWindow!

  // always assume you are not there.
  // can't mix them.

  // first set up the line size stuff.

  // here is where the texture stuff has to be done.
  // first setup the texture matrix.

  double model_mat[16]; // this is the modelview matrix

  glGetDoublev(GL_MODELVIEW_MATRIX,model_mat);
  glMatrixMode(GL_TEXTURE);
  glPushMatrix();

  // this is what you rip the view vector from
  // just use the "Z" axis, normalized

  view = Vector(model_mat[0*4+2],model_mat[1*4+2],model_mat[2*4+2]);

  view.normalize();

  for (int q=0;q<15;q++)
    model_mat[q] = 0.0;

  model_mat[0*4+0] = view.x()*0.5;
  model_mat[1*4+0] = view.y()*0.5;
  model_mat[2*4+0] = view.z()*0.5;
  model_mat[3*4+0] = 0.5;

  model_mat[15] = 1.0;

  // you might want to zero out the rest, but id doesn't matter for 1D
  // texture maps

  glLoadMatrixd(model_mat); // loads the matrix.

  if (!tmapid_)
  { // has the texture been created?
    tmap1d_.resize(256*3); // that is the size of the 1D texture.
    for (int i=0; i<256; i++)
    {
      double r,ks,LdotT;

      LdotT = i*2/(255.0)-1.0;
      ks =  0.3*(pow(2*LdotT*LdotT - 1,30));

      r = 0.05 + 0.6*(1-LdotT*LdotT) + ks;

      if (r>1.0)
        r = 1.0;

      if (r < 0.0)
        std::cerr << "Negative r!\n";

      if (r>1.0 || ks>1.0)
        std::cerr << r << " " << ks << " Error - out of range.\n";

      tmap1d_[i*3+0] = (unsigned char)(r*255);
      tmap1d_[i*3+1] = (unsigned char)(r*255);  // just have them be red for now.
      tmap1d_[i*3+2] = (unsigned char)(r*255);
    }

    // Now set the end conditions.
    tmapid_ = glGenLists(1);
    glNewList(tmapid_, GL_COMPILE_AND_EXECUTE);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    GLfloat brder[4];

    brder[3] = 1.0; // this is just the alpha component.

    brder[0] = (tmap1d_[0] + tmap1d_[255*3 + 0])/510.0;
    brder[0] = (tmap1d_[1] + tmap1d_[255*3 + 1])/510.0;
    brder[0] = (tmap1d_[2] + tmap1d_[255*3 + 2])/510.0;

    glTexParameterfv(GL_TEXTURE_1D,GL_TEXTURE_BORDER_COLOR,brder);

    glEnable(GL_TEXTURE_1D);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glTexImage1D(GL_TEXTURE_1D,0,3,
                 256,0,GL_RGB,GL_UNSIGNED_BYTE,
                 &tmap1d_[0]);
    glEndList();
  }
  else
  {
    glCallList(tmapid_);
  }

  glEnable(GL_TEXTURE_1D);

  if (!colors_.size()) // set if you don't have colors.
    glColor4f(1.0,0.0,0.0,1.0);  // this state always needs to be set.

  mutex_.lock();
  if (alpha_ == 1.0)
  {

    glBegin(GL_LINES);
    if (tex_per_seg_)
    {
      if (colors_.size())
      {
        for (int i=0;i<points_.size()/2;i++)
        {
          Point& pt1=points_[i*2];
          Point& pt2=points_[i*2+1];

          Vector& tan = tangents_[i];

          glColor3ubv(colors_[i].ptr());
          glTexCoord3d(tan.x(), tan.y(), tan.z());
          glVertex3d(pt1.x(), pt1.y(), pt1.z());
          glVertex3d(pt2.x(), pt2.y(), pt2.z());
        }
      }
      else
      {
        for (int i=0;i<points_.size()/2;i++)
        {
          Point& pt1=points_[i*2];
          Point& pt2=points_[i*2+1];
          Vector& tan = tangents_[i];

          glTexCoord3d(tan.x(), tan.y(), tan.z());
          glVertex3d(pt1.x(), pt1.y(), pt1.z());
          glVertex3d(pt2.x(), pt2.y(), pt2.z());
        }
      }
    }
    else
    {
      if (colors_.size())
      {
        for (int i=0;i<points_.size()/2;i++)
        {
          Point& pt1=points_[i*2];
          Point& pt2=points_[i*2+1];

          Vector& tan1 = tangents_[i*2];
          Vector& tan2 = tangents_[i*2+1];

          glColor3ubv(colors_[i*2].ptr());
          glTexCoord3d(tan1.x(), tan1.y(), tan1.z());
          glVertex3d(pt1.x(), pt1.y(), pt1.z());

          glColor3ubv(colors_[i*2+1].ptr());
          glTexCoord3d(tan2.x(), tan2.y(), tan2.z());
          glVertex3d(pt2.x(), pt2.y(), pt2.z());
        }
      }
      else
      {
        for (int i=0;i<points_.size()/2;i++)
        {
          Point& pt1=points_[i*2];
          Point& pt2=points_[i*2+1];

          Vector& tan1 = tangents_[i*2];
          Vector& tan2 = tangents_[i*2+1];

          glTexCoord3d(tan1.x(), tan1.y(), tan1.z());
          glVertex3d(pt1.x(), pt1.y(), pt1.z());

          glTexCoord3d(tan2.x(), tan2.y(), tan2.z());
          glVertex3d(pt2.x(), pt2.y(), pt2.z());
        }
      }
    }
    glEnd();
  }
  else
  {
    Sort(); // creates sorted lists.

    // render with transparency.
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

    if (!colors_.size())
      glColor4f(1,0,0,alpha_); // make sure it is used.

    int sort_start=0;
    int sort_dir=1; // positive direction

    //char which;

    if (Abs(view.x()) > Abs(view.y()))
    {
      if (Abs(view.x()) > Abs(view.z()))
      { // use x dir
        if (view.x() < 0)
        {
          sort_dir=-1; sort_start=points_.size()/2-1;
        }
        else
        {
          sort_start = 0;
        }
      }
      else
      { // use z dir
        if (view.z() < 0)
        {
          sort_dir=-1;sort_start = 2*points_.size()/2-1;
        } else
          sort_start = points_.size()/2;
      }
    }
    else if (Abs(view.y()) > Abs(view.z()))
    { // y greates
      if (view.y() < 0)
      {
        sort_dir=-1;sort_start = 3*(points_.size()/2)-1;
      }
      else
      {
        sort_start = 2*points_.size()/2-1;
      }
    }
    else
    { // z is the one
      if (view.z() < 0)
      {
        sort_dir=-1;sort_start = 2*points_.size()/2-1;
      }
      else
      {
        sort_start = points_.size()/2;
      }
    }

    glBegin(GL_LINES);
    int i = sort_start;
    if (tex_per_seg_)
    {
      for (int p=0;p<points_.size()/2;p++)
      {
        Point& pt1=points_[sindex_[i]];
        Point& pt2=points_[sindex_[i]+1]; // already times2.

        Vector& tan = tangents_[sindex_[i]/2];

        glTexCoord3d(tan.x(), tan.y(), tan.z());
        glVertex3d(pt1.x(), pt1.y(), pt1.z());
        glVertex3d(pt2.x(), pt2.y(), pt2.z());

        i += sort_dir; // increment i.
      }
    }
    else
    { // this is from the stream line data.
      if (colors_.size())
      {
        unsigned char aval = (unsigned char)(alpha_*255); // quantize this.
        for (int p=0; p<points_.size()/2; p++)
        {
          Point& pt1=points_[sindex_[i]];
          Point& pt2=points_[sindex_[i]+1]; // already times2.

          Vector& tan1 = tangents_[sindex_[i]];
          Vector& tan2 = tangents_[sindex_[i]+1];

          Colorub& col1 = colors_[sindex_[i]];
          Colorub& col2 = colors_[sindex_[i]+1];

          glColor4ub(col1.r(), col1.g(), col1.b(), aval);
          glTexCoord3d(tan1.x(), tan1.y(), tan1.z());
          glVertex3d(pt1.x(), pt1.y(), pt1.z());

          glColor4ub(col2.r(), col2.g(), col2.b(), aval);
          glTexCoord3d(tan2.x(), tan2.y(), tan2.z());
          glVertex3d(pt2.x(), pt2.y(), pt2.z());

          i += sort_dir; // increment i.
        }
      }
      else
      {
        for (int p=0; p<points_.size()/2; p++)
        {
          Point& pt1 = points_[sindex_[i]];
          Point& pt2 = points_[sindex_[i]+1]; // already times2.

          Vector& tan1 = tangents_[sindex_[i]];
          Vector& tan2 = tangents_[sindex_[i]+1];

          glTexCoord3d(tan1.x(),tan1.y(),tan1.z());
          glVertex3d(pt1.x(), pt1.y(), pt1.z());

          glTexCoord3d(tan2.x(),tan2.y(),tan2.z());
          glVertex3d(pt2.x(), pt2.y(), pt2.z());

          i += sort_dir; // increment i.
        }
      }
    }
    glEnd();
    glDisable(GL_BLEND);
  }
  mutex_.unlock();

  glDisable(GL_TEXTURE_1D);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  post_draw(di);
}

// this is used by the hedgehog...

void 
TexGeomLines::add(const Point& p1, const Point& p2,double scale)
{
  points_.add(p1);
  points_.add(p2);
  
  tangents_.add((p2-p1)*scale);
} 

void 
TexGeomLines::add(const Point& p1, const Vector& dir, const Colorub& c) {
  points_.add(p1);
  points_.add(p1+dir);

  Vector v(dir);
  v.normalize();
  tangents_.add(v);
  colors_.add(c);
}

// this is used by the streamline module...

void TexGeomLines::batch_add(Array1<double>&, Array1<Point>& ps)
{
  tex_per_seg_ = 0;  // this is not the hedgehog...
  int pstart = points_.size();
  int tstart = tangents_.size();

  points_.grow(2*(ps.size()-1));
  tangents_.grow(2*(ps.size()-1));  // ignore times for now...

  int i;
  for(i=0;i<ps.size()-1;i++) 
  {// forward differences to get tangents...
    Vector v = ps[i+1]-ps[i];
    v.normalize();

    tangents_[tstart++] = v; // vector is set...
    points_[pstart++] = ps[i];
    if (i) 
    { // only store it once...
      tangents_[tstart++] = v; // duplicate it otherwise
      points_[pstart++] = ps[i];
    }
  }
  tangents_[tstart] = tangents_[tstart-1]; // duplicate last guy...
  points_[pstart] = ps[i]; // last point...

}

void TexGeomLines::batch_add(Array1<double>&, Array1<Point>& ps,
			     Array1<Color>& cs)
{
  tex_per_seg_ = 0;  // this is not the hedgehog...
  int pstart = points_.size();
  int tstart = tangents_.size();
  int cstart = colors_.size();

  points_.grow(2*(ps.size()-1));
  tangents_.grow(2*(ps.size()-1));
  colors_.grow(2*(ps.size()-1));

  int i;
  for(i=0;i<ps.size()-1;i++)
  {// forward differences to get tangents...
    Vector v = ps[i+1]-ps[i];
    v.normalize();

    tangents_[tstart++] = v; // vector is set...
    points_[pstart++] = ps[i];
    colors_[cstart++] = Colorub(cs[i]);
    if (i)
    { // only store it once...
      tangents_[tstart++] = v; // duplicate it otherwise
      points_[pstart++] = ps[i];
      colors_[cstart++] = Colorub(cs[i]);
    }
  }
  tangents_[tstart] = tangents_[tstart-1]; // duplicate last guy...
  points_[pstart] = ps[i]; // last point...
  colors_[cstart] = Colorub(cs[i]);
}


// this code sorts in three axis...
struct SortHelper 
{
  static Point* points__array;
  int id; // id for this guy...
};

Point* SortHelper::points__array=0;

int CompX(const void* e1, const void* e2)
{
  SortHelper *a = (SortHelper*)e1;
  SortHelper *b = (SortHelper*)e2;

  if (SortHelper::points__array[a->id].x() >
      SortHelper::points__array[b->id].x())
    return 1;
  if (SortHelper::points__array[a->id].x() <
      SortHelper::points__array[b->id].x())
    return -1;

  return 0; // they are equal...
}

int CompY(const void* e1, const void* e2)
{
  SortHelper *a = (SortHelper*)e1;
  SortHelper *b = (SortHelper*)e2;

  if (SortHelper::points__array[a->id].y() >
      SortHelper::points__array[b->id].y())
    return 1;
  if (SortHelper::points__array[a->id].y() <
      SortHelper::points__array[b->id].y())
    return -1;

  return 0; // they are equal...
}

int CompZ(const void* e1, const void* e2)
{
  SortHelper *a = (SortHelper*)e1;
  SortHelper *b = (SortHelper*)e2;

  if (SortHelper::points__array[a->id].z() >
      SortHelper::points__array[b->id].z())
    return 1;
  if (SortHelper::points__array[a->id].z() <
      SortHelper::points__array[b->id].z())
    return -1;

  return 0; // they are equal...
}


void TexGeomLines::Sort()
{
  const unsigned int vsize = points_.size()/2;
  if ((unsigned int)sindex_.size() == vsize*3) return;

  SortHelper::points__array = &points_[0];
  
  // TODO: REPLACE WITH STL
  Array1<SortHelper> help; // list for help stuff...
  
  unsigned int realsize = points_.size()/2;
  unsigned int imul = 2;

  sindex_.resize(3*realsize); // resize the array...

  help.resize(realsize);

  unsigned int i;
  for(i=0; i<realsize; i++)
  {
    help[i].id = imul*i;  // start it in order...
  }


  // TODO: REPLACE WITH STL AS THAT ONE IS FASTER
  qsort(&help[0],help.size(),sizeof(SortHelper),CompX);
  //	int (*) (const void*,const void*)CompX);

  // now dump these ids..

  for(i=0; i<realsize; i++)
  {
    sindex_[i] = help[i].id;
    help[i].id = imul*i;  // reset for next list...
  }
  
  // TODO: REPLACE WITH STL AS THAT ONE IS FASTER
  qsort(&help[0],help.size(),sizeof(SortHelper),CompZ);

  unsigned int j;
  for(j=0; j<realsize; j++,i++)
  {
    sindex_[i] = help[j].id;
    help[j].id=imul*j;
  }

  qsort(&help[0],help.size(),sizeof(SortHelper),CompY);

  for(j=0; j<realsize; j++,i++)
  {
    sindex_[i] = help[j].id;
  }
}

GeomCLineStrips::GeomCLineStrips()
  : line_width_(0.0)
{
  DEBUG_CONSTRUCTOR("GeomCLineStrips")
}

GeomCLineStrips::GeomCLineStrips(const GeomCLineStrips& copy)
  : GeomObj(copy),
    line_width_(copy.line_width_),
    points_(copy.points_),
    colors_(copy.colors_)
{
  DEBUG_CONSTRUCTOR("GeomCLineStrips")
}

GeomCLineStrips::~GeomCLineStrips()
{
  DEBUG_DESTRUCTOR("GeomCLineStrips")
}

GeomObj* GeomCLineStrips::clone()
{
  return new GeomCLineStrips(*this);
}

void 
GeomCLineStrips::get_bounds(BBox& bb)
{
  for(size_t s = 0; s < points_.size(); s++) 
  {
    for (size_t i = 0; i < points_[s].size(); i+=3) 
    {
      bb.extend(Point(points_[s][i+0], points_[s][i+1], points_[s][i+2]));
    }
  }
}

void
GeomCLineStrips::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di, matl, 0)) return;

  glLineWidth(line_width_);

  const int n_strips = points_.size();
  for (int i = 0; i < n_strips; i++)
  {
    const int n_points = points_[i].size()/3;
    di->polycount_ += n_points-1;
    if (points_[i].size()) {
      glVertexPointer(3, GL_FLOAT, 0, &(points_[i].front()));
      glEnableClientState(GL_VERTEX_ARRAY);
    }
    else {
      glDisableClientState(GL_VERTEX_ARRAY);
    }
    if (colors_[i].size())
    {
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, &(colors_[i].front()));
      glEnableClientState(GL_COLOR_ARRAY);
    }
    else
    {
      glDisableClientState(GL_COLOR_ARRAY);
    }

    glDrawArrays(GL_LINE_STRIP, 0, n_points);
  }

  glLineWidth(di->line_width_);

  post_draw(di);
}

void
GeomCLineStrips::add(const std::vector<Point> &p, 
		     const std::vector<MaterialHandle> &c)
{
  points_.push_back(std::vector<float>());
  colors_.push_back(std::vector<unsigned char>());

  ASSERT(p.size() == c.size());
  for (unsigned int i = 0; i < p.size(); i++)
  {
    add(p[i], c[i]);
  }
}

void
GeomCLineStrips::add(const Point &p, 
		     const MaterialHandle c)
{
  if (points_.empty()) points_.push_back(std::vector<float>());
  if (colors_.empty()) colors_.push_back(std::vector<unsigned char>());
  
  points_.back().push_back(p.x());
  points_.back().push_back(p.y());
  points_.back().push_back(p.z());

  colors_.back().push_back(COLOR_FTOB(c->diffuse.r()));
  colors_.back().push_back(COLOR_FTOB(c->diffuse.g()));
  colors_.back().push_back(COLOR_FTOB(c->diffuse.b()));
  colors_.back().push_back(COLOR_FTOB(c->transparency));
}

void
GeomCLineStrips::newline()
{
  points_.push_back(std::vector<float>());
  colors_.push_back(std::vector<unsigned char>());
}

} // End namespace SCIRun

