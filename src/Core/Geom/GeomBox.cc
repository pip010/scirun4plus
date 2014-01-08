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
 *  GeomBox.cc:  A box object
 *
 *  Written by:
 *   Yarden Livnat
 *   Department of Computer Science
 *   University of Utah
 *   Feb. 1997
 *
 */

#include <Core/Util/Debug.h>
#include <Core/Geom/GeomBox.h>
#include <Core/Geom/DrawInfoOpenGL.h>

#include <Core/Math/MiscMath.h>

namespace SCIRun {

GeomBox::GeomBox(const Point& p, const Point& q, int op) : GeomObj()
{
  DEBUG_CONSTRUCTOR("GeomBox")

  min = Min( p, q );
  max = Max( p, q );

  for (int i=0; i<6; i++ )
    opacity[i] = op;
}

GeomBox::GeomBox(const GeomBox& copy)
: GeomObj(copy)
{
  DEBUG_CONSTRUCTOR("GeomBox")

  min = copy.min;
  max = copy.max;
  for (int s=0; s<6; s++)
    opacity[s] = copy.opacity[s];
}

GeomBox::~GeomBox()
{
  DEBUG_DESTRUCTOR("GeomBox")
}

GeomObj* GeomBox::clone()
{
    return new GeomBox(*this);
}

void GeomBox::get_bounds(BBox& bb)
{
  bb.extend(min);
  bb.extend(max);
}

void
GeomBox::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di,matl,1)) return;

  di->polycount_ += 6;

  if (di->currently_lit_)
  {
    glEnable(GL_NORMALIZE);

    switch(di->get_drawtype())
    {
    case DrawInfoOpenGL::WireFrame:
    case DrawInfoOpenGL::Flat:
    case DrawInfoOpenGL::Gouraud:
      { 
        glBegin(GL_QUADS);

        // top
        glNormal3f(0,0,1);
        glColor4f(0.0, 0.0, 1.0, 0.0);
        glVertex3d(min.x(),min.y(),max.z());
        glVertex3d(max.x(),min.y(),max.z());
        glVertex3d(max.x(),max.y(),max.z());
        glVertex3d(min.x(),max.y(),max.z());

        // bottom
        glNormal3f(0,0,-1);
        glColor4f(0.0, 0.0, 0.5, 0.0);
        glVertex3d(min.x(),min.y(),min.z());
        glVertex3d(min.x(),max.y(),min.z());
        glVertex3d(max.x(),max.y(),min.z());
        glVertex3d(max.x(),min.y(),min.z());
        
        // left
        glNormal3f(-1.0,0,0);
        glColor4f(0.5, 0.0, 0.0, 0.0);
        glVertex3d(min.x(),min.y(),min.z());
        glVertex3d(min.x(),min.y(),max.z());
        glVertex3d(min.x(),max.y(),max.z());
        glVertex3d(min.x(),max.y(),min.z());

        // right
        glNormal3f(1,0,0);
        glColor4f(1.0, 0.0, 0.0, 0.0);
        glVertex3d(max.x(),min.y(),min.z());
        glVertex3d(max.x(),max.y(),min.z());
        glVertex3d(max.x(),max.y(),max.z());
        glVertex3d(max.x(),min.y(),max.z());
                
        // top
        glNormal3f(0,1.0,0);
        glColor4f(0.0, 1.0, 0.0, 0.0);
        glVertex3d(min.x(),max.y(),min.z());
        glVertex3d(min.x(),max.y(),max.z());
        glVertex3d(max.x(),max.y(),max.z());
        glVertex3d(max.x(),max.y(),min.z());

        // back
        glNormal3f(0,-1,0);
        glColor4f(0.0, 0.5, 0.0, 0.0);
        glVertex3d(min.x(),min.y(),min.z());
        glVertex3d(max.x(),min.y(),min.z());
        glVertex3d(max.x(),min.y(),max.z());
        glVertex3d(min.x(),min.y(),max.z());

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
        glBegin(GL_QUADS);
        //front
        glVertex3d(max.x(),min.y(),max.z());
        glVertex3d(max.x(),max.y(),max.z());
        glColor4f(0.0,1.0,0.0,0.2);
        glVertex3d(min.x(),max.y(),max.z());
        glVertex3d(min.x(),min.y(),max.z());

        //back
        glVertex3d(max.x(),max.y(),min.z());
        glVertex3d(max.x(),min.y(),min.z());
        glVertex3d(min.x(),min.y(),min.z());
        glColor4f(0.0,1.0,0.0,0.2);
        glVertex3d(min.x(),max.y(),min.z());
        
        glColor4f(1.0,0.0,0.0,0.2);
        
        //left
        glVertex3d(min.x(),min.y(),max.z());
        glVertex3d(min.x(),max.y(),max.z());
        glVertex3d(min.x(),max.y(),min.z());
        glVertex3d(min.x(),min.y(),min.z());
        glColor4f(1.0,0.0,0.0,0.2);
        
        //right
        glVertex3d(max.x(),min.y(),min.z());
        glVertex3d(max.x(),max.y(),min.z());
        glVertex3d(max.x(),max.y(),max.z());
        glColor4f(1.0,0.0,0.0,0.2);
        glVertex3d(max.x(),min.y(),max.z());
        
        
        glColor4f(0.0,0.0,1.0,0.2);
        
        //top
        glVertex3d(min.x(),max.y(),max.z());
        glVertex3d(max.x(),max.y(),max.z());
        glColor4f(0.0,0.0,1.0,0.2);
        glVertex3d(max.x(),max.y(),min.z());
        glVertex3d(min.x(),max.y(),min.z());

        //bottom
        glVertex3d(min.x(),min.y(),min.z());
        glColor4f(0.0,0.0,1.0,0.2);
        glVertex3d(max.x(),min.y(),min.z());
        glVertex3d(max.x(),min.y(),max.z());
        glVertex3d(min.x(),min.y(),max.z());

        glEnd();
      }
      break;
    }
  }
  post_draw(di);
}

GeomSimpleBox::GeomSimpleBox(const Point& p, const Point& q) : GeomObj()
{
  min = Min( p, q );
  max = Max( p, q );
}

GeomSimpleBox::GeomSimpleBox(const GeomSimpleBox& copy)
  : GeomObj(copy), min(copy.min), max(copy.max)
{
}

GeomSimpleBox::~GeomSimpleBox()
{
}

GeomObj* GeomSimpleBox::clone()
{
  return new GeomSimpleBox(*this);
}

void
GeomSimpleBox::get_bounds(BBox& bb)
{
  bb.extend(min);
  bb.extend(max);
}

void
GeomSimpleBox::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di,matl,1)) return;

  di->polycount_ += 6;

  glEnable(GL_NORMALIZE);

  glBegin(GL_QUADS);

  // top
  glNormal3f(0,0,1);
  glVertex3d(min.x(),min.y(),max.z());
  glVertex3d(max.x(),min.y(),max.z());
  glVertex3d(max.x(),max.y(),max.z());
  glVertex3d(min.x(),max.y(),max.z());

  // bottom
  glNormal3f(0,0,-1);
  glVertex3d(min.x(),min.y(),min.z());
  glVertex3d(min.x(),max.y(),min.z());
  glVertex3d(max.x(),max.y(),min.z());
  glVertex3d(max.x(),min.y(),min.z());
        
  // left
  glNormal3f(-1.0,0,0);
  glVertex3d(min.x(),min.y(),min.z());
  glVertex3d(min.x(),min.y(),max.z());
  glVertex3d(min.x(),max.y(),max.z());
  glVertex3d(min.x(),max.y(),min.z());

  // right
  glNormal3f(1,0,0);
  glVertex3d(max.x(),min.y(),min.z());
  glVertex3d(max.x(),max.y(),min.z());
  glVertex3d(max.x(),max.y(),max.z());
  glVertex3d(max.x(),min.y(),max.z());
                
  // top
  glNormal3f(0,1.0,0);
  glVertex3d(min.x(),max.y(),min.z());
  glVertex3d(min.x(),max.y(),max.z());
  glVertex3d(max.x(),max.y(),max.z());
  glVertex3d(max.x(),max.y(),min.z());

  // back
  glNormal3f(0,-1,0);
  glVertex3d(min.x(),min.y(),min.z());
  glVertex3d(max.x(),min.y(),min.z());
  glVertex3d(max.x(),min.y(),max.z());
  glVertex3d(min.x(),min.y(),max.z());

  glEnd();

  post_draw(di);
}

GeomCBox::GeomCBox(const Point& p, const Point& q) : GeomSimpleBox(p, q)
{
}

GeomCBox::GeomCBox(const GeomCBox& copy)
  : GeomSimpleBox(copy)
{
}

GeomCBox::~GeomCBox()
{
}

GeomObj* 
GeomCBox::clone()
{
  return new GeomCBox(*this);
}

void
GeomCBox::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di,matl,1)) return;

  di->polycount_ += 6;

  glEnable(GL_NORMALIZE);

  glBegin(GL_QUADS);

  // top
  glNormal3f(0,0,1);
  glColor4f(0.0, 0.0, 1.0, 0.0);
  glVertex3d(min.x(),min.y(),max.z());
  glVertex3d(max.x(),min.y(),max.z());
  glVertex3d(max.x(),max.y(),max.z());
  glVertex3d(min.x(),max.y(),max.z());

  // bottom
  glNormal3f(0,0,-1);
  glColor4f(0.0, 0.0, 0.5, 0.0);
  glVertex3d(min.x(),min.y(),min.z());
  glVertex3d(min.x(),max.y(),min.z());
  glVertex3d(max.x(),max.y(),min.z());
  glVertex3d(max.x(),min.y(),min.z());
        
  // left
  glNormal3f(-1.0,0,0);
  glColor4f(0.5, 0.0, 0.0, 0.0);
  glVertex3d(min.x(),min.y(),min.z());
  glVertex3d(min.x(),min.y(),max.z());
  glVertex3d(min.x(),max.y(),max.z());
  glVertex3d(min.x(),max.y(),min.z());

  // right
  glNormal3f(1,0,0);
  glColor4f(1.0, 0.0, 0.0, 0.0);
  glVertex3d(max.x(),min.y(),min.z());
  glVertex3d(max.x(),max.y(),min.z());
  glVertex3d(max.x(),max.y(),max.z());
  glVertex3d(max.x(),min.y(),max.z());
                
  // top
  glNormal3f(0,1.0,0);
  glColor4f(0.0, 1.0, 0.0, 0.0);
  glVertex3d(min.x(),max.y(),min.z());
  glVertex3d(min.x(),max.y(),max.z());
  glVertex3d(max.x(),max.y(),max.z());
  glVertex3d(max.x(),max.y(),min.z());

  // back
  glNormal3f(0,-1,0);
  glColor4f(0.0, 0.5, 0.0, 0.0);
  glVertex3d(min.x(),min.y(),min.z());
  glVertex3d(max.x(),min.y(),min.z());
  glVertex3d(max.x(),min.y(),max.z());
  glVertex3d(min.x(),min.y(),max.z());

  glEnd();

  post_draw(di);
}

GeomBoxes::GeomBoxes(double edge, int nu, int nv)
  : GeomObj(),
    nu_(nu),
    nv_(nv),
    global_edge_(edge)
{
}

GeomBoxes::GeomBoxes(const GeomBoxes& copy)
  : GeomObj(copy),
    centers_(copy.centers_),
    edges_(copy.edges_),
    colors_(copy.colors_),
    indices_(copy.indices_),
    nu_(copy.nu_),
    nv_(copy.nv_),
    global_edge_(copy.global_edge_)
{
}

GeomBoxes::~GeomBoxes()
{
}


GeomObj *
GeomBoxes::clone()
{
  return new GeomBoxes(*this);
}


void
GeomBoxes::get_bounds(BBox& bb)
{
  const bool ugr = !(edges_.size() == centers_.size());
  for (size_t i=0; i < centers_.size(); i++)
  {
    bb.extend(centers_[i], ugr?global_edge_:edges_[i]);
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
GeomBoxes::add(const Point &center)
{
  centers_.push_back(center);
}


void
GeomBoxes::add(const Point &center, const MaterialHandle &mat)
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
}

void
GeomBoxes::add(const Point &center, float index)
{
  add(center);
  indices_.push_back(index);
}

bool
GeomBoxes::add_edge(const Point &c, double r)
{
  if (r < 1.0e-6) { return false; }
  centers_.push_back(c);
  edges_.push_back(r);
  return true;
}

bool
GeomBoxes::add_edge(const Point &c, double r, const MaterialHandle &mat)
{
  if (r < 1.0e-6) { return false; }
  add_edge(c, r);
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
GeomBoxes::add_edge(const Point &c, double r, float index)
{
  if (r < 1.0e-6) { return false; }
  add_edge(c, r);
  indices_.push_back(index);
  return true;
}

void
GeomBoxes::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  const bool ulr = edges_.size() == centers_.size();
  if (!ulr && global_edge_ < 1.0e-6) { return; }

  if (!pre_draw(di, matl, 1)) return;

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

    double edge = ulr ? edges_[i] : global_edge_;

    glTranslated(centers_[i].x()-edge/2.0,
                 centers_[i].y()-edge/2.0,
                 centers_[i].z()-edge/2.0);

    di->polycount_ += 6;

    glEnable(GL_NORMALIZE);

    glBegin(GL_QUADS);

    // top
    glNormal3f(0,0,1);
    glVertex3d(0,0,edge);
    glVertex3d(edge,0,edge);
    glVertex3d(edge,edge,edge);
    glVertex3d(0,edge,edge);

    // bottom
    glNormal3f(0,0,-1);
    glVertex3d(0,0,0);
    glVertex3d(0,edge,0);
    glVertex3d(edge,edge,0);
    glVertex3d(edge,0,0);
        
    // left
    glNormal3f(-1.0,0,0);
    glVertex3d(0,0,0);
    glVertex3d(0,0,edge);
    glVertex3d(0,edge,edge);
    glVertex3d(0,edge,0);

    // right
    glNormal3f(1,0,0);
    glVertex3d(edge,0,0);
    glVertex3d(edge,edge,0);
    glVertex3d(edge,edge,edge);
    glVertex3d(edge,0,edge);
                
    // top
    glNormal3f(0,1.0,0);
    glVertex3d(0,edge,0);
    glVertex3d(0,edge,edge);
    glVertex3d(edge,edge,edge);
    glVertex3d(edge,edge,0);

    // back
    glNormal3f(0,-1,0);
    glVertex3d(0,0,0);
    glVertex3d(edge,0,0);
    glVertex3d(edge,0,edge);
    glVertex3d(0,0,edge);

    glEnd();

    glPopMatrix();
  }

  glDisable(GL_TEXTURE_1D);

  post_draw(di);
}

} // End namespace SCIRun
