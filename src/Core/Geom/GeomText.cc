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

#include <Core/Util/Debug.h>

#include <Core/Datatypes/Color.h>
#include <Core/Geometry/Point.h>
#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomText.h>
#include <Core/Geom/GeomMaterial.h>

#include <Core/Persistent/PersistentSTL.h>
#include <Core/Thread/MutexPool.h>

namespace SCIRun {

GeomTexts::GeomTexts() : 
  GeomObj(),
  fontindex_(2),
  disable_depth_test_(false),
  is_2d_(false)
{
  DEBUG_CONSTRUCTOR("GeomTexts")
  renderer_ = new TextRenderer;
}

GeomTexts::GeomTexts(TextRendererHandle& renderer) : 
  GeomObj(),
  fontindex_(2),
  renderer_(renderer),
  disable_depth_test_(false),
  is_2d_(false)
{
  DEBUG_CONSTRUCTOR("GeomTexts")
}

GeomTexts::GeomTexts(const GeomTexts& copy) :
  GeomObj(copy),
  fontindex_(copy.fontindex_),
  text_(copy.text_),
  location_(copy.location_),
  color_(copy.color_),
  renderer_(copy.renderer_.get_rep()),
  disable_depth_test_(copy.disable_depth_test_),
  is_2d_(copy.is_2d_),
  offset_(0.05),
  scale_(0.0)
{
  DEBUG_CONSTRUCTOR("GeomTexts")
}


GeomTexts::~GeomTexts()
{
  DEBUG_DESTRUCTOR("GeomTexts")
}


GeomObj* GeomTexts::clone()
{
  return new GeomTexts(*this);
}


void
GeomTexts::get_bounds(BBox& in_bb)
{
  for (size_t i = 0; i < location_.size(); i++)
  {
    in_bb.extend( location_[i] );
  }
  
  in_bb.extend(in_bb.diagonal().length()*0.3);
  
  scale_ = in_bb.diagonal().length();
}

void
GeomTexts::reset_bbox()
{
}


void
GeomTexts::set_font_index(int a)
{
  if (a >= 0 && a < 50)
  {
    fontindex_ = a;
  }
}

void
GeomTexts::set_offset(double o)
{
  offset_ = o;
}


void
GeomTexts::add(const std::string &t, const Point &p)
{
  text_.push_back(t);
  location_.push_back(p);
  
}

void
GeomTexts::add(const std::string &t, const Point &p, const Color &c)
{
  text_.push_back(t);
  location_.push_back(p);
  color_.push_back(c);
}

void
GeomTexts::add(const std::string &t, const Point &p, float index)
{
  text_.push_back(t);
  location_.push_back(p);
  index_.push_back(index);
}


void
GeomTexts::add(const std::string &t, const Vector &norm, const Point &p)
{
  text_.push_back(t);
  location_.push_back(p);
  normal_.push_back(norm);
}

void
GeomTexts::add(const std::string &t, const Vector &norm, const Point &p, const Color &c)
{
  text_.push_back(t);
  location_.push_back(p);
  color_.push_back(c);
  normal_.push_back(norm);
}

void
GeomTexts::add(const std::string &t, const Vector &norm, const Point &p, float index)
{
  text_.push_back(t);
  location_.push_back(p);
  index_.push_back(index);
  normal_.push_back(norm);
}


void
GeomTexts::clear()
{
  text_.clear();
  location_.clear();
  color_.clear();
  index_.clear();
}

void
GeomTexts::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di,matl,0)) return;
  glDisable(GL_TEXTURE_1D);
  if (disable_depth_test_) glDisable(GL_DEPTH_TEST);
  glPushAttrib(GL_LIST_BIT);

  switch(fontindex_)
  {
    case 0:
      renderer_->set_size(14); break;
    case 1:
      renderer_->set_size(18); break;
    case 2:
      renderer_->set_size(24); break;
    case 3:
      renderer_->set_size(30); break;
    case 4:
      renderer_->set_size(36); break;
  }

  renderer_->set_offset(di->text_offset_*scale_);

  if (matl) 
  {
    renderer_->set_color(matl->diffuse.r(), matl->diffuse.g(), 
			 matl->diffuse.b(), 1.0); 
  }

  double sf = renderer_->find_scale_factor(location_, di->view_);

  const bool coloring = color_.size() == location_.size();
  const bool use_normal = normal_.size() == location_.size();
  
  for (unsigned int i = 0; i < text_.size(); i++)
  {
    if (coloring) 
    {
      renderer_->set_color(color_[i].r(), color_[i].g(), 
			   color_[i].b(), 1.0); 
    }
    if (is_2d_p()) 
    {
      GLint vmat[4];
      glGetIntegerv(GL_VIEWPORT, vmat);
      int w = vmat[2];
      int h = vmat[3];
      
      float sx = w * ((location_[i].x() + 1.) / 2.);
      float sy = h * ((location_[i].y() + 1.) / 2.);
      renderer_->render(text_[i], sx, sy, TextRenderer::SW | TextRenderer::SHADOW);
      glDisable(GL_BLEND);
      glShadeModel(GL_SMOOTH);
    } 
    else 
    {
      if (use_normal)
      {
        renderer_->render(text_[i], location_[i],normal_[i], sf,
          TextRenderer::SW);
      }
      else
      {
        renderer_->render(text_[i], location_[i].x(), 
          location_[i].y(), location_[i].z(), 
          di->view_, sf,
          TextRenderer::SW);
      }
    }
  }

  glPopAttrib();
  if (disable_depth_test_) glEnable(GL_DEPTH_TEST);
  post_draw(di);
}

GeomTextsCulled::GeomTextsCulled()
  : GeomTexts()
{
  DEBUG_CONSTRUCTOR("GeomTextsCulled")
}

GeomTextsCulled::GeomTextsCulled(const GeomTextsCulled& copy) :
  GeomTexts(copy),
  normal_(copy.normal_)
{
  DEBUG_CONSTRUCTOR("GeomTextsCulled")
}

GeomTextsCulled::~GeomTextsCulled()
{
  DEBUG_DESTRUCTOR("GeomTextsCulled")
}

GeomObj* GeomTextsCulled::clone()
{
  return new GeomTextsCulled(*this);
}

void
GeomTextsCulled::add(const std::string &t, const Point &p, const Vector &n)
{
  GeomTexts::add(t, p);
  normal_.push_back(n);
}

void
GeomTextsCulled::add(const std::string &t, const Point &p,
		     const Vector &n, const Color &c)
{
  GeomTexts::add(t, p, c);
  normal_.push_back(n);
}

void
GeomTextsCulled::add(const std::string &t, const Point &p,
		     const Vector &n, float index)
{
  GeomTexts::add(t, p, index);
  normal_.push_back(n);
}

void
GeomTextsCulled::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di,matl,0)) return;

  glDisable(GL_DEPTH_TEST);
  glPushAttrib(GL_LIST_BIT);

  switch(fontindex_)
  {
    case 0:
      renderer_->set_size(18); break;
    case 1:
      renderer_->set_size(24); break;
    case 2:
      renderer_->set_size(32); break;
    case 3:
      renderer_->set_size(40); break;
    case 4:
      renderer_->set_size(48); break;
  }
  
  if (matl) 
  {
    renderer_->set_color(matl->diffuse.r(), matl->diffuse.g(), 
			 matl->diffuse.b(), 1.0); 
  }

  double sf = renderer_->find_scale_factor(location_, di->view_);

  const bool coloring = color_.size() == location_.size();
  double mat[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, mat);
  const Vector view (mat[2], mat[6], mat[10]);

  for (size_t i = 0; i < text_.size(); i++)
  {
    if (Dot(view, normal_[i]) < 0)
    {
      if (coloring) 
      {
        renderer_->set_color(color_[i].r(), color_[i].g(), 
			     color_[i].b(), 1.0); 
      }
      renderer_->render(text_[i], location_[i].x(), 
			location_[i].y(), location_[i].z(), 
			di->view_, sf,
			TextRenderer::SW);
    }
  }

  glPopAttrib();
  glEnable(GL_DEPTH_TEST);
  post_draw(di);
}

GeomTextTexture::GeomTextTexture(const std::string& )
  : GeomObj(),
    text_(),
    origin_(0.,0.,0.),
    up_(0.,1.,0.),
    left_(0.,1.,0.),
    color_(Color(1,1,1)),
    anchor_(sw),
    up_hack_(false)
{
  DEBUG_CONSTRUCTOR("GeomTextTexture")
}

GeomTextTexture::GeomTextTexture( const std::string& /*fontfile*/,
				  const std::string &text, 
				  const Point &at,
				  const Vector &up,
				  const Vector &left,
				  const Color &c)
  : GeomObj(), 
    text_(text),
    origin_(at),
    up_(up),
    left_(left),
    color_(c),
    anchor_(sw),
    own_font_(true),
    up_hack_(false)
{
  DEBUG_CONSTRUCTOR("GeomTextTexture")
}

GeomTextTexture::GeomTextTexture(const GeomTextTexture& copy)
  : GeomObj(copy)
{
  DEBUG_CONSTRUCTOR("GeomTextTexture")
}

GeomObj* 
GeomTextTexture::clone()
{
  return new GeomTextTexture(*this);
}

void 
GeomTextTexture::get_bounds(BBox& in_bb)
{
  in_bb.extend(origin_);
}

GeomTextTexture::~GeomTextTexture()
{
  DEBUG_DESTRUCTOR("GeomTextTexture")
}

void GeomTextTexture::render() 
{
}

void
GeomTextTexture::build_transform(Transform& /*modelview*/) 
{
}

void
GeomTextTexture::draw(DrawInfoOpenGL* di, Material* matl, double)
{
  if (!pre_draw(di,matl,0)) return;

  glColor4d(1.0,1.0,1.0,1.0);
  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_TEXTURE_2D);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  double trans[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, trans);
  Transform modelview;
  modelview.set_trans(trans);
  build_transform(modelview);
  transform_.get_trans(trans);
  glMultMatrixd((GLdouble *)trans);
  render();
  glPopMatrix();
  glDepthMask(GL_TRUE);
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_BLEND);
  post_draw(di);
}

} // End namespace SCIRun
