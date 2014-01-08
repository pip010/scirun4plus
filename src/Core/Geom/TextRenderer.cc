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
 *  TextRenderer.cc
 *
 *  Written by:
 *   McKay Davis
 *   School of Computing
 *   University of Utah
 *   November, 2005
 *
 */

#include <Core/Util/StringUtil.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Geom/TextRenderer.h>
#include <Core/Geom/TextureObj.h>
#include <Core/Util/MemoryUtil.h>

#include <Core/Math/MiscMath.h>
#include <Core/Util/Environment.h>
#include <Core/Util/FileUtils.h>

#include <sci_gl.h>
#include <sci_glx.h>


#include <freetype/fterrors.h>

#ifdef _WIN32
# undef min
# undef max
#endif

namespace SCIRun {  

TextRenderer::TextRenderer() :
  ref_cnt(0),
  glyphs_(),
  shadow_offset_(std::make_pair(1,-1)),
  cursor_position_(0),
  textures_(),
  texture_(0),
  x_(0),
  y_(0),
  height_(0),
  offset_(0.0)
{
  face_ = new FreeTypeFace(get_default_font_name());
  set_color(1.0, 1.0, 1.0, 1.0);
  set_shadow_color(0.0, 0.0, 0.0, 1.0);
}

TextRenderer::TextRenderer(const std::string& fontname) :
  ref_cnt(0),
  glyphs_(),
  shadow_offset_(std::make_pair(1,-1)),
  cursor_position_(0),
  textures_(),
  texture_(0),
  x_(0),
  y_(0),
  height_(0),
  offset_(0.0)
{
  face_ = new FreeTypeFace(fontname);
  set_color(1.0, 1.0, 1.0, 1.0);
  set_shadow_color(0.0, 0.0, 0.0, 1.0);
}


TextRenderer::TextRenderer(FreeTypeFace *face) :
  ref_cnt(0),
  face_(face),
  glyphs_(),
  shadow_offset_(std::make_pair(1,-1)),
  cursor_position_(0),
  textures_(),
  texture_(0),
  x_(0),
  y_(0),
  height_(0),
  offset_(0.0)  
{
  set_color(1.0, 1.0, 1.0, 1.0);
  set_shadow_color(0.0, 0.0, 0.0, 1.0);
}

TextRenderer::~TextRenderer() 
{
  delete_all_items(textures_);
  delete_all_values(glyphs_);
}

int
TextRenderer::width(const std::string &text, int /*flags*/) 
{
  render_string_glyphs_to_texture(text);
  int wid = 0;
  for (unsigned int c = 0; c < text.size(); ++c)
    if (glyphs_[text[c]])
      wid += glyphs_[text[c]]->ft_metrics_.horiAdvance;
  return wid / 64;
}

int
TextRenderer::height(const std::string &text, int /*flags*/) 
{
  render_string_glyphs_to_texture(text);
  int hei = 0;
  for (unsigned int c = 0; c < text.size(); ++c)
    if (glyphs_[text[c]]) 
      hei = Max(hei, (int)(glyphs_[text[c]]->ft_metrics_.horiBearingY));
  return hei / 64;
}


void
TextRenderer::set_color(float r, float g, float b, float a)
{
  color_[0] = r;
  color_[1] = g;
  color_[2] = b;
  color_[3] = a;
}


void
TextRenderer::set_color(float rgba[4]) 
{
  set_color(rgba[0], rgba[1], rgba[2], rgba[3]);
}


void
TextRenderer::set_shadow_color(float r, float g, float b, float a)
{
  shadow_color_[0] = r;
  shadow_color_[1] = g;
  shadow_color_[2] = b;
  shadow_color_[3] = a;
}


void
TextRenderer::set_shadow_color(float rgba[4]) 
{
  set_shadow_color(rgba[0], rgba[1], rgba[2], rgba[3]);
}

void
TextRenderer::set_shadow_offset(int x, int y)
{
  shadow_offset_.first = x;
  shadow_offset_.second = y;
}



void
TextRenderer::layout_text(TextLayoutInfo &layout_info,
			  bool clamp_vals)
{
  const std::string &text = layout_info.text_;
  const Vector &left = layout_info.left_;
  const Vector &up = layout_info.up_;
  int &flags = layout_info.flags_;
  GlyphLayoutVector &layoutvec = layout_info.layout_vector_;
  unsigned int &total_glyphs = layout_info.total_glyphs_;
  
  unsigned int passes = 1;

  if (flags & SHADOW) {
    passes = 2;
  }

  if (flags & EXTRUDED) {
    passes = Max(Abs(shadow_offset_.first), 
                 Abs(shadow_offset_.second))+1;
  }

  total_glyphs = text.size()*passes; 
  layoutvec.reserve(total_glyphs+1);
  int num = total_glyphs + 1 - layoutvec.size();
  for (int i = 0; i < num; ++i) {
    layoutvec.push_back(new GlyphLayoutInfo);
  }
  //  ASSERT(layoutvec.size() == total_glyphs+1);

  Vector left_pt = left / 64.0;
  Vector up_pt = up / 64.0;

  Vector nextline = (face_->ft_face_->size->metrics.height + passes*64)*up_pt;
  Vector cursor_up = (face_->ft_face_->size->metrics.ascender)*up_pt;
  Vector cursor_down = (face_->ft_face_->size->metrics.descender)*up_pt;

  float delta_color[4] = {0.0, 0.0, 0.0, 0.0};
  Vector delta_position(0.0, 0.0, 0.0);
  if (passes > 1) {
    delta_position = (shadow_offset_.first * left + 
                      shadow_offset_.second * up) / (passes - 1);
    for (int rgba = 0; rgba < 4; ++rgba) 
    {
      delta_color[rgba] = (color_[rgba] - shadow_color_[rgba]) / (passes - 1);
    }
  }

  double bbox_x1 = AIR_POS_INF;
  double bbox_y1 = AIR_POS_INF;
  double bbox_z1 = AIR_POS_INF;

  double bbox_x2 = AIR_NEG_INF;
  double bbox_y2 = AIR_NEG_INF;
  double bbox_z2 = AIR_NEG_INF;

#define EXTEND_X1(exp) (bbox_x1 = exp < bbox_x1 ? exp : bbox_x1)
#define EXTEND_Y1(exp) (bbox_y1 = exp < bbox_y1 ? exp : bbox_y1)
#define EXTEND_Z1(exp) (bbox_z1 = exp < bbox_z1 ? exp : bbox_z1)
#define EXTEND_X2(exp) (bbox_x2 = exp > bbox_x2 ? exp : bbox_x2)
#define EXTEND_Y2(exp) (bbox_y2 = exp > bbox_y2 ? exp : bbox_y2)
#define EXTEND_Z2(exp) (bbox_z2 = exp > bbox_z2 ? exp : bbox_z2)
#define EXTEND(exp) tmppt = exp; \
  EXTEND_X1(tmppt.x()); EXTEND_Y1(tmppt.y()); EXTEND_Z1(tmppt.z()); \
  EXTEND_X2(tmppt.x()); EXTEND_Y2(tmppt.y()); EXTEND_Z2(tmppt.z());
  
  bool disable_cursor_flag = true;
  bool do_cursor_layout = false;
  Point cursor_ll(0,0,0);
  Point tmppt;
  
  GlyphInfo *last_glyph = 0;
  GlyphLayoutVector::iterator layout_iter = layoutvec.begin();
  total_glyphs = 0;
  for (unsigned int pass = 0; pass < passes; ++pass) 
  {
    Point position = (pass * delta_position).asPoint();
    int linenum = 0;
    for (unsigned int c = 0; c < text.size(); ++c) 
    {
      if (!(flags & VERTICAL) && text[c] == '\n') 
      {
        position = (pass * delta_position - nextline * (++linenum)).asPoint();
        continue;
      } 

      GlyphLayoutInfo &layout = **layout_iter;
      if (!glyphs_[text[c]]) { 
        layout.glyph_info_ = glyphs_[' '];
      } else {
        layout.glyph_info_ = glyphs_[text[c]];
      }

      if (!layout.glyph_info_) 
      {
        continue;
      }

      total_glyphs++;
      layout_iter++;


      for (int rgba = 0; rgba < 4; ++rgba) 
      {
        if (bool(flags & REVERSE) != bool(passes == 1)) 
        {
          layout.color_[rgba] = color_[rgba] - pass * delta_color[rgba];
        } 
        else 
        {
          layout.color_[rgba] = shadow_color_[rgba] + pass * delta_color[rgba];
        }
      }

      FT_Glyph_Metrics &metrics = layout.glyph_info_->ft_metrics_;

      Point ll;
      if (!(flags & VERTICAL)) 
      {
        if (c && face_->has_kerning_p() && last_glyph) 
        {
          FT_Vector kerning; 
          FT_Get_Kerning(face_->ft_face_, 
                         last_glyph->index_,
                         layout.glyph_info_->index_, 
                         FT_KERNING_DEFAULT, &kerning); 
          position += left_pt * kerning.x;
          // position += up_pt * kerning.y;
        }

        ll = (position + 
              up_pt * (metrics.horiBearingY - metrics.height) +
              left_pt * metrics.horiBearingX);

        if (pass == (passes - 1)) {
          // expand bbox by a a phantom cursor 
          // This fixes the layout from jumping around when 
          // using different strings and anchor combinations
          // ie, the string 'gpq' with its descenders
          // will have a different bounding box than 'ABC'... 
          // creating a phantom cursor that has a descender and ascender
          // at least as large as any glyph in the font fixes this.
          
          if (c == (text.size()-1)) {
            Point temp = position + 
              left_pt * (metrics.horiAdvance);

            EXTEND(temp + cursor_down);
            EXTEND(temp + 2*left + cursor_down);
            EXTEND(temp + 2*left + cursor_up);
            EXTEND(temp + cursor_up);

          }
            
          // Cursor is in current layout position
          if (c == cursor_position_) {
            do_cursor_layout = true;
            cursor_ll = position + left_pt * metrics.horiBearingX;
            if (c) 
              cursor_ll = cursor_ll - left;
          }

          // Cursor is at end of string
          if (cursor_position_ == text.size() && (c == text.size()-1)) {
            do_cursor_layout = true;
            cursor_ll = position + 
              left_pt * (metrics.horiAdvance);
          }
                         
        }

        position = position + left_pt * metrics.horiAdvance;
      } 
      else 
      {
        int halfwidth = (metrics.width) >> 1;
        int advance = (64 + (int)metrics.height);
        ll = (position - up_pt * advance - left_pt * halfwidth);
        position = position - up_pt * advance;
      }

      layout.vertices_[0] = ll;
      layout.vertices_[1] = ll + left_pt * metrics.width;
      layout.vertices_[2] = (ll + left_pt * metrics.width + 
                             up_pt * metrics.height);
      layout.vertices_[3] = ll + up_pt * metrics.height;

      for (int v = 0; v < 4; ++v) 
      {
        EXTEND(layout.vertices_[v]);
      }

      last_glyph = layout.glyph_info_;
    }
  }

  if (!(flags & VERTICAL) && 
      !do_cursor_layout && !cursor_position_) 
  {
    cursor_ll = Point(0.0, 0.0, 0.0);
    do_cursor_layout = true;
  }

  if (do_cursor_layout) 
  {
    // Cursor layout is now valid, dont turn off the cursor flag
    disable_cursor_flag = false; 
    GlyphLayoutInfo *lptr = layoutvec.back();
    lptr->vertices_[0] = cursor_ll + cursor_down;
    lptr->vertices_[1] = cursor_ll + 2 * left + cursor_down ;
    lptr->vertices_[2] = cursor_ll + 2 * left + cursor_up;
    lptr->vertices_[3] = cursor_ll + cursor_up;
    for (int v = 0; v < 4; ++v) 
    {
      lptr->color_[v] = color_[v];
      EXTEND(lptr->vertices_[v]);
    }
  }


  if (disable_cursor_flag) {
    flags &= ~CURSOR;
  }
  BBox &bbox = layout_info.bbox_;
  bbox = BBox(Point(bbox_x1, bbox_y1, bbox_z1),
              Point(bbox_x2, bbox_y2, bbox_z2));
            
  if (!bbox.valid()) {
    total_glyphs = 0;
    return;
  }

  if (clamp_vals) 
  {
    GlyphLayoutVector::iterator iter = layoutvec.begin();
    for (unsigned int c = 0; c < total_glyphs; ++c, ++iter) 
    {
      for (int v = 0; v < 4; ++v) 
      {
        for (int o = 0; o < 3; ++o) 
        {
          (*iter)->vertices_[v](o) = Floor((*iter)->vertices_[v](o));
        }
      }
    }
  }
}

void
TextRenderer::move_text_to_anchor(TextLayoutInfo &layout_info, int flags)
{
  const Vector &left = layout_info.left_;
  const Vector &up = layout_info.up_;
  BBox &bbox = layout_info.bbox_;

  Vector diagonal = bbox.diagonal();
  Vector width = diagonal * left;
  Vector height = diagonal * up;

  Vector offset = Vector(0,0,0);
  switch (flags & ANCHOR_MASK) {
  case N:  offset = - width/2.0;                break;
  case E:  offset = - width       + height/2.0; break;
  case S:  offset = - width/2.0   + height;     break;
  case W:  offset =   height/2.0;               break;
  case NE: offset = - width;                    break;
  case SE: offset = - width       + height;     break;
  case SW: offset =   height;                   break;
  case C:  offset = - width/2.0   + height/2.0; break;
  default:
  case NW: break;
  }
  
  // The fix vector moves the text to the NW anchor
  Vector fix = bbox.min().x()*left + bbox.max().y()*up;
  
  layout_info.offset_ = offset + layout_info.anchor_.asVector() - fix;
}



double
TextRenderer::find_scale_factor(const std::vector<Point> &pnts, 
				const View &view) const
{
  //loop through and find the closest point to the eye
  Vector e2look = view.lookat() - view.eyep();
  e2look.normalize();
  std::vector<Point>::const_iterator iter = pnts.begin();
  double dist = DBL_MAX;
  Point close(0, 0, 0);
  while (iter != pnts.end()) 
  {
    const Point& p = *iter++;
    Vector e2p = p - view.eyep();
    double d = e2p.normalize();
    
    // make sure the point is in front of us.
    double dp = Dot(e2p, e2look);
    double angle  = acos(dp);
    // reject any point more than a 15 degrees off.
    if (dp < 0.0 ||  angle > .261799) { continue; }
    // cache closest point.
    if (d < dist) {
      dist = d;
      close = p;
    }
  }
  // calculate the scale factor for text at the closest point
  GLint vmat[4];
  glGetIntegerv(GL_VIEWPORT, vmat);

  double mmat[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, mmat);
  double pmat[16];
  glGetDoublev(GL_PROJECTION_MATRIX, pmat);
  
  double wx, wy, wz;
  Vector scene_up(mmat[1], mmat[5], mmat[9]);


  Point end = close + (scene_up * 20.0);

  gluProject(end.x(), end.y(), end.z(),
	     mmat, pmat, vmat,
	     &wx, &wy, &wz);

  Point send(wx, wy, wz);

  gluProject(close.x(), close.y(), close.z(),
	     mmat, pmat, vmat,
	     &wx, &wy, &wz);
  
  Point sbeg(wx, wy, wz);
  Vector slen = send - sbeg;

  return  20. / slen.length();
}


void
TextRenderer::render(const std::string &text, float x, float y, float z, 
		     const View &view, double sf, int flags)
{
  render_string_glyphs_to_texture(text);
  if (!texture_) return;

  glMatrixMode(GL_TEXTURE);
  glPushMatrix();
  glLoadIdentity();

  CHECK_OPENGL_ERROR();
  glColor4d(1.0, 0.0, 0.0, 1.0);


  double mmat[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, mmat);
  double pmat[16];
  glGetDoublev(GL_PROJECTION_MATRIX, pmat);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  CHECK_OPENGL_ERROR();

  glEnable(GL_TEXTURE_2D);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

//  glDepthMask(GL_FALSE); // no zbuffering for now.

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(0.0,50.0);

  Vector scene_up(mmat[1], mmat[5], mmat[9]);
  Point anchor(x, y, z);
  Vector bbv = view.eyep() - anchor;
  bbv.normalize();

  Vector scene_right = Cross(scene_up, bbv);
  scene_right.normalize();

  scene_up  *= sf;
  scene_right *= sf;

  TextLayoutInfo layout;
  layout.text_ = text;
  layout.left_ = scene_right;
  layout.up_ = scene_up;
  layout.flags_ = flags;
  layout_text(layout, false);

  layout.anchor_ = anchor;
  move_text_to_anchor(layout, flags);

  x = x + offset_*bbv.x();
  y = y + offset_*bbv.y();
  z = z + offset_*bbv.z();

  glTranslated(x, y, z);
  texture_->bind();
  for (unsigned int c = 0; c < layout.total_glyphs_; ++c) 
  {
    GlyphLayoutInfo &glyph_layout = *layout.layout_vector_[c];
    GlyphInfo *glyph = glyph_layout.glyph_info_;
    if (!glyph) {
      continue;
    }
    glyph->texture_->set_color(glyph_layout.color_);
    glyph->texture_->draw(4, glyph_layout.vertices_, glyph->tex_coords_);
  }
  
  glDisable(GL_TEXTURE_2D);
  glDepthMask(GL_TRUE); // turn zbuff back on.
  glDisable(GL_BLEND);

  glMatrixMode(GL_TEXTURE);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);               
  glPopMatrix();
  CHECK_OPENGL_ERROR();
}



void
TextRenderer::render(const std::string &text, Point& anchor, 
		    Vector& normal, double sf, int flags)
{
  render_string_glyphs_to_texture(text);
  if (!texture_) return;

  glMatrixMode(GL_TEXTURE);
  glPushMatrix();
  glLoadIdentity();

  CHECK_OPENGL_ERROR();
  glColor4d(1.0, 0.0, 0.0, 1.0);


  double mmat[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, mmat);
  double pmat[16];
  glGetDoublev(GL_PROJECTION_MATRIX, pmat);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  CHECK_OPENGL_ERROR();

  glEnable(GL_TEXTURE_2D);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

//  glDepthMask(GL_FALSE); // no zbuffering for now.

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(0.0,50.0);

  Vector scene_up(mmat[1], mmat[5], mmat[9]);
  normal.normalize();

  Vector scene_right = Cross(scene_up, normal);
  scene_right.normalize();

  scene_up  *= sf;
  scene_right *= sf;

  TextLayoutInfo layout;
  layout.text_ = text;
  layout.left_ = scene_right;
  layout.up_ = scene_up;
  layout.flags_ = flags;
  layout_text(layout, false);

  layout.anchor_ = anchor;
  move_text_to_anchor(layout, flags);

  anchor = (anchor+offset_*normal);

  glTranslated(anchor.x(), anchor.y(), anchor.z());
  texture_->bind();
  
  for (unsigned int c = 0; c < layout.total_glyphs_; ++c) 
  {
    GlyphLayoutInfo &glyph_layout = *layout.layout_vector_[c];
    GlyphInfo *glyph = glyph_layout.glyph_info_;
    if (!glyph) {
      continue;
    }
    glyph->texture_->set_color(glyph_layout.color_);
    glyph->texture_->draw(4, glyph_layout.vertices_, glyph->tex_coords_);
  }
  
  glDisable(GL_TEXTURE_2D);
  glDepthMask(GL_TRUE); // turn zbuff back on.
  glDisable(GL_BLEND);

  glMatrixMode(GL_TEXTURE);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);               
  glPopMatrix();
  CHECK_OPENGL_ERROR();
}


void
TextRenderer::render(const std::string &text, float x, float y, int flags)
{
  render_string_glyphs_to_texture(text);
  CHECK_OPENGL_ERROR();
  if (!texture_) return;

  CHECK_OPENGL_ERROR();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  CHECK_OPENGL_ERROR();
  glLoadIdentity();
  glScaled(2.0, 2.0, 2.0);
  glTranslated(-.5, -.5, -.5);
  glColor4d(1.0, 0.0, 0.0, 1.0);
  
  GLint gl_viewport[4];
  glGetIntegerv(GL_VIEWPORT, gl_viewport);
  float vw = gl_viewport[2];
  float vh = gl_viewport[3];
  glScaled(1/vw, 1/vh, 1.0);

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glShadeModel(GL_FLAT);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  CHECK_OPENGL_ERROR();

  int hash_flags = flags & ~ANCHOR_MASK;
  TextLayoutInfo &layout = layouts_[text][hash_flags];

  if (layout.text_ == text && layout.flags_ == hash_flags) {

  } else {
    //    std::cerr << "Layout for: " << text 
    //    << " hash   flags: " << hash_flags
    //    << " layout flags : " << layout.flags_ << endl;
    layout.text_ = text;
    layout.left_ = Vector(1,0,0);
    layout.up_ = Vector(0,1,0);
    layout.flags_ = hash_flags;
    layout_text(layout, true);
  }
  
  layout.anchor_ = Point(x,y,0.0);
  move_text_to_anchor(layout, flags);

  glTranslated(Round(layout.offset_.x()), 
               Round(layout.offset_.y()), 
               Round(layout.offset_.z()));

  texture_->bind();
  for (unsigned int c = 0; c < layout.total_glyphs_; ++c) 
  {
    GlyphLayoutInfo &glyph_layout = *layout.layout_vector_[c];
    GlyphInfo *glyph = glyph_layout.glyph_info_;
    if (!glyph) {
      continue;
    }
    glyph->texture_->set_color(glyph_layout.color_);
    glyph->texture_->draw(4, glyph_layout.vertices_, glyph->tex_coords_);
  }

  glDisable(GL_TEXTURE_2D);

  if ((flags & CURSOR) && layout.layout_vector_.size()) 
  {
    GlyphLayoutInfo &cursor_layout = *layout.layout_vector_.back();
    glColor4fv(cursor_layout.color_);
    glBegin(GL_QUADS);
    for (int v = 0; v < 4; ++v) {
      glVertex3d(cursor_layout.vertices_[v](0),
                 cursor_layout.vertices_[v](1),
                 cursor_layout.vertices_[v](2));
    }
    glEnd();
  }
               
               
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  CHECK_OPENGL_ERROR();
}


void 
TextRenderer::render_string_glyphs_to_texture(const std::string &text) 
{
  for (unsigned int c = 0; c < text.size(); ++c)
    if (glyphs_.find(text[c]) == glyphs_.end()) 
      glyphs_[text[c]] = render_glyph_to_texture(text[c]);
}



TextRenderer::GlyphInfo *
TextRenderer::render_glyph_to_texture(const wchar_t &character) 
{
  std::string str = " ";
  str[0] = character;
  GlyphInfo *glyph_info = new GlyphInfo();
  glyph_info->index_ = FT_Get_Char_Index(face_->ft_face_, character);
  if (glyph_info->index_ == 0) 
  {
    return 0;
  }


  FT_Error err;
  err = FT_Load_Glyph(face_->ft_face_, glyph_info->index_, FT_LOAD_DEFAULT);

  if (err) {
    std::cerr << "Freetype error: " << err << " \n";
    std::cerr << "FreeType Unable to Load Glyph: ";
    std::cerr << character;
    std::cerr <<"  " << to_string(int(character));
    std::cerr <<"  " <<  std::string(__FILE__) << to_string(__LINE__) << std::endl;
    return 0;
  }

  FT_Glyph glyph;
  err = FT_Get_Glyph(face_->ft_face_->glyph, &glyph);
  if (err) {
    throw ("FreeType Unable to Get Glyph: "+character+
           std::string(__FILE__)+to_string(__LINE__));
  }

  glyph_info->ft_metrics_ = face_->ft_face_->glyph->metrics;
  err = FT_Glyph_To_Bitmap(&glyph, FT_RENDER_MODE_NORMAL, 0, 1);
  if (err) {
    FT_Done_Glyph(glyph);
    std::cerr << "Freetype error: " << err << std::endl;
    return 0;
  }

  FT_BitmapGlyph bitmap_glyph = (FT_BitmapGlyph)(glyph);
  
  ASSERT(bitmap_glyph->bitmap.num_grays == 256);
  ASSERT(bitmap_glyph->bitmap.pixel_mode == FT_PIXEL_MODE_GRAY);

  int width = bitmap_glyph->bitmap.width;
  int height = bitmap_glyph->bitmap.rows;

  if (!texture_ || x_ + width > texture_->width()) {
    x_ = 0;
    y_ += height_ + 1;
    height_ = 0;
  }

  if (!texture_ || y_ + height > texture_->height()) {
    texture_ = new TextureObj(1, 256, 256);
    textures_.insert(texture_);
    y_ = 0;
    x_ = 0;
    height_ = 0;
  }

  int tex_width = texture_->width();
  int tex_height = texture_->height();

  height_ = Max(height_, height);
  glyph_info->texture_ = texture_;

  int v = 0;
  glyph_info->tex_coords_[v++] = x_/float(tex_width);
  glyph_info->tex_coords_[v++] = (y_+height)/float(tex_height);
  glyph_info->tex_coords_[v++] = (x_+width)/float(tex_width);
  glyph_info->tex_coords_[v++] = (y_+height)/float(tex_height);
  glyph_info->tex_coords_[v++] = (x_+width)/float(tex_width);
  glyph_info->tex_coords_[v++] = y_/float(tex_height);
  glyph_info->tex_coords_[v++] = x_/float(tex_width);
  glyph_info->tex_coords_[v++] = y_/float(tex_height);

  Nrrd *nrrd = texture_->nrrd_handle_->nrrd_;
  unsigned char *data = (unsigned char *)nrrd->data;

  // render glyph to texture data
  //  int pos;
  for (int y = 0; y < height; ++y) {
    int Y = y_+y;
    if (Y < 0 || Y >= tex_height) continue;
    for (int x = 0; x < width; ++x) {
      int X = x_+x;
      if (X < 0 || X >= tex_width) continue;
      if (!(X>=0 && X < int(nrrd->axis[1].size) &&
           Y>=0 && Y < int(nrrd->axis[2].size))) {
        std::cerr << "X: " << X
             << "Y: " << Y
             << "A1: " << nrrd->axis[1].size 
             << "A2: " << nrrd->axis[2].size;
      }
          
      data[Y*tex_width+X] =
        bitmap_glyph->bitmap.buffer[y*Abs(bitmap_glyph->bitmap.pitch)+x];
    }
  }

  x_ += width + 1;
  FT_Done_Glyph(glyph);
  texture_->set_dirty();
  return glyph_info;
}


void
TextRenderer::set_cursor_position(unsigned int pos) 
{
  cursor_position_ = pos;
}

void
TextRenderer::set_size(double size) 
{
  face_->set_points(size);
}


void
TextRenderer::set_offset(double offset) 
{
  offset_ = offset;
}

std::string
TextRenderer::get_default_font_name()
{
  static std::string fontname;
  std::string filename = "scirun.ttf";
  
  if (fontname.size() == 0)
  {
    const char *path_ptr = sci_getenv("SCIRUN_FONT_PATH");
    std::string path = (path_ptr ? std::string(path_ptr) : "");
    
    // Always search SCIRUN_SRCDIR/pixmaps for the font
    path = std::string(sci_getenv("SCIRUN_SRCDIR")) + "/pixmaps" + ":" +  path;

    // Search for the font 'filename' in the font 'path'
    std::string font_dir = findFileInPath(filename, path);

    // If is not found in the path, print error and return
    if (font_dir.empty()) 
    {
      std::cerr << "FontManager::load_face(" << filename << ") Error\n";
      std::cerr << "SCIRUN_FONT_PATH=" << path << std::endl;
      std::cerr << "Does not contain a file named " << filename << std::endl;
      return 0;
    }

    fontname = font_dir+"/"+filename;
  }
  
  return (fontname);
}


}
