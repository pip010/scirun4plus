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
 *  ColorMappedNrrdTextureObj.cc
 *
 *  Written by:
 *   McKay Davis
 *   School of Computing
 *   University of Utah
 *   November, 2005
 *
 */

#include <sci_gl.h>
#include <sci_glx.h>

#include <Core/Util/StringUtil.h>
#include <Core/Datatypes/NrrdData.h>
#include <Core/Geom/ColorMappedNrrdTextureObj.h>
#include <Core/Math/MiscMath.h>
#include <Core/Geom/GeomResourceManager.h>

#include <teem/air.h>
#include <limits.h>

#include <Core/Thread/MutexPool.h>

namespace SCIRun {


ColorMappedNrrdTextureObj::ColorMappedNrrdTextureObj(NrrdDataHandle &nrrdh,
                                                     ColorMapHandle &cmap,
                                                     unsigned int label) :
  UsedWithLockingHandle<Mutex>("ColorMappedNrrdTextureObj"),
  colormap_(0),
  colormap_rgba_(0),
  nrrd_handle_(nrrdh),
  dirty_(true),
  clut_min_(0.0), 
  clut_max_(1.0),
  data_(0),
  label_(label),
  opacity_(1.0),
  min_(0.0, 0.0, 0.0),
  xdir_(1.0, 0.0, 0.0),
  ydir_(0.0, 1.0, 0.0),
  thresholding_(false)
{
  ASSERT(nrrdh.get_rep());
  ASSERT(nrrdh->nrrd_);
  ASSERT(nrrdh->nrrd_->axis[0].size == 1);
  size_t dim = nrrdh->nrrd_->axis[1].size * nrrdh->nrrd_->axis[2].size;
  if (label_)
  {
    data_ = new unsigned char[dim];
    ASSERT(data_);
  }
  else
  {
    data_ = new unsigned char[dim * 4];
    set_colormap(cmap);
    ASSERT(data_);
  }

  // 64 and false are least common denominator graphics cards.
  const unsigned int max_texture_size = 64;
  // MacOS and maybe other machines are reporting this to GLEW incorrrectly.
#ifdef __APPLE__
  const bool non_power_2_textures = false;
#else
  const bool non_power_2_textures = GLEW_ARB_texture_non_power_of_two;
#endif

  unsigned int s, e;
  s = 0; e = 0;
  unsigned int mw = nrrdh->nrrd_->axis[1].size;
  while (e < mw)
  {
    const unsigned int left = mw - e;
    if (max_texture_size < left)
    {
      // Max texture sized blocks
      s = e;
      e += max_texture_size;
    }
    else if (LargestPowerOf2(left) == left || non_power_2_textures)
    {
      // The rest if we don't have to worry about power of 2s.
      s = e;
      e = mw;
    }
    else if (left > Min((unsigned int)32, LargestPowerOf2(mw)))
    {
      // Dice up the remainder into large power of 2 blocks.
      s = e;
      e += LargestPowerOf2(left);
    }
    else
    {
      // The non-power-of-two leftover.
      s = mw - Pow2(left);
      e = mw;
    }
    if (e != s)
    {
      xdiv_.push_back(std::pair<int, int>(s, e));
    }
    else break;
  }

  s = 0; e = 0;
  mw = nrrdh->nrrd_->axis[2].size;
  while (e < mw)
  {
    const unsigned int left = mw - e;
    if (max_texture_size < left)
    {
      // Max texture sized blocks
      s = e;
      e += max_texture_size;
    }
    else if (LargestPowerOf2(left) == left || non_power_2_textures)
    {
      // The rest if we don't have to worry about power of 2s.
      s = e;
      e = mw;
    }
    else if (left > Min((unsigned int)32, LargestPowerOf2(mw)))
    {
      // Dice up the remainder into large power of 2 blocks.
      s = e;
      e += LargestPowerOf2(left);
    }
    else
    {
      // The non-power-of-two leftover.
      s = mw - Pow2(left);
      e = mw;
    }
    if (e != s)
    {
      ydiv_.push_back(std::pair<int, int>(s, e));
    }
    else break;
  }

  texture_id_.resize(xdiv_.size() * ydiv_.size(), 0);
  texture_id_dirty_.resize(xdiv_.size() * ydiv_.size(), true);
}


ColorMappedNrrdTextureObj::~ColorMappedNrrdTextureObj()
{
  if (data_) { delete[] data_; }
  if (colormap_rgba_) { delete[] colormap_rgba_; }
  for (unsigned int t = 0; t < texture_id_.size(); ++t)
  {
    GeomResourceManager::delete_texture_id(texture_id_[t]);
  }
}


void
ColorMappedNrrdTextureObj::set_clut_minmax(float min, float max)
{
  clut_min_ = min;
  clut_max_ = max;
  dirty_ = true;
}


void
ColorMappedNrrdTextureObj::set_colormap(ColorMapHandle &cmap)
{
  if (label_ || thresholding_) return; // Uses label_color_ instead.

  if (colormap_.get_rep() != cmap.get_rep())
  {
    colormap_ = cmap;
    dirty_ = true;

    if (colormap_rgba_) { delete [] colormap_rgba_; }
    const int ncolors = colormap_->resolution();
    const float *rgba = colormap_->get_rgba();
    colormap_rgba_ = new unsigned char[ncolors * 4];
    for (int i = 0; i < ncolors * 4; i++)
    {
      colormap_rgba_[i] = (unsigned char)(rgba[i] * 255.0);
    }
  }
}


template <class T>
void
ColorMappedNrrdTextureObj::apply_colormap_to_data(unsigned char *dst, 
                                                  T *src, 
                                                  int row_width,
                                                  int region_start,
                                                  int region_width,
                                                  int region_height,
                                                  const unsigned char *rgba, 
                                                  int ncolors,
                                                  float scale, 
                                                  float bias)
{
  const int maxcolor = ncolors-1;
  for (int row = 0; row < region_height; ++row) {
    int start_pos = region_start + row*row_width;
    int end_pos = start_pos + region_width;
    for (int pos = start_pos; pos < end_pos; ++pos) {
      float val = float(src[pos])*scale+bias;
      int color = Round(Clamp(val, 0.0f, 1.0f)*maxcolor);
      memcpy(dst + pos * 4, rgba + color * 4, 4);
    }
  }
}


template <class T>
void
ColorMappedNrrdTextureObj::apply_label_to_data(unsigned char *dst, 
                                               T *src, 
                                               int row_width,
                                               int region_start,
                                               int region_width,
                                               int region_height,
                                               unsigned int mask)
{
  for (int row = 0; row < region_height; ++row)
  {
    int start_pos = region_start + row*row_width;
    int end_pos = start_pos + region_width;
    for (int pos = start_pos; pos < end_pos; ++pos)
    {
      dst[pos] = (src[pos] & mask) ? 0xff : 0x00;
    }
  }
}


template <class T> 
void
ColorMappedNrrdTextureObj::apply_threshold_to_data(unsigned char *dst,
                                                   T *src,
                                                   int row_width,
                                                   int region_start,
                                                   int region_width,
                                                   int region_height)
{
  for (int row = 0; row < region_height; ++row)
  {
    int start_pos = region_start + row*row_width;
    int end_pos = start_pos + region_width;
    for (int pos = start_pos; pos < end_pos; ++pos)
    {
      dst[pos] = (src[pos] >= lower_threshold_ && src[pos] <= upper_threshold_) ? 0xff : 0x00;
    }
  }
}


void
ColorMappedNrrdTextureObj::apply_colormap(int x1, int y1, int x2, int y2,
                                          int border)
{
  Nrrd *nrrd = nrrd_handle_->nrrd_;

  if (x1 > x2) std::swap(x1,x2);
  if (y1 > y2) std::swap(y1,y2);
  x1 = Clamp(x1-border,   0, nrrd->axis[1].size);
  x2 = Clamp(x2+border+1, 0, nrrd->axis[1].size);
  y1 = Clamp(y1-border,   0, nrrd->axis[2].size);
  y2 = Clamp(y2+border+1, 0, nrrd->axis[2].size);

  const int row_width = nrrd->axis[1].size;
  const int region_start = row_width * y1 + x1;
  const int region_wid = x2 - x1;
  const int region_hei = y2 - y1;

  if (thresholding_)
  {
    switch (nrrd->type)
    {
    case nrrdTypeChar:
      apply_threshold_to_data(data_, (char *)nrrd->data,
                              row_width, region_start,
                              region_wid, region_hei);
      break;

    case nrrdTypeUChar:
      apply_threshold_to_data(data_, (unsigned char *)nrrd->data,
                              row_width, region_start,
                              region_wid, region_hei);
      break;

    case nrrdTypeShort:
      apply_threshold_to_data(data_, (short *)nrrd->data,
                              row_width, region_start,
                              region_wid, region_hei);
      break;

    case nrrdTypeUShort:
      apply_threshold_to_data(data_, (unsigned short *)nrrd->data,
                              row_width, region_start,
                              region_wid, region_hei);
      break;

    case nrrdTypeInt:
      apply_threshold_to_data(data_, (int *)nrrd->data,
                              row_width, region_start,
                              region_wid, region_hei);
      break;
    
    case nrrdTypeUInt:
      apply_threshold_to_data(data_, (unsigned int *)nrrd->data,
                              row_width, region_start,
                              region_wid, region_hei);
      break;

    case nrrdTypeFloat:
      apply_threshold_to_data(data_, (float *)nrrd->data,
                              row_width, region_start,
                              region_wid, region_hei);
      break;

    case nrrdTypeDouble:
      apply_threshold_to_data(data_, (double *)nrrd->data,
                              row_width, region_start,
                              region_wid, region_hei);
      break;

    default: throw "Unsupported data type for thresholding: " + 
               to_string(nrrd->type);
    }
  }
  else if (label_)
  {
    switch (nrrd->type) {
    case nrrdTypeChar:
      apply_label_to_data(data_, (char *)nrrd->data,
                          row_width, region_start,
                          region_wid, region_hei,
                          label_);
      break;

    case nrrdTypeUChar:
      apply_label_to_data(data_, (unsigned char *)nrrd->data,
                          row_width, region_start, 
                          region_wid, region_hei,
                          label_);
      break;

    case nrrdTypeShort:
      apply_label_to_data(data_, (short *)nrrd->data,
                          row_width, region_start, 
                          region_wid, region_hei,
                          label_);
      break;

    case nrrdTypeUShort:
      apply_label_to_data(data_, (unsigned short *)nrrd->data,
                          row_width, region_start, 
                          region_wid, region_hei,
                          label_);
      break;
    
    case nrrdTypeInt:
      apply_label_to_data(data_, (int *)nrrd->data,
                          row_width, region_start, 
                          region_wid, region_hei,
                          label_);
      break;
    
    case nrrdTypeUInt:
      apply_label_to_data(data_, (unsigned int *)nrrd->data,
                          row_width, region_start, 
                          region_wid, region_hei,
                          label_);
      break;

    case nrrdTypeLLong:
      apply_label_to_data(data_, (signed long long *)nrrd->data,
                          row_width, region_start, 
                          region_wid, region_hei,
                          label_);
      break;

    case nrrdTypeULLong:
      apply_label_to_data(data_, (unsigned long long *)nrrd->data,
                          row_width, region_start, 
                          region_wid, region_hei,
                          label_);
      break;

    default: throw "Unsupported data type for label: " + 
               to_string(nrrd->type);
    }
  }
  else
  {
    const int ncolors = colormap_->resolution();

    const float range = (clut_max_ - clut_min_);
    const float scale = (range > 0.00001) ? (1.0 / range) : 1.0;
    const float bias =  (range > 0.00001) ? -clut_min_*scale : 0.0;

    switch (nrrd->type) {
    case nrrdTypeChar:
      apply_colormap_to_data(data_, (char *)nrrd->data,
                             row_width, region_start, region_wid, region_hei,
                             colormap_rgba_, ncolors, scale, bias);
      break;

    case nrrdTypeUChar:
      apply_colormap_to_data(data_, (unsigned char *)nrrd->data,
                             row_width, region_start, region_wid, region_hei,
                             colormap_rgba_, ncolors, scale, bias);
      break;

    case nrrdTypeShort:
      apply_colormap_to_data(data_, (short *)nrrd->data,
                             row_width, region_start, region_wid, region_hei,
                             colormap_rgba_, ncolors, scale, bias);
      break;

    case nrrdTypeUShort:
      apply_colormap_to_data(data_, (unsigned short *)nrrd->data,
                             row_width, region_start, region_wid, region_hei,
                             colormap_rgba_, ncolors, scale, bias);
      break;
    
    case nrrdTypeInt:
      apply_colormap_to_data(data_, (int *)nrrd->data,
                             row_width, region_start, region_wid, region_hei,
                             colormap_rgba_, ncolors, scale, bias);
      break;
    
    case nrrdTypeUInt:
      apply_colormap_to_data(data_, (unsigned int *)nrrd->data,
                             row_width, region_start, region_wid, 
                             region_hei, colormap_rgba_, ncolors, scale, bias);
      break;

    case nrrdTypeLLong:
      apply_colormap_to_data(data_, 
                             (signed long long *)nrrd->data,
                             row_width, region_start, region_wid, region_hei,
                             colormap_rgba_, ncolors, scale, bias);
      break;

    case nrrdTypeULLong:
      apply_colormap_to_data(data_, 
                             (unsigned long long *)nrrd->data,
                             row_width, region_start, region_wid, region_hei,
                             colormap_rgba_, ncolors, scale, bias);
      break;
    
    case nrrdTypeFloat:
      apply_colormap_to_data(data_, (float *)nrrd->data,
                             row_width, region_start,
                             region_wid, region_hei,
                             colormap_rgba_, ncolors, scale, bias);
      break;

    case nrrdTypeDouble:
      apply_colormap_to_data(data_, (double *)nrrd->data,
                             row_width, region_start, region_wid, region_hei,
                             colormap_rgba_, ncolors, scale, bias);
      break;

    default: throw "Unsupported data type: " + 
               to_string(nrrd->type);
    }
  }
  
  Point min(x1, y1, 0), max(x2, y2, 1);
  BBox bbox(min,max);
  for (unsigned int i = 0; i < xdiv_.size(); i++)
  {
    for (unsigned int j = 0; j < ydiv_.size(); j++)
    {
      Point min2(xdiv_[i].first,  ydiv_[j].first,  0);
      Point max2(xdiv_[i].second, ydiv_[j].second, 1);
      BBox bbox2(min2, max2);
      if (bbox.overlaps(bbox2)) {
        texture_id_dirty_[j * xdiv_.size() + i] = true;
      }
    }
  }

  dirty_ = false;
}


void
ColorMappedNrrdTextureObj::set_threshold(double lower, double upper)
{
  thresholding_ = true;
  if (lower_threshold_ == lower && upper_threshold_ == upper) return;
  lower_threshold_ = lower;
  upper_threshold_ = upper;
  dirty_ = true;
}

  
bool
ColorMappedNrrdTextureObj::bind(int x, int y)
{
  const int pos = y * xdiv_.size() + x;
  const bool bound = glIsTexture(texture_id_[pos]);

  if (!bound) {
    glGenTextures(1, (GLuint *)&texture_id_[pos]);
    CHECK_OPENGL_ERROR();
  }

  glBindTexture(GL_TEXTURE_2D, texture_id_[pos]);
  CHECK_OPENGL_ERROR();  

  if (bound && !dirty_ && !texture_id_dirty_[pos]) {
    return true;
  }

  Nrrd *nrrd = nrrd_handle_->nrrd_;

  if (dirty_) {
    apply_colormap(0, 0, nrrd->axis[1].size, nrrd->axis[2].size);
  }

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, nrrd->axis[1].size);  
  glPixelStorei(GL_UNPACK_IMAGE_HEIGHT, nrrd->axis[2].size);  
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, xdiv_[x].first);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, ydiv_[y].first);
  CHECK_OPENGL_ERROR();

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  CHECK_OPENGL_ERROR();  
  
  if (label_ || thresholding_)
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA,
                 xdiv_[x].second - xdiv_[x].first, 
                 ydiv_[y].second - ydiv_[y].first,
                 0, GL_ALPHA, GL_UNSIGNED_BYTE, data_);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
                 xdiv_[x].second - xdiv_[x].first, 
                 ydiv_[y].second - ydiv_[y].first,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, data_);
  }
  CHECK_OPENGL_ERROR();

  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(GL_UNPACK_IMAGE_HEIGHT, 0); 
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0); 

  CHECK_OPENGL_ERROR();  
  texture_id_dirty_[pos] = false;

  return true;
}


void
ColorMappedNrrdTextureObj::set_coords(const Point &min, 
                                      const Vector &xdir, 
                                      const Vector &ydir)
{
  min_ = min;
  xdir_ = xdir;
  ydir_ = ydir;
}


void
ColorMappedNrrdTextureObj::get_bounds(BBox &bbox)
{
  bbox.extend(min_);
  bbox.extend(min_ + xdir_);
  bbox.extend(min_ + ydir_);
}

  
void
ColorMappedNrrdTextureObj::draw_quad(bool noalpha)
{
  double opacity = opacity_;
  if (noalpha) opacity = 1.0;

  glEnable(GL_TEXTURE_2D);
  glColor4d(opacity, opacity, opacity, opacity);
  if (label_ || thresholding_)
  {
    glDisable(GL_LIGHTING);
    glColor4d(label_color_.r() * opacity,
              label_color_.g() * opacity,
              label_color_.b() * opacity,
              opacity);
  }

  unsigned int yoff = 0;
  for (unsigned int y = 0; y < ydiv_.size(); ++y)
  {
    unsigned int xoff = 0;
    for (unsigned int x = 0; x < xdiv_.size(); ++x)
    {
      if (!bind(x,y)) continue;

      float x1 = xoff / float(nrrd_handle_->nrrd_->axis[1].size);
      float y1 = yoff / float(nrrd_handle_->nrrd_->axis[2].size);
      float x2 = xdiv_[x].second / float(nrrd_handle_->nrrd_->axis[1].size);
      float y2 = ydiv_[y].second / float(nrrd_handle_->nrrd_->axis[2].size);

      float tx = ((xoff - xdiv_[x].first) / 
                  float(xdiv_[x].second - xdiv_[x].first));
      float ty = ((yoff - ydiv_[y].first) / 
                  float(ydiv_[y].second - ydiv_[y].first));
      xoff = xdiv_[x].second;

      glBegin(GL_QUADS);
      
      glTexCoord2d(tx, ty); 
      Point p = min_ + x1*xdir_ + y1*ydir_;
      glVertex3f(p.x(), p.y(), p.z());
      
      glTexCoord2d(1.0, ty);
      p = min_ + x2*xdir_ + y1*ydir_;
      glVertex3f(p.x(), p.y(), p.z());
      
      glTexCoord2d(1.0, 1.0);
      p = min_ + x2*xdir_ + y2*ydir_;
      glVertex3f(p.x(), p.y(), p.z());
      
      glTexCoord2d(tx, 1.0);
      p = min_ + x1*xdir_ + y2*ydir_;
      glVertex3f(p.x(), p.y(), p.z());
      glEnd();

      CHECK_OPENGL_ERROR(); 
    }
    yoff = ydiv_[y].second;
  }
  glDisable(GL_TEXTURE_2D);
}


} // end namespace SCIRun
