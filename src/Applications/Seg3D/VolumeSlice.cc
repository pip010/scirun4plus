//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : VolumeSlice.cc
//    Author : McKay Davis
//    Date   : Fri Oct 13 15:35:55 2006

#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/VolumeOps.h>
#include <sci_comp_warn_fixes.h>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <typeinfo>
#include <iostream>
#include <sci_gl.h>
#include <algorithm>
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/ColorMappedNrrdTextureObj.h>
#include <Core/Geom/GeomSwitch.h>
#include <Core/Skinner/GeomSkinnerVarSwitch.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/TexSquare.h>
#include <Core/Math/MiscMath.h>
#include <Core/Math/MinMax.h>
#include <Core/Geom/GeomColorMappedNrrdTextureObj.h>
#include <Core/Skinner/Variables.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/SceneGraphEvent.h>
#include <Core/Util/FileUtils.h>


namespace SCIRun {

VolumeSlice::VolumeSlice(NrrdVolume *volume,
                         const SLIVR::Plane &plane,
                         NrrdDataHandle nrrd,
                         unsigned int label) :
  lock("Volume Slice"),
  ref_cnt(0),
  volume_(volume),
  nrrd_handle_(nrrd),
  outline_(0),
  texture_(0),
  geom_texture_(0),
  plane_(plane),
  label_(label),
  tex_dirty_(true),
  pos_(),
  xdir_(),
  ydir_()
{
  vector<int> sindex = volume_->world_to_index(plane_.project(SLIVR::Point(0,0,0)));
  unsigned int ax = get_axis();

  // Lower-Left origin corner for slice quad,
  // project into plane of window, to ensure volumes with different sample
  // spacings share same plane in view space.
  Point p = volume_->min(ax, sindex[ax]);
  SLIVR::Point p2 = plane_.project(SLIVR::Point(p.x(),p.y(),p.z()));

  pos_ = Point(p2.x(),p2.y(),p2.z());
  vector<int> index(volume_->nrrd_handle_->nrrd_->dim,0);

  Point origin = volume_->min(ax, 0);
  int primary = (ax == 1) ? 2 : 1;
  index[primary] = volume_->nrrd_handle_->nrrd_->axis[primary].size;
  xdir_ = volume_->index_to_world(index) - origin;
  index[primary] = 0;

  int secondary = (ax == 3) ? 2 : 3;
  index[secondary] = volume_->nrrd_handle_->nrrd_->axis[secondary].size;
  ydir_ = volume_->index_to_world(index) - origin;


  if (!nrrd_handle_.get_rep()) {
    extract_nrrd_slice_from_volume();
  }

  if (!nrrd_handle_.get_rep()) {
    return;
  }

  ColorMapHandle cmap = volume_->get_colormap();
  texture_ =
    new ColorMappedNrrdTextureObj(nrrd_handle_, cmap, volume_->label_);
  texture_->set_label_color(volume_->get_label_color());
  outline_ = new NrrdBitmaskOutline(nrrd_handle_);
  geom_texture_ = new GeomColorMappedNrrdTextureObj(texture_);
  tex_dirty_ = true;
}


void
VolumeSlice::set_tex_dirty()
{
  tex_dirty_ = true;
}


void
VolumeSlice::set_region_dirty(int x1, int y1, int x2, int y2, int border)
{
  texture_->apply_colormap(x1, y1, x2, y2, border);
  outline_->set_dirty();
}


unsigned int 
VolumeSlice::get_axis()
{
  ASSERT(volume_);
  return max_vector_magnitude_index(volume_->vector_to_index(plane_.normal()));
}


#define NRRD_EXEC(__nrrd_command__) \
  if (__nrrd_command__) { \
    char *err = biffGetDone(NRRD); \
    cerr << string("Error on line #")+to_string(__LINE__) + \
	    string(" executing nrrd command: \n")+ \
            string(#__nrrd_command__)+string("\n")+ \
            string("Message: ")+string(err); \
    free(err); \
    return; \
  }


void
VolumeSlice::extract_nrrd_slice_from_volume()
{
  vector<int> sindex = volume_->world_to_index(plane_.project(SLIVR::Point(0,0,0)));
  unsigned int ax = get_axis();

  int min_slice = sindex[ax];
  int max_slice = sindex[ax];
  int slice = sindex[ax];  
  if (slice < 0 || slice >= (int)volume_->nrrd_handle_->nrrd_->axis[ax].size)
  {
    return;
  }

  volume_->nrrd_handle_->lock.lock();
  nrrd_handle_ = new NrrdData();
  Nrrd *dst = nrrd_handle_->nrrd_;
  Nrrd *src = volume_->nrrd_handle_->nrrd_;
  ASSERT(src);

  if (min_slice != max_slice) {
    size_t *min = new size_t[src->dim];
    size_t *max = new size_t[src->dim];
    for (unsigned int i = 0; i < src->dim; i++) {
      min[i] = 0;
      max[i] = src->axis[i].size-1;
    }
    min[ax] = Min(min_slice, max_slice);
    max[ax] = Max(min_slice, max_slice);
    NrrdDataHandle tmp1_handle = new NrrdData();
    NRRD_EXEC(nrrdCrop(tmp1_handle->nrrd_, src, min, max));
    NRRD_EXEC(nrrdProject(dst, tmp1_handle->nrrd_, ax, 
                          nrrdMeasureMax, nrrdTypeDefault));
  } else {
    NRRD_EXEC(nrrdSlice(dst, src, ax, min_slice));
  }
  volume_->nrrd_handle_->lock.unlock();

  if (label_ && dst->type == nrrdTypeFloat) {
    nrrd_handle_ = VolumeOps::float_to_bit(nrrd_handle_, 0, label_);
  }
}


#if 0
static GLubyte stripe[4*32] = { 
  0x33, 0x33, 0x33, 0x33,
  0x66, 0x66, 0x66, 0x66,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x99, 0x99, 0x99, 0x99,
  
  0x33, 0x33, 0x33, 0x33,
  0x66, 0x66, 0x66, 0x66,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x99, 0x99, 0x99, 0x99,
  
  0x33, 0x33, 0x33, 0x33,
  0x66, 0x66, 0x66, 0x66,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x99, 0x99, 0x99, 0x99,
  
  0x33, 0x33, 0x33, 0x33,
  0x66, 0x66, 0x66, 0x66,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x99, 0x99, 0x99, 0x99,
  
  0x33, 0x33, 0x33, 0x33,
  0x66, 0x66, 0x66, 0x66,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x99, 0x99, 0x99, 0x99,
  
  0x33, 0x33, 0x33, 0x33,
  0x66, 0x66, 0x66, 0x66,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x99, 0x99, 0x99, 0x99,
  
  0x33, 0x33, 0x33, 0x33,
  0x66, 0x66, 0x66, 0x66,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x99, 0x99, 0x99, 0x99,
  
  0x33, 0x33, 0x33, 0x33,
  0x66, 0x66, 0x66, 0x66,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x99, 0x99, 0x99, 0x99,
};


static GLubyte stripe2[4*36] = { 
  0x33, 0x33, 0x33, 0x33,
  0x99, 0x99, 0x99, 0x99,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x66, 0x66, 0x66, 0x66,

  0x33, 0x33, 0x33, 0x33,
  0x99, 0x99, 0x99, 0x99,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x66, 0x66, 0x66, 0x66,

  0x33, 0x33, 0x33, 0x33,
  0x99, 0x99, 0x99, 0x99,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x66, 0x66, 0x66, 0x66,

  0x33, 0x33, 0x33, 0x33,
  0x99, 0x99, 0x99, 0x99,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x66, 0x66, 0x66, 0x66,

  0x33, 0x33, 0x33, 0x33,
  0x99, 0x99, 0x99, 0x99,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x66, 0x66, 0x66, 0x66,

  0x33, 0x33, 0x33, 0x33,
  0x99, 0x99, 0x99, 0x99,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x66, 0x66, 0x66, 0x66,

  0x33, 0x33, 0x33, 0x33,
  0x99, 0x99, 0x99, 0x99,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x66, 0x66, 0x66, 0x66,

  0x33, 0x33, 0x33, 0x33,
  0x99, 0x99, 0x99, 0x99,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x66, 0x66, 0x66, 0x66,

  0x33, 0x33, 0x33, 0x33,
  0x99, 0x99, 0x99, 0x99,
  0xCC, 0xCC, 0xCC, 0xCC,
  0x66, 0x66, 0x66, 0x66,
};
#endif

static GLubyte stripe3[4*36] = { 
  0x11, 0x11, 0x11, 0x11,
  0x22, 0x22, 0x22, 0x22,
  0x44, 0x44, 0x44, 0x44,
  0x88, 0x88, 0x88, 0x88,

  0x11, 0x11, 0x11, 0x11,
  0x22, 0x22, 0x22, 0x22,
  0x44, 0x44, 0x44, 0x44,
  0x88, 0x88, 0x88, 0x88,

  0x11, 0x11, 0x11, 0x11,
  0x22, 0x22, 0x22, 0x22,
  0x44, 0x44, 0x44, 0x44,
  0x88, 0x88, 0x88, 0x88,

  0x11, 0x11, 0x11, 0x11,
  0x22, 0x22, 0x22, 0x22,
  0x44, 0x44, 0x44, 0x44,
  0x88, 0x88, 0x88, 0x88,

  0x11, 0x11, 0x11, 0x11,
  0x22, 0x22, 0x22, 0x22,
  0x44, 0x44, 0x44, 0x44,
  0x88, 0x88, 0x88, 0x88,

  0x11, 0x11, 0x11, 0x11,
  0x22, 0x22, 0x22, 0x22,
  0x44, 0x44, 0x44, 0x44,
  0x88, 0x88, 0x88, 0x88,

  0x11, 0x11, 0x11, 0x11,
  0x22, 0x22, 0x22, 0x22,
  0x44, 0x44, 0x44, 0x44,
  0x88, 0x88, 0x88, 0x88,

  0x11, 0x11, 0x11, 0x11,
  0x22, 0x22, 0x22, 0x22,
  0x44, 0x44, 0x44, 0x44,
  0x88, 0x88, 0x88, 0x88,

  0x11, 0x11, 0x11, 0x11,
  0x22, 0x22, 0x22, 0x22,
  0x44, 0x44, 0x44, 0x44,
  0x88, 0x88, 0x88, 0x88,

};
    

static GLubyte stripe4[4*36] = { 
  0x11, 0x11, 0x11, 0x11,
  0x88, 0x88, 0x88, 0x88,
  0x44, 0x44, 0x44, 0x44,
  0x22, 0x22, 0x22, 0x22,

  0x11, 0x11, 0x11, 0x11,
  0x88, 0x88, 0x88, 0x88,
  0x44, 0x44, 0x44, 0x44,
  0x22, 0x22, 0x22, 0x22,

  0x11, 0x11, 0x11, 0x11,
  0x88, 0x88, 0x88, 0x88,
  0x44, 0x44, 0x44, 0x44,
  0x22, 0x22, 0x22, 0x22,

  0x11, 0x11, 0x11, 0x11,
  0x88, 0x88, 0x88, 0x88,
  0x44, 0x44, 0x44, 0x44,
  0x22, 0x22, 0x22, 0x22,

  0x11, 0x11, 0x11, 0x11,
  0x88, 0x88, 0x88, 0x88,
  0x44, 0x44, 0x44, 0x44,
  0x22, 0x22, 0x22, 0x22,

  0x11, 0x11, 0x11, 0x11,
  0x88, 0x88, 0x88, 0x88,
  0x44, 0x44, 0x44, 0x44,
  0x22, 0x22, 0x22, 0x22,

  0x11, 0x11, 0x11, 0x11,
  0x88, 0x88, 0x88, 0x88,
  0x44, 0x44, 0x44, 0x44,
  0x22, 0x22, 0x22, 0x22,

  0x11, 0x11, 0x11, 0x11,
  0x88, 0x88, 0x88, 0x88,
  0x44, 0x44, 0x44, 0x44,
  0x22, 0x22, 0x22, 0x22,

  0x11, 0x11, 0x11, 0x11,
  0x88, 0x88, 0x88, 0x88,
  0x44, 0x44, 0x44, 0x44,
  0x22, 0x22, 0x22, 0x22,

};


void
VolumeSlice::draw()
{
  if (!volume_->visible()) return;
  if (!nrrd_handle_.get_rep()) return;
  if (!texture_.get_rep() || !outline_.get_rep()) return;

  const bool use_stipple_ =
    volume_->painter_->get_vars()->get_bool("Painter::appearance::stipple");
  const bool fat_outlines_ = 
    volume_->painter_->get_vars()->get_bool("Painter::appearance::fatlines");

  glDisable(GL_CULL_FACE);
  glEnable(GL_BLEND);
  glDisable(GL_LIGHTING);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glShadeModel(GL_FLAT);
  CHECK_OPENGL_ERROR();

  if (tex_dirty_) {
    tex_dirty_ = false;
    texture_->set_label_color(volume_->get_label_color());
    texture_->set_clut_minmax(volume_->clut_min_, volume_->clut_max_);
    ColorMapHandle cmap = volume_->get_colormap();
    texture_->set_colormap(cmap);
    outline_->set_colormap(cmap);
  }
   
  if (texture_.get_rep())
  {
    double opacity = volume_->get_opacity();

    if (volume_->label_)
    {
      if (use_stipple_)
      {
        const int depth = volume_->bit();
        GLubyte *pattern = stripe4;
        switch (depth % 4) {
        case 0: pattern = stripe4; break;
        case 1: pattern = stripe3+8; break;
        case 2: pattern = stripe4+8; break;
        default:
        case 3: pattern = stripe3; break;
        }

        glPolygonStipple(pattern);
        glEnable(GL_POLYGON_STIPPLE);
      }
      else
      {
        opacity = 0.5;
      }
    }

    texture_->set_opacity(opacity);
    texture_->set_coords(pos_, xdir_, ydir_);
    texture_->draw_quad();

    if (use_stipple_ && volume_->label_) {
      glDisable(GL_POLYGON_STIPPLE);
    }    
  }
  
  if (volume_->label_ && outline_.get_rep()) {
    outline_->set_line_style(fat_outlines_);
    outline_->set_coords(pos_, xdir_, ydir_);
    outline_->draw_lines(3.0, volume_->label_);
  }
}


}
