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
 *  TextureObj.cc
 *
 *  Written by:
 *   McKay Davis
 *   SCI Institute
 *   University of Utah
 *
 */

#include <sci_gl.h>
#include <sci_glx.h>
#include <string.h>

#include <Core/Geom/TextureObj.h>
#include <Core/Math/MiscMath.h>
#include <Core/Util/StringUtil.h>
#include <slivr/ShaderProgramARB.h>
#include <Core/Geom/GeomResourceManager.h>

namespace SCIRun {

TextureObj::TextureObj(NrrdDataHandle &nrrd_handle) :
  nrrd_handle_(0),
  width_(-1),
  height_(-1),
  dirty_(true),
  texture_id_(0)
{
  set_color(1.0, 1.0, 1.0, 1.0);
  set_nrrd(nrrd_handle);
}


TextureObj::TextureObj(int components, int x, int y) :
  nrrd_handle_(0),
  width_(-1),
  height_(-1),
  dirty_(true),
  texture_id_(0)
{
  set_color(1.0, 1.0, 1.0, 1.0);

  size_t size[NRRD_DIM_MAX];
  size[0] = components;
  size[1] = x;
  size[2] = y;
  NrrdDataHandle nrrd_handle = new NrrdData();
  nrrdAlloc_nva(nrrd_handle->nrrd_, nrrdTypeUChar, 3, size);
    
  set_nrrd(nrrd_handle);
}

int
TextureObj::width() { return width_; }

int             
TextureObj::height() { return height_; }

void            
TextureObj::set_dirty() { dirty_ = true; }

void
TextureObj::set_nrrd(NrrdDataHandle &nrrd_handle) {
  if (!nrrd_handle.get_rep() ||
      !nrrd_handle->nrrd_ ||
       nrrd_handle->nrrd_->dim != 3) 
    throw "TextureObj::set_nrrd(NrrdDataHandle &nrrd_handle): nrrd not valid";
  nrrd_handle_ = nrrd_handle;
  width_  = nrrd_handle_->nrrd_->axis[1].size;
  height_ = nrrd_handle_->nrrd_->axis[2].size;

  pad_to_power_of_2();
}


TextureObj::~TextureObj()
{
  GeomResourceManager::delete_texture_id(texture_id_);
  nrrd_handle_ = 0;
}


void
TextureObj::set_color(float r, float g, float b, float a)
{
  color_[0] = r;
  color_[1] = g;
  color_[2] = b;
  color_[3] = a;
}

void
TextureObj::set_color(float rgba[4])
{
  set_color(rgba[0], rgba[1], rgba[2], rgba[3]);
}



void
TextureObj::resample_to_power_of_2()
{
  //  if (ShaderProgramARB::texture_non_power_of_two()) 
  //    return;
  if (!nrrd_handle_.get_rep() || !nrrd_handle_->nrrd_)
    return;
  if (IsPowerOf2(nrrd_handle_->nrrd_->axis[1].size) &&
      IsPowerOf2(nrrd_handle_->nrrd_->axis[2].size)) 
    return;
  NrrdDataHandle nout_handle = new NrrdData();

  NrrdResampleInfo *info = nrrdResampleInfoNew();
  double p[NRRD_KERNEL_PARMS_NUM];
  memset(p, 0, NRRD_KERNEL_PARMS_NUM * sizeof(double));
  p[0] = 1.0;
  info->boundary = nrrdBoundaryBleed;
  info->type = nrrd_handle_->nrrd_->type;
  info->renormalize = AIR_FALSE;

  info->kernel[0] = 0;
  info->kernel[1] = nrrdKernelBox;
  info->kernel[2] = nrrdKernelBox;
  info->samples[0] = nrrd_handle_->nrrd_->axis[0].size;
  info->samples[1] = width_ = Pow2(nrrd_handle_->nrrd_->axis[1].size);
  info->samples[2] = height_ = Pow2(nrrd_handle_->nrrd_->axis[2].size);
  memcpy(info->parm[0], p, NRRD_KERNEL_PARMS_NUM * sizeof(double));
  memcpy(info->parm[1], p, NRRD_KERNEL_PARMS_NUM * sizeof(double));
  memcpy(info->parm[2], p, NRRD_KERNEL_PARMS_NUM * sizeof(double));

  Nrrd * nin = nrrd_handle_->nrrd_;

  for (unsigned int a = 0; a < nin->dim; ++a) {
    nrrdAxisInfoMinMaxSet(nin, a, nin->axis[a].center ? 
                          nin->axis[a].center : nrrdCenterCell);
    
    info->min[a] = nin->axis[a].min;
    info->max[a] = nin->axis[a].max;
  }


  if (nrrdSpatialResample(nout_handle->nrrd_, nrrd_handle_->nrrd_, info)) {
    char *err = biffGetDone(NRRD);
    std::string error = std::string("Trouble resampling: ") + err;
    free (err);
    throw error;
  }
  nrrd_handle_ = nout_handle;
  set_dirty();
}



void
TextureObj::pad_to_power_of_2()
{
  //  if (ShaderProgramARB::texture_non_power_of_two()) 
  //    return;
  if (!nrrd_handle_.get_rep() || !nrrd_handle_->nrrd_)
    return;
  if (IsPowerOf2(nrrd_handle_->nrrd_->axis[1].size) &&
      IsPowerOf2(nrrd_handle_->nrrd_->axis[2].size)) 
    return;
  NrrdDataHandle nout_handle = new NrrdData();

  const int pady = (Pow2(nrrd_handle_->nrrd_->axis[2].size) - 
                    nrrd_handle_->nrrd_->axis[2].size);

  
  ptrdiff_t minp[3] = { 0, 0, -pady };
  ptrdiff_t maxp[3] = { nrrd_handle_->nrrd_->axis[0].size -1,
                        Pow2(nrrd_handle_->nrrd_->axis[1].size) - 1, 
                        nrrd_handle_->nrrd_->axis[2].size -1 };
    
  if (nrrdPad_nva(nout_handle->nrrd_, nrrd_handle_->nrrd_,
                  minp, maxp, nrrdBoundaryBleed, 1)) {
    char *err = biffGetDone(NRRD);
    std::string error = std::string("Trouble resampling: ") + err;
    free (err);
    throw error;
  }

  nrrd_handle_ = nout_handle;

  set_dirty();
}




bool
TextureObj::bind()
{
  CHECK_OPENGL_ERROR();

  if (!nrrd_handle_.get_rep() || !nrrd_handle_->nrrd_) return false;

  const bool bound = glIsTexture(texture_id_);

  if (!bound)
    glGenTextures(1, (GLuint*)&texture_id_);

  glBindTexture(GL_TEXTURE_2D, texture_id_);
  CHECK_OPENGL_ERROR();
  if (bound && !dirty_) return true;
  dirty_ = false;
  Nrrd nrrd = *nrrd_handle_->nrrd_;
  int prim = 1;
  GLenum pixtype;
  if (nrrd.axis[0].size == 1) 
    pixtype = GL_ALPHA;  
  else if (nrrd.axis[0].size == 2) 
    pixtype = GL_LUMINANCE_ALPHA;
  else if (nrrd.axis[0].size == 3) 
    pixtype = GL_RGB;
  else if (nrrd.axis[0].size == 4) 
    pixtype = GL_RGBA;
  else {
    prim = 0;
    pixtype = GL_ALPHA;
  }
  GLenum type = 0;
  switch (nrrd.type) {
  case nrrdTypeChar: type = GL_BYTE; break;
  case nrrdTypeUChar: type = GL_UNSIGNED_BYTE; break;
  case nrrdTypeShort: type = GL_SHORT; break;
  case nrrdTypeUShort: type = GL_UNSIGNED_SHORT; break;	
  case nrrdTypeInt: type = GL_INT; break;
  case nrrdTypeUInt: type = GL_UNSIGNED_INT; break;
  case nrrdTypeFloat: type = GL_FLOAT; break;
  default: throw std::string("Cannot bind nrrd of unknown type"); break;
  }
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  CHECK_OPENGL_ERROR();

  glTexImage2D(GL_TEXTURE_2D, 0, pixtype,
	       nrrd.axis[prim].size, nrrd.axis[prim+1].size, 
	       0, pixtype, type, nrrd.data);
  CHECK_OPENGL_ERROR();
  return true;
}


void
TextureObj::draw(int n, Point *vertices, float *tex_coords)
{
  CHECK_OPENGL_ERROR();
  ASSERT(n == 4);
  if (bind()) {
    glEnable(GL_TEXTURE_2D);
  } else {
    glDisable(GL_TEXTURE_2D);
    return;
  }

  glColor4fv(color_);            
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

  if (n == 4) {
    float x_scale = float(width_ )/float(nrrd_handle_->nrrd_->axis[1].size);
    float y_scale = float(height_)/float(nrrd_handle_->nrrd_->axis[2].size);

    glBegin(GL_QUADS);    
    for (int v = 0; v < 4; ++v) {
      glTexCoord2f(tex_coords[v*2+0] * x_scale, 
                   tex_coords[v*2+1] * y_scale);
      glVertex3dv(&vertices[v](0));
    }
    glEnd();
  }

  glDisable(GL_TEXTURE_2D);
  CHECK_OPENGL_ERROR();
}

}
