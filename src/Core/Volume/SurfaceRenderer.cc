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
//    File   : SurfaceRenderer.cc
//    Author : Milan Ikits
//    Date   : Wed Jul  7 23:37:16 2004


#include <Core/Volume/SurfaceRenderer.h>
#include <Core/Geometry/BBox.h>
#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomResourceManager.h>


#include   <iostream>
#include   <string>


using std::cerr;
using std::endl;
using std::string;

namespace SCIRun {

#ifndef GL_TEXTURE_3D
#  define GL_TEXTURE_3D 0x806F
#endif

#ifdef _WIN32
#  undef min
#  undef max
#endif

SurfaceRenderer::SurfaceRenderer(SLIVR::Texture             *tex,
                                 SLIVR::ColorMap            *cmap1, 
                                 std::vector<SLIVR::ColorMap2*>  &cmap2,
                                 int                         tex_mem) : 
  SLIVR::SurfaceRenderer(tex, cmap1, cmap2, tex_mem),
  mutex_("SCIRun::SurfaceRenderer mutex")
{
  set_pending_delete_texture_callback(GeomResourceManager::delete_texture_id);
}


SurfaceRenderer::SurfaceRenderer(const SurfaceRenderer& copy) : 
  SLIVR::SurfaceRenderer(copy),
  mutex_("SCIRun::SurfaceRenderer mutex")
{
}


SurfaceRenderer::~SurfaceRenderer()
{
}


GeomObj*
SurfaceRenderer::clone()
{
  return new SurfaceRenderer(*this);
}


void
SurfaceRenderer::get_bounds(BBox& bb) 
{ 
  if (!tex_) return;
  double xmin;
  double ymin;
  double zmin;  
  double xmax;
  double ymax;
  double zmax;
  tex_->get_bounds(xmin, ymin, zmin, xmax, ymax, zmax); 

  bb.extend(Point(xmin, ymin, zmin));
  bb.extend(Point(xmax, ymax, zmax));
}


void
SurfaceRenderer::draw(DrawInfoOpenGL* di, Material* mat, double)
{
  if (!tex_) return;
  if(!pre_draw(di, mat, lighting_)) return;
  mutex_.lock();
  if (di->get_drawtype() == DrawInfoOpenGL::WireFrame)
  {
    draw_wireframe();
  }
  else
  {
    draw_surface();
  }
  mutex_.unlock();
}


} // namespace SCIRun
