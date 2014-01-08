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
//    File   : VolumeRenderer.h
//    Author : Milan Ikits
//    Date   : Sat Jul 10 11:26:26 2004

#ifndef VolumeRenderer_h
#define VolumeRenderer_h


#include <slivr/Texture.h>
#include <slivr/VolumeRenderer.h>

#include <Core/Containers/LockingHandle.h>

#include <Core/Volume/Texture.h>
#include <Core/Volume/ColorMap2.h>
#include <Core/Datatypes/ColorMap.h>

#include <Core/Geom/GeomObj.h>
#include <Core/Volume/share.h>

namespace SCIRun {

class VolShaderFactory;

class SCISHARE VolumeRenderer : public SLIVR::VolumeRenderer, public GeomObj
{
  public:
    VolumeRenderer(TextureHandle tex, 
                   ColorMapHandle cmap1, 
                   std::vector<ColorMap2Handle>  &cmap2,
                   std::vector<Plane> &planes,
                   int tex_mem);

    VolumeRenderer(const VolumeRenderer&);
    virtual ~VolumeRenderer();

    virtual void draw(DrawInfoOpenGL*, Material*, double time);
    virtual void get_bounds(BBox& bb);

    virtual GeomObj* clone();

    virtual void io(Piostream&);
    static PersistentTypeID type_id;

    // Protect functions for raw point access
    // These functions store the handle to the objects inside the
    // object, which will force the deletion of the objects after
    // the VolumeRenderer is deleted.
    // These functions overload the functions that change these objects
    void set_texture(TextureHandle tex)
    {
      tex_handle_ = tex;
      SLIVR::TextureRenderer::set_texture(tex.get_rep());
    }
    
    void set_colormap1(ColorMapHandle cmap)
    {
      colormap_handle_ = cmap;
      SLIVR::TextureRenderer::set_colormap1(cmap.get_rep());
    }
    
    void set_colormap2(const std::vector<ColorMap2Handle> cmap2)
    {
      colormap2_handles_ = cmap2;
      std::vector<SLIVR::ColorMap2*> cmap2_ptrs(cmap2.size());
      for (size_t j=0;j<cmap2_ptrs.size();j++)
        cmap2_ptrs[j] = cmap2[j].get_rep();
      SLIVR::TextureRenderer::set_colormap2(cmap2_ptrs);
    }
    
    void set_planes(const std::vector<Plane>& planes)
    {
      planes_.clear();
      planes_.resize(planes.size());
        
      std::vector<SLIVR::Plane*> plane_ptrs(planes.size());
      for (size_t j=0; j<planes.size(); j++)
      {
        Vector n = planes[j].normal();
        planes_[j] = SLIVR::Plane(n.x(),n.y(),n.z(),planes[j].distance());
        plane_ptrs[j] = &(planes_[j]);
      }
      SLIVR::VolumeRenderer::set_planes(plane_ptrs);    
    }
    
  private:
    TextureHandle tex_handle_;
    ColorMapHandle colormap_handle_;
    std::vector<ColorMap2Handle> colormap2_handles_;
    std::vector<SLIVR::Plane> planes_;
};

typedef LockingHandle<VolumeRenderer> VolumeRendererHandle;


} // End namespace SCIRun

#endif // VolumeRenderer_h
