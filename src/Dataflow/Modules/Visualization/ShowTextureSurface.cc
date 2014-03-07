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
//    File   : ShowTextureSurface.cc
//    Author : Michael Callahan
//    Date   : April 2008

#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/ColorMapPort.h>
#include <Dataflow/Network/Ports/ColorMap2Port.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/TexturePort.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/Field.h>
#include <Core/Util/StringUtil.h>
#include <Core/Volume/SurfaceRenderer.h>
#include <Core/Geometry/Transform.h>
#include <Core/Geometry/SLIVRConvert.h>
#include <Core/Geom/View.h>
#include <slivr/VideoCardInfo.h>
#include <slivr/ShaderProgramARB.h>

namespace SCIRun {

using namespace SLIVR;

class ShowTextureSurface : public Module {
public:
  ShowTextureSurface(GuiContext*);
  virtual ~ShowTextureSurface() {}
  virtual void execute();

private:
  TextureHandle tex_;
  TextureHandle tex_handle_copy_;

  GeometryOPortHandle ogeom_;

  int cmap1_prevgen_;
  int cmap2_prevgen_;
  int field_prevgen_;
  int card_mem_;
  
  GuiInt interp_mode_;

  TextureHandle old_tex_;
  ColorMapHandle old_cmap1_;
  ColorMap2Handle old_cmap2_;
  Point old_min_, old_max_;
  GeomID geom_id_;
  SurfaceRenderer* renderer_;
  std::vector<SLIVR::ColorMap2*> cmap2_array_;
};


DECLARE_MAKER(ShowTextureSurface)

ShowTextureSurface::ShowTextureSurface(GuiContext* ctx) :
  Module("ShowTextureSurface", ctx, Source, "Visualization", "SCIRun"),
  tex_(0),
  cmap1_prevgen_(0),
  cmap2_prevgen_(0),
  field_prevgen_(0),
  card_mem_(video_card_memory_size()),
  interp_mode_(get_ctx()->subVar("interp_mode"), 1),
  old_tex_(0),
  old_cmap1_(0),
  old_cmap2_(0),
  old_min_(Point(0,0,0)),
  old_max_(Point(0,0,0)),
  geom_id_(-1),
  renderer_(0)
{
  get_oport_handle("Geometry",ogeom_);
}

void
ShowTextureSurface::execute()
{
  bool geom_changed = false;

  get_input_handle("Texture", tex_,true);
  
  ColorMapHandle cmap1;
  ColorMap2Handle cmap2;
  bool c1 = get_input_handle("ColorMap", cmap1, false);
  bool c2 = get_input_handle("ColorMap2", cmap2, false);

  FieldHandle ifieldhandle;
  if (!get_input_handle("Surface", ifieldhandle)) return;
  
  VMesh *trimesh = ifieldhandle->vmesh();
  if (!(trimesh->is_trisurfmesh()))
  {
    error("Surface did not contain a triangle mesh.");
    return;
  }

  // Inform module that execution started
  update_state(Executing);

  if (c2 && !ShaderProgramARB::shaders_supported())
  {
    warning("ColorMap2 usage is not supported by this machine.");
    cmap2 = 0;
    c2 = false;
  }

  reset_vars();

  cmap2_array_.clear(); 
  if (c2) { cmap2_array_.push_back(cmap2.get_rep()); }
    
  if (!c1 && !c2)
  {
    error("No colormap available to render.  Nothing to draw.");
    return;
  }

  if (!renderer_)
  {
    tex_handle_copy_ = tex_; //assure positive ref_cnt
    renderer_ = new SurfaceRenderer((SLIVR::Texture*)tex_.get_rep(), 
                                    (SLIVR::ColorMap*)cmap1.get_rep(), 
                                    cmap2_array_,
                                    int(card_mem_*1024*1024*0.8));
    geom_id_ = ogeom_->addObj(renderer_, "Surface Embedded in Volume");
    old_tex_ = tex_;
    old_cmap1_ = cmap1;
    old_cmap2_ = cmap2;
    SLIVR::BBox b;
    tex_->get_bounds(b);
    old_min_ = FROM_SLIVR_POINT(b.min());
    old_max_ = FROM_SLIVR_POINT(b.max());
    geom_changed = true;
  }
  else
  {
    SLIVR::BBox b;
    tex_->get_bounds(b);
    if (tex_.get_rep() != old_tex_.get_rep() ||
        FROM_SLIVR_POINT(b.min()) != old_min_ || 
        FROM_SLIVR_POINT(b.max()) != old_max_)
    {
      old_tex_ = tex_;
      old_min_ = FROM_SLIVR_POINT(b.min());
      old_max_ = FROM_SLIVR_POINT(b.max());
      if(geom_id_ != -1)
      {
        ogeom_->delObj(geom_id_);
        geom_id_ = ogeom_->addObj(renderer_, "Volume Slicer");
      }
      tex_handle_copy_ = tex_; // assure positive ref_cnt
      renderer_->lock(); // dont render while changing state.
      renderer_->set_texture(tex_.get_rep());
      renderer_->unlock();
      geom_changed = true;
    }

    renderer_->lock(); // dont render while changing state.
    if (cmap1 != old_cmap1_)
    {
      renderer_->set_colormap1(cmap1.get_rep());
      old_cmap1_ = cmap1;
      geom_changed = true;
    }
    if (cmap2 != old_cmap2_)
    {
      renderer_->set_colormap2(cmap2_array_);
      old_cmap2_ = cmap2;
      geom_changed = true;
    }
    renderer_->unlock();
  }

  renderer_->set_interp(bool(interp_mode_.get()));

  if (field_prevgen_ != ifieldhandle->generation)
  {
    field_prevgen_ = ifieldhandle->generation;

    // Push the tris to the volume renderer.
    renderer_->lock();
    renderer_->clear_triangles();
    VMesh::Elem::iterator abi, aei;
    VMesh::Node::array_type nodes;
    trimesh->begin(abi);
    trimesh->end(aei);
    std::vector<Point> pts(3);
    while (abi != aei)
    {
      trimesh->get_nodes(nodes, *abi);
      Point pts[3];
      for (int i = 0; i < 3; i++)
      {
        trimesh->get_point(pts[i], nodes[i]);
      }
      renderer_->add_triangle(TO_SLIVR_POINT(pts[0]), 
                              TO_SLIVR_POINT(pts[1]), 
                              TO_SLIVR_POINT(pts[2]));
      ++abi;
    }
    renderer_->unlock();
    geom_changed = true;
  }

  if (geom_changed) { ogeom_->flushViews(); }

  if (cmap1.get_rep())
  {
    ColorMapHandle outcmap;
    outcmap = new ColorMap(*cmap1.get_rep()); 
    double vmin = tex_->vmin();
    double vmax = tex_->vmax();
    outcmap->Scale(vmin, vmax);
    send_output_handle("ColorMap", outcmap);
  }
}


} // namespace SCIRun
