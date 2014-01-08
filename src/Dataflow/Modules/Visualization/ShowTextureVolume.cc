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
//    File   : ShowTextureVolume.cc
//    Author : Milan Ikits
//    Author : Kurt Zimmerman
//    Date   : Sat Jul 10 21:55:34 2004

#include <slivr/VideoCardInfo.h>
#include <slivr/ShaderProgramARB.h>
#include <slivr/Utils.h>
#include <Core/Geom/GeomSticky.h>
#include <Core/Volume/VolumeRenderer.h>
#include <Core/Volume/CM2Widget.h>

#include <Core/Geometry/SLIVRConvert.h>
#include <Core/Util/StringUtil.h>
#include <Core/Datatypes/ColorMap.h>
#include <Core/Thread/CrowdMonitor.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/ColorMapPort.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/GuiInterface/GuiVar.h>

#include <Dataflow/Widgets/PointWidget.h>

#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Widgets/PointWidget.h>
#include <Dataflow/Network/Ports/ColorMapPort.h>
#include <Dataflow/Network/Ports/ColorMap2Port.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/TexturePort.h>
#include <Dataflow/Widgets/ArrowWidget.h>

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>


namespace SCIRun {
using SLIVR::ShaderProgramARB;

class ShowTextureVolume : public Module {
  public:
    ShowTextureVolume(GuiContext*);
    virtual ~ShowTextureVolume();
    
    virtual void execute();
    virtual void widget_moved(bool last, BaseWidget*);
    virtual void post_read(); // get the widget state...
  private:

    void push_back_new_widget(const BBox &bbox);
    void remove_last_widget();

    TextureHandle texture_;
    TextureHandle texture_handle_copy_; // make sure the ref_cnt is maintained
                                        // until we can lock the volume renderer
                                        // draw.  get_input_handle sets to 0.

    GeometryOPortHandle geom_oport_;

    int cmap1_prevgen_;
    std::vector<int> cmap2_prevgen_;
    int tex_prevgen_;
    int card_mem_;

    CrowdMonitor              widget_lock_;
    std::vector<int>               widget_id_;
    std::vector<ArrowWidget*>      widget_;
    std::vector<GeomHandle>        widget_switch_;
    std::vector<Plane>             clipping_planes_;
    std::map<ArrowWidget *, int>   widget_index_map_;

    GuiDouble gui_sampling_rate_hi_;
    GuiDouble gui_sampling_rate_lo_;
    GuiDouble gui_gradient_min_;
    GuiDouble gui_gradient_max_;
    GuiInt gui_adaptive_;
    GuiInt gui_cmap_size_;
    GuiInt gui_sw_raster_;
    GuiInt gui_render_style_;
    GuiDouble gui_alpha_scale_;
    GuiInt gui_interp_mode_;
    GuiInt gui_shading_;
    GuiDouble gui_ambient_;
    GuiDouble gui_diffuse_;
    GuiDouble gui_specular_;
    GuiDouble gui_shine_;
    GuiInt gui_light_;
    GuiInt gui_blend_res_;
    GuiInt gui_multi_level_;
    GuiInt gui_use_stencil_;
    GuiInt gui_invert_opacity_;
    GuiInt gui_level_flag_;
    GuiInt gui_num_slices_;  // unused except for backwards compatability
    GuiInt gui_num_clipping_planes_;
    std::vector<GuiDouble*> gui_widget_state_;
    GuiInt gui_show_clipping_widgets_;
    GuiString gui_level_on_; // used for saving to net
    GuiString gui_level_vals_; // used for saving to net

    VolumeRendererHandle volren_;
    
    
    Point  oldmin_;
    Point  oldmax_;
    int    oldni_;
    int    oldnj_;
    int    oldnk_;
    GeomID geomid_;  
};


DECLARE_MAKER(ShowTextureVolume)
ShowTextureVolume::ShowTextureVolume(GuiContext* ctx) :
  Module("ShowTextureVolume", ctx, Source, "Visualization", "SCIRun"),
  texture_(0),
  cmap1_prevgen_(0),
  cmap2_prevgen_(0),
  tex_prevgen_(0),
  card_mem_(video_card_memory_size()),
  widget_lock_("Clipping planes widget lock"),
  widget_id_(),
  widget_(),
  widget_switch_(),
  clipping_planes_(),
  widget_index_map_(),
  gui_sampling_rate_hi_(get_ctx()->subVar("sampling_rate_hi"), 4.0),
  gui_sampling_rate_lo_(get_ctx()->subVar("sampling_rate_lo"), 1.0),
  gui_gradient_min_(get_ctx()->subVar("gradient_min"), 0.0),
  gui_gradient_max_(get_ctx()->subVar("gradient_max"), 0.0),
  gui_adaptive_(get_ctx()->subVar("adaptive"), 1),
  gui_cmap_size_(get_ctx()->subVar("cmap_size"), 8),
  gui_sw_raster_(get_ctx()->subVar("sw_raster"), 0),
  gui_render_style_(get_ctx()->subVar("render_style"), 0),
  gui_alpha_scale_(get_ctx()->subVar("alpha_scale"), 0),
  gui_interp_mode_(get_ctx()->subVar("interp_mode"), 1),
  gui_shading_(get_ctx()->subVar("shading"), 0),
  gui_ambient_(get_ctx()->subVar("ambient"), 0.5),
  gui_diffuse_(get_ctx()->subVar("diffuse"), 0.5),
  gui_specular_(get_ctx()->subVar("specular"), 0.0),
  gui_shine_(get_ctx()->subVar("shine"), 30.0),
  gui_light_(get_ctx()->subVar("light"), 0),
  gui_blend_res_(get_ctx()->subVar("blend_res"), 8),
  gui_multi_level_(get_ctx()->subVar("multi_level"), 1),
  gui_use_stencil_(get_ctx()->subVar("use_stencil"), 0),
  gui_invert_opacity_(get_ctx()->subVar("invert_opacity"), 0),
  gui_level_flag_(get_ctx()->subVar("show_level_flag", false), 1),
  gui_num_slices_(get_ctx()->subVar("num_slices", false), -1), // dont save
  gui_num_clipping_planes_(get_ctx()->subVar("num_clipping_planes"), 0),
  gui_show_clipping_widgets_(get_ctx()->subVar("show_clipping_widgets"), 1),
  gui_level_on_(get_ctx()->subVar("level_on"), ""),
  gui_level_vals_(get_ctx()->subVar("level_vals"), "")
{
  get_oport_handle("Geometry",geom_oport_);
}

ShowTextureVolume::~ShowTextureVolume()
{
  for (size_t j=0;j<widget_.size(); j++)
    if (widget_[j]) delete widget_[j];
  widget_.clear();
}

void
ShowTextureVolume::widget_moved(bool release, BaseWidget *widget)
{
  ArrowWidget *arrow = dynamic_cast<ArrowWidget*>(widget);
  if (!arrow) return;
  
  std::map<ArrowWidget *, int>::iterator pos = widget_index_map_.find(arrow);
  if (pos == widget_index_map_.end()) return;
  
  const int idx = pos->second;
  Point up = FROM_SLIVR_TRANSFORM(texture_->transform()).unproject(arrow->GetPosition());
  clipping_planes_[idx] = Plane(up, -arrow->GetDirection());
  //cache the widget state

  if (release) 
  {
    int gw_idx = idx * 7;
    //scale, direction, and point get cached per widget. 7 doubles
    gui_widget_state_[gw_idx  ]->set(arrow->GetScale());
    Vector dir = arrow->GetDirection();
    gui_widget_state_[gw_idx+1]->set(dir.x());
    gui_widget_state_[gw_idx+2]->set(dir.y());
    gui_widget_state_[gw_idx+3]->set(dir.z());
    Point wpos = arrow->GetPosition();
    gui_widget_state_[gw_idx+4]->set(wpos.x());
    gui_widget_state_[gw_idx+5]->set(wpos.y());
    gui_widget_state_[gw_idx+6]->set(wpos.z());
  }

  //update volume renderers planes.
  volren_->lock.lock();
  volren_->set_planes(clipping_planes_);
  volren_->lock.unlock();
}


void 
ShowTextureVolume::post_read()
{
  // restore the widgets.
  gui_num_clipping_planes_.reset();
  // Deal with bad values from old .srns..
  if (gui_num_clipping_planes_.get() < 0) 
  {
    gui_num_clipping_planes_.set(0);
    gui_num_clipping_planes_.reset();
  }
  
  std::string ps("widget-state-");
  for (int i = 0; i < gui_num_clipping_planes_.get(); ++i) 
  {
    std::string istr = to_string(i);
    std::string id = ps + istr;
    int b = i*7;
    gui_widget_state_.insert(gui_widget_state_.end(), 7, (GuiDouble*)0);
    //scale, direction, and point get cached per widget. 7 doubles
    gui_widget_state_[b  ] = new GuiDouble(get_ctx()->subVar(id + "-scale"));
    gui_widget_state_[b+1] = new GuiDouble(get_ctx()->subVar(id + "-dir-x"));
    gui_widget_state_[b+2] = new GuiDouble(get_ctx()->subVar(id + "-dir-y"));
    gui_widget_state_[b+3] = new GuiDouble(get_ctx()->subVar(id + "-dir-z"));
    gui_widget_state_[b+4] = new GuiDouble(get_ctx()->subVar(id + "-pnt-x"));
    gui_widget_state_[b+5] = new GuiDouble(get_ctx()->subVar(id + "-pnt-y"));
    gui_widget_state_[b+6] = new GuiDouble(get_ctx()->subVar(id + "-pnt-z"));

    for (int off = 0; off < 7; ++off) 
    {
      gui_widget_state_[b + off]->reset();
    }
  }  
}

// Return the plane the new widget defines.
void
ShowTextureVolume::push_back_new_widget(const BBox &bbox)
{
  Vector dir = bbox.diagonal();
  double scale = dir.length() / 60.0;
  Point p = bbox.min() + dir / 2.0;
  
  int ib = widget_.size() * 7;
  
  if (gui_num_clipping_planes_.get() > (int)gui_widget_state_.size() / 7)
  {
    gui_widget_state_.insert(gui_widget_state_.end(), 7, (GuiDouble*)0);
    std::string istr = to_string(widget_.size());
    std::string ps("widget-state-");
    std::string id = ps + istr;
    
    gui_widget_state_[ib  ] = new GuiDouble(get_ctx()->subVar(id + "-scale"), 
					    scale);
    gui_widget_state_[ib+1] = new GuiDouble(get_ctx()->subVar(id + "-dir-x"),
					    dir.x());
    gui_widget_state_[ib+2] = new GuiDouble(get_ctx()->subVar(id + "-dir-y"),
					    dir.y());
    gui_widget_state_[ib+3] = new GuiDouble(get_ctx()->subVar(id + "-dir-z"),
					    dir.z());
    gui_widget_state_[ib+4] = new GuiDouble(get_ctx()->subVar(id + "-pnt-x"),
					    p.x());
    gui_widget_state_[ib+5] = new GuiDouble(get_ctx()->subVar(id + "-pnt-y"),
					    p.y());
    gui_widget_state_[ib+6] = new GuiDouble(get_ctx()->subVar(id + "-pnt-z"),
					    p.z());
  }

  double gscale = gui_widget_state_[ib]->get();
  Vector gdir(gui_widget_state_[ib+1]->get(), gui_widget_state_[ib+2]->get(), 
              gui_widget_state_[ib+3]->get());
  Point gp(gui_widget_state_[ib+4]->get(), gui_widget_state_[ib+5]->get(), 
           gui_widget_state_[ib+6]->get());

  ArrowWidget *widget = new ArrowWidget(this, &widget_lock_, gscale, true);
  widget_.push_back(widget);
  widget->Connect(geom_oport_.get_rep());
  widget->SetCurrentMode(0);
  GeomSwitch *swtch = (GeomSwitch *)widget->GetWidget().get_rep();
  swtch->set_state(1);

  widget_switch_.push_back(swtch);
  ((GeomSwitch *)(widget_switch_.back().get_rep()))->set_state(1);
  widget_id_.push_back(geom_oport_->addObj(widget_switch_.back(),
					   "Clipping plane " +
					   to_string(widget_.size()),
					   &widget_lock_));

  widget->SetDirection(gdir);
  widget->SetScale(gscale);
  widget->SetLength(gscale*2);
  widget->Move(gp);
  widget->redraw();
  int idx = widget_.size() - 1;
  widget_index_map_[widget] = idx;
  
  // add the plane.
  Point pnt = FROM_SLIVR_TRANSFORM(texture_->transform()).unproject(widget->GetPosition());
  clipping_planes_.push_back(Plane(pnt,-widget->GetDirection()));

}


void
ShowTextureVolume::remove_last_widget()
{
  if (widget_id_.size()) 
  {
    int id = widget_id_.back();
    geom_oport_->delObj(id);
    widget_id_.pop_back();
  }

  if (widget_.size()) 
  {
    widget_.pop_back();
  }

  // remap the widget indeces.
  widget_index_map_.clear();
  for (size_t i = 0; i < widget_.size(); ++i) 
  {
    widget_index_map_[widget_[i]] = i;
  }
  
  if (gui_widget_state_.size()) 
  {
    int i = gui_widget_state_.size() - 7;
    delete gui_widget_state_[i  ];
    delete gui_widget_state_[i+1];
    delete gui_widget_state_[i+2];
    delete gui_widget_state_[i+3];
    delete gui_widget_state_[i+4];
    delete gui_widget_state_[i+5];
    delete gui_widget_state_[i+6];

    std::vector<GuiDouble*>::iterator iter = gui_widget_state_.end() - 7;
    gui_widget_state_.erase(iter, gui_widget_state_.end());
  }
  
  if (clipping_planes_.size()) 
  {
    clipping_planes_.pop_back();
  }
}

void
ShowTextureVolume::execute()
{
  get_input_handle("Texture", texture_);

  bool shading_state = false;
  if (ShaderProgramARB::shaders_supported())
  {
    shading_state = (texture_->nb(0) == 1);
  }

  TCLInterface::execute(get_id() + " change_shading_state " + (shading_state?"0":"1"));
  
  ColorMapHandle cmap1;
  
  const bool c1 = get_input_handle("ColorMap", cmap1, false);
  std::vector<int> cmap2_generation;

  std::vector<ColorMap2Handle> cmap2;
  get_dynamic_input_handles("ColorMap2",cmap2,false);
  
  // Inform module that execution started
  update_state(Executing);
  
  cmap2_generation.clear();
  for (size_t j=0; j<cmap2.size(); j++)
  {
    cmap2_generation.push_back(cmap2[j]->generation);
  }

    
  bool c2 = !cmap2.empty();

  if (c2)
  {
    if (!ShaderProgramARB::shaders_supported())
    {
      warning("ColorMap2 usage is not supported by this machine.");
      c2 = false;
    }
  }

  if (!c1 && !c2)
  {
    error("No colormap available to render.  Nothing drawn.");
    return;
  }

  bool cmap1_dirty = false;
  if(c1 && (cmap1->generation != cmap1_prevgen_)) 
  {
    cmap1_prevgen_ = cmap1->generation;
    cmap1_dirty = true;
  }

  bool cmap2_dirty = cmap2_generation.size() != cmap2_prevgen_.size();
  if(c2 && !cmap2_dirty) 
  {
    for (size_t g = 0; g < cmap2_generation.size(); ++g) 
    {
      cmap2_dirty = cmap2_generation[g] != cmap2_prevgen_[g];
      if (cmap2_dirty) break;
    }
  }
  cmap2_prevgen_ = cmap2_generation;

  bool tex_dirty = false;
  if (texture_.get_rep() && texture_->generation != tex_prevgen_) 
  {
    tex_prevgen_ = texture_->generation;
    tex_dirty = true;
  }
   
  if (!cmap1_dirty && !cmap2_dirty && !tex_dirty && 
      !gui_sampling_rate_hi_.changed() && !gui_sampling_rate_lo_.changed() &&
      !gui_adaptive_.changed() && !gui_cmap_size_.changed() && 
      !gui_sw_raster_.changed() && !gui_render_style_.changed() &&
      !gui_alpha_scale_.changed() && !gui_interp_mode_.changed() &&
      !gui_shading_.changed() && !gui_ambient_.changed() &&
      !gui_diffuse_.changed() && !gui_specular_.changed() &&
      !gui_shine_.changed() && !gui_light_.changed() &&
      !gui_blend_res_.changed() && !gui_multi_level_.changed() &&
      !gui_use_stencil_.changed() && !gui_invert_opacity_.changed() &&
      !gui_level_flag_.changed() && (gui_multi_level_.get() == 0) )
  {
    if (texture_.get_rep())
    {
      for (size_t i = 0; i < texture_->bricks().size(); i++)
      {
        if (texture_->bricks()[i]->dirty())
        {
          geom_oport_->flushViews();
          break;
        }
      }
    }
    return;
  }


  std::string s;
  TCLInterface::eval(get_id() + " hasUI", s);
  if( s == "0" )
    TCLInterface::execute(get_id() + " buildTopLevel");

  if( texture_->nlevels() > 1 && gui_multi_level_.get() == 1)
  {
    gui_multi_level_.set(texture_->nlevels());
    TCLInterface::execute(get_id() + " build_multi_level");
  } 
  else if(texture_->nlevels() == 1 && gui_multi_level_.get() > 1)
  {
    gui_multi_level_.set(1);
    TCLInterface::execute(get_id() + " destroy_multi_level");
  }

  // cache gui vars for use while we have the volren_ lock.
  double gradient_min = gui_gradient_min_.get();
  double gradient_max = gui_gradient_max_.get();
  double sampling_rate_hi = gui_sampling_rate_hi_.get();
  double sampling_rate_lo = gui_sampling_rate_lo_.get();
  int colormap_size = 1 << gui_cmap_size_.get();
  int adaptive = gui_adaptive_.get();
  double alpha_scale = gui_alpha_scale_.get();
  bool use_stencil = gui_use_stencil_.get();
  bool invert_opacity = gui_invert_opacity_.get();
  int nslices = gui_num_slices_.get();
  bool shading = static_cast<bool>(gui_shading_.get());
  double ambient = gui_ambient_.get();
  double diffuse = gui_diffuse_.get();
  double specular = gui_specular_.get();
  double shine = gui_shine_.get();
  int light = gui_light_.get();
  int rstyle = gui_render_style_.get();
  bool interp_mode = gui_interp_mode_.get();
  int blend_num_bits = gui_blend_res_.get();


  if(!(volren_.get_rep())) 
  {
    texture_handle_copy_ = texture_; // maintain the ref_cnt.
    volren_ = new VolumeRenderer(texture_, cmap1, 
				 cmap2, clipping_planes_, 
				 int(card_mem_*1024*1024*0.8));
         
    oldmin_ = FROM_SLIVR_POINT(texture_->bbox().min());
    oldmax_ = FROM_SLIVR_POINT(texture_->bbox().max());
    oldni_ = texture_->nx();
    oldnj_ = texture_->ny();
    oldnk_ = texture_->nz();
    // geom_oport_->delAll();
    geomid_ = geom_oport_->addObj(volren_.get_rep(), "VolumeRenderer Transparent");
    volren_->lock.lock(); //stop rendering while state changes.
  } 
  else 
  {
    volren_->lock.lock(); //stop rendering while state changes.
                     // don't trigger a gui lock while locked.

    texture_handle_copy_ = texture_; // maintain the ref_cnt.
    volren_->set_texture(texture_);
    if(c1 && cmap1_dirty) 
    {
      volren_->set_colormap1(cmap1);
    } 
    else if (!c1) 
    {
      volren_->set_colormap1(0);
    }
    if(c2 && cmap2_dirty) 
    {
      volren_->set_colormap2(cmap2);
    } 
    else if (!c2) 
    {
      volren_->set_colormap2(cmap2);
    }

    volren_->lock.unlock(); 

    if(oldmin_ != FROM_SLIVR_POINT(texture_->bbox().min()) || 
       oldmax_ !=  FROM_SLIVR_POINT(texture_->bbox().max()) ||
       oldni_  != texture_->nx() || 
       oldnj_  != texture_->ny() || 
       oldnk_  != texture_->nz()) 
    {
      geom_oport_->delObj(geomid_);
      geomid_ = geom_oport_->addObj(volren_.get_rep(), "VolumeRenderer Transparent");
      oldni_  = texture_->nx(); 
      oldnj_  = texture_->ny(); 
      oldnk_  = texture_->nz();
      oldmin_ = FROM_SLIVR_POINT(texture_->bbox().min());
      oldmax_ = FROM_SLIVR_POINT(texture_->bbox().max());
    }
    
    volren_->lock.lock(); 

  }
  
  volren_->set_interp(interp_mode);
  switch(rstyle) 
  {
  case 0:
    volren_->set_mode(VolumeRenderer::MODE_OVER);
    break;
  case 1:
    volren_->set_mode(VolumeRenderer::MODE_MIP);
    break;
  default:
    warning("Unsupported blend mode.  Using over.");
    volren_->set_mode(VolumeRenderer::MODE_OVER);
    break;
  }
  double rate = 0.0;
  if (nslices > 0) 
  {
    rate = volren_->num_slices_to_rate(nslices);
  }
  
  volren_->set_gradient_range(gradient_min,gradient_max);
  volren_->set_sampling_rate(sampling_rate_hi);
  volren_->set_interactive_rate(sampling_rate_lo);
  volren_->set_adaptive(adaptive);
  volren_->set_colormap2_width(colormap_size);
  volren_->set_slice_alpha(alpha_scale);
  volren_->set_stencil(use_stencil );
  volren_->invert_opacity(invert_opacity);
  volren_->set_blend_num_bits(blend_num_bits);
  volren_->set_shading(shading);
  volren_->set_material(ambient, diffuse, specular, shine);
  volren_->set_light(light);

  volren_->lock.unlock(); //done changing state, allow rendering.

  if (nslices > 0) 
  {
    gui_sampling_rate_hi_.set(rate);
    gui_num_slices_.set(-1);
  }

  // Always one fewer widget than colormap2Ds.
  int ncmps = cmap2.size();
  gui_num_clipping_planes_.set(ncmps ? ncmps - 1 : 0);

  if (gui_num_clipping_planes_.get() != static_cast<int>(widget_.size())) 
  {
    SLIVR::BBox bbox;
    texture_->get_bounds(bbox);
    
    while (gui_num_clipping_planes_.get() != static_cast<int>(widget_.size()))
    {
      if (static_cast<int>(widget_.size()) < gui_num_clipping_planes_.get()) 
      {
        // add a widget.
        push_back_new_widget(FROM_SLIVR_BBOX(bbox));
      } 
      else 
      {
        // remove a widget.
        remove_last_widget();
      }
    }
  }

  //Sync the planes in the Volume Renderer.
  volren_->lock.lock();
  volren_->set_planes(clipping_planes_);
  volren_->set_colormap2_dirty();
  volren_->lock.unlock();

  // Sync all of the widget switches with the checkbox.
  for (unsigned int w = 0; w < widget_switch_.size(); ++w) 
    ((GeomSwitch*)widget_switch_[w].get_rep())->
      set_state(gui_show_clipping_widgets_.get());

  geom_oport_->flushViews();				  

  if (c1)
  {
    ColorMapHandle outcmap;
    outcmap = new ColorMap(*cmap1.get_rep()); 
    outcmap->Scale(texture_->vmin(), texture_->vmax());
    send_output_handle("ColorMap", outcmap);
  }
}

} // End namespace SCIRun
