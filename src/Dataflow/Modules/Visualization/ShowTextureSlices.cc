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
//    File   : ShowTextureSlices.cc
//    Author : Milan Ikits
//    Author : Kurt Zimmerman
//    Date   : Sat Jul 10 21:55:08 2004


#include <Core/Util/StringUtil.h>
#include <Core/Datatypes/ColorMap.h>
#include <Core/Geom/View.h>
#include <Core/Thread/CrowdMonitor.h>

#include <Core/Geometry/Transform.h>
#include <Core/Geometry/SLIVRConvert.h>

#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/ColorMapPort.h>
#include <Dataflow/Network/Ports/ColorMap2Port.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/TexturePort.h>
#include <Dataflow/Widgets/PointWidget.h>


#include <iostream>
#include <algorithm>

#include <slivr/VideoCardInfo.h>
#include <Core/Volume/SliceRenderer.h>
#include <slivr/ShaderProgramARB.h>

namespace SCIRun {

using namespace SLIVR;

class ShowTextureSlices : public Module {
public:
  ShowTextureSlices(GuiContext*);
  virtual ~ShowTextureSlices();
  virtual void execute();
  virtual void widget_moved(bool last, BaseWidget*);
  virtual void tcl_command(GuiArgs&, void*);

private:
  TextureHandle tex_;
  TextureHandle tex_handle_copy_;

  GeometryOPortHandle ogeom_;

  int cmap1_prevgen_;
  int cmap2_prevgen_;
  int card_mem_;
  
  CrowdMonitor control_lock_; 
  PointWidget *control_widget_;
  GeomID control_id_;

  GuiInt control_pos_saved_;
  GuiDouble control_x_;
  GuiDouble control_y_;
  GuiDouble control_z_;
  GuiInt draw_x_;
  GuiInt draw_y_;
  GuiInt draw_z_;
  GuiInt draw_view_;
  GuiInt interp_mode_;
  GuiInt draw_phi0_;
  GuiInt draw_phi1_;
  GuiDouble phi0_;
  GuiDouble phi1_;
  GuiInt cyl_active_;
  GuiInt gui_multi_level_;
  GuiInt gui_color_changed_;
  GuiString gui_colors_;
  GuiString gui_level_on_;
  GuiInt gui_outline_levels_;
  GuiInt gui_use_stencil_;

  TextureHandle old_tex_;
  ColorMapHandle old_cmap1_;
  ColorMap2Handle old_cmap2_;
  Point old_min_, old_max_;
  GeomID geom_id_;
  SliceRenderer* slice_ren_;
  Point dmin_;
  Vector ddx_;
  Vector ddy_;
  Vector ddz_;
  double ddview_;

  int color_changed_;
};


static std::string control_name("Control Widget");
			 
DECLARE_MAKER(ShowTextureSlices)

ShowTextureSlices::ShowTextureSlices(GuiContext* ctx) :
  Module("ShowTextureSlices", ctx, Source, "Visualization", "SCIRun"),
  tex_(0),
  cmap1_prevgen_(0),
  cmap2_prevgen_(0),
  card_mem_(video_card_memory_size()),
  control_lock_("ShowTextureSlices resolution lock"),
  control_widget_(0),
  control_id_(-1),
  control_pos_saved_(get_ctx()->subVar("control_pos_saved"), 0),
  control_x_(get_ctx()->subVar("control_x")),
  control_y_(get_ctx()->subVar("control_y")),
  control_z_(get_ctx()->subVar("control_z")),
  draw_x_(get_ctx()->subVar("drawX"), 0),
  draw_y_(get_ctx()->subVar("drawY"), 0),
  draw_z_(get_ctx()->subVar("drawZ"), 0),
  draw_view_(get_ctx()->subVar("drawView"), 0),
  interp_mode_(get_ctx()->subVar("interp_mode"), 1),
  draw_phi0_(get_ctx()->subVar("draw_phi_0"), 0),
  draw_phi1_(get_ctx()->subVar("draw_phi_1"), 0),
  phi0_(get_ctx()->subVar("phi_0"), 30.0),
  phi1_(get_ctx()->subVar("phi_1"), 60.0),
  cyl_active_(get_ctx()->subVar("cyl_active")), 
  gui_multi_level_(get_ctx()->subVar("multi_level"), 1),
  gui_color_changed_(get_ctx()->subVar("color_changed"), 1),
  gui_colors_(get_ctx()->subVar("colors"), ""),
  gui_level_on_(get_ctx()->subVar("level_on"), ""),
  gui_outline_levels_(get_ctx()->subVar("outline_levels"), 0),
  gui_use_stencil_(get_ctx()->subVar("use_stencil"), 0),
  old_tex_(0),
  old_cmap1_(0),
  old_cmap2_(0),
  old_min_(Point(0,0,0)), old_max_(Point(0,0,0)),
  geom_id_(-1),
  slice_ren_(0),
  color_changed_(true)
{
  get_oport_handle("Geometry",ogeom_);
}

ShowTextureSlices::~ShowTextureSlices()
{
  if (control_widget_) delete control_widget_;
}

void
ShowTextureSlices::execute()
{
  ColorMapHandle cmap1;
  ColorMap2Handle cmap2;

  get_input_handle("Texture", tex_,true);
  bool c1 = get_input_handle("ColorMap", cmap1, false);
  bool c2 = get_input_handle("ColorMap2", cmap2, false);

  // Inform module that execution started
  update_state(Executing);

  if (c2)
  {
    if (!ShaderProgramARB::shaders_supported())
    {
      warning("ColorMap2 usage is not supported by this machine.");
      cmap2 = 0;
      c2 = false;
    }
    else
    {
      if (tex_->nc() == 1)
      {
        warning("ColorMap2 requires gradient magnitude in the texture.");
        cmap2 = 0;
        c2 = false;
      }
    }
  }
  reset_vars();

  std::vector<ColorMap2Handle> cmap2_array(0);
  if (c2) 
  {
    cmap2_array.push_back(cmap2);
  }
    

  if (!c1 && !c2)
  {
    error("No colormap available to render.  Nothing drawn.");
    return;
  }

  if(!control_widget_)
  {
    control_widget_= new PointWidget(this, &control_lock_, 0.2);
    control_widget_->Connect(ogeom_.get_rep());
    
    double minx,miny,minz;
    double maxx,maxy,maxz;
    tex_->get_bounds(minx, miny, minz, maxx, maxy, maxz);
    BBox b;
    b.extend(Point(minx, miny, minz));
    b.extend(Point(maxx, maxy, maxz));

    Vector dv(b.diagonal());
    if(control_pos_saved_.get()) 
    {
      control_widget_->SetPosition(Point(control_x_.get(),
					 control_y_.get(),
					 control_z_.get()));
      control_widget_->SetScale(dv.length()/80.0);
    } 
    else 
    {
      control_widget_->SetPosition(Interpolate(b.min().asVector(), 
					       b.max().asVector(), 
					       0.5).asPoint());
      control_widget_->SetScale(dv.length()/80.0);
    }
    int nx = tex_->nx();
    int ny = tex_->ny();
    int nz = tex_->nz();
    Transform t(FROM_SLIVR_TRANSFORM(tex_->transform()));
    dmin_=t.project(Point(0,0,0));
    ddx_ = t.project(Vector(1.0/(nx-1), 0, 0));
    ddy_ = t.project(Vector(0, 1.0/(ny-1), 0));
    ddz_ = t.project(Vector(0, 0, 1.0/(nz-1)));
    ddview_ = (dv.length()/(std::max(nx, std::max(ny,nz)) -1));
  }

  std::string s;
  TCLInterface::eval(get_id() + " hasUI", s);
  if( s == "0" )
    TCLInterface::execute(get_id() + " buildTopLevel");

  if( tex_->nlevels() > 1 && gui_multi_level_.get() == 1)
  {
    gui_multi_level_.set(tex_->nlevels());
    TCLInterface::execute(get_id() + " build_multi_level");
  } 
  else if(tex_->nlevels() == 1 && gui_multi_level_.get() > 1)
  {
    gui_multi_level_.set(1);
    TCLInterface::execute(get_id() + " destroy_multi_level");
  }

  if( !slice_ren_ && gui_multi_level_.get() >= 1 ) 
  { 
    gui_color_changed_.set(1);
  }
  
  int cyl_active = cyl_active_.get();
  int dphi0 = draw_phi0_.get();
  double phi0 = phi0_.get();
  int dphi1 = draw_phi1_.get();
  double phi1 = phi1_.get();

  if(!slice_ren_) 
  {
    tex_handle_copy_ = tex_; //assure positive ref_cnt
    slice_ren_ = new SliceRenderer(tex_, cmap1, cmap2_array,
 				   int(card_mem_*1024*1024*0.8));
    
    slice_ren_->
      set_control_point(tex_->transform().unproject(TO_SLIVR_POINT(control_widget_->
						    ReferencePoint())));
    //    ogeom->delAll();
    geom_id_ = ogeom_->addObj(slice_ren_, "Volume Slicer");
    slice_ren_->set_cylindrical(cyl_active, dphi0, phi0, dphi1, phi1);
    old_tex_ = tex_;
    old_cmap1_ = cmap1;
    old_cmap2_ = cmap2;
    SLIVR::BBox b;
    tex_->get_bounds(b);
    old_min_ = FROM_SLIVR_POINT(b.min());
    old_max_ = FROM_SLIVR_POINT(b.max());
  } 
  else 
  {
    SLIVR::BBox b;
    tex_->get_bounds(b);
    if(tex_.get_rep() != old_tex_.get_rep() ||
       FROM_SLIVR_POINT(b.min()) != old_min_ || 
       FROM_SLIVR_POINT(b.max()) != old_max_) 
    {
      old_tex_ = tex_;
      old_min_ = FROM_SLIVR_POINT(b.min());
      old_max_ = FROM_SLIVR_POINT(b.max());
      if(geom_id_ != -1) 
      {
        ogeom_->delObj(geom_id_);
        geom_id_ = ogeom_->addObj(slice_ren_, "Volume Slicer");
      }
      Vector dv(FROM_SLIVR_VECTOR(b.diagonal()));
      int nx = tex_->nx();
      int ny = tex_->ny();
      int nz = tex_->nz();
      Transform t = FROM_SLIVR_TRANSFORM(tex_->transform());
      dmin_=t.project(Point(0,0,0));
      ddx_ = t.project(Vector(1.0/(nx-1), 0, 0));
      ddy_ = t.project(Vector(0, 1.0/(ny-1), 0));
      ddz_ = t.project(Vector(0, 0, 1.0/(nz-1)));
      ddview_ = (dv.length()/(std::max(nx, std::max(ny,nz)) -1));
      if(!b.inside(TO_SLIVR_POINT(control_widget_->GetPosition()))) 
      {
        control_widget_->SetPosition(FROM_SLIVR_POINT(Interpolate(b.min().asVector(), 
						 b.max().asVector(), 
						 0.5).asPoint()));
      }
      control_widget_->SetScale(dv.length() / 80.0);
      tex_handle_copy_ = tex_; //assure positive ref_cnt

      slice_ren_->lock.lock(); // dont render while changing state.
      slice_ren_->set_texture(tex_);
      slice_ren_->set_control_point(tex_->transform().unproject(
					TO_SLIVR_POINT(control_widget_->ReferencePoint())));
      slice_ren_->lock.unlock();
    }

    slice_ren_->lock.lock(); // dont render while changing state.
    if (cmap1 != old_cmap1_)
    {
      slice_ren_->set_colormap1(cmap1);
      old_cmap1_ = cmap1;
    }
    if (cmap2 != old_cmap2_)
    {
      slice_ren_->set_colormap2(cmap2_array);
      old_cmap2_ = cmap2;
    }
    slice_ren_->lock.unlock();
  }

  if(gui_multi_level_.get() > 1 && gui_color_changed_.get() == 1){
    gui_color_changed_.set(0);
    std::string outline_colors;
    TCLInterface::eval(get_id()+" getOutlineColors", outline_colors);

    std::istringstream is( outline_colors );
    // Slurp in the rgb values.
    unsigned int rgbsize;
    is >> rgbsize;
    std::vector< Color > rgbs(rgbsize);
    for (unsigned int i = 0; i < rgbsize; i++)
    {
      double r, g, b;
      is >> r >> g >> b;
      rgbs[i] = Color(r, g, b);
    }
    slice_ren_->set_outline_colors( rgbs );
  }

  slice_ren_->set_interp(bool(interp_mode_.get()));
  if(draw_x_.get() || draw_y_.get() || draw_z_.get()){
    if(control_id_ == -1) {
      GeomHandle w = control_widget_->GetWidget();
      control_id_ = ogeom_->addObj( w, control_name, &control_lock_);
    }
  } else {
    if(control_id_ != -1) {
      ogeom_->delObj(control_id_, 0);
      control_id_ = -1;
    }
  }
  // Separate gui locking from draw locking to avoid deadlocks.
  int dx = draw_x_.get();
  int dy = draw_y_.get();
  int dz = draw_z_.get();
  int dv = draw_view_.get();
  bool us = gui_use_stencil_.get();
  bool ol = gui_outline_levels_.get();

  //\todo {This should be stored in GuiVars, we shouldn't ever use eval}
  std::vector<bool> onvals;
  if( tex_->nlevels() > 1 )
  {
    onvals.resize(tex_->nlevels());
    for(int i = 0; i < tex_->nlevels(); i++)
    {
      std::string result;
      TCLInterface::eval(get_id() + " isOn l" + to_string(i), result);
      if ( result == "0")
      {
        onvals[i] = false;
      } 
      else 
      {
        onvals[i] = true;
      }
    }
  }

  slice_ren_->lock.lock();
  slice_ren_->set_x(dx);
  slice_ren_->set_y(dy);
  slice_ren_->set_z(dz);
  slice_ren_->set_view(dv);
  slice_ren_->set_stencil(us);
  slice_ren_->set_level_outline(ol);
  slice_ren_->set_cylindrical(cyl_active, dphi0, phi0, dphi1, phi1);
  
  if( tex_->nlevels() > 1 )
  {
    for(int i = 0; i < tex_->nlevels(); i++)
    {
      slice_ren_->set_draw_level(tex_->nlevels()-1 -i, onvals[i]);
    }
  }

  slice_ren_->lock.unlock(); // done changing state, allow rendering.
  ogeom_->flushViews();		  
  
  if(cmap1.get_rep())
  {
    ColorMapHandle outcmap;
    outcmap = new ColorMap(*cmap1.get_rep()); 
    double vmin = tex_->vmin();
    double vmax = tex_->vmax();
    outcmap->Scale(vmin, vmax);
    send_output_handle("ColorMap", outcmap);
  }
}

void
ShowTextureSlices::tcl_command(GuiArgs& args, void* userdata)
{
  if (args[1] == "MoveWidget") {
    if (!control_widget_) return;
    Point w(control_widget_->ReferencePoint());
    if (args[2] == "xplus") 
    {
      double val; from_string(args[3],val);
      w+=ddx_*val;
    } 
    else if (args[2] == "xat") 
    {
      double val; from_string(args[3],val);
      w=dmin_+ddx_*val;
    } 
    else if (args[2] == "yplus") 
    {
      double val; from_string(args[3],val);
      w+=ddy_*val;
    } 
    else if (args[2] == "yat") 
    {
      double val; from_string(args[3],val);
      w=dmin_+ddy_*val;
    } 
    else if (args[2] == "zplus") 
    {
      double val; from_string(args[3],val);
      w+=ddz_*val;
    } 
    else if (args[2] == "zat") 
    {
      double val; from_string(args[3],val);
      w=dmin_+ddz_*val;
    } 
    else if (args[2] == "vplus")
    {
      GeometryData* data = ogeom_->getData(0, 0, 1);
      Vector view = data->view->lookat() - data->view->eyep();
      view.normalize();
      double val; from_string(args[3],val);
      w += view*ddview_*val;
    }
    control_widget_->SetPosition(w);
    widget_moved(true, 0);
    control_x_.set(w.x());
    control_y_.set(w.y());
    control_z_.set(w.z());
    control_pos_saved_.set( 1 );
    ogeom_->flushViews();
  } 
  else if (args[1] == "color_changed") 
  {
    color_changed_ = true;
  } 
  else 
  {
    Module::tcl_command(args, userdata);
  }
}

void
ShowTextureSlices::widget_moved(bool,BaseWidget*)
{
  if(slice_ren_) 
  {
    SLIVR::Point p = tex_->transform().unproject(TO_SLIVR_POINT(control_widget_->ReferencePoint()));
    slice_ren_->lock.lock();
    slice_ren_->set_control_point(p);
    slice_ren_->lock.unlock();
  }
  Point w(control_widget_->ReferencePoint());
  control_x_.set(w.x());
  control_y_.set(w.y());
  control_z_.set(w.z());
  control_pos_saved_.set(1);
}

} // namespace SCIRun
