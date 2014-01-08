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
  1 , 0 = Uniform cubic b-spline
  0 , 0.5 = Catmull-Rom
  1/3, 1/3 = sweet spot
*/
  
TODO: review for usefulness in SCIRun 5 (not currently included in build)

#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/PathPort.h>
#include <Dataflow/Network/Module.h>

#include <Dataflow/Widgets/CrosshairWidget.h>
#include <Dataflow/Widgets/ViewWidget.h>

#include <Core/Geom/View.h>

#include <Core/Util/StringUtil.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Core/Geom/Path.h>

#include <Core/Util/Timer.h>
#include <Core/Thread/Time.h>
#include <Core/Thread/CrowdMonitor.h>
#include <Core/Thread/Semaphore.h>
#include <Core/Thread/Mutex.h>

#include <Core/Geometry/Quaternion.h>
#include <Core/Geometry/Transform.h>

#include <teem/limn.h>

#ifndef _WIN32
#include <unistd.h>
#endif

#include <math.h>

#include <vector>
#include <map>
#include <iostream>



namespace SCIRun {

class ShowAndEditCameraWidget : public Module {
  GuiInt		frame_;
  GuiInt		num_frames_;
  GuiDouble		time_;
  GuiString		playmode_;
  GuiString		execmode_;
  GuiInt		track_;
  GuiDouble		B_;
  GuiDouble		C_;

  std::vector<ViewWidget*>		keyframe_widgets_;
  std::vector<std::pair<double,double> >	keyframe_znear_zfar_;
  std::vector<GeomID>		keyframe_widget_geom_ids_;
  std::vector<int>			keyframe_ms_gap_delay_;
  CrowdMonitor          keyframe_widget_lock_;

  std::vector<View>		path_views_;
  std::vector<GeomID>	path_geom_ids_;
  Array1<limnCamera>	path_limnCameras_;
  
public:
  ShowAndEditCameraWidget(GuiContext* ctx);
  virtual ~ShowAndEditCameraWidget();
  virtual void execute();
  virtual void tcl_command(GuiArgs&, void*);
  virtual void widget_moved(bool last, BaseWidget *widget);

  void add_keyframe();
  void insert_keyframe(unsigned int after_frame_num);
  void delete_keyframe(unsigned int frame_to_delete);

  void draw_keyframe_widgets();
  void hide_keyframe_widgets();
  void show_keyframe_widgets();

  void draw_path();
  void hide_path();
  void show_path();
  
};

DECLARE_MAKER(ShowAndEditCameraWidget)
ShowAndEditCameraWidget::ShowAndEditCameraWidget(GuiContext* ctx)
: Module("ShowAndEditCameraWidget", ctx, Filter, "Render", "SCIRun"),
  frame_(get_ctx()->subVar("frame")),
  num_frames_(get_ctx()->subVar("num_frames")),
  time_(get_ctx()->subVar("time")),
  playmode_(get_ctx()->subVar("playmode")),
  execmode_(get_ctx()->subVar("execmode")),
  track_(get_ctx()->subVar("track")),
  B_(get_ctx()->subVar("B")),
  C_(get_ctx()->subVar("C")),
  keyframe_widgets_(),
  keyframe_znear_zfar_(),
  keyframe_widget_geom_ids_(),
  keyframe_ms_gap_delay_(),
  keyframe_widget_lock_("ShowAndEditCameraWidget module widget lock"),
  path_views_(),
  path_geom_ids_(),
  path_limnCameras_()
{
}

ShowAndEditCameraWidget::~ShowAndEditCameraWidget()
{
}

void ShowAndEditCameraWidget::execute()
{  
  GuiInt frameReady(get_ctx()->subVar("frameReady"),false);
  do frameReady.reset(); while (!frameReady.get());
  frame_.reset();
  std::cerr << "displaying frame " << frame_.get() << std::endl;

  GeometryOPortHandle ogeom;
  get_oport_handle("Geometry",ogeom);
  
  ogeom->setView(0, 0, path_views_[frame_.get()-1]);
  ogeom->flushViewsAndWait();
  TCLInterface::eval(get_id()+" post_next_frame");
}


void ShowAndEditCameraWidget::tcl_command(GuiArgs& args, void* userdata)
{
  int frame_num;
  if (args[1] == "add_frame") {
    add_keyframe();
  } else if (args[1] == "insert_keyframe") {
    if (string_to_int(args[2],frame_num)) {
      insert_keyframe(frame_num);
    }
  } else if (args[1] == "delete_keyframe") {
    if (string_to_int(args[2],frame_num)) {
      delete_keyframe(frame_num);
    }
  } else if (args[1] == "hide_keyframe_widgets") {
    hide_keyframe_widgets();
  } else if (args[1] == "show_keyframe_widgets") {
    show_keyframe_widgets();
  } else if (args[1] == "sync_to_keyframe_gui") {
  } else if (args[1] == "sync_from_keyframe_gui") {
  } else if (args[1] == "create_frames"){
    draw_path();
  } else if (args[1] == "stop") {
    want_to_execute();
  } else {
    Module::tcl_command(args, userdata);
  }
}



void
ShowAndEditCameraWidget::widget_moved(bool /*last*/, BaseWidget *widget) 
{  
  for (unsigned int w = 0; w < keyframe_widgets_.size(); w++) {
    //    
    if (keyframe_widgets_[w] == widget) {
      if (keyframe_widgets_[w]->GetMode() != 0) {
	keyframe_widgets_[w]->SetCurrentMode(0);
      }
    } else {
      if (keyframe_widgets_[w]->GetMode() != 1) {
        keyframe_widgets_[w]->SetCurrentMode(1);
      }
    }
  }
}


void
ShowAndEditCameraWidget::add_keyframe() 
{
  insert_keyframe(keyframe_widgets_.size());
}


void
ShowAndEditCameraWidget::insert_keyframe(unsigned int frame_number)
{
  ASSERT(frame_number <= keyframe_widgets_.size());

  GeometryOPortHandle ogeom;
  get_oport_handle("Geometry",ogeom);

  GeometryData* geometry=ogeom->getData(0, 0, GEOM_VIEW);
  ASSERT(geometry && geometry->view);

  
  const View *view = geometry->view;
  const double scale = (view->lookat() - view->eyep()).length()/25.0;
  const double aspect = double(geometry->xres)/double(geometry->yres);
  std::cerr << "aspect: " << aspect << "  FOV: " << view->fov() << std::endl 
       << view->eyep() << std::endl 
       << view->lookat() << std::endl 
       << view->up() << std::endl;

  std::vector<ViewWidget *>::iterator widget_iter = keyframe_widgets_.begin();
  std::vector<GeomID>::iterator geomid_iter = keyframe_widget_geom_ids_.begin();
  std::vector<std::pair<double,double> >::iterator z_iter = keyframe_znear_zfar_.begin();

  advance(widget_iter, frame_number);
  advance(geomid_iter, frame_number);
  advance(z_iter, frame_number);

  widget_iter = keyframe_widgets_.insert
    (widget_iter,new ViewWidget(this,&keyframe_widget_lock_,scale,aspect));
  ViewWidget &widget = *(*widget_iter);

  widget.Connect(ogeom.get_rep());
  widget.SetView(*view);
  widget.SetState(1);

  geomid_iter = keyframe_widget_geom_ids_.insert
    (geomid_iter, ogeom->addObj(widget.GetWidget(), 
				 "Keyframe" + to_string(frame_number),
				 &keyframe_widget_lock_));
  z_iter = keyframe_znear_zfar_.insert
    (z_iter, std::make_pair(geometry->znear, geometry->zfar));
     
  ogeom->flush();
}

void
ShowAndEditCameraWidget::delete_keyframe(unsigned int frame_number)
{
  ASSERT(frame_number < keyframe_widgets_.size());
  std::vector<ViewWidget *>::iterator widget_iter = keyframe_widgets_.begin();
  std::vector<GeomID>::iterator geomid_iter = keyframe_widget_geom_ids_.begin();
  std::vector<std::pair<double,double> >::iterator z_iter = keyframe_znear_zfar_.begin();
  advance(widget_iter, frame_number);
  advance(geomid_iter, frame_number);
  advance(z_iter, frame_number);
  keyframe_widgets_.erase(widget_iter);
  keyframe_widget_geom_ids_.erase(geomid_iter);
  keyframe_znear_zfar_.erase(z_iter);
}


void
ShowAndEditCameraWidget::hide_keyframe_widgets() 
{
}

void
ShowAndEditCameraWidget::show_keyframe_widgets() 
{
}


void
ShowAndEditCameraWidget::draw_path() 
{
  num_frames_.reset();
  const int unsigned num = num_frames_.get();
  const unsigned int num_keyframes = keyframe_widgets_.size();

  Array1<limnCamera> keycameras(num_keyframes);
  Array1<double> times(num_keyframes);

  B_.reset();
  C_.reset();

  limnSplineTypeSpec all;
  all.type = limnSplineTypeBC;
  all.B = B_.get();
  all.C = C_.get();

  limnSplineTypeSpec quatType = all, posType = all, 
    distType = all, viewType = all;
  distType.type = limnSplineTypeLinear;
  int current_time_ms = 0;
  for (unsigned int k = 0; k < num_keyframes; ++k) 
  {

    times[k] = current_time_ms;
    current_time_ms +=1;
    //    if (k < keyframe_ms_gap_delay_[k]) {
    //  current_time_ms += keyframe_ms_gap_delay_[k];
      //}
    //    current_time_ms = 1000;

    ViewWidget* widget = keyframe_widgets_[k];
    keycameras[k].fov = widget->GetFOV();
    keycameras[k].aspect = widget->GetAspectRatio();
    
    const View view = widget->GetView();

    std::cerr << k << ": aspect: " << keycameras[k].aspect << "  FOV: " 
         << view.fov() << std::endl 
         << view.eyep() << std::endl 
         << view.lookat() << std::endl 
         << view.up() << std::endl;


    keycameras[k].from[0] = view.eyep().x();
    keycameras[k].from[1] = view.eyep().y();
    keycameras[k].from[2] = view.eyep().z();

    keycameras[k].at[0] = view.lookat().x();
    keycameras[k].at[1] = view.lookat().y();
    keycameras[k].at[2] = view.lookat().z();
    
    keycameras[k].up[0] = view.up().x();
    keycameras[k].up[1] = view.up().y();
    keycameras[k].up[2] = view.up().z();

    keycameras[k].dist = (view.lookat() - view.eyep()).length();

    keycameras[k].neer = keyframe_znear_zfar_[k].first;
    keycameras[k].faar = keyframe_znear_zfar_[k].second;

    keycameras[k].atRelative = 0;
    keycameras[k].orthographic = 0;
    keycameras[k].rightHanded = 1;
  }

  track_.reset();
  
  std::cerr << "\n\nGENERATING\n\n";

  path_limnCameras_.resize(num);
  if (limnCameraPathMake(path_limnCameras_.get_objs(), 
			 num, keycameras.get_objs(), 
			 times.get_objs(), num_keyframes, track_.get(),
			 &quatType, &posType, &distType, &viewType)) 
  {
    std::cerr << "ERROR:" << biffGetDone(LIMN) << std::endl;
    //airMopAdd(mop, err = biffGetDone(LIMN), airFree, airMopAlways);
    //fprintf(stderr, "%s: trouble making camera path:\n%s\n", me, err);
    //airMopError(mop); return 1;
  }

  path_views_.clear();
  path_views_.resize(num);
  for (unsigned int k = 0; k < num; ++k) 
  {
    Point eyep(path_limnCameras_[k].from[0], path_limnCameras_[k].from[1], path_limnCameras_[k].from[2]);
    Point lookat(path_limnCameras_[k].at[0], path_limnCameras_[k].at[1], path_limnCameras_[k].at[2]);
    Vector up(path_limnCameras_[k].up[0], path_limnCameras_[k].up[1], path_limnCameras_[k].up[2]);
    path_views_[k] = View(eyep, lookat, up, path_limnCameras_[k].fov);
    View *view = &path_views_[k];
    std::cerr << "aspect: " << path_limnCameras_[k].aspect << "  FOV: " << view->fov() << std::endl 
             << view->eyep() << std::endl 
             << view->lookat() << std::endl 
             << view->up() << std::endl;
    

    // pathviews_.push_back();
    
  }
}


void
ShowAndEditCameraWidget::hide_path()
{
}

void
ShowAndEditCameraWidget::show_path()
{
}

} // End namespace SCIRun

