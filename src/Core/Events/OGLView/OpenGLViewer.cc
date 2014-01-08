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
//    File   : OpenGLViewer.cc
//    Author : Martin Cole
//    Date   : Sat May 27 08:51:31 2006

#include <sci_gl.h>
#include <sci_values.h>
#include <string.h>

#include <sci_defs/bits_defs.h>

#include <png.h>

#include <Core/Thread/CrowdMonitor.h>
#include <Core/Events/OGLView/OpenGLViewer.h>
#include <Core/Events/SceneGraphEvent.h>
#include <Core/Events/BaseTool.h>
#include <Core/Events/OGLView/Ball.h>
#include <Core/Events/OGLView/ViewRotateTool.h>
#include <Core/Events/OGLView/ViewScaleTool.h>
#include <Core/Events/OGLView/ViewTranslateTool.h>
#include <Core/Events/OGLView/AutoviewTool.h>
#include <Core/Events/OGLView/FrameBufferPickTool.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/Environment.h>
#include <Core/Util/FileUtils.h>
#include <Core/Geom/OpenGLContext.h>
#include <Core/Geom/GeomObj.h>
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/GeomPick.h>
#include <Core/Geom/PointLight.h>
#include <Core/Math/MiscMath.h>
#include <Core/Geom/GeomScene.h>
#include <Core/Geom/GeomCone.h>
#include <Core/Geom/GeomCylinder.h>
#include <Core/Geom/GeomTri.h>
#include <Core/Geom/GeomText.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/GeomSticky.h>
#include <Core/Geom/GeomRenderMode.h>
#include <Core/Geom/HeadLight.h>
#include <Core/Geom/DirectionalLight.h>
#include <slivr/ShaderProgramARB.h>



#include   <iostream>
#include   <fstream>

#ifdef _WIN32
#  include <windows.h>
#  include <winbase.h>
#  include <Core/Thread/Time.h>
#  undef near
#  undef far
#  undef min
#  undef max
#  if !defined(BUILD_CORE_STATIC)
#    define SCISHARE __declspec(dllimport)
#  else
#    define SCISHARE
#  endif
#else
#  define SCISHARE
#  include <sys/time.h>
#endif

namespace SCIRun {

#define DO_REDRAW     0
#define DO_PICK       1
#define DO_GETDATA    2
#define REDRAW_DONE   4
#define PICK_DONE     5
#define DO_IMAGE      6
#define IMAGE_DONE    7

static const int pick_buffer_size = 512;
static const double pick_window = 10.0;
using std::cerr;
using std::cout;
using std::endl;

class ToolInterface :
    public ViewToolInterface, public AutoviewToolInterface,
    public SSTInterface
{
public:
  ToolInterface(View &v, OpenGLViewer *ov) :
    ViewToolInterface(v),
    AutoviewToolInterface(v),
    SSTInterface(),
    viewer_(ov)
  {}

  // Viewing Tools interface
  virtual int width() const { return viewer_->width(); }
  virtual int height() const { return viewer_->height(); }
  virtual bool compute_depth(const View& view, double& near, double& far)
  { return viewer_->compute_depth(view, near, far); }
  virtual void update_mode_string(const std::string& s) const
  { viewer_->update_mode_string(s); }
  virtual void need_redraw() const { viewer_->need_redraw(); }

  // Autoview Interface
  virtual bool get_bounds(BBox &b) const {
    viewer_->get_bounds(b);
    return b.valid();
  }

  // Selection Set Interface
  virtual void set_selection_set_visible(bool b) {
    viewer_->set_selection_set_visible(b);
  }
  virtual OpenGLViewer::sel_set_t &get_selection_set() {
    return viewer_->get_selection_set();
  }
  virtual void set_selection_geom(GeomHandle geom) {
    viewer_->set_selection_geom(geom);
  }
  virtual void add_selection(unsigned int idx) {
    viewer_->add_selection(idx);
  }
  virtual void remove_selection(unsigned int idx) {
    viewer_->remove_selection(idx);
  }
  virtual void clear_selection_set() {
    viewer_->clear_selection_set();
  }

  OpenGLViewer *viewer_;
};


class SGTool : public BaseTool
{
public:
  SGTool(const std::string& name,
	 GeomIndexedGroup *sg,
	 std::map<std::string, bool> &v) :
    BaseTool(name),
    ids_(0),
    sg_(sg),
    visible_(v),
    scene_events_()
  {}

  struct EventFilter
  {
    std::string gname;
    bool visibility_only;

    bool operator()(const event_handle_t &event)
    {
      SceneGraphEvent *ev = dynamic_cast<SceneGraphEvent*>(event.get_rep());
      if (ev->get_geom_obj_name() == gname)
      {
        if (visibility_only && !ev->toggle_visibility_p()) return false;
        return true;
      }
      return false;
    }
  };

  virtual
  propagation_state_e process_event(event_handle_t event)
  {
    SceneGraphEvent *ev = 0;
    if ((ev = dynamic_cast<SceneGraphEvent*>(event.get_rep())))
    {
      // Check for events that replace other events and clear them.
      // Otherwise the events just pile up indefinitely waiting on
      // redraws on windows that may not even be open/visible.
      // Replace any geometry set event because the prior contents
      // will never be drawn.  Also replace any visibility event if
      // we're forcing the visibility to be a particular value later.
      if (!ev->toggle_visibility_p() || ev->force_visibility_p())
      {
        // This is a replacement event, we can safely remove all of the
        // other events with this identifier.
        const std::string &gname = ev->get_geom_obj_name();
        EventFilter filter;
        filter.gname = gname;
        filter.visibility_only = ev->toggle_visibility_p();
        scene_events_.erase(std::remove_if(scene_events_.begin(),
                                           scene_events_.end(), filter),
                            scene_events_.end());
      }

      // Cache for current context.
      scene_events_.push_back(ev);
      return STOP_E;
    }

    return CONTINUE_E;
  }

  //! Call when context is current!
  void handle_waiting_events()
  {
    for (unsigned int e = 0; e < scene_events_.size(); ++e)
    {
      SceneGraphEvent *ev =
	dynamic_cast<SceneGraphEvent*>(scene_events_[e].get_rep());
      ASSERT(ev);

      const std::string &gname = ev->get_geom_obj_name();

      if (ev->toggle_visibility_p()) {
	// Toggle all sg objects containing the geometry id string.
	std::map<std::string, bool>::iterator iter = visible_.begin();
	while (iter != visible_.end())
        {
	  std::string name = (*iter).first;
	  bool val = ! ((*iter).second);
          if (ev->force_visibility_p()) val = ev->forced_visibility_flag();
	  ++iter;
	  if (name.find(gname) != std::string::npos)
          {
	    visible_[name] = val;
	  }
	}
      } else {

	std::map<std::string, int>::iterator iter;
	iter = obj_ids_.find(gname);
	if (iter != obj_ids_.end())
        {
	  // Already showing this geometry.  Delete it first.
	  int id = (*iter).second;
	  visible_[gname] = false;
	  sg_->delObj(id);
	  obj_ids_.erase(gname);
	}
        GeomHandle geom_obj = ev->get_geom_obj();
        if (geom_obj.get_rep())
        {
          GeomViewerItem* si =
            new GeomViewerItem(geom_obj, gname, 0);
          bool vis = true;
          if (ev->force_visibility_p()) vis = ev->forced_visibility_flag();
          visible_[gname] = vis;
          obj_ids_[gname] = ids_;
          sg_->addObj(si, ids_++);
        }
      }
    }
    scene_events_.clear();
  }

private:
  int                     ids_;
  GeomIndexedGroup       *sg_;
  std::map<std::string, bool>	 &visible_;
  std::map<std::string, int>        obj_ids_;
  std::vector<event_handle_t>  scene_events_;
};


class FBI : public FBInterface
{
public:
  FBI(OpenGLViewer *v) :
    oglv_(v)
  {}

  virtual int width() { return oglv_->width(); }
  virtual int height() { return oglv_->height(); }
  virtual std::vector<unsigned char> &img() { return oglv_->get_fbpick_image(); }
  virtual void do_pick_draw() { oglv_->do_fbpick_draw(); }
  virtual void add_selection(unsigned int id)
  {
    oglv_->add_selection(id);
  }
  virtual void remove_selection(unsigned int id)
  {
    oglv_->remove_selection(id);
  }

private:
  OpenGLViewer *oglv_;
};


class ViewerCommandTool : public CommandTool
{
public:

  ViewerCommandTool(const std::string& name, OpenGLViewer *oglv) :
    CommandTool(name),
    oglv_(oglv)
  {}

  virtual propagation_state_e issue_command(const std::string &cmmd,
					    unsigned int /*time*/)
  {
    if (cmmd == "delete faces") {

      oglv_->delete_selected_faces();
      return STOP_E;

    } else if (cmmd == "add face") {

      oglv_->add_selected_face();
      return STOP_E;
    }
    return CONTINUE_E;
  }
private:
  OpenGLViewer           *oglv_;
};


class ToolManip : public TMNotificationTool
{
public:
  ToolManip(OpenGLViewer *oglv) :
    TMNotificationTool("OpenGLViewer:ToolManager:ToolManip"),
    oglv_(oglv),
    fbpt_(0)
  {}

  virtual propagation_state_e start_tool(const std::string& id, unsigned int /*time*/,
					 const std::string& mode)
  {
    const std::string fbpick_name("FBPickTool");
    const std::string rmfaces_name("RMFacesTool");
    cerr << "start_tool: " << id << ", " << mode << endl;
    if (id == fbpick_name) {
      if (! fbpt_) {
	fbpt_ = new FrameBufferPickTool(fbpick_name,
					new FBI(oglv_));
      }
      ToolManager &tm = oglv_->get_tm();
      if (tm.query_tool_id(OpenGLViewer::SELECTION_TOOL_E) == fbpick_name)
      {
	//Already active, tell the tool to reset.
	fbpt_->reset();
      } else {
	tm.add_tool(fbpt_, OpenGLViewer::SELECTION_TOOL_E);
      }

      if (mode == "nodes") {
	cerr << "sel mode nodes" << endl;
	oglv_->set_selection_mode(SelectionSetTool::NODES_E);
      } else {
	cerr << "sel modes faces" << endl;
	oglv_->set_selection_mode(SelectionSetTool::FACES_E);
      }
      return STOP_E;
    }

    return CONTINUE_E;
  }

  virtual propagation_state_e stop_tool(const std::string& /*id*/, unsigned int /*time*/)
  {
    return CONTINUE_E;
  }

  virtual propagation_state_e suspend_tool(const std::string& /*id*/, unsigned int /*time*/)
  {
    if (! fbpt_) { return STOP_E; }
    ToolManager &tm = oglv_->get_tm();
    cerr << "suspend_tool(string id, unsigned int time)" << endl;

    tool_handle_t th = new BaseTool("suspend selection");

    tm.add_tool(th, OpenGLViewer::SELECTION_TOOL_E);
    return STOP_E;
  }

  virtual propagation_state_e resume_tool(const std::string& id, unsigned int time,
					  const std::string& mode)
  {
    ToolManager &tm = oglv_->get_tm();
    cerr << "resume_tool(...)" << id << ", " << mode << endl;
    if (! fbpt_) {
      start_tool(id, time, mode);
    } else if (tm.query_tool_id(OpenGLViewer::SELECTION_TOOL_E) ==
	       "suspend selection")
    {
      if (mode == "nodes") {
	cerr << "sel mode nodes" << endl;
	oglv_->set_selection_mode(SelectionSetTool::NODES_E);
      } else {
	cerr << "sel modes faces" << endl;
	oglv_->set_selection_mode(SelectionSetTool::FACES_E);
      }
      tm.rm_tool(OpenGLViewer::SELECTION_TOOL_E);
    }
    return STOP_E;
  }

private:
  OpenGLViewer                    *oglv_;
  FrameBufferPickTool             *fbpt_;
};


OpenGLViewer::OpenGLViewer(OpenGLContext *oglc, const std::string& id) :
  OpenGLBase(oglc, id),
  do_autoview_(true),
  doing_image_p_(false),
  doing_movie_p_(false),
  make_MPEG_p_(false),
  current_movie_frame_(0),
  movie_name_("./movie.%04d"),
  // private member variables
  drawinfo_(0),
  do_hi_res_(false),
  encoding_mpeg_(false),
  max_gl_lights_(0),
  animate_num_frames_(0),
  animate_time_begin_(0.0),
  animate_time_end_(0.0),
  animate_framerate_(0.0),
  znear_(0.0),
  zfar_(0.0),
  current_time_(0.0),
  frame_count_(1),
  cached_view_(),
  send_pick_x_(0),
  send_pick_y_(0),
  ret_pick_index_(0),
  ret_pick_obj_(0),
  ret_pick_pick_(0),
  homeview_(Point(2.1, 1.6, 11.5), Point(.0, .0, .0), Vector(0,1,0), 20),
  view_(homeview_),
  do_stereo_(false),
  ambient_scale_(1.0),	
  diffuse_scale_(1.0),	
  specular_scale_(0.4),	
  shininess_scale_(1.0),	
  emission_scale_(1.0),	
  line_width_(1.0),	
  point_size_(1.0),	
  polygon_offset_factor_(1.0),
  polygon_offset_units_(0.0),
  eye_sep_base_(0.4),
  ortho_view_(false),
  fogusebg_(true),
  fogcolor_(Color(0.0,0.0,1.0)),
  fog_start_(0.0),
  fog_end_(0.714265),
  fog_visibleonly_(true),
  fbpick_(false),
  focus_sphere_(new GeomSphere),
  scene_graph_(new GeomIndexedGroup()),
  visible_(),
  obj_tag_(),
  draw_type_(GOURAUD_E),
  capture_z_data_(false),
  need_redraw_(true),
  selection_set_(),
  selection_set_visible_(false),
  selection_geom_(0),
  selection_set_tool_(0),
  scene_graph_tool_(0)
{
  fps_timer_.start();

  // Add a headlight
  lighting_.lights.add(new HeadLight("Headlight", Color(1,1,1)));
  for(int i = 1; i < 4; i++){ // only set up 3 more lights
    std::ostringstream str;
    str << "Light" << i;
    lighting_.lights.add(new DirectionalLight(str.str(),
						 Vector(0,0,1),
						 Color(1,1,1), false, false));
  }

  // 0 - Axes, visible
  internal_objs_.push_back(new GeomViewerItem(create_viewer_axes(),
						 "Axis",0));
  internal_objs_visible_p_.push_back(true);

  // 1 - Unicam control sphere, not visible by default.
  MaterialHandle focus_color = new Material(Color(0.0, 0.0, 1.0));
  internal_objs_.push_back(new GeomMaterial(focus_sphere_, focus_color));
  internal_objs_visible_p_.push_back(false);

  default_material_ =
    new Material(Color(.1,.1,.1), Color(.6,0,0), Color(.7,.7,.7), 50);

  // Setup the tools in the right slots.
  init_tool_manager();
}


OpenGLViewer::~OpenGLViewer()
{
  // Finish up the mpeg that was in progress.
  if (encoding_mpeg_)
  {
    encoding_mpeg_ = false;
    EndMpeg();
  }
  fps_timer_.stop();

  delete drawinfo_;
  drawinfo_ = 0;
}


void
OpenGLViewer::set_draw_axis(bool drawit)
{
  if (internal_objs_visible_p_.size()) internal_objs_visible_p_[0] = drawit;
}


OpenGLViewer::sel_set_t&
OpenGLViewer::get_selection_set()
{
  return selection_set_;
}


void
OpenGLViewer::delete_selected_faces()
{
  SelectionSetTool* sst =
    dynamic_cast<SelectionSetTool*>(selection_set_tool_.get_rep());
  if (sst) { sst->delete_faces(); }
}


void
OpenGLViewer::add_selected_face()
{
  SelectionSetTool* sst =
    dynamic_cast<SelectionSetTool*>(selection_set_tool_.get_rep());
  if (sst) { sst->add_face(); }
}


void
OpenGLViewer::add_selection(unsigned int idx)
{
  // Only add an index once.
  sel_set_t::iterator it = find(selection_set_.begin(),
                                selection_set_.end(), idx);
  if (it == selection_set_.end()) {
    selection_set_.push_back(idx);
    selection_set_visible_ = true;
  }
}


void
OpenGLViewer::remove_selection(unsigned int idx)
{
  sel_set_t::iterator it = find(selection_set_.begin(),
                                selection_set_.end(), idx);
  if (it != selection_set_.end()) {
    selection_set_.erase(it);
  }
  if (selection_set_.size() == 0) selection_set_visible_ = false;
}


void
OpenGLViewer::clear_selection_set()
{
  selection_set_.clear();
  selection_set_visible_ = false;
}


void
OpenGLViewer::set_selection_geom(GeomHandle geom)
{
  GeomViewerItem* si = new GeomViewerItem(geom, "Selection Geom", 0);
  selection_geom_ = si;
}


void
OpenGLViewer::set_selection_set_visible(bool b)
{
  selection_set_visible_ = b;
}


void
OpenGLViewer::set_selection_mode(SelectionSetTool::selection_mode_e m)
{
  SelectionSetTool *sst =
    dynamic_cast<SelectionSetTool*>(selection_set_tool_.get_rep());
  if (sst) { sst->set_selection_mode(m); }
}


void
OpenGLViewer::init_tool_manager()
{
  // viewer tools (always on the priority queue)
  unsigned int view_pri = VIEWER_TOOLS_E;
  ToolInterface* ti = new ToolInterface(view_, this);
  tool_handle_t rot = new ViewRotateTool("OpenGLViewer Rotate Tool", ti);
  tm_.add_tool(rot, view_pri++);

  tool_handle_t scale = new ViewScaleTool("OpenGLViewer Scale Tool", ti);
  tm_.add_tool(scale, view_pri++);

  tool_handle_t trans = new ViewTranslateTool("OpenGLViewer Translate Tool",
					      ti);
  tm_.add_tool(trans, view_pri++);

  tool_handle_t av = new AutoviewTool("OpenGLViewer Autoview Tool", ti);
  tm_.add_tool(av, view_pri++);


  // data tools (always on the queue)
  unsigned int data_pri = DATA_TOOLS_E;
  scene_graph_tool_ = new SGTool("OpenGLViewer Scene Graph Tool",
			    scene_graph_, visible_);
  tm_.add_tool(scene_graph_tool_, data_pri++);

  selection_set_tool_ = new SelectionSetTool("OpenGLViewer SelectionSet Tool",
					     ti);
  tm_.add_tool(selection_set_tool_, data_pri++);

  tool_handle_t cmmdtool;
  cmmdtool = new ViewerCommandTool("OpenGLViewer Command Tool",
				   this);
  tm_.add_tool(cmmdtool, data_pri++);

  // tool modifiers
  unsigned int tm_pri = TOOL_MODIFIERS_E;
  ToolManip *tm = new ToolManip(this);
  tool_handle_t tmnotify(tm);
  tm_.add_tool(tmnotify, tm_pri++);
}


GeomHandle
OpenGLViewer::create_viewer_axes()
{
  const Color black(0,0,0), grey(.5,.5,.5);
  MaterialHandle dk_red =   new Material(black, Color(.2,0,0), grey, 20);
  MaterialHandle dk_green = new Material(black, Color(0,.2,0), grey, 20);
  MaterialHandle dk_blue =  new Material(black, Color(0,0,.2), grey, 20);
  MaterialHandle lt_red =   new Material(black, Color(.8,0,0), grey, 20);
  MaterialHandle lt_green = new Material(black, Color(0,.8,0), grey, 20);
  MaterialHandle lt_blue =  new Material(black, Color(0,0,.8), grey, 20);

  GeomGroup* xp = new GeomGroup;
  GeomGroup* yp = new GeomGroup;
  GeomGroup* zp = new GeomGroup;
  GeomGroup* xn = new GeomGroup;
  GeomGroup* yn = new GeomGroup;
  GeomGroup* zn = new GeomGroup;

  const double sz = 1.0;
  xp->add(new GeomCylinder(Point(0,0,0), Point(sz, 0, 0), sz/20));
  xp->add(new GeomCone(Point(sz, 0, 0), Point(sz+sz/5, 0, 0), sz/10, 0));
  yp->add(new GeomCylinder(Point(0,0,0), Point(0, sz, 0), sz/20));
  yp->add(new GeomCone(Point(0, sz, 0), Point(0, sz+sz/5, 0), sz/10, 0));
  zp->add(new GeomCylinder(Point(0,0,0), Point(0, 0, sz), sz/20));
  zp->add(new GeomCone(Point(0, 0, sz), Point(0, 0, sz+sz/5), sz/10, 0));
  xn->add(new GeomCylinder(Point(0,0,0), Point(-sz, 0, 0), sz/20));
  xn->add(new GeomCone(Point(-sz, 0, 0), Point(-sz-sz/5, 0, 0), sz/10, 0));
  yn->add(new GeomCylinder(Point(0,0,0), Point(0, -sz, 0), sz/20));
  yn->add(new GeomCone(Point(0, -sz, 0), Point(0, -sz-sz/5, 0), sz/10, 0));
  zn->add(new GeomCylinder(Point(0,0,0), Point(0, 0, -sz), sz/20));
  zn->add(new GeomCone(Point(0, 0, -sz), Point(0, 0, -sz-sz/5), sz/10, 0));
  GeomGroup* all=new GeomGroup;
  all->add(new GeomMaterial(xp, lt_red));
  all->add(new GeomMaterial(yp, lt_green));
  all->add(new GeomMaterial(zp, lt_blue));
  all->add(new GeomMaterial(xn, dk_red));
  all->add(new GeomMaterial(yn, dk_green));
  all->add(new GeomMaterial(zn, dk_blue));

  return all;
}


void
OpenGLViewer::redraw(double tbeg, double tend, int nframes, double framerate)
{
  if (dead_) return;
  animate_time_begin_ = tbeg;
  animate_time_end_ = tend;
  animate_num_frames_ = nframes;
  animate_framerate_ = framerate;
  //FIX_ME  make sure a redraw happens...
}


void
OpenGLViewer::do_fbpick_draw()
{
  fbpick_ = true;
  redraw_frame();
}


void
OpenGLViewer::render_and_save_image()
// FIX_ME these parameters come from the tool via the event.
// 				    int x, int y,
// 				    const std::string & fname,
// 				    const std::string & ftype)

{
  ASSERTFAIL("render_and_save_image unimplimented");
} // end render_and_save_image()


void
OpenGLViewer::redraw_frame()
{
  if (dead_) return;

  do_autoview_ = do_autoview_ || (!xres_ && !yres_);

  // Get the window size
  xres_ = width();
  yres_ = height();

  if (gl_context_) {
    // Make sure our GL context is current
    if (! gl_context_->make_current()) {
      gl_context_->release();
      return; // wait for the context.
    }
    SLIVR::ShaderProgramARB::init_shaders_supported();

    glViewport(0, 0, xres_, yres_);
    // Clear the screen.
    glClearColor(bgcolor().r(), bgcolor().g(), bgcolor().b(), 0);
  }

  // Some GeomObj's can only be deleted while the context is current,
  // so processing the scene graph events here assures this.
  // These sort of delete issues should be managed by the context obj.
  if (scene_graph_tool_.get_rep()) {
    SGTool *sgt = dynamic_cast<SGTool*>(scene_graph_tool_.get_rep());
    sgt->handle_waiting_events();
  }

  if (do_autoview_) {
    tm_.propagate_event(new AutoviewEvent());
    do_autoview_ = false;
  }

  GLint data[1];
  glGetIntegerv(GL_MAX_LIGHTS, data);
  max_gl_lights_=data[0];

  // dump a frame if we are Saving an image, or Recording a movie.
  const bool dump_frame = doing_image_p_ || doing_movie_p_;

  // Setup the view...
  cached_view_ = view_;
  double fovy  = view_.fov();
  const double aspect = double(xres_)/double(yres_);

  if (aspect < 1)
    fovy = RtoD(2*atan(1.0/aspect*tan(DtoR(view_.fov()/2.))));

  if (!drawinfo_) {
	  drawinfo_ = new DrawInfoOpenGL();
  }
  drawinfo_->reset();

  if (do_stereo_p())
  {
    GLboolean supported;
    glGetBooleanv(GL_STEREO, &supported);
    if (!supported)
    {
      do_stereo_ = false;
      static bool warnonce = true;
      if (warnonce)
      {
        cout << "Stereo display selected but not supported.\n";
        warnonce = false;
      }
    }
  }

  drawinfo_->view_                  = view_;
  drawinfo_->ambient_scale_         = ambient_scale_;	
  drawinfo_->diffuse_scale_         = diffuse_scale_;	
  drawinfo_->specular_scale_        = specular_scale_;	
  drawinfo_->shininess_scale_       = shininess_scale_;	
  drawinfo_->emission_scale_        = emission_scale_;	
  drawinfo_->line_width_            = line_width_;	
  drawinfo_->point_size_            = point_size_;	
  drawinfo_->polygon_offset_factor_ = polygon_offset_factor_;
  drawinfo_->polygon_offset_units_  = polygon_offset_units_;

  if (compute_depth(view_, znear_, zfar_))
  {
    // Set up graphics state.
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    set_state(drawinfo_);
    drawinfo_->pickmode_=0;

    CHECK_OPENGL_ERROR();

    // Do the redraw loop for each time value.
    const double dt = (animate_time_end_ - animate_time_begin_)
      / animate_num_frames_;
    const double frametime = animate_framerate_==0?0:1.0/animate_framerate_;
    TimeThrottle throttle;
    throttle.start();
    Vector eyesep(0.0, 0.0, 0.0);
    if (do_stereo_p())
    {
      // gui_sr_ was always 1 as far as I could tell, no idea what it was for.
      // (gui_sr_.get() ? 0.048 : 0.0125);
      const double eye_sep_dist = eye_sep_base_ * 0.048;
      Vector u, v;
      view_.get_viewplane(aspect, 1.0, u, v);
      u.safe_normalize();
      const double zmid = (znear_ + zfar_) / 2.0;
      eyesep = u * eye_sep_dist * zmid;
    }

    if (need_redraw_ && animate_num_frames_ < 1) {
      animate_num_frames_ = 1;
    }
    for (int t = 0; t < animate_num_frames_; t++)
    {
      need_redraw_ = false;
      int n = 1;
      if ( do_stereo_p() ) n = 2;
      for (int i = 0; i < n; i++)
      {
        if ( do_stereo_p() )
        {
          glDrawBuffer( (i == 0) ? GL_BACK_LEFT : GL_BACK_RIGHT);
        }
        else
        {
          if (dump_frame)
          {
            glDrawBuffer(GL_FRONT);
          }
          else
          {
            glDrawBuffer(GL_BACK);
          }
        }
        if (gl_context_) {
          glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
          glViewport(0, 0, xres_, yres_);
        }

        const double modeltime = t * dt + animate_time_begin_;
        //set_current_time(modeltime);

        // render normal
	// Setup view.
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	if (do_ortho_view_p())
	{
	  const double len = (view_.lookat() - view_.eyep()).length();
	  const double yval = tan(fovy * M_PI / 360.0) * len;
	  const double xval = yval * aspect;
	  glOrtho(-xval, xval, -yval, yval, znear_, zfar_);
	}
	else
	{
	  gluPerspective(fovy, aspect, znear_, zfar_);
	}
	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
	Point eyep(view_.eyep());
	Point lookat(view_.lookat());

	if (do_stereo_p()) {
	  if (i == 0) {
	    eyep -= eyesep;
	  } else {
	    eyep += eyesep;
	  }
	}
	Vector up(view_.up());
	gluLookAt(eyep.x(), eyep.y(), eyep.z(),
		  lookat.x(), lookat.y(), lookat.z(),
		  up.x(), up.y(), up.z());
	if (do_hi_res_)
	{
	  setFrustumToWindowPortion();
	}


        // Set up Lighting
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
        const Lighting& lighting = lighting_;
        int idx=0;
        int ii;
        for (ii=0; ii<lighting.lights.size(); ii++)
        {
          LightHandle light = lighting.lights[ii];
          light->opengl_setup(view_, drawinfo_, idx);
        }
        for (ii=0; ii<idx && ii<max_gl_lights_; ii++)
        {
          glEnable((GLenum)(GL_LIGHT0 + ii));
        }
        for (;ii<max_gl_lights_;ii++)
        {
          glDisable((GLenum)(GL_LIGHT0 + ii));
        }

        // Now set up the fog stuff.
        double fognear, fogfar;
        compute_fog_depth(view_, fognear, fogfar, fog_visibleonly_p());
        glFogi(GL_FOG_MODE, GL_LINEAR);
        const float fnear =
          fognear + (fogfar - fognear) * fog_start_;
        glFogf(GL_FOG_START, fnear);
        const double ffar =
          fognear + (fogfar - fognear) /
          Max(fog_end_, 0.001);
        glFogf(GL_FOG_END, ffar);
        GLfloat bgArray[4];
        if (fogusebg_)
        {
          bgArray[0] = bgcolor().r();
          bgArray[1] = bgcolor().g();
          bgArray[2] = bgcolor().b();
        }
        else
        {
          Color fogcolor(fogcolor_);
          bgArray[0] = fogcolor.r();
          bgArray[1] = fogcolor.g();
          bgArray[2] = fogcolor.b();
        }
        bgArray[3] = 1.0;
        glFogfv(GL_FOG_COLOR, bgArray);

        // FIX_ME clip tool and mouse_action in drawinfo
        //setClip(drawinfo_);
        //setMouse(drawinfo_);

        // UNICAM addition
        glGetDoublev (GL_MODELVIEW_MATRIX, modelview_matrix_);
        glGetDoublev (GL_PROJECTION_MATRIX, projection_matrix_);
        glGetIntegerv(GL_VIEWPORT, viewport_matrix_);

        // set up point size, line size, and polygon offset
        glPointSize(drawinfo_->point_size_);
        glLineWidth(drawinfo_->line_width_);
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
        if (drawinfo_->polygon_offset_factor_ ||
            drawinfo_->polygon_offset_units_)
        {
          glPolygonOffset(drawinfo_->polygon_offset_factor_,
                          drawinfo_->polygon_offset_units_);
          glEnable(GL_POLYGON_OFFSET_FILL);
        }
        else
        {
          glDisable(GL_POLYGON_OFFSET_FILL);
        }

        // Draw it all.
        current_time_ = modeltime;
        draw_visible_scene_graph();

        if (do_rotation_axis_p())
        {
          render_rotation_axis(view_, do_stereo_p(), i, eyesep);
        }
      }

      // Save z-buffer data.
      if (capture_z_data_)
      {
	depth_buffer_.resize(xres_ * yres_);
	capture_z_data_ = false;
	glReadPixels(0, 0, xres_, yres_, GL_DEPTH_COMPONENT, GL_FLOAT,
		     &depth_buffer_[0]);
      }

      if (do_fbpick_p())
      {
	// Save frame buffer pick data.
	glReadBuffer(GL_BACK);
	fbpick_image_.resize(xres_ * yres_ * 4);
	glReadPixels(0, 0, xres_, yres_, GL_RGBA, GL_UNSIGNED_BYTE,
		     &fbpick_image_[0]);

	fbpick_ = false;  // just do the draw once.
	if (gl_context_) { gl_context_->release(); }
	return;
      }

      // Wait for the right time before swapping buffers
      const double realtime = t * frametime;
      if (animate_num_frames_>1)
      {
        throttle.wait_for_time(realtime);
      }

      //total_frames_.set(total_frames_+1);

      // Show the pretty picture.
      if (gl_context_ && !(dump_frame))
      {
        gl_context_->swap();
      }
    }
    throttle.stop();
    double fps;
    if (throttle.time() > 0)
    {
      fps = animate_num_frames_ / throttle.time();
    }
    else
    {
      fps = animate_num_frames_;
    }
    //set_current_time(animate_time_end_);
  }
  else
  {
    // Just show the cleared screen
    //set_current_time(animate_time_end_);
    if (gl_context_) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      if (! dump_frame)
      {
	gl_context_->swap();
      }
    }
  }
  if (gl_context_) {
    gl_context_->swap();
    gl_context_->release();
  }
}


void
OpenGLViewer::get_pick(int /*x*/, int /*y*/,
		       GeomHandle& /*pick_obj*/, GeomPickHandle& /*pick_pick*/,
		       int& /*pick_index*/)
{
#if 0
  send_pick_x_ = x;
  send_pick_y_ = y;
  send_mailbox_.send(DO_PICK);
  for (;;)
  {
    const int r = recv_mailbox_.receive();
    if (r != PICK_DONE)
    {
      cerr << "WANTED A PICK!!! (got back " << r << "\n";
    }
    else
    {
      pick_obj = ret_pick_obj_;
      pick_pick = ret_pick_pick_;
      pick_index = ret_pick_index_;
      break;
    }
  }
#endif
}


void
OpenGLViewer::real_get_pick(int /*x*/, int /*y*/,
                      GeomHandle& /*pick_obj*/, GeomPickHandle& /*pick_pick*/,
                      int& /*pick_index*/)
{
#if 0
  pick_obj = 0;
  pick_pick = 0;
  pick_index = 0x12345678;
  // Make ourselves current

  // Make sure our GL context is current
  if ((tk_gl_context_ != old_tk_gl_context_))
  {
    old_tk_gl_context_ = tk_gl_context_;
    gui_->lock();
    tk_gl_context_->make_current();
    gui_->unlock();
  }

  // Setup the view...
  View view(view_);
  viewer_->geomlock_.readLock();

  // Compute znear and zfar.
  double znear;
  double zfar;
  if (compute_depth(view, znear, zfar))
  {
    // Setup picking.
    gui_->lock();

    GLuint pick_buffer[pick_buffer_size];
    glSelectBuffer(pick_buffer_size, pick_buffer);
    glRenderMode(GL_SELECT);
    glInitNames();
#ifdef SCI_64BITS
    glPushName(0);
    glPushName(0);
    glPushName(0x12345678);
    glPushName(0x12345678);
    glPushName(0x12345678);
#else
    glPushName(0); //for the pick
    glPushName(0x12345678); //for the object
    glPushName(0x12345678); //for the object's face index
#endif

    // Picking
    const double aspect = double(xres_)/double(yres_);
    // XXX - UNICam change-- should be '1.0/aspect' not 'aspect' below
    const double fovy = RtoD(2*Atan(1.0/aspect*Tan(DtoR(view.fov()/2.))));
    glViewport(0, 0, xres_, yres_);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    gluPickMatrix(x, viewport[3]-y, pick_window, pick_window, viewport);
    if (do_ortho_view_p())
    {
      const double len = (view.lookat() - view.eyep()).length();
      const double yval = tan(fovy * M_PI / 360.0) * len;
      const double xval = yval * aspect;
      glOrtho(-xval, xval, -yval, yval, znear, zfar);
    }
    else
    {
      gluPerspective(fovy, aspect, znear, zfar);
    }

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    Point eyep(view.eyep());
    Point lookat(view.lookat());
    Vector up(view.up());
    gluLookAt(eyep.x(), eyep.y(), eyep.z(),
              lookat.x(), lookat.y(), lookat.z(),
              up.x(), up.y(), up.z());

    drawinfo_->lighting_=0;
    drawinfo_->set_drawtype(DrawInfoOpenGLViewer::Flat);
    drawinfo_->pickmode_=1;

    // Draw it all.
    do_for_visible(this, true);

#ifdef SCI_64BITS
    glPopName();
    glPopName();
    glPopName();
#else
    glPopName();
    glPopName();
#endif

    glFlush();
    const int hits = glRenderMode(GL_RENDER);

//    CHECK_OPENGL_ERROR("OpenGLViewer::real_get_pick");

    gui_->unlock();
    GLuint min_z;
#ifdef SCI_64BITS
    unsigned long hit_obj=0;
    //    GLuint hit_obj_index = 0x12345678;
    unsigned long hit_pick=0;
    //    GLuint hit_pick_index = 0x12345678;  // need for object indexing
#else
    GLuint hit_obj = 0;
    //GLuint hit_obj_index = 0x12345678;  // need for object indexing
    GLuint hit_pick = 0;
    //GLuint hit_pick_index = 0x12345678;  // need for object indexing
#endif
    if (hits >= 1)
    {
      int idx = 0;
      min_z = 0;
      bool have_one = false;
      for (int h=0; h<hits; h++)
      {
        int nnames = pick_buffer[idx++];
        GLuint z=pick_buffer[idx++];
        if (nnames > 1 && (!have_one || z < min_z))
        {
          min_z = z;
          have_one = true;
          idx++; // Skip Max Z
#ifdef SCI_64BITS
          idx += nnames - 5; // Skip to the last one.
          const unsigned int ho1 = pick_buffer[idx++];
          const unsigned int ho2 = pick_buffer[idx++];
          hit_pick = ((long)ho1<<32) | ho2;
          //hit_obj_index = pick_buffer[idx++];
          const unsigned int hp1 = pick_buffer[idx++];
          const unsigned int hp2 = pick_buffer[idx++];
          hit_obj = ((long)hp1<<32)|hp2;
          //hit_pick_index = pick_buffer[idx++];
          idx++;
#else
          // hit_obj=pick_buffer[idx++];
          // hit_obj_index=pick_buffer[idx++];
          //for (int i=idx; i<idx+nnames; ++i) cerr << pick_buffer[i] << "\n";
          idx += nnames - 3; // Skip to the last one.
          hit_pick = pick_buffer[idx++];
          hit_obj = pick_buffer[idx++];
          idx++;
          //hit_pick_index=pick_buffer[idx++];
#endif
        }
        else
        {
          idx += nnames + 1;
        }
      }

      pick_obj = (GeomObj*)hit_obj;
      pick_pick = (GeomPick*)hit_pick;
      pick_obj->getId(pick_index); //(int)hit_pick_index;
    }
  }
  viewer_->geomlock_.readUnlock();
#endif
}


void
OpenGLViewer::draw_visible_scene_graph()
{
  // Do internal objects first...
  unsigned int i;
  for (i = 0; i < internal_objs_.size(); i++){
    if (internal_objs_visible_p_[i] == 1) {
      if (do_picking_p()) {
	pick_draw_obj(default_material_, internal_objs_[i].get_rep());
      } else {
	redraw_obj(default_material_, internal_objs_[i].get_rep());
      }
    }
  }

  if (!scene_graph_) return;
  for (int pass=0; pass < 4; pass++)
  {
    GeomIndexedGroup::IterIntGeomObj iter = scene_graph_->getIter();
    for ( ; iter.first != iter.second; iter.first++) {
      GeomViewerItem *si = (GeomViewerItem*)((*iter.first).second.get_rep());
      // Look up the name to see if it should be drawn...

      if (item_visible_p(si))
      {
	const bool transparent =
	  strstr(si->getString().c_str(), "TransParent") ||
	  strstr(si->getString().c_str(), "Transparent");
	const bool culledtext = strstr(si->getString().c_str(), "Culled Text");
	const bool sticky = strstr(si->getString().c_str(), "Sticky");
	if ((pass == 0 && !transparent && !culledtext && !sticky) ||
	    (pass == 1 && transparent && !culledtext && !sticky) ||
	    (pass == 2 && culledtext && !sticky) ||
	    (pass == 3 && sticky))
	{
	  if(si->crowd_lock()){
	    si->crowd_lock()->readLock();
	  }
	
	  if (do_picking_p()) {
	    pick_draw_obj(default_material_, si);
	  } else {
	    redraw_obj(default_material_, si);
	  }
	
	  if(si->crowd_lock()) {
	    si->crowd_lock()->readUnlock();
	  }
	}
      }
    }
  }

  // Render the selected geometry last.
  if (selection_set_visible_)
  {
    // make sure the selection_geom_ is up to date.
    SelectionSetTool *sst =
      (SelectionSetTool*)(selection_set_tool_.get_rep());

    if (!sst) {
      cerr << "null selection set tool" << endl;
      return;
    }

    sst->render_selection_set();
    if (selection_geom_.get_rep()) {
      redraw_obj(default_material_, selection_geom_.get_rep());
    }
  }
}


// Dump a ppm image.
void
OpenGLViewer::dump_image(const std::string& fname, const std::string& ftype)
{
  GLint vp[4];
  glGetIntegerv(GL_VIEWPORT, vp);
  const int pix_size = 3;  // for RGB
  const int n = pix_size * vp[2] * vp[3];
  unsigned char* pixels = new unsigned char[n];
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadBuffer(GL_FRONT);
  glReadPixels(0, 0, vp[2], vp[3], GL_RGB, GL_UNSIGNED_BYTE, pixels);

  if (ftype == "png")
  {
    // Create the PNG struct.
    png_structp png =
      png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);

    if (png == 0) {
      //setMovieMessage( "ERROR - Failed to create PNG write struct", true );
      return;
    }

    // Create the PNG info struct.
    png_infop info = png_create_info_struct(png);

    if (info == 0) {
      //setMovieMessage( "ERROR - Failed to create PNG info struct", true );
      png_destroy_write_struct(&png, 0);
      return;
    }

    if (setjmp(png_jmpbuf(png))) {
      //setMovieMessage( "ERROR - Initializing PNG.", true );
      png_destroy_write_struct(&png, &info);
      return;
    }

    // Initialize the PNG IO.
    FILE *fp = fopen(fname.c_str(), "wb");
    if (0 == fp) {
      std::string errorMsg = "ERROR opening file: " + fname;
      cerr << errorMsg << "\n";
      return;
    }
    png_init_io(png, fp);

    png_set_IHDR(png, info, vp[2], vp[3],
		 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    // Write the PNG header.
    png_write_info(png, info);

    // Run loop to divide memory into "row" chunks
    png_bytep *rows = (png_bytep*)malloc(sizeof(png_bytep) * vp[3]);
    if (0 == rows) {
      std::string errorMsg = "ERROR allocating memory for png";
      cerr << errorMsg << "\n";
      return;
    }
    for (int hi = 0; hi < vp[3]; hi++) {
      rows[hi] = &((png_bytep)pixels)[(vp[3] - hi - 1) * vp[2] * pix_size];
    }

    png_set_rows(png, info, rows);

    png_write_image(png, rows);

    /* end write */
    if (setjmp(png_jmpbuf(png))) {
      //setMovieMessage( "Error during end of PNG write", true );
      png_destroy_write_struct(&png, &info);
      return;
    }

    // Finish writing.
    png_write_end(png, 0);

    // More clean up.
    png_destroy_write_struct(&png, &info);
    fclose(fp);
    free(rows);
  }
  else
  {
    std::ofstream dumpfile(fname.c_str());
    if ( !dumpfile )
    {
      std::string errorMsg = "ERROR opening file: " + fname;
      cerr << errorMsg << "\n";
      return;
    }

    // Print out the ppm  header.
    dumpfile << "P6" << std::endl;
    dumpfile << vp[2] << " " << vp[3] << std::endl;
    dumpfile << 255 << std::endl;

    // OpenGL renders upside-down to ppm_file writing.
    unsigned char *top_row, *bot_row;
    unsigned char *tmp_row = new unsigned char[ vp[2] * pix_size];
    int top, bot;
    for ( top = vp[3] - 1, bot = 0; bot < vp[3]/2; top --, bot++){
      top_row = pixels + vp[2] * top * pix_size;
      bot_row = pixels + vp[2] * bot * pix_size;
      memcpy(tmp_row, top_row, vp[2] * pix_size);
      memcpy(top_row, bot_row, vp[2] * pix_size);
      memcpy(bot_row, tmp_row, vp[2] * pix_size);
    }
    // Now dump the file.
    dumpfile.write((const char *)pixels, n);
    delete [] tmp_row;
  }

  delete [] pixels;
}


void
OpenGLViewer::put_scanline(int y, int width, Color* scanline, int repeat)
{
  float* pixels = new float[width*3];
  float* p = pixels;
  int i;
  for (i=0; i<width; i++)
  {
    *p++ = scanline[i].r();
    *p++ = scanline[i].g();
    *p++ = scanline[i].b();
  }
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslated(-1.0, -1.0, 0.0);
  glScaled(2.0 / xres_, 2.0 / yres_, 1.0);
  glDepthFunc(GL_ALWAYS);
  glDrawBuffer(GL_FRONT);
  for (i=0; i<repeat; i++)
  {
    glRasterPos2i(0, y + i);
    glDrawPixels(width, 1, GL_RGB, GL_FLOAT, pixels);
  }
  glDepthFunc(GL_LEQUAL);
  glDrawBuffer(GL_BACK);
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  delete[] pixels;
}


void
OpenGLViewer::pick_draw_obj(MaterialHandle def, GeomHandle obj)
{
#ifdef SCI_64BITS
   /// @todo SKETCHY CODE!  Fixme.
  unsigned long long o = reinterpret_cast<unsigned long long>(obj.get_rep());
  unsigned int o1 = (o>>32)&0xffffffff;
  unsigned int o2 = o&0xffffffff;
  glPopName();
  glPopName();
  glPopName();
  glPushName(o1);
  glPushName(o2);
  glPushName(0x12345678);
#else 
  glPopName();
  glPushName((GLuint)(obj.get_rep()));
  glPushName(0x12345678);
#endif
  obj->draw(drawinfo_, def.get_rep(), current_time_);
}


void
OpenGLViewer::redraw_obj(MaterialHandle def, GeomHandle obj)
{
  GeomViewerItem *gvi = dynamic_cast<GeomViewerItem *>(obj.get_rep());
  ASSERT(gvi);
  std::string name = gvi->getString();
  std::map<std::string,int>::iterator tag_iter = obj_tag_.find(name);
  if (tag_iter != obj_tag_.end()) {
    // if found
    set_state(drawinfo_);
  }
  if (do_fbpick_p()) {
    // for now only draw faces
    SelectionSetTool *sst =
      dynamic_cast<SelectionSetTool*>(selection_set_tool_.get_rep());
    if (! sst) return;
    SelectionSetTool::selection_mode_e mode;
    mode = sst->get_selection_mode();
    bool do_draw = false;

    if (mode == SelectionSetTool::NODES_E) {
      if (name.find("Node") != std::string::npos) {
	cerr << "fbpick draw for : " << name << endl;
	do_draw = true;
      }
    } else if (mode == SelectionSetTool::FACES_E) {
      if (name.find("Face") != std::string::npos) {
	do_draw = true;
      }
    }

    if (do_draw) {
      obj->fbpick_draw(drawinfo_, def.get_rep(), current_time_);
    }

  } else {
    obj->draw(drawinfo_, def.get_rep(), current_time_);
  }
}


void
OpenGLViewer::deriveFrustum()
{
  double pmat[16];
  glGetDoublev(GL_PROJECTION_MATRIX, pmat);
  const double G = (pmat[10]-1)/(pmat[10]+1);
  frustum_.znear = -(pmat[14]*(G-1))/(2*G);
  frustum_.zfar = frustum_.znear*G;

  if (do_ortho_view_p())
  {
    frustum_.left = (pmat[8]-1)/pmat[0];
    frustum_.right = (pmat[8]+1)/pmat[0];
    frustum_.bottom = (pmat[9]-1)/pmat[5];
    frustum_.top = (pmat[9]+1)/pmat[5];
    frustum_.width = frustum_.right - frustum_.left;
    frustum_.height = frustum_.top - frustum_.bottom;
  }
  else
  {
    frustum_.left = frustum_.znear*(pmat[8]-1)/pmat[0];
    frustum_.right = frustum_.znear*(pmat[8]+1)/pmat[0];
    frustum_.bottom = frustum_.znear*(pmat[9]-1)/pmat[5];
    frustum_.top = frustum_.znear*(pmat[9]+1)/pmat[5];
    frustum_.width = frustum_.right - frustum_.left;
    frustum_.height = frustum_.top - frustum_.bottom;
  }
}


void
OpenGLViewer::setFrustumToWindowPortion()
{
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  if (do_ortho_view_p())
  {
    glOrtho(frustum_.left + frustum_.width / hi_res_.ncols * hi_res_.col,
            frustum_.left + frustum_.width / hi_res_.ncols * (hi_res_.col+1),
            frustum_.bottom + frustum_.height / hi_res_.nrows * hi_res_.row,
            frustum_.bottom + frustum_.height / hi_res_.nrows *(hi_res_.row+1),
            znear_, zfar_);
  }
  else
  {
    glFrustum(frustum_.left + frustum_.width / hi_res_.ncols * hi_res_.col,
              frustum_.left + frustum_.width / hi_res_.ncols * (hi_res_.col+1),
              frustum_.bottom + frustum_.height / hi_res_.nrows* hi_res_.row,
              frustum_.bottom + frustum_.height /hi_res_.nrows*(hi_res_.row+1),
              frustum_.znear, frustum_.zfar);
  }
}

void
OpenGLViewer::StartMpeg(const std::string& fname)
{
#ifdef HAVE_MPEG
  // Get a file pointer pointing to the output file.
  mpeg_file_ = fopen(fname.c_str(), "w");
  if (!mpeg_file_)
  {
    cerr << "Failed to open file " << fname << " for writing\n";
    return;
  }
  // Get the default options.
  MPEGe_default_options( &mpeg_options_ );
  // Change a couple of the options.
  char *pattern = new char[4];
  strncpy(pattern, "II\0", 4);
  mpeg_options_.frame_pattern = pattern;
  mpeg_options_.search_range[1]=0;
  mpeg_options_.gop_size=1;
  mpeg_options_.IQscale=1;
  mpeg_options_.PQscale=1;
  mpeg_options_.BQscale=1;
  mpeg_options_.pixel_search=MPEGe_options::FULL;
  if ( !MPEGe_open(mpeg_file_, &mpeg_options_ ) )
  {
    cerr << "MPEGe library initialisation failure!:" <<
      mpeg_options_.error << "\n";
    return;
  }
#endif // HAVE_MPEG
}


void
OpenGLViewer::AddMpegFrame()
{
#ifdef HAVE_MPEG
  // Looks like we'll blow up if you try to make more than one movie
  // at a time, as the memory per frame is static.
  static ImVfb *image=0; /* Only alloc memory for these once. */
  int width, height;
  ImVfbPtr ptr;

  GLint vp[4];
  glGetIntegerv(GL_VIEWPORT, vp);

  width = vp[2];
  height = vp[3];

  // Set up the ImVfb used to store the image.
  if ( !image )
  {
    image=MPEGe_ImVfbAlloc( width, height, IMVFBRGB, true );
    if ( !image )
    {
      cerr<<"Couldn't allocate memory for frame buffer\n";
      exit(2);
    }
  }

  // Get to the first pixel in the image.
  ptr = ImVfbQPtr( image, 0, 0 );
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadBuffer(GL_FRONT);
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, ptr);


  const int r = 3 * width;
  unsigned char* row = new unsigned char[r];
  unsigned char* p0, *p1;

  int k, j;
  for (k = height -1, j = 0; j < height/2; k--, j++)
  {
    p0 = ptr + r * j;
    p1 = ptr + r * k;
    memcpy( row, p0, r);
    memcpy( p0, p1, r);
    memcpy( p1, row, r);
  }
  delete[] row;

  if ( !MPEGe_image(image, &mpeg_options_) )
  {
    //setMovieMessage( std::string("ERROR creating MPEG frame: ") + mpeg_options_.error, true );
  }
#endif // HAVE_MPEG
}



void
OpenGLViewer::EndMpeg()
{
#ifdef HAVE_MPEG
  if ( !MPEGe_close(&mpeg_options_) )
  {
    std::string errorMsg = std::string("ERROR closing MPEG file: ") + mpeg_options_.error;
    //setMovieMessage( errorMsg, true );
  }
  else
  {
    std::string message = "Ending Mpeg.";
    //setMovieMessage( message );
  }

  //setMovieStopped();
#endif // HAVE_MPEG
}


// Return world-space depth to point under pixel (x, y).
bool
OpenGLViewer::pick_scene( int x, int y, Point *p )
{
  // y = 0 is bottom of screen (not top of screen, which is what X
  // events reports)
  y = (yres_ - 1) - y;
  int index = x + (y * xres_);
  double z = depth_buffer_[index];
  if (p)
  {
    // Unproject the window point (x, y, z).
    GLdouble world_x, world_y, world_z;
    gluUnProject(x, y, z,
                 modelview_matrix_, projection_matrix_, viewport_matrix_,
                 &world_x, &world_y, &world_z);

    *p = Point(world_x, world_y, world_z);
  }

  // if z is close to 1, then assume no object was picked
  return (z < .999999);
}


bool
OpenGLViewer::item_visible_p(GeomViewerItem* si)
{
  std::map<std::string, bool>::iterator viter = visible_.find(si->getString());
  if (viter != visible_.end()) { return (*viter).second; }
  return false;
}


void
OpenGLViewer::get_bounds(BBox &bbox, bool check_visible)
{
  bbox.reset();
  if (scene_graph_) {
    GeomIndexedGroup::IterIntGeomObj iter = scene_graph_->getIter();
    for ( ; iter.first != iter.second; iter.first++) {
      GeomViewerItem *si = (GeomViewerItem*)((*iter.first).second.get_rep());
      // Look up the name to see if it should be drawn...
      if (!check_visible || item_visible_p(si)) {

	if(si->crowd_lock()) si->crowd_lock()->readLock();
	si->get_bounds(bbox);
	if(si->crowd_lock()) si->crowd_lock()->readUnlock();

      }
    }
  }
  const unsigned int objs_size = internal_objs_.size();
  for(unsigned int i = 0; i < objs_size; i++) {
    if (!check_visible || internal_objs_visible_p_[i])
      internal_objs_[i]->get_bounds(bbox);
  }

  // If the bounding box is empty, make it default to sane view.
  if (! bbox.valid()) {
    bbox.extend(Point(-1.0, -1.0, -1.0));
    bbox.extend(Point(1.0, 1.0, 1.0));
  }
}


bool
OpenGLViewer::compute_depth(const View& view, double& znear, double& zfar)
{
  znear = DBL_MAX;
  zfar =- DBL_MAX;
  BBox bb;
  get_bounds(bb);
  if (bb.valid())
  {
    // We have something to draw.
    Point min(bb.min());
    Point max(bb.max());
    Point eyep(view.eyep());
    Vector dir(view.lookat()-eyep);
    const double dirlen2 = dir.length2();
    if (dirlen2 < 1.0e-6 || dirlen2 != dirlen2)
    {
      return false;
    }
    dir.safe_normalize();
    const double d = -Dot(eyep, dir);
    for (int i=0;i<8;i++)
    {
      const Point p((i&1)?max.x():min.x(),
                    (i&2)?max.y():min.y(),
                    (i&4)?max.z():min.z());
      const double dist = Dot(p, dir) + d;
      znear = Min(znear, dist);
      zfar = Max(zfar, dist);
    }
    znear *= 0.99;
    zfar  *= 1.01;

    if (znear <= 0.0)
    {
      if (zfar <= 0.0)
      {
        // Everything is behind us - it doesn't matter what we do.
        znear = 1.0;
        zfar = 2.0;
      }
      else
      {
        znear = zfar * 0.001;
      }
    }
    return true;
  }
  else
  {
    return false;
  }
}


bool
OpenGLViewer::compute_fog_depth(const View &view, double &znear, double &zfar,
                                bool visible_only)
{
  znear = DBL_MAX;
  zfar = -DBL_MAX;
  BBox bb;
  if (visible_only)
  {
    get_bounds(bb);
  }
  else
  {
    get_bounds_all(bb);
  }
  if (bb.valid())
  {
    // We have something to draw.
    Point eyep(view.eyep());
    Vector dir(view.lookat()-eyep);
    const double dirlen2 = dir.length2();
    if (dirlen2 < 1.0e-6 || dirlen2 != dirlen2)
    {
      return false;
    }
    dir.safe_normalize();
    const double d = -Dot(eyep, dir);

    // Compute distance to center of bbox.
    const double dist = Dot(bb.center(), dir);
    // Compute bbox view radius.
    const double radius = bb.diagonal().length() * dir.length2() * 0.5;

    znear = d + dist - radius;
    zfar = d + dist + radius;

    return true;
  }
  else
  {
    return false;
  }
}


// i is the frame number, usually refers to left or right when do_stereo_p()
// is set.
void
OpenGLViewer::render_rotation_axis(const View &view,
				   bool /*do_stereo*/, int i,
				   const Vector &eyesep)
{
  static GeomHandle axis_obj = 0;
  if (axis_obj.get_rep() == 0) axis_obj = create_viewer_axes();

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  const int xysize = Min(viewport[2], viewport[3]) / 4;
  glViewport(viewport[0] + viewport[2] - xysize,
             viewport[1] + viewport[3] - xysize, xysize, xysize);
  const double aspect = 1.0;

  // fovy 16 eyedist 10 is approximately the default axis view.
  // fovy 32 eyedist 5 gives an exagerated perspective.
  const double fovy = 32.0;
  const double eyedist = 5.0;
  const double znear = eyedist - 2.0;
  const double zfar = eyedist + 2.0;

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluPerspective(fovy, aspect, znear, zfar);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  Vector oldeye(view.eyep().asVector() - view.lookat().asVector());
  oldeye.safe_normalize();
  Point eyep((oldeye * eyedist).asPoint());
  Point lookat(0.0, 0.0, 0.0);

  if (do_stereo_p()) {
    if (i == 0) {
      eyep -= eyesep;
    } else {
      eyep += eyesep;
    }
  }

  Vector up(view.up());
  gluLookAt(eyep.x(), eyep.y(), eyep.z(),
            lookat.x(), lookat.y(), lookat.z(),
            up.x(), up.y(), up.z());
  if (do_hi_res_)
  {
    // Draw in upper right hand corner of total image, not viewport image.
    const int xysize = Min(hi_res_.resx, hi_res_.resy) / 4;
    const int xoff = hi_res_.resx - hi_res_.col * viewport[2];
    const int yoff = hi_res_.resy - hi_res_.row * viewport[3];
    glViewport(xoff - xysize, yoff - xysize, xysize, xysize);
  }

  set_state(drawinfo_);

  // Disable fog for the orientation axis.
  const bool fog = drawinfo_->fog_;
  if (fog) { glDisable(GL_FOG); }
  drawinfo_->fog_ = false;

  // Set up Lighting
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  const Lighting& l = lighting_;
  int idx=0;
  int ii;
  for (ii=0;ii<l.lights.size();ii++)
  {
    LightHandle light=l.lights[ii];
    light->opengl_setup(view, drawinfo_, idx);
  }
  for (ii=0;ii<idx && ii<max_gl_lights_;ii++)
    glEnable((GLenum)(GL_LIGHT0+ii));
  for (;ii<max_gl_lights_;ii++)
    glDisable((GLenum)(GL_LIGHT0+ii));

  // Disable clipping planes for the orientation icon.
  std::vector<bool> cliplist(6, false);
  for (ii = 0; ii < 6; ii++)
  {
    if (glIsEnabled((GLenum)(GL_CLIP_PLANE0+ii)))
    {
      glDisable((GLenum)(GL_CLIP_PLANE0+ii));
      cliplist[ii] = true;
    }
  }

  // Use depthrange to force the icon to move forward.
  // Ideally the rest of the scene should be drawn at 0.05 1.0,
  // so there was no overlap at all, but that would require
  // mucking about in the picking code.
  glDepthRange(0.0, 0.05);
  axis_obj->draw(drawinfo_, 0, current_time_);
  glDepthRange(0.0, 1.0);

  drawinfo_->fog_ = fog;  // Restore fog state.
  if (fog) { glEnable(GL_FOG); }

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);

  // Reenable clipping planes.
  for (ii = 0; ii < 6; ii++)
  {
    if (cliplist[ii])
    {
      glEnable((GLenum)(GL_CLIP_PLANE0+ii));
    }
  }
}


void
OpenGLViewer::set_state(DrawInfoOpenGL* drawinfo)
{
  switch (draw_type_) {
  case WIRE_E :
    {
      drawinfo->set_drawtype(DrawInfoOpenGL::WireFrame);
      drawinfo->lighting_=0;
    }
    break;
  case FLAT_E :
    {
      drawinfo->set_drawtype(DrawInfoOpenGL::Flat);
      drawinfo->lighting_=0;
    }
    break;
  case DEFAULT_E  :
  case GOURAUD_E :
  default:
    {
      drawinfo->set_drawtype(DrawInfoOpenGL::Gouraud);
      drawinfo->lighting_=1;
    }
  };

  // Now see if they want a bounding box.
  drawinfo->show_bbox_ = do_bbox_p();

#if 0 // FIX_ME make a movie tool
  if (!doing_movie_p())
  {
    doing_movie_p_ = 0;
    make_MPEG_p_ = 0;
  } else {
    current_movie_frame_ = movie_tool_->frame();

    doing_movie_p_ = 1;
    if (movie == 1)
    {
      renderer_->make_MPEG_p_ = 0;
      renderer_->movie_frame_extension_ = "ppm";
    }
    if (movie == 3)
    {
      renderer_->make_MPEG_p_ = 0;
      renderer_->movie_frame_extension_ = "png";
    }
    else if (movie == 2)
    {
      renderer_->make_MPEG_p_ = 1;
    }
  }
#endif

#if 0 //FIX_ME make clipping tool.
  if (clip.valid())
  {
    drawinfo->check_clip_ = clip;
  }
  setClip(drawinfo);
#endif

  drawinfo->cull_            = do_backface_cull_p();
  drawinfo->display_list_p_  = do_display_list_p();
  drawinfo->fog_             = do_fog_p();
  drawinfo->lighting_        = do_lighting_p();
  drawinfo->currently_lit_   = drawinfo->lighting_;
  drawinfo->init_lighting(drawinfo->lighting_);
}

} // End namespace SCIRun
