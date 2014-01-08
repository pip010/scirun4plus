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
//    File   : OpenGLViewer.h
//    Author : Martin Cole
//    Date   : Sat May 27 08:51:31 2006
//    Much of this code was taken from Dataflow/Modules/Render/*
//    which was mostly written by Steve Parker.


#if !defined(OpenGLViewer_h)
#define OpenGLViewer_h

#include <string>
#include <map>
#include <vector>


#include <Core/Events/OpenGLBase.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/ToolManager.h>
#include <Core/Thread/Runnable.h>
#include <Core/Thread/Thread.h>
#include <Core/Geom/DrawInfoOpenGL.h>
#include <Core/Geom/Light.h>
#include <Core/Geom/Lighting.h>
#include <Core/Geom/GeomObj.h>
#include <Core/Geom/GeomPick.h>
#include <Core/Geom/GeomSphere.h>
#include <Core/Geom/IndexedGroup.h>
#include <Core/Geom/GeomViewerItem.h>
#include <Core/Geom/View.h>
#include <Core/Util/Timer.h>
#include <Core/Events/OGLView/SelectionSetTool.h>

#ifdef HAVE_MPEG
#  include <mpege.h>
#endif // HAVE_MPEG

#include <Core/Events/OGLView/share.h>
namespace SCIRun {

struct Frustum {
  double znear;
  double zfar;
  double left;
  double right;
  double bottom;
  double top;
  double width;
  double height;
};

struct HiRes {
  double nrows;
  double ncols;
  int row;
  int col;
  int resx;
  int resy;
};


class SCISHARE OpenGLViewer : public OpenGLBase
{
public:
  typedef std::vector<unsigned int> sel_set_t;

  enum draw_type_e {
    DEFAULT_E,
    WIRE_E,
    FLAT_E,
    GOURAUD_E
  };
  //! The id is the name the mailbox gets. setting event targets to the
  //! name constructed here, uniquely sends events to this object.
  OpenGLViewer(OpenGLContext*, const std::string& id = "OpenGLViewer");
  virtual ~OpenGLViewer();

  void                  get_pick(int, int, GeomHandle&,
				 GeomPickHandle&, int&);
  virtual void          need_redraw() { need_redraw_ = true; }
  void                  redraw(double tbeg, double tend,
                               int ntimesteps, double frametime);
  bool                  compute_depth(const View& view,
                                      double& near, double& far);
  bool                  compute_fog_depth(const View& view,
                                          double& near, double& far,
                                          bool visible_only);
  void                  save_image(const std::string& fname,
    const std::string& type = "ppm",
				   int x=640, int y=512);

  void                  set_draw_axis(bool drawit);

  void                   do_fbpick_draw();
  std::vector<unsigned char> &get_fbpick_image() { return fbpick_image_; }

  // Compute world space point under cursor (x,y).  If successful,
  // set 'p' to that value & return true.  Otherwise, return false.
  bool                  pick_scene(int x, int y, Point *p);

  bool                  do_stereo_p() { return do_stereo_; }
  bool                  do_hi_res_p() { return do_hi_res_; }

  //! \todo {hook these up to tools for setting state.}
  bool                do_backface_cull_p() { return false; }
  bool                do_display_list_p()  { return false; }
  bool                do_fog_p()           { return false; }
  bool                do_lighting_p()      { return true; }
  bool                do_ortho_view_p()    { return ortho_view_; }
  bool                fog_visibleonly_p()  { return fog_visibleonly_; }
  bool                do_rotation_axis_p() { return true; }
  bool                do_picking_p()       { return false; }
  bool                do_fbpick_p() const  { return fbpick_; }
  bool                do_bbox_p()          { return false; }

  void                update_mode_string(std::string) {}
  void                get_bounds(BBox &bbox, bool check_visible = true);

  // Selection tool functionality.
  sel_set_t&          get_selection_set();
  void                delete_selected_faces();
  void                add_selected_face();
  void                add_selection(unsigned int idx);
  void                remove_selection(unsigned int idx);
  void                clear_selection_set();
  void                set_selection_geom(GeomHandle geom);
  void                set_selection_set_visible(bool b);
  void                set_selection_mode(SelectionSetTool::selection_mode_e m);


  //! the following enum defines the stack priorities for the tool manager,
  //! and classifies the priority ranges. For now leave room for 100 stacks
  //! in each range.  See init_tool_manager() for use of these ranges.
  enum {
    EVENT_MODIFIERS_E = 0, // Tools to modify incoming events.
    TOOL_MODIFIERS_E = 100, // Tools to manipulate the set of active tools.
    DATA_TOOLS_E = 200, // Tools that handle data,
    SELECTION_TOOL_E = 298,  // Tool that gets pushed used and popped...
    ACTIVE_TOOL_E = 299,  // Tool that gets pushed used and popped...
    VIEWER_TOOLS_E = 300,  // Tools to manipulate the current view.
                           // always on the stack (so last)
    LAST_CHANCE_E = 500
  };

protected:
  void                  init_tool_manager();
  void                  redraw_frame();
  GeomHandle            create_viewer_axes() ;
  void                  draw_visible_scene_graph();
  bool                  item_visible_p(GeomViewerItem* si);
  void                  get_bounds_all(BBox &bbox) { get_bounds(bbox, false); }
  void                  set_state(DrawInfoOpenGL* drawinfo);
  void                  setFrustumToWindowPortion();
  void                  deriveFrustum();
  void                  redraw_obj(MaterialHandle def, GeomHandle obj);
  void                  pick_draw_obj(MaterialHandle def, GeomHandle obj);
  void                  dump_image(const std::string&, const std::string& type = "raw");
  void                  put_scanline(int, int, Color* scanline, int repeat=1);
  void                  StartMpeg(const std::string& fname);
  void                  AddMpegFrame();
  void                  EndMpeg();
  void                  real_get_pick(int, int, GeomHandle&,
                                      GeomPickHandle&, int&);

  void                  render_and_save_image();

  void                  render_rotation_axis(const View &view,
                                             bool do_stereo, int i,
                                             const Vector &eyesep);


  // Private Member Variables
  bool                  do_autoview_;
  bool                  doing_image_p_;
  bool                  doing_movie_p_;
  bool                  make_MPEG_p_;
  std::string           movie_frame_extension_;  // Currently "png" or "ppm".
  int                   current_movie_frame_;
  std::string           movie_name_;
  DrawInfoOpenGL*       drawinfo_;
  WallClockTimer        fps_timer_;
  Frustum               frustum_;
  HiRes                 hi_res_;
  bool                  do_hi_res_;
  bool                  encoding_mpeg_;
  int                   max_gl_lights_;
  int                   animate_num_frames_;
  double                animate_time_begin_;
  double                animate_time_end_;
  double                animate_framerate_;
  double                znear_;
  double                zfar_;
  double                current_time_;
  unsigned int          frame_count_;
  std::vector<float>    depth_buffer_;
  std::vector<unsigned char> fbpick_image_;
  GLdouble              modelview_matrix_[16];
  GLdouble              projection_matrix_[16];
  GLint                 viewport_matrix_[4];
  View                  cached_view_;

  // Mouse Picking variables
  int                   send_pick_x_;
  int                   send_pick_y_;
  int                   ret_pick_index_;
  GeomHandle            ret_pick_obj_;
  GeomPickHandle        ret_pick_pick_;

#ifdef HAVE_MPEG
  FILE *                mpeg_file_;
  MPEGe_options         mpeg_options_;
#endif // HAVE_MPEG

  Lighting               lighting_;
  View                   homeview_;
  View                   view_;
  bool                   do_stereo_;
  float                  ambient_scale_;	
  float                  diffuse_scale_;	
  float                  specular_scale_;	
  float                  shininess_scale_;	
  float                  emission_scale_;	
  float                  line_width_;	
  float                  point_size_;	
  float                  polygon_offset_factor_;
  float                  polygon_offset_units_;
  double                 eye_sep_base_; //! used in stereo rendering.
  bool                   ortho_view_;
  bool                   fogusebg_;
  Color                  fogcolor_;
  double                 fog_start_;
  double                 fog_end_;
  bool                   fog_visibleonly_;
  bool                   fbpick_;
  std::vector<GeomHandle>	 internal_objs_;
  std::vector<bool>		 internal_objs_visible_p_;
  GeomSphere            *focus_sphere_;
  GeomIndexedGroup      *scene_graph_;
  std::map<std::string, bool>	 visible_;
  std::map<std::string, int>	 obj_tag_;
  MaterialHandle         default_material_;
  draw_type_e            draw_type_;
  bool                   capture_z_data_;
  bool                   need_redraw_;
  sel_set_t              selection_set_;
  bool                   selection_set_visible_;
  GeomHandle             selection_geom_;
  tool_handle_t          selection_set_tool_;
  tool_handle_t          scene_graph_tool_;
};

} // End namespace SCIRun

#endif // OpenGLViewer_h
