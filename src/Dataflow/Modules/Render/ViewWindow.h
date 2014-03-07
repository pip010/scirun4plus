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
 *  ViewWindow.h: The Geometry Viewer Window
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 */

#ifndef SCI_project_module_ViewWindow_h
#define SCI_project_module_ViewWindow_h

#include <Dataflow/Modules/Render/BallAux.h>
#include <Dataflow/Comm/MessageBase.h>
#include <Dataflow/GuiInterface/GuiCallback.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Core/Containers/Array1.h>
#include <Core/Datatypes/Color.h>
#include <Core/Geometry/BBox.h>
#include <Core/Geometry/Transform.h>
#include <Core/Thread/FutureValue.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/GeomSphere.h>
#include <Core/Geom/GeomCylinder.h>
#include <Core/Geom/GeomPick.h>
#include <Core/Geom/TextRenderer.h>

#include <Dataflow/GuiInterface/GuiGeom.h>
#include <Dataflow/GuiInterface/GuiView.h>
#include <Core/Geom/View.h>
#include <map>

// TODO: REMOVE THIS AND REPLACE THIS WITH THE CLASS FROM CORE/THREAD
// THIS SHOULD NOT HAVE BEEN REINVENTED

// define 'the_time()' function for UniCam
#ifdef _WIN32
#  include <windows.h>
#  include <winbase.h>
   inline double the_time() {
      return double(GetTickCount())/1000.0;
   }
#else
#  include <sys/time.h>
   inline double the_time() {
      struct timeval ts; struct timezone tz;
      gettimeofday(&ts, &tz);
      return (double)(ts.tv_sec + (double)ts.tv_usec/1e6);
   }
#endif



namespace SCIRun {

struct DrawInfoOpenGL;
struct GeometryData;

class GeomObj;
class GeomSphere;
class Light;
class Transform;
class OpenGL;
class ViewScene;
class GeomViewerItem;
class BallData;
class OpenGL;
class ViewWindow;
class ViewWindowClipFrame;

typedef void (ViewWindow::*MouseHandler)(int, int x, int y, 
				  int state, int btn, int time);
typedef void (OpenGL::*ViewWindowVisPMF)(ViewScene*, ViewWindow*, GeomHandle);


// This class contains the variables that record state of the geometry that
// is rendered
class ViewWindowState {
  public:
    ViewWindowState(GuiContext* ctx, const std::string& tclID);
    void request();

    int         use_global()  { return(use_global_.get()); }
    std::string type()        { return(type_.get()); }
    int         show_bbox()   { return(show_bbox_.get()); }
    int         movie_frame() { return(movie_frame_.get()); }
    std::string movie_name()  { return(movie_name_.get()); }
    int         movie()       { return(movie_.get()); }
    int         sync()        { return(sync_.get()); }
    int         clip()        { return(clip_.get()); }
    int         cull()        { return(cull_.get()); }
    int         dl()          { return(dl_.get()); }
    int         fog()         { return(fog_.get()); }
    int         lighting()    { return(lighting_.get()); }

  private:
    GuiInt     use_global_;
    GuiString  type_;
    GuiInt     show_bbox_;
    GuiInt     movie_frame_;
    GuiString  movie_name_;
    GuiInt     movie_;
    GuiInt     sync_;
    GuiInt     clip_;
    GuiInt     cull_;
    GuiInt     dl_;
    GuiInt     fog_;
    GuiInt     lighting_;
};


class ViewWindow : public GuiCallback {
    friend class ViewScene;
  public:
    ViewWindow(ViewScene *s, GuiContext* ctx);
    ~ViewWindow();

    void			itemAdded(GeomViewerItem*);
    void			itemDeleted(GeomViewerItem*);
    void			itemRenamed(GeomViewerItem*, const std::string& newname);
    void			redraw_if_needed();
    // Mouse Callbacks
    void			mouse_dolly(int, int, int, int, int, int);
    void			mouse_translate(int, int, int, int, int, int);
    void			mouse_scale(int, int, int, int, int, int);
    void			mouse_unicam(int, int, int, int, int, int);
    void			mouse_rotate(int, int, int, int, int, int);
    void			mouse_rotate_eyep(int, int, int, int, int, int);
    void			mouse_pick(int, int, int, int, int, int);
    void			tcl_command(GuiArgs&, void*);
    void			get_bounds(BBox&);
    void			get_bounds_all(BBox&);
    void			autoview(const BBox&);
    void			scaled_autoview(const BBox&);
    // sets up the state (OGL) for a tool/viewwindow
    void      requestState(const std::string& tclID);
    void			setState(DrawInfoOpenGL*, const std::string& tclID);
    void			requestDI(const std::string& tclID); // setup DI for drawinfo
    void			setDI(DrawInfoOpenGL*,const std::string& tclID); // setup DI for drawinfo

    void      requestClip(); // force clip settings to be cached in the ViewWindow
    void			setClip(DrawInfoOpenGL*); // setup OGL clipping planes
    void			setMouse(DrawInfoOpenGL*); // setup mouse state
    void			do_for_visible(OpenGL*, ViewWindowVisPMF);
    void			set_current_time(double time);
    void			setView(View view);
    void			getData(int mask, FutureValue<GeometryData*>* result);

    void      requestVisible();

    // Sticky objects
    GeomHandle		createGenAxes();   
    GeomHandle    createScaleBar();
    
    // Request the GUI parameters to be cached for the ScaleBar
    void          requestScaleBar();
    
    ////////////////
    // Set movie recording values.
    //   - stops recording a movie
    void                  setMovieStopped();
    //   - set the frame at which to start numbering output files
    void                  setMovieFrame( int movieframe, const std::string& tclID );
    //   - Places a message on the movie dialog.  If there is an error, stops movie recording.
    void                  setMovieMessage( const std::string & message, bool error = false );

    // This is called by the Viewer when the module has executed.  If a
    // movie is being recorded, then a frame is created based on the
    // parameters of the movie GUI.
    void                  maybeSaveMovieFrame();

    // UNICAM START
    void			unicam_choose(int X, int Y);
    void			unicam_rot(int X, int Y);
    void			unicam_zoom(int X, int Y);
    void			unicam_pan(int X, int Y);
    void			ShowFocusSphere();
    void			HideFocusSphere();
    void			MyTranslateCamera(Vector offset);
    void			MyRotateCamera(Point center,Vector axis, double angle);
    Vector		CameraToWorld(Vector v);
    void			NormalizeMouseXY( int X, int Y, float *NX, float *NY);
    void			UnNormalizeMouseXY(float NX, float NY, int *X, int *Y);
    float			WindowAspect();
    // for 'film_dir' and 'film_pt', x & y should be in the range [-1, 1].
    Vector		film_dir(double x, double y); 
    Point			film_pt(double x, double y, double z=1.0);
    // UNICAM END

    // Functions for accessing gui and ctx so we can disable them
    GuiContext* get_ctx() { return ctx_; }
    
    // Public Variables, (public for OpenGL class access)
    std::string id_;

    BallData *	ball_;		// this is the ball for arc ball stuff
    int         loop_count_; // number of frames in a loop rotation
    View			  rot_view_;	// pre-rotation view
    Transform		prev_trans_;
    double		  eye_dist_;
    GuiView		  gui_view_;
    GuiView		  gui_homeview_;
    GuiInt	  	gui_sr_;
    GuiInt		  gui_do_stereo_;		// Stereo
    GuiInt		  gui_ortho_view_;
    GuiInt		  gui_track_view_window_0_;
    GuiInt      gui_lock_view_window_;
    GuiInt		  gui_raxes_;

    GuiDouble		gui_ambient_scale_;	// Scene material scales
    GuiDouble		gui_diffuse_scale_;
    GuiDouble		gui_specular_scale_;
    GuiDouble		gui_emission_scale_;
    GuiDouble		gui_shininess_scale_;
    GuiDouble		gui_polygon_offset_factor_;
    GuiDouble		gui_polygon_offset_units_;
    GuiDouble		gui_point_size_;
    GuiDouble		gui_line_width_;
    GuiDouble		gui_text_offset_;
    GuiDouble		gui_sbase_;
    GuiColor		gui_bgcolor_;		// Background Color
    GuiInt      gui_fogusebg_;
    GuiColor    gui_fogcolor_;
    GuiDouble   gui_fog_start_;
    GuiDouble   gui_fog_end_;
    GuiInt      gui_fog_visibleonly_;
    GuiInt      gui_inertia_loop_count_;
    GuiDouble		gui_inertia_x_;
    GuiDouble		gui_inertia_y_;
    GuiInt      gui_inertia_recalculate_;
    GuiInt      gui_inertia_mode_;
    GuiInt		  gui_current_time_;
    GuiString		gui_currentvisual_;
    GuiInt		  gui_autoav_;
    GuiInt		  gui_caxes_;
    GuiString		gui_pos_;

    // For ScaleBar
    GuiInt      gui_scalebar_;          // on/off
    GuiString   gui_scalebar_unit_;     // unit of the space dimensions
    GuiDouble   gui_scalebar_length_;   // length of the bar
    GuiDouble   gui_scalebar_height_;   // height of the bar
    GuiDouble   gui_scalebar_multiplier_;   // multiplier of the bar
    GuiColor    gui_scalebar_color_;    // color of the scalarbar
    GuiInt      gui_scalebar_nticks_;   // number of ticks along the bar
    GuiInt      gui_scalebar_fontsize_;   // size of label font
    GuiDouble   gui_scalebar_linesize_;   // size of label font
    
    std::vector<ViewWindowClipFrame*> viewwindow_clip_frames_;
    std::vector<bool>          viewwindow_clip_frames_draw_;

    GeomHandle  orientation_axis_;
    GeomHandle  scale_bar_;

  private:
    ViewWindow(const ViewWindow&); // Should not be called
    void			do_mouse(MouseHandler, GuiArgs&);
    void			animate_to_view(const View& v, double time);
    void			redraw();
    void			redraw(double tbeg, double tend, 
               int nframes, double framerate);

    // Private Member variables
    ViewScene*            viewer_;
    OpenGL*               renderer_;
    GuiContext*           ctx_;
    std::string                tclID_;
    std::map<std::string,GuiInt*>   visible_;   // Which of the objects do we draw?
    std::map<std::string,int>       obj_tag_;
    bool                  need_redraw_;
    int                   pick_n_;
    int                   maxtag_;
    int                   last_x_;
    int                   last_y_;
    int                   pick_x_;
    int                   pick_y_;
    int                   last_time_;
    bool                  mouse_action_;
    double                total_x_;
    double                total_y_;
    double                total_z_;
    double                total_scale_;
    Vector                translate_u_;
    Vector                translate_v_;
    View                  start_translate_view_;
    View                  homeview_;
    int                   rotate_valid_p_;
    GeomPickHandle        pick_pick_;
    GeomHandle            pick_obj_;
    std::vector<GeomHandle>	viewwindow_objs_;
    std::vector<bool>     viewwindow_objs_draw_;   

    // Variables for mouse_dolly method (dollying into/outof a view)
    double		dolly_total_;
    Vector		dolly_vector_;
    double		dolly_throttle_;
    double		dolly_throttle_scale_;
    // history for quaternions and time
    int       prev_time_[3]; // These come from tcl and tcl gets them from the X
                             // server timestamp (typically the time since the
                             // last server reset) in milliseconds.
    HVect			prev_quat_[3];
    // UNICAM variables
    enum {UNICAM_CHOOSE = 0, UNICAM_ROT, UNICAM_PAN, UNICAM_ZOOM};
    int         unicam_state_;
    int         down_x_;
    int         down_y_;
    Point       down_pt_;
    double      dtime_;
    double      uni_dist_;
    float       last_pos_[3];
    float       start_pix_[2];
    float       last_pix_[2];
    GeomSphere *focus_sphere_;
    // TCL GUI Variables
    bool				mpick_;
    
    std::vector<GuiInt*>    clip_visible_;
    std::vector<GuiInt*>    clip_frame_;
    std::vector<GuiInt*>    clip_normal_reverse_;
    std::vector<GuiDouble*> clip_normal_x_;
    std::vector<GuiDouble*> clip_normal_y_;
    std::vector<GuiDouble*> clip_normal_z_;
    std::vector<GuiDouble*> clip_normal_d_;
    
    std::map<std::string,ViewWindowState*> state_;
    
    TextRendererHandle      text_renderer_;
};

class ViewWindowMouseMessage : public MessageBase {
  public:
    std::string rid;
    MouseHandler handler;
    int action;
    int x, y;
    int state;
    int btn;
    int time;
    
    ViewWindowMouseMessage(const std::string& rid, MouseHandler handler,
                           int action, int x, int y, int state, int btn,
                           int time);
    virtual ~ViewWindowMouseMessage();
};


class ViewWindowClipFrame 
{
  public:
    friend class ViewWindow;
    friend class OpenGL;

    ViewWindowClipFrame();
    ViewWindowClipFrame(const ViewWindowClipFrame& copy);
    virtual ~ViewWindowClipFrame();
    
    void Set(const Point& c, const Vector& n,
             double width, double height, double scale);  
   
  private:
    Point center_;
    Vector normal_;
    double width_, height_, scale_;

    std::vector<Point>          verts_;
    std::vector<GeomSphere*>    corners_;
    std::vector<GeomCylinder*>  edges_;

    void adjust();
    void set_position(const Point& c, const Vector& n){ center_ = c; normal_ = n;}
    void set_size(double w, double h) { width_ = w;  height_ = h;}
    void set_scale(double s) { scale_ = s;}

};

} // End namespace SCIRun



#endif
