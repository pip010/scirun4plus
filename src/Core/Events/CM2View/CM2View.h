//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//  License for the specific language governing rights and limitations under
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
//    File   : CM2View.h
//    Author : Martin Cole
//    Date   : Wed Sep 19 08:55:30 2007


#if !defined(CM2View_h)
#define CM2View_h

#include <sci_gl.h>
#include <Core/Events/OpenGLBase.h>
#include <slivr/CM2Shader.h>
#include <slivr/ShaderProgramARB.h>
#include <Core/Volume/ColorMap2.h>
#include <Core/Volume/CM2Widget.h>
#include <Core/Events/CM2View/share.h>

namespace SCIRun {
class CM2View;
using SLIVR::deep_copy_widgets;

//! Event to synchronize widget state between threads.
class SCISHARE SetWidgetsEvent : public BaseEvent
{
public:
  SetWidgetsEvent(const std::vector <SLIVR::CM2Widget*> &w,
		  const std::vector <SLIVR::CM2Widget*> &add,
		  const std::string& target = "") :
    BaseEvent(target),
    widgets_(0),
    additional_(0),
    updating_(true)
  {
    deep_copy_widgets(w, widgets_);
    deep_copy_widgets(add, additional_);
  }

  ~SetWidgetsEvent()
  {}
  virtual SetWidgetsEvent *clone() { return new SetWidgetsEvent(*this); }

  // copy of widget set
  std::vector<SLIVR::CM2Widget*>       widgets_;
  std::vector<SLIVR::CM2Widget*>       additional_;
  bool                            updating_;
};

//! Event to synchronize widget state between threads.
class SetHistoEvent : public BaseEvent
{
public:
  SetHistoEvent(Nrrd* h, double opacity, const std::string& target = "") :
    BaseEvent(target),
    histo_(h),
    opacity_(opacity)
  {}
  ~SetHistoEvent()
  {}
  virtual SetHistoEvent *clone() { return new SetHistoEvent(*this); }

  Nrrd*                  histo_;
  double                 opacity_;
};

//! Tools for widget manipulation and event handling.
class SCISHARE WidgetTool : public BaseTool
{
public:
  WidgetTool(const std::string& name, CM2View *cv) :
    BaseTool(name),
    cm2view_(cv)
  {}
  virtual ~WidgetTool() {}

  virtual propagation_state_e process_event(event_handle_t event);
private:
  CM2View                *cm2view_;
};

//! Tools for widget manipulation and event handling.
class SCISHARE HistoTool : public BaseTool
{
public:
  HistoTool(const std::string& name, CM2View *cv) :
    BaseTool(name),
    cm2view_(cv)
  {}
  virtual ~HistoTool() {}

  virtual propagation_state_e process_event(event_handle_t event);
private:
  CM2View                *cm2view_;
};

//! Main widget manipulation opengl view.
class SCISHARE CM2View : public OpenGLBase
{
public:
  CM2View(OpenGLContext *oglc, const std::string& parent_id, const std::string& id = "CM2View");
  virtual ~CM2View();

  virtual void                redraw_frame();

  //! commands
  static const std::string Save_PPM_c;
  static const std::string Add_Triangle_Widget_c;
  static const std::string Add_Rectangle_Widget_c;
  static const std::string Add_Ellipsoid_Widget_c;
  static const std::string Add_Paraboloid_Widget_c;
  static const std::string Del_Selected_Widget_c;

  //! Update the histogram nrrd.
  void                        set_histo(Nrrd* h, double opacity);
  //! Update the widget set.
  void                        set_widgets(const std::vector<SLIVR::CM2Widget*> &w,
					  const std::vector<SLIVR::CM2Widget*> &a);

  //! access the widgets.
  std::vector<SLIVR::CM2Widget*> &widgets() { return widgets_; }

  //! Window width.
  int                          width();
  //! Window height.
  int                          height();

  double pan_x() const { return pan_x_; }
  double pan_y() const { return pan_y_; }
  double scale() const { return scale_; }

  void set_pan_x(double v) { pan_x_ = v; }
  void set_pan_y(double v) { pan_y_ = v; }
  void set_scale(double v) { scale_ = v; }

  //! The event target for parent.
  std::string                     parent_id() const { return parent_id_; }
  void                       set_notify_id(const std::string& id) { parent_id_ = id; }

  void                       set_ppm(const std::string& f) { ppm_ = f; save_ppm_ = true; }
  void                       add_triangle_widget();
  void                       add_rectangle_widget();
  void                       add_ellipsoid_widget();
  void                       add_paraboloid_widget();
  void                       delete_selected_widget();
  void                       widget_changed_notify(bool updating = false);
private:
  void select_widget(int widget, int object);
  void build_histogram_texture();
  void draw_texture(const GLuint &texture_id);
  void clear_widgets();
  //! saves the image, at size == w/4 and h/4
  //! where w,h are relative to the image in buf
  void save_ppm_thumb(const std::string& filename, int w, int h,
		      const unsigned char * buf) const;


  SLIVR::CM2ShaderFactory      *shader_factory_;

  //! Our own copy, no shared data.
  std::vector<SLIVR::CM2Widget*>	widgets_;
  //! Additional widgets are read only.
  std::vector<SLIVR::CM2Widget*>	additional_widgets_;

  Nrrd*				histo_; //! read only.
  bool				histo_dirty_;
  GLuint			histogram_texture_id_;
  std::pair<float,float>		value_range_;
  double                        hdraw_opacity_;

  HistoTool                    *histo_tool_;
  WidgetTool                   *widget_tool_;

  double                        pan_x_;
  double                        pan_y_;
  double                        scale_;
  //! the id of parent, used to target update events.
  std::string                        parent_id_;
  bool                          save_ppm_;
  std::string                        ppm_;

};

} // namespace SCIRun

#endif //CM2View_h
