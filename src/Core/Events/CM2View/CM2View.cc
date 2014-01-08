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
//    File   : CM2View.cc
//    Author : Martin Cole
//    Date   : Wed Sep 19 09:01:35 2007

#include <Core/Events/CM2View/CM2View.h>
#include <Core/Events/CM2View/CM2TranslateTool.h>
#include <Core/Events/CM2View/CM2ScaleTool.h>
#include <Core/Events/CM2View/CM2CursorTool.h>
#include <Core/Events/CM2View/CM2WidgetTransformTool.h>
#include <Core/Events/CM2View/CM2UndoTool.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/Endian.h>
#include <Core/Util/MemoryUtil.h>
#include <Core/Geom/OpenGLContext.h>
#include <fstream>

namespace SCIRun {



// command strings.
const std::string CM2View::Save_PPM_c = "SAVE_IMAGE:";
const std::string CM2View::Add_Triangle_Widget_c = "ADD_TRI:";
const std::string CM2View::Add_Rectangle_Widget_c = "ADD_RECT:";
const std::string CM2View::Add_Ellipsoid_Widget_c = "ADD_ELLI:";
const std::string CM2View::Add_Paraboloid_Widget_c = "ADD_PARA:";
const std::string CM2View::Del_Selected_Widget_c = "DEL_SEL:";



BaseTool::propagation_state_e
WidgetTool::process_event(event_handle_t event)
{
  SetWidgetsEvent *swe = dynamic_cast<SetWidgetsEvent*>(event.get_rep());
  if (! swe) {
    return CONTINUE_E;
  }
  cm2view_->set_widgets(swe->widgets_, swe->additional_);

  return STOP_E;
}

BaseTool::propagation_state_e
HistoTool::process_event(event_handle_t event)
{
  SetHistoEvent *she = dynamic_cast<SetHistoEvent*>(event.get_rep());
  if (! she) {
    return CONTINUE_E;
  }
  cm2view_->set_histo(she->histo_, she->opacity_);

  return STOP_E;
}

//! A CommandTool waits for various user defined commands, and executes
//! the command when the command is issued.
class CM2CommandTool :  public CommandTool
{
public:
  CM2CommandTool(std::string name, CM2View* v) :
    CommandTool(name),
    view_(v)
  {}
  virtual ~CM2CommandTool() {}

  virtual propagation_state_e issue_command(const std::string &cmmd,
					    unsigned int /*time*/)
  {

    size_t idx;
    if ((idx = cmmd.find(CM2View::Save_PPM_c)) != std::string::npos) {
      std::string fname = cmmd.substr(idx + CM2View::Save_PPM_c.size());
      view_->set_ppm(fname);
      return STOP_E;
    } else if ((idx = cmmd.find(CM2View::Add_Triangle_Widget_c)) !=
	       std::string::npos) {
      view_->add_triangle_widget();
      return STOP_E;
    } else if ((idx = cmmd.find(CM2View::Add_Rectangle_Widget_c)) !=
	       std::string::npos) {
      view_->add_rectangle_widget();
      return STOP_E;
    } else if ((idx = cmmd.find(CM2View::Add_Ellipsoid_Widget_c) !=
		std::string::npos)) {
      view_->add_ellipsoid_widget();
      return STOP_E;
    } else if ((idx = cmmd.find(CM2View::Add_Paraboloid_Widget_c)) !=
	       std::string::npos) {
      view_->add_paraboloid_widget();
      return STOP_E;
    } else if ((idx = cmmd.find(CM2View::Del_Selected_Widget_c)) !=
	       std::string::npos) {
      view_->delete_selected_widget();
      return STOP_E;
    }
    return CONTINUE_E;
  }

private:
  CM2View *view_;
};


CM2View::CM2View(OpenGLContext *oglc, const std::string& parent_id, const std::string& id) :
  OpenGLBase(oglc, id),
  shader_factory_(new SLIVR::CM2ShaderFactory()),
  widgets_(),
  additional_widgets_(),
  histo_(0),
  histo_dirty_(false),
  histogram_texture_id_(0),
  value_range_(0.0, -1.0),
  hdraw_opacity_(0.5),
  pan_x_(0.0),
  pan_y_(0.0),
  scale_(1.0),
  parent_id_(parent_id),
  save_ppm_(false),
  ppm_("")
{
  unsigned int view_pri = VIEWER_TOOLS_E;
  tool_handle_t cur = new CM2CursorTool("CM2View Cursor Tool", this);
  tm_.add_tool(cur, view_pri++);

  tool_handle_t wtt = new CM2WidgetTransformTool("CM2View Widget Transf Tool",
						 this);
  tm_.add_tool(wtt, view_pri++);

  tool_handle_t trans = new CM2TranslateTool("CM2View Translate Tool", this);
  tm_.add_tool(trans, view_pri++);

  tool_handle_t sc = new CM2ScaleTool("CM2View Scale Tool", this);
  tm_.add_tool(sc, view_pri++);

  unsigned int data_pri = DATA_TOOLS_E;
  widget_tool_ = new WidgetTool("CM2View WidgetTool", this);
  tm_.add_tool(widget_tool_, data_pri++);

  tm_.add_tool(new CM2CommandTool("CM2View CommandTool", this), data_pri++);

  histo_tool_ = new HistoTool("CM2View HistoTool", this);
  tm_.add_tool(histo_tool_, data_pri++);

  tm_.add_tool(new CM2UndoTool("CM2View UndoTool", this), data_pri++);
}

CM2View::~CM2View()
{
  clear_widgets();
}

void
CM2View::clear_widgets()
{
  delete_all_items(widgets_);
  widgets_.clear();

  delete_all_items(additional_widgets_);
  additional_widgets_.clear();
}


void
CM2View::set_widgets(const std::vector<SLIVR::CM2Widget*> &w,
		     const std::vector<SLIVR::CM2Widget*> &add)
{
  clear_widgets();
  widgets_ = w;
  additional_widgets_ = add;
}

void
CM2View::set_histo(Nrrd *h, double opacity)
{
  // Could just be an opacity update, only build a new texture if the nrrd is
  // different.
  if (h == histo_) {
    hdraw_opacity_ = opacity;
  } else {
    histo_ = h;
    histo_dirty_ = true;
    if (histo_ && histo_->kvp) {
      const char *min = nrrdKeyValueGet(histo_, "jhisto_nrrd0_min");
      const char *max = nrrdKeyValueGet(histo_, "jhisto_nrrd0_max");
      double dmin, dmax;
      if (min && max && string_to_double(min, dmin) &&
	  string_to_double(max, dmax))
      {
	      value_range_ = std::make_pair(float(dmin), float(dmax));
      }
    }
  }
}
void
CM2View::build_histogram_texture()
{
  if (!histo_dirty_) return;
  histo_dirty_ = false;

  if(glIsTexture(histogram_texture_id_)) {
    glDeleteTextures(1, &histogram_texture_id_);
    histogram_texture_id_ = 0;
  }
  if (!histo_) return;

  glGenTextures(1, &histogram_texture_id_);
  glBindTexture(GL_TEXTURE_2D, histogram_texture_id_);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
#ifndef GL_CLAMP_TO_EDGE
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
#else
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
#endif
  size_t axis_size[3];

  nrrdAxisInfoGet_nva(histo_, nrrdAxisInfoSize, axis_size);
  glTexImage2D(GL_TEXTURE_2D, 0, 1, axis_size[histo_->dim-2],
	       axis_size[histo_->dim-1], 0, GL_LUMINANCE, GL_UNSIGNED_BYTE,
	       histo_->data);
  glBindTexture(GL_TEXTURE_2D, 0);
}


int
CM2View::width()
{
  if (!gl_context_) return -1;
  return gl_context_->width();
}

int
CM2View::height()
{
  if (!gl_context_) return -1;
  return gl_context_->height();
}

void
CM2View::draw_texture(const GLuint &texture_id)
{
  if (texture_id <= 0) return;
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, texture_id);

  glBegin(GL_QUADS);
  {
    glTexCoord2f( 0.0,  0.0);
    glVertex2f( 0.0,  0.0);
    glTexCoord2f( 1.0,  0.0);
    glVertex2f( 1.0,  0.0);
    glTexCoord2f( 1.0,  1.0);
    glVertex2f( 1.0,  1.0);
    glTexCoord2f( 0.0,  1.0);
    glVertex2f( 0.0,  1.0);
  }
  glEnd();

  glBindTexture(GL_TEXTURE_2D, 0);
  glDisable(GL_TEXTURE_2D);
}

void
CM2View::redraw_frame()
{
  if (dead_) return;
   
  if (!gl_context_) return;
  gl_context_->lock();
  if (dead_) 
  { 
    gl_context_->unlock();
    return;
  }
  
  if(! gl_context_->make_current())
  {
    std::cerr << "make_current failed." << std::endl;
    gl_context_->unlock();
    return;
  }

  CHECK_OPENGL_ERROR();
  glViewport(0, 0, width(), height());

  glDrawBuffer(GL_BACK);
  CHECK_OPENGL_ERROR();
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);


  build_histogram_texture();
  CHECK_OPENGL_ERROR();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();


  const float scale_factor = 2.0 * scale_;
  glScalef(scale_factor, scale_factor, scale_factor);
  glTranslatef(-0.5 - pan_x_, -0.5 - pan_y_ , 0.0);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glColor4f(hdraw_opacity_, hdraw_opacity_, hdraw_opacity_, 1.0);
  draw_texture(histogram_texture_id_);

  // Draw Colormap
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

  // draw the widgets.
  glEnable(GL_BLEND);
  glDrawBuffer(GL_BACK);
  glReadBuffer(GL_BACK);
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  for(unsigned int i = 0; i < widgets_.size(); i++) {
    widgets_[i]->set_value_range(value_range_);
    widgets_[i]->rasterize(*shader_factory_);
  }
  for(unsigned int i = 0; i < additional_widgets_.size(); i++) {
    additional_widgets_[i]->set_value_range(value_range_);
    additional_widgets_[i]->rasterize(*shader_factory_);
  }
  CHECK_OPENGL_ERROR();
  glDisable(GL_BLEND);
  // Draw Colormap Widget Frames
  for(unsigned int i=0; i<widgets_.size(); i++)
    widgets_[i]->draw();

  CHECK_OPENGL_ERROR();


  if (save_ppm_) {
    save_ppm_ = false;
    unsigned char* img = new unsigned char[width() * height() * 3];
    glReadBuffer(GL_BACK);
    glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, width());
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

    glReadPixels(0, 0, width(), height(),
		 GL_RGB, GL_UNSIGNED_BYTE, img);
    std::string fn = ppm_ + ".ppm";
    save_ppm_thumb(fn, width(), height(), img);
    delete[] img;
  }

  gl_context_->swap();
  // check for errors before the release, as opengl will barf (at least in
  // windows) if you do it after
  CHECK_OPENGL_ERROR();
  gl_context_->unlock();
}

void
CM2View::save_ppm_thumb(const std::string& filename, int sx, int sy,
			const unsigned char * buf) const
{
  int R = 0;
  int G = 1;
  int B = 2;

  //  int A = 0;
  std::ofstream output(filename.c_str(), std::ios::out);
  if (!output) {
    std::cerr << "ERROR: can't open file: " << filename << std::endl;
    return;
  }

  output << "P6\n# CREATOR: " << std::endl << std::endl;
  output << sx / 4 << " " << sy / 4 << std::endl;
  output << 255 << std::endl;

  for (int row = sy - 1; row >= 0; row -= 4) {
    for (int col = 0; col < sx; col += 4) {
      int p = 3 * (row * sx + col);
      output << buf[p + R] << buf[p + G] << buf[p + B];
    }
  }

  output.close();
  std::cerr << filename << " closed" << std::endl;
}

void
CM2View::select_widget(int widget, int object)
{
  for (size_t i = 0; i < widgets_.size(); ++i) {
    widgets_[i]->select(i == (size_t)widget ? object : 0);
  }
}

void
CM2View::widget_changed_notify(bool updating)
{
  static std::vector<SLIVR::CM2Widget*> add;
  SetWidgetsEvent *e = new SetWidgetsEvent(widgets(), add, parent_id());
  e->updating_ = updating;
  event_handle_t event = e;
  EventManager::add_event(event);
}


void
CM2View::add_rectangle_widget()
{
  // record the widget state.
  event_handle_t e = new CM2UndoActionEvent(CM2UndoActionEvent::ADD_E,
					    widgets_.size() - 1,
					    0, id());
  EventManager::add_event(e);
  widgets_.push_back(new RectangleCM2Widget());
  EventManager::add_event(e);
  select_widget(widgets_.size() - 1, 1);
  widget_changed_notify();
}

void
CM2View::add_triangle_widget()
{
  // record the widget state.
  event_handle_t e = new CM2UndoActionEvent(CM2UndoActionEvent::ADD_E,
					    widgets_.size() - 1,
					    0, id());
  EventManager::add_event(e);
  widgets_.push_back(new TriangleCM2Widget());
  EventManager::add_event(e);
  select_widget(widgets_.size() - 1, 1);
  widget_changed_notify();
}

void
CM2View::add_ellipsoid_widget()
{
  // record the widget state.
  event_handle_t e = new CM2UndoActionEvent(CM2UndoActionEvent::ADD_E,
					    widgets_.size() - 1,
					    0, id());
  EventManager::add_event(e);
  widgets_.push_back(new EllipsoidCM2Widget());
  EventManager::add_event(e);
  select_widget(widgets_.size() - 1, 1);
  widget_changed_notify();
}

void
CM2View::add_paraboloid_widget()
{
  // record the widget state.
  event_handle_t e = new CM2UndoActionEvent(CM2UndoActionEvent::ADD_E,
					    widgets_.size() - 1,
					    0, id());
  EventManager::add_event(e);

  widgets_.push_back(new ParaboloidCM2Widget());
  EventManager::add_event(e);
  select_widget(widgets_.size() - 1, 1);
  widget_changed_notify();
}


void
CM2View::delete_selected_widget()
{
  size_t widget = 0;
  for (size_t i = 0; i < widgets_.size(); ++i) {
    if (widgets_[i]->selected()) {
      widget = i;
      break;
    }
  }


  event_handle_t e = new CM2UndoActionEvent(CM2UndoActionEvent::DELETE_E,
					    widget,
					    widgets_[widget]->duplicate(),
					    id());
  EventManager::add_event(e);
  delete widgets_[widget];
  widgets_.erase(widgets_.begin() + widget);
  if (widget >= widgets_.size()) {
    select_widget(widgets_.size() - 1, 1);
  } else {
    select_widget(widget, 1);
  }
  widget_changed_notify();
}

} //namespace SCIRun
