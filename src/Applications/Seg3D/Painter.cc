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
//    File   : Painter.cc
//    Author : McKay Davis
//    Date   : Nov 2005

#include <sci_comp_warn_fixes.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <math.h>
#include <map>
#include <typeinfo>
#include <iostream>
#include <sci_gl.h>
#include <algorithm>
#include <Core/Containers/Array3.h>
#include <Core/Datatypes/Field.h> 
#include <Core/Exceptions/GuiException.h>
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/ColorMappedNrrdTextureObj.h>
#include <Core/Geom/GeomSwitch.h>
#include <Core/Skinner/GeomSkinnerVarSwitch.h>
#include <Core/Geom/GeomCull.h>
#include <Core/Geom/GeomGroup.h>
#include <Core/Geom/TexSquare.h>
#include <Core/Geom/OpenGLViewport.h>
#include <Core/Geom/FreeType.h>

#include <Core/Math/MiscMath.h>
#include <Core/Math/MinMax.h>
#include <Core/Thread/CleanupManager.h>
#include <Core/Thread/Runnable.h>
#include <Core/Thread/Mutex.h>
#include <Core/Util/Environment.h>
#include <Core/Volume/CM2Widget.h>
#include <Core/Geom/TextRenderer.h>
#include <Core/Geom/FontManager.h>
#include <Core/Geom/GeomColorMappedNrrdTextureObj.h>
#include <Core/Skinner/Variables.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/SceneGraphEvent.h>
#include <Core/Util/FileUtils.h>
#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/PointerToolSelectorTool.h>
#include <Applications/Seg3D/KeyToolSelectorTool.h>
#include <Applications/Seg3D/Seg3DFrame.h>
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>
#include <itkGDCMSeriesFileNames.h>

#include <Applications/Seg3D/GuiCode/cursorinformation.h>
#include <Applications/Seg3D/GuiCode/seg3devents.h>

// needed for Painter::maker
#if defined(_WIN32) && !defined(BUILD_SCIRUN_STATIC)
#  define SCISHARE __declspec(dllexport)
#else
#  define SCISHARE
#endif


namespace SCIRun {

Seg3DFrame* Painter::global_seg3dframe_pointer_ = NULL;
static Painter *global_painter_pointer = NULL;


Painter::Painter(Skinner::Variables *variables, VarContext* ctx) :
  Parent(variables),
  cur_window_(0),
  tm_("Painter"),
  pointer_pos_(),
  windows_(),
  volumes_(),
  current_volume_(0),
  mask_volume_(0),
  current_vrender_target_(0),
  current_vrender_target_deferred_(false),
  volume_lock_("Volume")
{
  undo_manager_ = new UndoManager();

  tm_.add_tool(new PointerToolSelectorTool(this), 50);
  tm_.add_tool(new SliceWindowKeyToolSelectorTool(this), 51);
  global_key_tool_ = new GlobalKeyToolSelectorTool(this);

  event_handle_t event;
  InitializeSignalCatcherTargets(event);

  global_painter_pointer = this;
}


Painter::~Painter()
{
}


void
Painter::redraw_all()
{
  for (SliceWindows::iterator i = windows_.begin(); i != windows_.end(); ++i) {
    (*i)->mark_redraw();
  }
}


void
Painter::extract_all_window_slices()
{
  for (SliceWindows::iterator i = windows_.begin(); i != windows_.end(); ++i)
  {
    (*i)->extract_slices();
  }
}



void
Painter::get_data_from_layer_buttons()  
{
  volume_lock_.lock();
  event_handle_t event = 0;
  for (unsigned int i = 0; i < layer_buttons_.size(); ++i) {
    layer_buttons_[i]->update_from_gui(event);
  }
  volume_lock_.unlock();
}


void
Painter::rebuild_layer_buttons()  
{
  get_data_from_layer_buttons();
  unsigned int bpos = 0;  
  volume_lock_.lock();
  for (int i = volumes_.size()-1; i >= 0 ; --i) {
    build_layer_button(bpos, volumes_[i]);
  }
  for (; bpos < layer_buttons_.size(); ++bpos) {
    layer_buttons_[bpos]->visible_ = false;
    layer_buttons_[bpos]->volume_ = 0;
  }
  volume_lock_.unlock();
  EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
}


void
Painter::build_layer_button(unsigned int &bpos, NrrdVolumeHandle &volume)  
{
  // TODO: Allocate the layer buttons dynamically instead of
  // statically declaring 25 of them in the .skin file (when the
  // layer UI moves to wxwidgets).
  if (bpos >= layer_buttons_.size())
  {
    set_status("Volume button overflow.");
    return;
  }

  volume->lock.lock();
  LayerButton *button = layer_buttons_[bpos];
  button->visible_ = true;
  button->layer_name_ = volume->name_;
  button->layer_visible_ = volume->tmp_visible_;
  button->volume_ = volume;
  volume->button_ = button;
  volume->set_geom_switch(button->layer_visible_);
  
  button->indent_ = volume->depth() * 10 + 5;

  if (volume->label_) {
    Color c = volume->get_label_color();
    button->color_ = Skinner::Color(c.r(), c.g(), c.b(), 1.0);
    button->color_is_selectable_ = true;
  }
  else
  {
    button->color_ = Skinner::Color(0.0, 0.0, 0.0, 0.0);
    button->color_is_selectable_ = false;
  }

  if (volume == current_volume_ && volume == mask_volume_)
  {
    button->background_color_ = Skinner::Color(0.6, 0.3, 0.6, 0.75);
  }
  else if (volume == current_volume_)
  {
    button->background_color_ = Skinner::Color(0.3, 0.3, 0.6, 0.75);
  }
  else if (volume == mask_volume_)
  {
    button->background_color_ = Skinner::Color(0.6, 0.3, 0.3, 0.75);
  }
  else
  {
    button->background_color_ = Skinner::Color(0.0, 0.0, 0.0, 0.0);
  }

  button->active_render_target_ = (volume == current_vrender_target_);

  bpos++;  

  volume->lock.unlock();
}


void
Painter::build_volume_list(NrrdVolumes &volumes)
{
  volume_lock_.lock();
  for (unsigned int i = 0; i < volumes_.size(); ++i) {
    volumes.push_back(volumes_[i]);
  }
  volume_lock_.unlock();
}


void
Painter::opacity_down()
{
  if (current_volume_.get_rep()) {
    current_volume_->
      set_opacity( Clamp(current_volume_->get_opacity()-0.05, 0.0, 1.0) );
    redraw_all();
  }
}


void
Painter::opacity_up()
{
  if (current_volume_.get_rep()) {
    current_volume_->
      set_opacity( Clamp(current_volume_->get_opacity()+0.05, 0.0, 1.0) );
    redraw_all();
  }
}


void
Painter::current_layer_down()
{
  if (!current_volume_.get_rep()) return;

  volume_lock_.lock();
  int pos = 0;
  int i;
  for (i = 0; i < (int)layer_buttons_.size(); i++)
  {
    if (layer_buttons_[i]->volume_ == current_volume_)
    {
      pos = i+1;
    }
    if (layer_buttons_[i]->volume_.get_rep() == 0) break;
  }
  if (pos == i) pos = 0;
  current_volume_ = layer_buttons_[pos]->volume_;
  volume_lock_.unlock();

  rebuild_layer_buttons();
}


void
Painter::current_layer_up()
{
  if (!current_volume_.get_rep()) return;

  volume_lock_.lock();
  int pos = 0;
  int i;
  for (i = 0; i < (int)layer_buttons_.size(); i++)
  {
    if (layer_buttons_[i]->volume_ == current_volume_)
    {
      pos = i-1;
    }
    if (layer_buttons_[i]->volume_.get_rep() == 0) break;
  }
  if (pos == -1) pos = i-1;
  current_volume_ = layer_buttons_[pos]->volume_;
  volume_lock_.unlock();

  rebuild_layer_buttons();
}


void
Painter::move_layer_down()
{
  if (!current_volume_.get_rep()) return;

  bool redraw = false;
  volume_lock_.lock();
  for (int i = 1; i < (int)volumes_.size(); i++)
  {
    if (volumes_[i] == current_volume_)
    {
      NrrdVolumeHandle tmp = volumes_[i];
      volumes_[i] = volumes_[i-1];
      volumes_[i-1] = tmp;
      redraw = true;
      break;
    }
  }
  volume_lock_.unlock();

  if (redraw)
  {
    rebuild_layer_buttons();
    extract_all_window_slices();
    redraw_all();
  }
}


void
Painter::move_layer_up()
{
  if (!current_volume_.get_rep()) return;

  bool redraw = false;
  volume_lock_.lock();
  for (int i = 0; i < (int)volumes_.size()-1; i++)
  {
    if (volumes_[i] == current_volume_)
    {
      NrrdVolumeHandle tmp = volumes_[i];
      volumes_[i] = volumes_[i+1];
      volumes_[i+1] = tmp;
      redraw = true;
      break;
    }
  }
  volume_lock_.unlock();
  
  if (redraw)
  {
    rebuild_layer_buttons();
    extract_all_window_slices();
    redraw_all();
  }
}


void
Painter::set_probe()
{
  for (SliceWindows::iterator i = windows_.begin(); i != windows_.end(); ++i) {
    (*i)->set_probe();
  }
  redraw_all();
}


void
Painter::set_all_slices_tex_dirty()
{
  volume_lock_.lock();
  for (NrrdVolumes::iterator iter = volumes_.begin(); 
       iter != volumes_.end(); ++iter) {
    (*iter)->set_slices_dirty();
  }
  volume_lock_.unlock();
}


NrrdVolumeHandle
Painter::find_volume_by_name(const string &name, NrrdVolumeHandle parent)
{
  NrrdVolumeHandle retval = 0;
  volume_lock_.lock();
  for (size_t i = 0; i < volumes_.size(); ++i) {
    if (volumes_[i]->name_ == name) {
      retval = volumes_[i]; break;
    }
  }
  volume_lock_.unlock();
  return retval;
}


string
Painter::unique_layer_name(string base)
{
  if (!find_volume_by_name(base).get_rep())
  {
    return base;
  }
  string::size_type pos = base.find_last_not_of(" 0123456789");
  base = base.substr(0, pos+1);
  int i = 0;
  string name = base;
  while (find_volume_by_name(name).get_rep())
    name = base + " "+to_string(++i);
  return name;
}


string
Painter::unique_layer_name_from_filename(string filename)
{
  string tmp = changeExtension(filename, "");
  string base = tmp.substr(0, tmp.size()-1);
 
  int i = 1;
  string name = base;
  while (find_volume_by_name(name).get_rep())
    name = base + "(" + to_string(++i) + ")";
  return name;
}


Skinner::Drawable *
Painter::maker(Skinner::Variables *vars) 
{
  return new Painter(vars, 0);
}


NrrdVolumeHandle
Painter::make_layer(string name, NrrdDataHandle &nrrdh, unsigned int mask)
{
  volume_lock_.lock();
  NrrdVolume *vol = new NrrdVolume(this, unique_layer_name(name), nrrdh, mask);
  volumes_.push_back(vol);
  current_volume_ = vol;
  rebuild_layer_buttons();
  volume_lock_.unlock();  
  return vol;
}


NrrdVolumeHandle
Painter::copy_current_layer(string suff)
{
  if (!current_volume_.get_rep()) return 0;
  NrrdDataHandle nrrdh = current_volume_->nrrd_handle_;
  nrrdh.detach(); // Copies the layer memory to the new layer
  return make_layer(current_volume_->name_+suff, nrrdh, current_volume_->label_);
}


void
Painter::show_visible_item(const string &id, const string &group)
{
  // Find the top of the drawable tree.
  Drawable *ptmp = this;
  while (ptmp->get_parent()) ptmp = ptmp->get_parent();

  // Grab all the callbacks with the name we care about in the tree.
  Skinner::Callbacks_t callbacks;
  Parent *parent = dynamic_cast<Parent *>(ptmp);
  if (parent)
    parent->find_callbacks(callbacks, "VisibilityGroup::show_VisibleItem");

  // Set the callback's variables.
  Skinner::Callbacks_t::iterator itr = callbacks.begin();
  while (itr != callbacks.end())
  {
    (*itr)->variables_->insert("id", id);
    (*itr)->variables_->insert("group", group);

    ++itr;
  }

  // Call all the callbacks.  Only the ones with the correct variables
  // will do anything, the others will just return.
  Skinner::Signal *s =
    new Skinner::Signal("VisibilityGroup::show_VisibleItem", this, NULL);
  event_handle_t event(s);

  SignalThrower::throw_signal(callbacks, event);
}


vector<string>
Painter::get_filename_series(string filename)
{
  filename = substituteTilde(filename);
  convertToUnixPath(filename);

  pair<string,string> dfile = split_filename(filename);
  if (validDir(filename)) {
    dfile.first = filename;
    dfile.second = "";
  }

  vector<string> files;
  try {
    typedef itk::GDCMSeriesFileNames series_t;
    series_t::Pointer names = series_t::New();
    names->SetUseSeriesDetails(true);
    names->SetInputDirectory(dfile.first);
    vector<string> series = names->GetSeriesUIDs();
    bool found = false;
    int selected = 0;
    for (size_t i = 0; i < series.size(); i++)
    {
      files = names->GetFileNames(series[i]);
      if (dfile.second == "")
      {
        // If the user picked the directory, open the first available series.
        found = true;
        break;
      }
      for (size_t j = 0; j < files.size(); j++)
      {
        // If the user picked a file, open the series that contains
        // that particular file.
        pair<string, string> dfile2 = split_filename(files[j]);
        if (dfile2.second == dfile.second)
        {
          found = true;
          break;
        }
      }
      if (found) { selected = i; break; }
    }

    if (!found) { files.clear(); }

#if 0
    // TODO:  Fix dicom series selection so the user doesn't have to guess
    // from the license plate style names which one they really want to open.
    else if (series.size() > 1)
    {
      wxString *choices = new wxString[series.size()];
      for (size_t i = 0; i < series.size(); i++)
      {
        choices[i] = std2wx(series[i]);
      }

      wxSingleChoiceDialog dialog(global_seg3dframe_pointer_,
                                  _T("There are multiple series dicoms in this directory.\nPlease select the one you would like to open."),
                                  _T("Dicom series selection"),
                                  series.size(), choices);
      dialog.SetSelection(selected);

      int newselected = -1;
      if (dialog.ShowModal() == wxID_OK)
      {
        newselected = dialog.GetSelection();
      }
      else
      {
        files.clear();
      }
      if (newselected != -1 && selected != newselected)
      {
        files = names->GetFileNames(series[newselected]);
      }
    }
#endif
  }
  catch (...)
  {
  }

  if (files.empty())
  {
    files = GetFilenamesInSequence(dfile.first, dfile.second);
    sort (files.begin(), files.end());
    vector<string> files2;
    for (vector<string>::size_type i = 0; i < files.size(); ++i)
    {
      const string filename = dfile.first + "/" + files[i];
      if (validFile(filename)) {
        files2.push_back(filename);
      }
    }
    files = files2;
  }

  return files;
}


BaseTool::propagation_state_e
Painter::process_event(event_handle_t &event)
{
  KeyEvent *ke = dynamic_cast<KeyEvent *>(event.get_rep());
  if (ke)
  {
    unsigned int s = ke->get_key_state();

    if (s & KeyEvent::KEY_PRESS_E)
    {
      if (global_key_tool_->key_press(ke->get_key_string(),
                                      ke->get_keyval(),
                                      ke->get_modifiers(),
                                      ke->get_time()) == STOP_E)
      {
        return STOP_E;
      }
    }
  }

  ThrowSkinnerSignalEvent *tsse =
    dynamic_cast<ThrowSkinnerSignalEvent *>(event.get_rep());
  if (tsse)
  {
    tsse->set_vars(get_vars());
    throw_signal(tsse->get_name());
#if defined(__WX_GTK__) || defined(_WIN32)
    // Finish the event for threaded synchronous throw.
    tsse->up();
#endif
    return STOP_E;
  }

  return Parent::process_event(event);
}


void
Painter::ThrowSkinnerSignal(const string &name, bool sync)
{
  ThrowSkinnerSignalEvent *tsse = new ThrowSkinnerSignalEvent(name);
  ThrowSkinnerSignal(tsse, sync);
}


void
Painter::ThrowSkinnerSignal(ThrowSkinnerSignalEvent *tsse, bool sync)
{
  if (sync)
  {
#if defined(__WX_GTK__)
    event_handle_t event(tsse);
    EventManager::add_event(event);
    // Wait for the event to finish.
    wxMutexGuiLeave();
    tsse->down();
    wxMutexGuiEnter();
#elif defined(_WIN32)
    event_handle_t event(tsse);
    EventManager::add_event(event);
    // Wait for the event to finish.
    tsse->down();
#else
    // Throw in this thread, backwards compatable (thread unsafe)
    tsse->set_vars(global_painter_pointer->get_vars());
    global_painter_pointer->throw_signal(tsse->get_name());
#endif
  }
  else
  {
    // Asynchronous throw.
    event_handle_t event(tsse);
    EventManager::add_event(event);
  }
}


string
Painter::get_current_layer_name()
{
  if (global_painter_pointer->current_volume_.get_rep())
  {
    return global_painter_pointer->current_volume_->name_;
  }
  return "";
}


void
Painter::update_progress(int progress)
{
  wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
  if (panel)
  {
    wxCommandEvent wxevent(wxEVT_SET_PROGRESS, wxID_ANY);
    wxevent.SetInt(progress);

    wxPostEvent(panel, wxevent);
  }
}


void
Painter::start_progress()
{
  update_progress(-1);
}


void
Painter::finish_progress()
{
  update_progress(101);
}


NrrdVolumeHandle
Painter::create_new_label(NrrdVolumeHandle &likethis, const string &name)
{
  NrrdVolumeHandle result;
  volume_lock_.lock();

  // See if it fits in the current label.
  if (likethis->label_)
  {
    result = likethis->create_child_label_volume();
    if (result.get_rep() && name != "")
    {
      result->name_ = unique_layer_name(name);
    }
  }
  
  // Search for a space amongst the current volumes.
  if (!result.get_rep())
  {
    for (size_t i = 0; i < volumes_.size(); i++)
    {
      if (volumes_[i]->label_ &&
          VolumeOps::compare_nrrd_info(likethis->nrrd_handle_,
                                       volumes_[i]->nrrd_handle_))
      {
        result = volumes_[i]->create_child_label_volume();
      }
      if (result.get_rep())
      {
        if (name != "")
        {
          result->name_ = unique_layer_name(name);
        }
        break;
      }
    }
  }

  // If we reused space make sure it is clear.
  if (result.get_rep())
  {
    VolumeOps::bit_clear(result->nrrd_handle_, result->label_);
  }

  // If no space was found just make a new one.
  if (!result.get_rep())
  {
    NrrdDataHandle nrrdh = 
      VolumeOps::create_clear_nrrd(likethis->nrrd_handle_, LabelNrrdType);
    string newname = unique_layer_name(name);
    if (newname == "")
    {
      newname = unique_layer_name(likethis->name_ + " Label 1");
    }
    result = new NrrdVolume(this, newname, nrrdh, 1);
  }
  volumes_.push_back(result);
  current_volume_ = result;

  volume_lock_.unlock();
  
  return result;
}


void
Painter::set_pointer_position(const Point &p)
{
  pointer_pos_ = p;

  CursorInformationStruct *CI = new CursorInformationStruct();
  CI->x = pointer_pos_.x();
  CI->y = pointer_pos_.y();
  CI->z = pointer_pos_.z();
  CI->value = 0.0;
  CI->valuer = 0.0;

  if (current_volume_.get_rep() && current_volume_->inside_p(pointer_pos_))
  {
    float val = 0.0;
    const std::vector<int> index = current_volume_->world_to_index(pointer_pos_);
    current_volume_->get_value(index, val);
    CI->value = val;
    CI->valuer = 255.0 * (val - current_volume_->clut_min_) /
      (current_volume_->clut_max_ - current_volume_->clut_min_) ;

    CI->xi = index[1];
    CI->yi = index[2];
    CI->zi = index[3];
    CI->value_is_valid = true;
  }
  else
  {
    CI->value_is_valid = false;
  }
	
  wxCommandEvent event(wxEVT_CURSOR_INFORMATION_CHANGE, wxID_ANY);
  event.SetClientData((void *)CI);
  wxPostEvent(global_seg3dframe_pointer_->cursorInformation_, event);
}


void
Painter::set_session_appearance_vars()
{
  // Set the active volume index.
  int avi = -1;
  for (int i = 0; i < (int)volumes_.size(); i++)
  {
    if (current_volume_ == volumes_[i]) avi = i;
  }
  get_vars()->insert("Painter::active_volume_index", to_string(avi), "int");
  
  // Grab the point from a random slice window, they should all be the same.
  // TODO:  This doesn't work, fixme.
  const Point &probe = windows_[0]->center_;
  get_vars()->insert("Painter::probe::x", to_string(probe.x()), "double");
  get_vars()->insert("Painter::probe::y", to_string(probe.y()), "double");
  get_vars()->insert("Painter::probe::z", to_string(probe.z()), "double");
}


void
Painter::get_session_appearance_vars(Skinner::Variables *vars)
{
  if (vars->exists("Painter::volume_visible"))
  {
    const bool vis = vars->get_bool("Painter::volume_visible");
    get_vars()->set_by_string("Painter::volume_visible", to_string(vis));
  }

  if (vars->exists("Painter::active_volume_index"))
  {
    const int avi = vars->get_int("Painter::active_volume_index");
    if (avi >= 0 && avi < (int)volumes_.size())
    {
      current_volume_ = volumes_[avi];
    }
  }

  if (vars->exists("Painter::probe::x") &&
      vars->exists("Painter::probe::y") &&
      vars->exists("Painter::probe::z"))
  {
    const double px = vars->get_double("Painter::probe::x");
    const double py = vars->get_double("Painter::probe::y");
    const double pz = vars->get_double("Painter::probe::z");
    const Point probe(px, py, pz);
    for (unsigned int i = 0; i < windows_.size(); i++)
    {
      windows_[i]->center_ = probe;
    }
  }
  else
  {
    if (current_volume_.get_rep())
    {
      for (unsigned int i = 0; i < windows_.size(); i++)
      {
        windows_[i]->center_ = current_volume_->center();
      }
    }
  }
}


void
Painter::toggle_current_volume_visibility()
{
  if (current_volume_.get_rep() && current_volume_->button_)
  {
    const bool vis = current_volume_->visible();
    if (vis)
    {
      set_status("Toggling current volume visibility off.");
    }
    else
    {
      set_status("Toggling current volume visibility on.");
    }
    current_volume_->button_->layer_visible_ = !vis;
  }
  redraw_all();    
}


double
Painter::scene_scale()
{
  if (volumes_.empty()) return 1.0;
  
  double mscale = AIR_POS_INF;
  for (size_t i = 0; i < volumes_.size(); i++)
  {
    const Vector vscale = volumes_[i]->scale();
    mscale = Min(mscale, vscale.x());
    mscale = Min(mscale, vscale.y());
    mscale = Min(mscale, vscale.z());
  }
  return mscale;
}


} // end namespace SCIRun
