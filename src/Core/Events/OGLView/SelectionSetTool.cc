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
//    File   : SelectionSetTool.cc
//    Author : Martin Cole
//    Date   : Mon Sep 18 10:04:12 2006

#include <Core/Events/EventManager.h>
#include <Core/Events/DataManager.h>
#include <Core/Geom/GeomSwitch.h>
#include <Core/Geom/GeomColorMap.h>
#include <Core/Events/OGLView/SelectionSetTool.h>
#include <Core/Events/OGLView/SelectionTargetEvent.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <Core/Datatypes/TriSurfMesh.h>

namespace SCIRun {

  using std::cerr;
  using std::endl;
  using std::cout;

SelectionSetTool::SelectionSetTool(const std::string& name, SSTInterface *ssti) :
  BaseTool(name),
  mode_(FACES_E),
  sel_fld_(0),
  sel_fld_id_(0),
  ssti_(ssti),
  def_material_(new Material(Color(0.3, 0.7, 0.3))),
  text_material_(new Material(Color(0.25, 0.25, 0.45)))
{
  float ramp[1024];
  float m = 1./1024.;
  for (int i = 0; i < 1024; i += 4) {
    ramp[i] = m * i;
    ramp[i+1] = 0.0;
    ramp[i+2] = 0.0;
    ramp[i+3] = 0.9;
  }
  color_map_ = new ColorMap(ramp);
}


SelectionSetTool::~SelectionSetTool()
{
}


BaseTool::propagation_state_e
SelectionSetTool::process_event(event_handle_t e)
{
  SelectionTargetEvent *st =
    dynamic_cast<SelectionTargetEvent*>(e.get_rep());
  if (st) {
    sel_fld_ = st->get_selection_target();
    sel_fld_id_ = st->get_selection_id();
    return STOP_E;
  }
  return CONTINUE_E;
}


void
SelectionSetTool::delete_faces()
{
  if (! sel_fld_.get_rep()) return;
  sel_fld_->lock.lock();

  // Turn this call into a general algorithm, but for now assume trisurf.
  MeshHandle mb = sel_fld_->mesh();
  VMesh *tsm = sel_fld_->vmesh();
  if (!tsm->is_trisurfmesh()) 
  {
    cerr << "Error:: not a TriSurf in SelectionSetTool::delete_faces"
	 << endl;
  }
  std::vector<unsigned int> &sfaces = ssti_->get_selection_set();

  std::vector<int> faces;
  std::vector<unsigned int>::iterator si = sfaces.begin();
  while (si != sfaces.end()) {
    faces.push_back(*si++);
  }

  bool altered = false;
  // Remove last index first.
  sort(faces.begin(), faces.end());
  std::vector<int>::reverse_iterator iter  = faces.rbegin();
  while (iter != faces.rend()) {
    int face = *iter++;
//    altered |= tsm->remove_face(face);
    cout << "removed face " << face << endl;
  }

  // Clear the selection set.
  sfaces.clear();
  ssti_->set_selection_set_visible(false);
  sel_fld_->lock.unlock();

  // Notify the data manager that this model has changed.
  CommandEvent *c = new CommandEvent();
  c->set_command("selection field modified");
  event_handle_t event = c;
  EventManager::add_event(event);
}


void
SelectionSetTool::add_face()
{
  if (! sel_fld_.get_rep()) return;
  sel_fld_->lock.lock();

  // Turn this call into a general algorithm, but for now assume trisurf.
  MeshHandle mb = sel_fld_->mesh();
  VMesh* tsm = sel_fld_->vmesh();
  if (!tsm) 
  {
    cerr << "Error:: not a TriSurf in SelectionSetTool::delete_faces"
	 << endl;
  }
  
  std::vector<unsigned int> &snodes = ssti_->get_selection_set();
  if (snodes.size() < 3) 
  {
    cerr << "Error: must select 3 nodes to add a face" << endl;
  }

  std::vector<int> nodes;
  std::vector<unsigned int>::iterator si = snodes.begin();
  VMesh::Node::array_type array(3);
  array[0] = *si++;
  array[1] = *si++;
  array[2] = *si++;

  tsm->add_elem(array);

  //clear the selection set.
  snodes.erase(snodes.begin(), si);
  if (snodes.size() == 0) {
    ssti_->set_selection_set_visible(false);
  }
  sel_fld_->lock.unlock();

  //notify the data manager that this model has changed.
  CommandEvent *c = new CommandEvent();
  c->set_command("selection field modified");
  event_handle_t event = c;
  EventManager::add_event(event);
}


void
SelectionSetTool::render_selection_set()
{
  //create a new field with the selection items in it;
  if (! sel_fld_.get_rep()) { return; }
  if (! ssti_->get_selection_set().size()) {
    ssti_->set_selection_geom(GeomHandle(0));
    return;
  }
}


} // namespace SCIRun
