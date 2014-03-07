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
//    File   : LineTool.cc
//    Author : Michael Callahan
//    Date   : March 2008

#include <Applications/Seg3D/Painter.h>

#include <Applications/Seg3D/LineTool.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/SceneGraphEvent.h>
#include <Core/Geom/GeomCylinder.h>
#include <Core/Geom/GeomMaterial.h>
#include <sci_gl.h>

#include <Core/Events/keysyms.h>


namespace SCIRun {

SeedTool::seeds_t LineTool::seed_cache_;

LineTool::LineTool(Painter *painter) :
  SeedTool("LineTool", painter)
{
  seeds_ = seed_cache_;

  // Clear the text on tool startup.
  wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
  if (panel)
  {
    vector<Point> *seedcopy = new vector<Point>(seeds_);
    wxCommandEvent wxevent(wxEVT_MEASUREMENT_UPDATE, wxID_ANY);
    wxevent.SetClientData((void *)seedcopy);
    wxPostEvent(panel, wxevent);
  }
}


LineTool::~LineTool()
{
  // Clear the 3D seed lines..
  event_handle_t scene_event = new SceneGraphEvent(0, "Seedlines");
  EventManager::add_event(scene_event);

  seed_cache_ = seeds_;
}


void
LineTool::run_filter()
{
}


BaseTool::propagation_state_e 
LineTool::process_event(event_handle_t event)
{
  return SeedTool::process_event(event);
}


void
LineTool::seed_change_callback()
{
  SeedTool::seed_change_callback();

  if (seeds_.size() > 1)
  {
    GeomCylinders *lines = new GeomCylinders();
    lines->set_radius(2.0 * painter_->scene_scale());
    for (unsigned int i = 0; i < seeds_.size()-1; i++)
    {
      lines->add(seeds_[i], seeds_[i+1]);
    }             

    MaterialHandle green = new Material(Color(0.0, 0.8, 0.0));
    GeomHandle seedgeom = new GeomMaterial(lines, green);

    event_handle_t scene_event = new SceneGraphEvent(seedgeom, "Seedlines");
    EventManager::add_event(scene_event);
  }
  wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
  if (panel)
  {
    vector<Point> *seedcopy = new vector<Point>(seeds_);
    wxCommandEvent wxevent(wxEVT_MEASUREMENT_UPDATE, wxID_ANY);
    wxevent.SetClientData((void *)seedcopy);
    wxPostEvent(panel, wxevent);
  }
}


void
LineTool::draw_gl(SliceWindow &window)
{
  if (seeds_.size() > 1)
  {
    glColor4d(0.0, 0.8, 0.0, 1.0);

    glLineWidth(2.0);

    glBegin(GL_LINES);
    for (unsigned int i = 0; i < seeds_.size()-1; i++)
    {
      glVertex3d(seeds_[i].x(), seeds_[i].y(), seeds_[i].z());
      glVertex3d(seeds_[i+1].x(), seeds_[i+1].y(), seeds_[i+1].z());
    }
    glEnd();

    glLineWidth(1.0);
  }

  SeedTool::draw_gl(window);
}


}
