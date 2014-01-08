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
//    File   : SceneGraphEvent.h
//    Author : Martin Cole
//    Date   : Tue Jun  6 10:08:25 2006

#if !defined(SceneGraphEvent_h)
#define SceneGraphEvent_h

#include <Core/Events/BaseEvent.h>
#include <Core/Geom/GeomObj.h>
#include <string>


#include <Core/Events/share.h>
namespace SCIRun {

class SCISHARE SceneGraphEvent : public BaseEvent
{
public:
  SceneGraphEvent();
  SceneGraphEvent(GeomHandle o, const std::string& n,
		  const std::string &target = "",
		  unsigned int time = 0);
  virtual ~SceneGraphEvent();

  virtual bool          is_scene_graph_event() { return true; }

  SceneGraphEvent *     clone() { return new SceneGraphEvent(*this); }

  //! Accessors
  GeomHandle          get_geom_obj() const { return obj_; }
  const std::string &      get_geom_obj_name() const { return name_; }
  int                 get_scene_graph_id() const { return sg_id_; }
  bool                toggle_visibility_p() const { return toggle_visibility_;}
  bool                force_visibility_p() const { return force_visibility_; }
  bool                forced_visibility_flag() const {return visibility_flag_;}

  //! Mutators
  void                  set_geom_obj(GeomHandle obj) { obj_ = obj; }
  void                  set_geom_obj_name(GeomHandle obj) { obj_ = obj; }
  void                  set_scene_graph_id(int id) { sg_id_ = id; }
  void                  set_toggle_visibility() { toggle_visibility_ = true; }
  void                  set_force_visibility(bool vis)
                        {
                          force_visibility_ = true;
                          visibility_flag_ = vis;
                        }

private:
  GeomHandle          obj_;
  std::string              name_;
  int                 sg_id_;
  bool                toggle_visibility_;
  bool                force_visibility_;
  bool                visibility_flag_;
};


} // namespace SCIRun

#endif // SceneGraphEvent_h
