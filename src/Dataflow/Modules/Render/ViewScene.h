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
 *  ViewScene.h: The Geometry ViewScene!
 *
 *  Written by:
 *   Steven G. Parker & Dave Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 */

#ifndef SCIRun_src_Dataflow_Modules_Render_ViewScene_h
#define SCIRun_src_Dataflow_Modules_Render_ViewScene_h

#include <Dataflow/Network/Module.h>
#include <Dataflow/Comm/MessageBase.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/GeometryComm.h>
#include <Core/Geom/GeomObj.h>
#include <Core/Geom/GeomMaterial.h>
#include <Core/Geom/Lighting.h>
#include <Core/Geom/IndexedGroup.h>
#include <Dataflow/Modules/Render/ViewGeom.h>
#include <Core/Thread/RecursiveMutex.h>
#include <Core/Thread/CrowdMonitor.h>

#include <map>
#include <vector>
#include <string>

namespace SCIRun {

class ViewWindow;

class ViewScene : public Module {
public:
  ViewScene(GuiContext*);
  virtual ~ViewScene();
  virtual void			do_execute();
  virtual void			execute();

  MaterialHandle      default_material_;
  Lighting            lighting_;
  CrowdMonitor        geomlock_;
  GeomIndexedGroup		ports_;

private:
  void				tcl_command(GuiArgs&, void*);  
  void				delete_viewwindow(const std::string &);
  void				initPort(Mailbox<GeomReply>*);
  void				detachPort(int portno);
  int         real_portno(int portid);
  void				delete_patch_portnos(int portid);
  void				append_port_msg(GeometryComm*);
  void				addObj(GeomViewerPort* port, 
				       GeomID serial, GeomHandle obj,
				       const std::string&, CrowdMonitor* lock);
  void				delObj(GeomViewerPort* port, GeomID serial);
  void				delAll(GeomViewerPort* port);
  void				flushPort(int portid);
  void				finishPort(int portid);
  void				flushViews();
  int         process_event();

  virtual void set_context(Network* network);

  virtual void prepare_cleanup();

  std::map<int, std::map<LightID, int> >	pli_;  // port->light->index
  Mutex                         view_window_lock_;
  std::vector<ViewWindow*>           view_window_;
  int                           max_portno_;
  bool                          stop_rendering_;
  std::vector<int>                   portno_map_;
  std::vector<int>                   synchronized_map_;
  std::list<unsigned int>            synchronized_serials_;
  int                           synchronized_debt_;

  static bool regression_save_image_callback(void *stuff);
  static bool save_image_callback(void *stuff);
  
  bool delete_check_autoview_on_load_;
};



class ViewSceneMessage : public MessageBase {
public:
  std::string			rid;
  std::string			filename;
  std::string			format;
  int				resx;
  int				resy;
  double			tbeg;
  double			tend;
  int				nframes;
  double			framerate;
  Point                         loc;
  Vector			dir;
  Color				color;
  int				index;
  bool				on;

  ViewSceneMessage(const std::string& rid);
  ViewSceneMessage(const std::string& rid, double tbeg, double tend,
		int nframes, double framerate);
  ViewSceneMessage(MessageTypes::MessageType,
		const std::string& rid, int lightNo, bool on, 
		const Vector& lightDir, const Color& lightColor);
  ViewSceneMessage(MessageTypes::MessageType,
                const std::string& rid, int clip_no,
                const Point& center, const Vector& normal,
                double width, double height,double scale);
  ViewSceneMessage(MessageTypes::MessageType,
		const std::string& rid, const std::string& filename);
  ViewSceneMessage(MessageTypes::MessageType,
		const std::string& rid, const std::string& filename,
		const std::string& format, const std::string &resx_string,
		const std::string& resy_string);
  virtual ~ViewSceneMessage();
};

} // End namespace SCIRun

#endif // of #ifndef SCIRun_src_Dataflow_Modules_Render_ViewScene_h
