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
 *  GeometryPort.h: Handle to the Geometry Data type
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 */

#ifndef DATAFLOW_NETWORK_PORTS_GEOMETRYPORT_H
#define DATAFLOW_NETWORK_PORTS_GEOMETRYPORT_H 1

#include <Core/Thread/Mailbox.h>
#include <Core/Geom/GeomObj.h>
#include <Core/Geom/Light.h>
#include <Core/Geometry/BBox.h>
#include <Core/Containers/LockingHandle.h>

#include <Dataflow/Network/Ports/Port.h>
#include <Dataflow/Comm/MessageBase.h>

#include <string>
#include <vector>
#include <list>

#include <Dataflow/Network/share.h>

namespace SCIRun {

class View;
class CrowdMonitor;
class Mutex;
class ColorImage;
class DepthImage;
class GeometryComm;
struct GeomReply;

typedef int GeomID;
typedef short LightID;

//! This class is a virtual dummy class as it is only a place holder.
//! all the data is send to the mailbox of the module instead, this way
//! all messages end up in the same mailbox and hence that thread has only
//! to wait on one mailbox
class SCISHARE GeometryIPort : public IPort {
public:
  GeometryIPort(Module*, const std::string& name);
  virtual ~GeometryIPort();
};


 
struct SCISHARE GeometryData {

  GeometryData();
  ~GeometryData();

  ColorImage* colorbuffer;
  DepthImage* depthbuffer;

  View* view;
  int xres, yres;
  double znear, zfar;
  BBox   view_bounds_;
  
  void Print();
};



#define GEOM_VIEW		1
#define GEOM_COLORBUFFER	2
#define GEOM_DEPTHBUFFER	4
#define GEOM_VIEW_BOUNDS	8
#define GEOM_ALLDATA		(GEOM_VIEW|GEOM_COLORBUFFER|GEOM_DEPTHBUFFER)
#define GEOM_TRIANGLES		16
// CollabVis code begin
#define GEOM_MATRICES           32
// CollabVis code end

  
class SCISHARE GeometryOPort : public OPort {
private:

  GeomID    serial_;
  LightID   lserial_;
  bool      dirty_;

  std::list<GeometryComm* > saved_msgs_;

  std::vector<int> portid_;
  std::vector<ModuleHandle> modules_;

  void save_msg(GeometryComm*);

  Mailbox<GeomReply> reply_mailbox_;

public:
  GeometryOPort(Module*, const std::string& name);
  virtual ~GeometryOPort();

  GeomID addObj(GeomHandle, const std::string& name, CrowdMonitor* lock=0);
  LightID addLight(LightHandle, const std::string& name, CrowdMonitor* lock=0);
  void delObj(GeomID, int del=1);
  void delLight(LightID, int del = 1);
  void delAll();
  void flush();
  void flushViews();
  bool get_view_bounds(BBox &bbox);
  void flushViewsAndWait();

  void forward(GeometryComm* msg);
  bool direct_forward(GeometryComm* msg);

  virtual void reset();
  virtual void finish();
  virtual void synchronize();

  virtual void attach(Connection* c);
  virtual void detach(Connection* c, bool blocked);

  virtual bool have_data();
  virtual void resend(ConnectionHandle);

  GeometryData* getData(int which_viewer, int which_viewwindow, int mask);
  void setView(int which_viewer, int which_viewwindow, View view);
};

typedef LockingHandle<GeometryIPort> GeometryIPortHandle;
typedef LockingHandle<GeometryOPort> GeometryOPortHandle;


} // End namespace SCIRun


#endif /* SCI_project_GeometryPort_h */
