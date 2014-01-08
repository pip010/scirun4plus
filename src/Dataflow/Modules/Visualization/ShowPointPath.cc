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
 *  ShowPointPath.cc:
 *
 *  Written by:
 *   burak
 *   TODAY'S DATE HERE
 *
 */

#include <Core/Util/StringUtil.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Widgets/PointWidget.h>
#include <Core/Util/MemoryUtil.h>

#include <Core/Thread/CrowdMonitor.h>

#include <vector>
#include <string>


namespace SCIRun {

using namespace SCIRun;

class ShowPointPath : public Module {
  public:
    ShowPointPath(GuiContext*);

    CrowdMonitor control_lock_;
    std::vector<PointWidget*> traj;

    virtual ~ShowPointPath();
    virtual void execute();
};


DECLARE_MAKER(ShowPointPath)

ShowPointPath::ShowPointPath(GuiContext* ctx) :
  Module("ShowPointPath", ctx, Source, "Visualization", "SCIRun"), control_lock_("ShowPointPath resolution lock")
{
}

ShowPointPath::~ShowPointPath()
{
  delete_all_items(traj);
}

void
ShowPointPath::execute()
{
  // DO TO: fix memory leak in wdiget allocation/deletion. It should be done with handles

  // Get handle to output port
  GeomHandle currpoint;

  GeometryOPortHandle geomport;
  get_oport_handle("GeomOut",geomport);
  
  MatrixHandle inmatrix;
  get_input_handle("PointPath",inmatrix,true);

  // Inform module that execution started
  update_state(Executing);

  geomport->delAll();

  delete_all_items(traj);
  traj.clear();

  if(inmatrix->nrows() == 3)
  {
    for(index_type m=0;m<inmatrix->ncols();m++)
    {
      traj.push_back(new PointWidget(this, &control_lock_, 1.0));
      traj[m]->SetPosition(Point(inmatrix->get(0,m), inmatrix->get(1,m), inmatrix->get(2,m)));

      currpoint = traj[m]->GetWidget();
      std::string str = "ShowPointPathPoint " + to_string(m); 
      geomport->addObj(currpoint,str,&control_lock_);
    }
  }
  else if(inmatrix->ncols() == 3)
  for(index_type m=0;m<inmatrix->nrows();m++)
  {
    traj.push_back(new PointWidget(this, &control_lock_, 1.0));
    traj[m]->SetPosition(Point(inmatrix->get(m,0), inmatrix->get(m,1), inmatrix->get(m,2)));

    currpoint = traj[m]->GetWidget();
    std::string str = "ShowPointPathPoint " + to_string(m); 
    geomport->addObj(currpoint,str,&control_lock_);
  }
  else
	error("PointPath must be a 3xN or Nx3 matrix with N>1");

  geomport->flushViews();
}


} // End namespace SCIRun


