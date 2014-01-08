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
 *  ShowLinePath.cc:
 *
 *  Written by:
 *   burak
 *   TODAY'S DATE HERE
 *
 */

#include <Core/Geom/GeomLine.h>
#include <Core/Datatypes/ColorMap.h>
#include <Core/Util/StringUtil.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Ports/ColorMapPort.h>

#include <Core/Thread/CrowdMonitor.h>

#include <vector>
#include <string>

namespace SCIRun {

using namespace SCIRun;

class ShowLinePath : public Module {
  public:
    ShowLinePath(GuiContext*);

    CrowdMonitor control_lock_;
    GeomLines *lines_;
    GeomID     lineID_;

    virtual ~ShowLinePath();
    virtual void execute();
};


DECLARE_MAKER(ShowLinePath)

ShowLinePath::ShowLinePath(GuiContext* ctx) :
  Module("ShowLinePath", ctx, Source, "Visualization", "SCIRun"),
  control_lock_("ShowLinePath resolution lock"),
  lines_(0), lineID_(0)
{
}

ShowLinePath::~ShowLinePath()
{
}


void
ShowLinePath::execute()
{
  GeometryOPortHandle geomport;
  get_oport_handle("GeomOut",geomport);
  geomport->delAll();

  MatrixHandle inmatrix;
  ColorMapHandle incolors;
  
  // Handle checking input as well
  get_input_handle("PointPath",inmatrix,true);
  get_input_handle("TransitionColor",incolors,true);

  // Inform module that execution started
  update_state(Executing);

  lines_ = new GeomLines();

  if(inmatrix->nrows() == 3 && inmatrix->ncols() >= 2)
  {
    for(index_type m=0;m<inmatrix->ncols()-1;m++)
    {
      lines_->add(Point(inmatrix->get(0,m), inmatrix->get(1,m), inmatrix->get(2,m)),
                   incolors->lookup(2*(float(m)/float(inmatrix->ncols()))-1.0),
                   Point(inmatrix->get(0,m+1), inmatrix->get(1,m+1), inmatrix->get(2,m+1)),
                   incolors->lookup(2*(float(m+1)/float(inmatrix->ncols()))-1.0)
                   );
    }
  }
  else
  {
	  if(inmatrix->ncols() == 3 && inmatrix->nrows() >= 2)
	  {
		  for(index_type m=0;m<inmatrix->nrows()-1;m++)
		  {
			  lines_->add(Point(inmatrix->get(m,0), inmatrix->get(m,1), inmatrix->get(m,2)),
						  incolors->lookup(2*(float(m)/float(inmatrix->nrows()))-1.0),
						  Point(inmatrix->get(m+1,0), inmatrix->get(m+1,1), inmatrix->get(m+1,2)),
						  incolors->lookup(2*(float(m+1)/float(inmatrix->nrows()))-1.0)
						  );
			  
		  }
	  }
	  else
	  {
	      error("PointPath must be a 3xN or Nx3 matrix with N>1");
	  }
  }

  lineID_ = geomport->addObj(lines_,"ShowLinePathLines",&control_lock_);
  geomport->flushViews();
}

} // End namespace SCIRun


