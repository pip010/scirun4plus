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
 *  ManageSeriesModule.h: Manage a series of data
 *
 *  Written by:
 *   Allen R. Sanderson
 *   SCI Institute
 *   University of Utah
 *   April 2007
 *
 */

#include <Dataflow/Network/Module.h>
#include <Core/Util/StringUtil.h>
#include <Core/Containers/Handle.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <iostream>

namespace SCIRun {

template< class HANDLE_TYPE >
class ManageSeriesModule : public Module {

public:
  ManageSeriesModule(const std::string& module_name,
		    GuiContext* ctx, SchedClass,
		    const std::string& catagory  = "unknown",
		    const std::string& package   = "unknown",
		    const std::string& port_name = "unknown" );

  virtual ~ManageSeriesModule();
  virtual void execute();
  virtual void tcl_command(GuiArgs& args, void* userdata);

private:
  std::string port_name_;

  GuiInt gui_nports_;

  std::vector< HANDLE_TYPE > handles_;
};

template< class HANDLE_TYPE >
ManageSeriesModule< HANDLE_TYPE >::ManageSeriesModule(const std::string& module_name,
						    GuiContext* context,
						    SchedClass,
						    const std::string& catagory,
						    const std::string& package,
						    const std::string& port_name )
  : Module(module_name, context, Filter, catagory, package),
    port_name_(port_name),
    gui_nports_(context->subVar("num-ports"), 2)
{
}


template< class HANDLE_TYPE >
ManageSeriesModule< HANDLE_TYPE >::~ManageSeriesModule()
{
}

template< class HANDLE_TYPE >
void
ManageSeriesModule< HANDLE_TYPE >::execute()
{
  HANDLE_TYPE input_handle;

  if( !get_input_handle( "Input", input_handle, true ) )
    return;

  // Add the new handle to the stack.
  if( inputs_changed_ || handles_.size() == 0 ) {
    handles_.insert( handles_.begin(), input_handle );
  }

  // Make sure there is a handle for each port.
  while( handles_.size() < (unsigned int)gui_nports_.get() ) {
    handles_.push_back( handles_[ handles_.size()-1] );
  }

  // Remove any excess handles.
  while( handles_.size() > (unsigned int)gui_nports_.get() ) {
    handles_.pop_back();
  }

  // Send the data downstream - NOTE: manage the handles within the module
  for( unsigned int i=0; i<handles_.size(); i++ ) {
    std::string port_name = "Output " + to_string(i);
    send_output_handle( port_name, handles_[i], true );
  }
}

template< class HANDLE_TYPE >
void
ManageSeriesModule< HANDLE_TYPE >::tcl_command(GuiArgs& args, void* userdata)
{
  if(args.count() < 2){
    args.error("ManageSeriesModule needs a minor command");
    return;
  }

  if (args[1] == "clear") {
    handles_.clear();
  } else {
    Module::tcl_command(args, userdata);
  }
}

} // End namespace SCIRun
