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
 *  NetworkEditor.cc: The network editor...
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 *  Distributed Dataflow changes:
 *   Michelle Miller
 *   Nov. 1997
 *
 */

#ifdef _WIN32
#  pragma warning(disable:4786)
#  include <objbase.h>
#  include <shellapi.h>
#endif

#include <Dataflow/Network/NetworkEditor.h>

#include <Dataflow/Comm/MessageBase.h>
#include <Dataflow/Network/Connection.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Network.h>
#include <Dataflow/Network/NetworkIO.h>
#include <Dataflow/Network/PackageDB.h>
#include <Dataflow/Network/Ports/Port.h>
#include <Dataflow/Network/Scheduler.h>
#include <Dataflow/Network/ComponentNode.h>
#include <Dataflow/Network/GenFiles.h>
#include <Core/XMLUtil/XMLUtil.h>

#include <Core/Math/MiscMath.h>
#include <Dataflow/GuiInterface/GuiCallback.h>
#include <Dataflow/GuiInterface/TCLInterface.h>
#include <Core/Util/StringUtil.h>
#include <Core/Thread/Thread.h>
#include <Core/Util/sci_system.h>
#include <Core/Util/Environment.h>
#include <Core/Exceptions/GuiException.h>
#include <Core/Events/EventManager.h>

#include <sci_debug.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
  
namespace SCIRun {

static bool
scheduling_starting_callback(void *data)
{
  const std::string proc("set_network_executing 1");
  TCLInterface::execute(proc);
  return true;
}

static bool
scheduling_done_callback(void *data)
{  
  if (sci_getenv("SCI_REGRESSION_TESTING"))
  {
#if DEBUG
    NetworkEditor* ne = (NetworkEditor*)data;

    //! Cleanup the network
    TCLInterface::execute("ClearCanvas 0");

    ne->get_network()->wait_all_deleted();
#endif
    //! Give modules some time to clean up
  }

  if (sci_getenv("SCI_REGRESSION_TESTING")||sci_getenv("SCIRUN_QUIT_WHEN_DONE"))
  {
    //! Quit SCIRun
    TCLInterface::execute("netedit quit");   
  }

  const std::string proc("set_network_executing 0");
  TCLInterface::execute(proc);

  return true;
}


NetworkEditor::NetworkEditor(Network* net) :
  net_(net)
{
  netio_ = new NetworkIO(net);
  //! We add these with highest priority so it is called first
  net->get_scheduler()->add_start_callback(scheduling_starting_callback, this, 255);
  //! We add these with lowest priority so it is called last
  net->get_scheduler()->add_callback(scheduling_done_callback, this, -255);

  // Create User interface...
  TCLInterface::add_command("netedit", this, 0);
}

NetworkEditor::~NetworkEditor()
{
  delete netio_;
}


void
NetworkEditor::tcl_command(GuiArgs& args, void*)
{
  if (args.count() < 2) 
  {
    throw "netedit needs a minor command";
  }
  if(args[1] == "quit") 
  {
    if (sci_getenv("SCI_REGRESSION_TESTING"))
    {
      std::cout<<"SCIRUN FINISHED NETWORK EXECUTION WITH "<< 
        net_->get_network_error_code() << " ERROR(S) AND " << 
        debug_number_of_objects() << " MEMORY LEAK(S)\n";
      std::ostringstream oss;
      oss <<"SCIRUN FINISHED NETWORK EXECUTION WITH "<< 
        net_->get_network_error_code() << " ERROR(S) AND " << 
        debug_number_of_objects() << " MEMORY LEAK(S)\n";
        net_->add_log(oss.str());

#ifdef DEBUG
// Print a list of objects that wasn't freed
// Check if ay objects remain in memory space
      if (debug_number_of_objects() != 0)
      {
        std::cerr << "LIST OF MEMORY LEAKS:" << std::endl;
        debug_print_objects();
      }
#endif
    }

    Thread::exitAll(net_->get_network_error_code()||debug_number_of_objects());
  } 
  else if (args[1] == "addmodule") 
  {
    if(args.count() < 5)
      throw "netedit addmodule needs a package name, category name and module name";

    ModuleHandle handle = net_->add_module(args[2],args[3],args[4]);
    if(!(handle.get_rep()))
      throw "netedit addmodule cannot add module "+args[2]+" "+args[3]+" "+args[4];
    
    TCLInterface::add_command(handle->id_+"-c", handle.get_rep(), 0);
    args.result(handle->id_);
  } 
  else if (args[1] == "deletemodule") 
  {
    if(args.count() < 3)
      throw "netedit deletemodule needs a module name";

    if(!(net_->delete_module(args[2]))) 
      throw GuiException("Cannot delete module "+args[2]);
  } 
  else if (args[1] == "deletemodule_warn") 
  {
    if(args.count() < 3)
      throw "netedit deletemodule_warn needs a module name";

    ModuleHandle handle=net_->get_module_by_id(args[2]);
    
    if (!(handle.get_rep()))
      throw "get_module_by_name failed for "+args[2];
    
    handle->delete_warn();
  } 
  else if(args[1] == "addconnection") 
  {
    if(args.count() < 6)
      throw "netedit addconnection needs 4 args";
    
    ModuleHandle ohandle = net_->get_module_by_id(args[2]);
    if(!(ohandle.get_rep()))
      throw "netedit addconnection can't find output module";
  
    int owhich = args.get_int(3);
    
    ModuleHandle ihandle = net_->get_module_by_id(args[4]);
    if(!(ihandle.get_rep()))
      throw "netedit addconnection can't find input module";
    int iwhich = args.get_int(5);
    args.result(net_->connect(ohandle, owhich, ihandle, iwhich));
  } 
  else if(args[1] == "deleteconnection") 
  {
    if (args.count() < 3)
      throw "netedit deleteconnection needs 1 arg";
  
    if (args.count() == 4 && args[3] == "1")
      net_->disable_connection(args[2]);
    
    if (!(net_->disconnect(args[2])))
      throw "Cannot find connection "+args[2]+" for deletion";
  } 
  else if(args[1] == "supportsPortCaching") 
  {
    if(args.count() < 4)
      throw "netedit supportsPortCaching needs 2 args";
    ModuleHandle ohandle = net_->get_module_by_id(args[2]);
    if(!ohandle.get_rep())
    {
      args.result("0");
    }
    else
    {
      const int owhich = args.get_int(3);
      if (owhich >= ohandle->num_output_ports())
        throw "netedit supportsPortCaching can't find output port";
      
      args.result(ohandle->oport_supports_cache_flag(owhich)?"1":"0");
    }
  } 
  else if(args[1] == "isPortCaching") 
  {
    if(args.count() < 4)
      throw "netedit isPortCaching needs 4 args";
    
    ModuleHandle ohandle = net_->get_module_by_id(args[2]);
    if(!(ohandle.get_rep()))
      throw "netedit isPortCaching can't find output module";
    
    const int owhich = args.get_int(3);
    
    if (owhich >= ohandle->num_output_ports())
      throw "netedit isPortCaching can't find output port";
  
    args.result(ohandle->get_oport_cache_flag(owhich)?"1":"0");
  } 
  else if(args[1] == "setPortCaching") 
  {
    if(args.count() < 5)
      throw "netedit setPortCaching needs 5 args";
    
    ModuleHandle ohandle = net_->get_module_by_id(args[2]);
    if(!ohandle.get_rep())
      throw "netedit setPortCaching can't find output module";
      
    const int owhich = args.get_int(3);
    if (owhich >= ohandle->num_output_ports())
      throw "netedit setPortCaching can't find output port";
    const int cache =  args.get_int(4);
    ohandle->set_oport_cache_flag(owhich, cache);
  } 
  else if(args[1] == "packageNames") 
  {
    args.result(args.make_list(packageDB->packageNames()));
  } 
  else if(args[1] == "categoryNames") 
  {
    if(args.count() != 3)
      throw "Usage: netedit categoryNames <packageName>";
    args.result(args.make_list(packageDB->categoryNames(args[2])));
  } 
  else if(args[1] == "categoryVNames") 
  {
    if(args.count() != 3)
      throw "Usage: netedit categoryVNames <packageName>";
    args.result(args.make_list(packageDB->categoryVNames(args[2])));
  }   
  else if(args[1] == "moduleNames") 
  {
    if(args.count() != 4)
      throw "Usage: netedit moduleNames <packageName> <categoryName>";
    args.result(args.make_list(packageDB->moduleNames(args[2],args[3])));
  } 
  else if(args[1] == "moduleVNames") 
  {
    if(args.count() != 4)
      throw "Usage: netedit moduleNames <packageName> <categoryName>";
    args.result(args.make_list(packageDB->moduleVNames(args[2],args[3])));
  } 
  else if(args[1] == "getCategoryName") 
  {
    if(args.count() != 5)
      throw "Usage: netedit getCategoryName <packageName> <categoryName> <moduleName>";
    args.result(packageDB->getCategoryName(args[2], args[3], args[4]));
  } 
  else if(args[1] == "dontschedule")
  {
  } 
  else if(args[1] == "scheduleok")
  {
    net_->schedule();
  } 
  else if(args[1] == "scheduleall")
  {
    net_->schedule_all();
  } 
  else if(args[1] == "reset_scheduler")
  {
    net_->reset_scheduler();
  } 
  else if(args[1] == "packageName")
  {
    if(args.count() != 3)
      throw "packageName needs a module id";
    
    ModuleHandle handle=net_->get_module_by_id(args[2]);
    if(!(handle.get_rep()))
      throw "cannot find module "+args[2];
    
    args.result(handle->package_name_);
  } 
  else if(args[1] == "moduleVersion")
  {
    if(args.count() != 3)
      throw "moduleVersion needs a module id";
    
    ModuleHandle handle=net_->get_module_by_id(args[2]);
    if(!(handle.get_rep()))
      throw "cannot find module "+args[2];
    
    args.result(handle->get_version());
  } 
  else if(args[1] == "categoryName")
  {
    if(args.count() != 3)
      throw "categoryName needs a module id";
    
    ModuleHandle handle=net_->get_module_by_id(args[2]);
    if(!(handle.get_rep()))
      throw "cannot find module "+args[2];
    
    args.result(handle->category_name_);
  } 
  else if(args[1] == "moduleName")
  {
    if(args.count() != 3)
      throw "moduleName needs a module id";
  
    ModuleHandle handle=net_->get_module_by_id(args[2]);
    if(!(handle.get_rep()))
      throw "cannot find module "+args[2];
    
    args.result(handle->module_name_);
  } 
  else if (args[1] == "create_pac_cat_mod") 
  {
    if (args.count()!=7)
      throw "create_pac_cat_mod needs 5 arguments";
    ModuleInfo mi;
    bool success = read_component_file(mi, args[6].c_str());
    if (! success)
      throw "NetworkEditor: 0) XML file did not pass validation: " + 
        args[2] + ".  Please see the messages window for details.";

    if (!(GenPackage((char*)args[3].c_str(),(char*)args[2].c_str()) &&
          GenCategory((char*)args[4].c_str(),(char*)args[3].c_str(),
                      (char*)args[2].c_str()) &&
          GenComponent(mi, (char*)args[3].c_str(),(char*)args[2].c_str())))
      throw "Unable to create new package, category or module. Check your paths and names and try again.";

  } 
  else if (args[1] == "create_cat_mod") 
  {
    if (args.count()!=7)
      throw "create_cat_mod needs 3 arguments";
    
    ModuleInfo mi;
    bool success = read_component_file(mi, args[6].c_str());
    if (!success)
      throw "NetworkEditor: 1) XML file did not pass validation: " + 
        args[2] + ".  Please see the messages window for details.";

    if (!(GenCategory((char*)args[4].c_str(),(char*)args[3].c_str(),
                      (char*)args[2].c_str()) &&
          GenComponent(mi, (char*)args[3].c_str(),(char*)args[2].c_str())))
      throw "Unable to create new category or module. Check your paths and names and try again.";
  } 
  else if (args[1] == "create_mod")
  {
    if (args.count()!=7) 
      throw "create_mod needs 6 arguments";
    ModuleInfo mi;
    bool success = read_component_file(mi, args[6].c_str());
    if (! success)
      throw "NetworkEditor: 2) XML file did not pass validation: " + 
          args[2] + ".  Please see the messages window for details.";
    if (!(GenComponent(mi, (char*)args[3].c_str(),(char*)args[2].c_str())))
      throw "Unable to create new module. Check your paths and names and try again.";
  } 
  else if (args[1] == "sci_system" && args.count() > 2) 
  {
    std::string command = args[2];
    for (int i = 3; i < args.count(); i++) {
      command = command + " " + args[i];
    }
    args.result(to_string(sci_system(command.c_str())));
  } 
  else if (args[1] == "run_windows_browser" && args.count() >= 2)
  {
#if defined (_WIN32)
    //std::cout << "runbrowser " << args[2] << std::endl;
    std::string url = args[2];
    CoInitializeEx(0, COINIT_APARTMENTTHREADED | COINIT_DISABLE_OLE1DDE);
    // using ANSI version
    HINSTANCE r = ShellExecuteA(0, LPCTSTR("open"), LPCTSTR(url.c_str()), 0, 0, SW_SHOWDEFAULT);
    // According to the MSDN documentation (http://msdn.microsoft.com/en-us/library/bb762153(VS.85).aspx)
    // the returned HINSTANCE must be cast to an int to check for errors.
    // A successful call to ShellExecute will return a value > 32.
    if ( (int) r > 32) // success!
      return;
    else if ( (int) r == ERROR_FILE_NOT_FOUND || (int) r == ERROR_PATH_NOT_FOUND || (int) r == SE_ERR_PNF )
      throw url + " could not be found";
    else if ( (int) r == SE_ERR_ACCESSDENIED )
      throw "Access to web browser was denied.";
    else
      throw "Opening " + url + " in a web browser failed";
#else
    throw "Invoking a web browser using this callback on this platform is not currently supported.";
#endif
  }
  else if (args[1] == "getenv" && args.count() == 3)
  {
    const char *result = sci_getenv( args[2] );
    if (result) 
    {
      args.result(std::string(result));
    }
  } 
  else if (args[1] == "setenv" && args.count() == 4)
  {
    sci_putenv(args[2], args[3]);
  } 
  else if (args[1] == "update_rcfile" && args.count() == 4)
  {
    update_rcfile(args[2], args[3]);
  } 
  else if (args[1] == "update_env" && args.count() == 4)
  {
    sci_putenv(args[2],args[3]);
    update_rcfile(args[2], args[3]);
  } 
  else if (args[1] == "module_oport_datatypes") 
  {
    if (args.count() != 5)
      throw "netedit module_oport_datatypes expects a package, category, and module";
    
    const ModuleInfo* info = packageDB->GetModuleInfo(args[4],args[3],args[2]);
    if (!info)
      throw "netedit module_oports cant find "+ args[2]+"->"+args[3]+"->"+args[4];
    std::string result;

    std::vector<OPortInfo*>::const_iterator i2 = info->oports_.begin();
    while (i2 < info->oports_.end())
    {
      OPortInfo* op = *i2++;
      result += op->datatype + " ";
    }
    args.result(result);
  } 
  else if (args[1] == "module_iport_datatypes") 
  {
    if (args.count() != 5)
      throw "netedit module_iport_datatypes expects a package, category, and module";
    
    const ModuleInfo* info = packageDB->GetModuleInfo(args[4],args[3],args[2]);
    if (!info)
      throw "netedit module_oports cant find "+args[2]+"->"+args[3]+"->"+args[4];
    
    std::string result;
    std::vector<IPortInfo*>::const_iterator i1 = info->iports_.begin();
    while (i1 < info->iports_.end())
    {
      IPortInfo* ip = *i1++;
      result += ip->datatype + " ";
    }
    args.result(result);
  } 
  else if (args[1] == "presave") 
  {
    net_->pre_save_network();
  } 
  else if (args[1] == "start-net-doc") 
  {
    netio_->start_net_doc(args[2], args[3], args[4]);
  } 
  else if (args[1] == "write-net-doc") 
  {
    netio_->write_net_doc();
  } 
  else if (args[1] == "network-variable") 
  {
    netio_->add_net_var(args[2], args[3]);
  } 
  else if (args[1] == "net-add-env-var") 
  {
    netio_->add_environment_sub(args[2], args[3]);
  } 
  else if (args[1] == "network-note") 
  {
    netio_->add_net_note(args[2]);
  } 
  else if (args[1] == "add-module") 
  {
    if (args.count() == 6)
    {
      std::string version = "1.0";
      netio_->add_module_node(args[2], args[3], args[4], args[5],version);
    }
    else 
      netio_->add_module_node(args[2], args[3], args[4], args[5], args[6]);
  } 
  else if (args[1] == "module-position") 
  {
    netio_->add_module_position(args[2], args[3], args[4]);
  } 
  else if (args[1] == "mod-note") 
  {
    netio_->add_module_note(args[2], args[3]);
  } 
  else if (args[1] == "mod-note-pos") 
  {
    netio_->add_module_note_position(args[2], args[3]);
  } 
  else if (args[1] == "mod-note-col") 
  {
    netio_->add_module_note_color(args[2], args[3]);
  } 
  else if (args[1] == "mod-connection") 
  {
    netio_->add_connection_node(args[2], args[3], args[4], args[5], args[6]);
  } 
  else if (args[1] == "conn-disabled") 
  {
    netio_->set_disabled_connection(args[2]);
  } 
  else if (args[1] == "conn-route") 
  {
    netio_->add_connection_route(args[2], args[3]);
  } 
  else if (args[1] == "conn-note") 
  {
    netio_->add_connection_note(args[2], args[3]);
  } 
  else if (args[1] == "conn-note-pos") 
  {
    netio_->add_connection_note_position(args[2], args[3]);
  } 
  else if (args[1] == "conn-note-col") 
  {
    netio_->add_connection_note_color(args[2], args[3]);
  } 
  else if (args[1] == "set-port-caching") 
  {
    netio_->set_port_caching(args[2], args[3], args[4]);
  } 
  else if (args[1] == "set-modgui-visible") 
  {
    netio_->set_module_gui_visible(args[2]);
  } 
  else if (args[1] == "add-mod-var") 
  {
    netio_->add_module_variable(args[2], args[3], args[4],false,false);
  } 
  else if (args[1] == "add-mod-substvar") 
  {
    netio_->add_module_variable(args[2], args[3], args[4],false,true);
  } 
  else if (args[1] == "add-mod-filevar") 
  {
		if (args[5] == "1" || args[5] == "true" || args[5] == "yes" || args[5] == "on")
		{
		  netio_->add_module_variable(args[2], args[3], args[4],true,true,true);
		}
		else
		{
		  netio_->add_module_variable(args[2], args[3], args[4],true,true,false);
		}
  } 
  else if (args[1] == "add-modgui-callback") 
  {
    netio_->add_module_gui_callback(args[2], args[3]);
  } 
  else if (args[1] == "subnet-start") 
  {
    netio_->push_subnet_scope(args[2], args[3]);
  } 
  else if (args[1] == "subnet-end") 
  {
    netio_->pop_subnet_scope();
  } 
  else if (args[1] == "load_srn") 
  {
    net_->add_log("NETWORK EDITOR: LOADING NETWORK "+args[2]);
    netio_->load_net(args[2]);
  }
  else if (args[1] == "reload_rcfile")
  {
    const char* rcfile = sci_getenv("SCIRUN_RCFILE");
    std::cout << "reloading rcfile....\n";
    parse_rcfile(rcfile);
  }
  else  
  {
    throw "Unknown minor command for netedit: "+args[1];
  }
} // end tcl_command()
  
} // End namespace SCIRun
