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
 *  main.cc: 
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   Feb. 1994
 *
 */

#include <sci_defs/version_defs.h>
#include <sci_defs/opengl_defs.h>

#include <slivr/ShaderProgramARB.h>

#include <Core/Util/StringUtil.h>
#include <Core/Init/init.h>
#include <Core/Util/Environment.h>
#include <Core/Util/FileUtils.h>
#include <Core/Util/sci_system.h>
#include <Core/Thread/Thread.h>
#include <Core/Thread/Time.h>
#include <Core/Events/EventManager.h>

#include <Core/Services/ServiceLog.h>
#include <Core/Services/ServiceDB.h>
#include <Core/Services/ServiceManager.h>
#include <Core/SystemCall/SystemCallManager.h>

#ifdef _WIN32
#  include <Core/Geom/Win32OpenGLContext.h>
#else
#  include <X11/Xlib.h>
#  include <Core/Geom/X11OpenGLContext.h>
#endif

#include <Dataflow/Network/Network.h>
#include <Dataflow/Network/NetworkEditor.h>
#include <Dataflow/Network/PackageDB.h>
#include <Dataflow/Network/NetworkIO.h>
#include <Dataflow/Network/Scheduler.h>
#include <Dataflow/GuiInterface/TCLThread.h>
#include <Dataflow/GuiInterface/TCLInterface.h>

#include <sci_debug.h>

#include <string>
#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <boost/asio/io_service.hpp>

namespace SCIRun {

#ifdef _WIN32
#  include <windows.h>

void InvalidParameterHandler(const wchar_t* expression, 
                             const wchar_t* function, 
                             const wchar_t* file, 
                             unsigned int line,	
                             uintptr_t pReserved ) 
{
	fprintf(
		stderr,
		"Invalid parameter in function %s. File: %s Line: %d Expression: %s\n",
		function,
		file,
		line,
    expression
	);
  
	// Cause a Debug Breakpoint.
	DebugBreak();
}

#endif
class MainThread : public TCLThread
{
  public:
    MainThread(int argc, char* argv[], Network* net, int startnetno) : 
      TCLThread(argc,argv), net_(net), startnetno_(startnetno) {}
    virtual ~MainThread() {}
    virtual int run_tcl();

  private:
    Network*  net_;
    int       startnetno_;
};


void
usage()
{
  std::cout << "Usage: scirun [args] [net_file] [session_file]\n";
  std::cout << "    [-]-d[atadir]            : scirun data directory\n";
  std::cout << "    [-]-r[egression]         : regression test a network\n";
  std::cout << "    [-]-R[egressionimagedir] : output directory for regression tests\n";
  std::cout << "    [-]-s[erver] [PORT]      : start a TCL server on port number PORT\n";
  std::cout << "    [-]-e[xecute]            : executes the given network on startup\n";
  std::cout << "    [-]-E[xecute]            : executes the given network on startup and quits when done\n";
  std::cout << "    [-]-c[onvert]            : converts a .net to a .srn network and exits\n";
  std::cout << "    [-]-v[ersion]            : prints out version information\n";
  std::cout << "    [-]-h[elp]               : prints usage information\n";
  std::cout << "    [-]-p[ort] [PORT]        : start remote services port on port number PORT\n";
  std::cout << "    [-]-l[ogfile] file       : add output messages to a logfile\n";
  std::cout << "    [-]-t[imeout] N          : kill scirun after N seconds\n";
  std::cout << "    [-]-I[mage] file         : Save resulting images in file\n";
  std::cout << "    [--nosplash]             : disable the splash screen\n";
  std::cout << "    net_file                 : SCIRun Network Input File\n";
  std::cout << "    session_file             : PowerApp Session File\n";
  exit( 0 );
}

static bool doing_convert_ = false;
// Parse the supported command-line arugments.
// Returns the argument # of the .net file
int
parse_args( int argc, char *argv[] )
{
  int found = 0;
  bool powerapp = false;
  int cnt = 1;
  
  while (cnt < argc)
  {
    std::string arg( argv[ cnt ] );
    if( ( arg == "--version" ) || ( arg == "-version" ) || ( arg == "-v" ) || ( arg == "--v" ) )
    {
      std::cout << "Version: " << SCIRUN_VERSION << "\n";
      std::cout << "SVN Revision: " << SCIRUN_SVN_REVISION << "\n";
      std::cout << "SVN URL: " << SCIRUN_SVN_URL << "\n";
      exit( 0 );
    }
    else if ( ( arg == "--help" ) || ( arg == "-help" ) || ( arg == "-h" ) ||  ( arg == "--h" ) )
    {
      usage();
    }
    else if ( ( arg == "--execute" ) || ( arg == "-execute" ) || ( arg == "-e" ) ||  ( arg == "--e" ) )
    {
      sci_putenv("SCIRUN_EXECUTE_ON_STARTUP","1");
    }
    else if ( ( arg == "--Execute" ) || ( arg == "-Execute" ) || ( arg == "-E" ) ||  ( arg == "--E" ) )
    {
      sci_putenv("SCIRUN_EXECUTE_ON_STARTUP","1");
      sci_putenv("SCIRUN_QUIT_WHEN_DONE","1");
    }
    else if ( ( arg == "--Image" ) || ( arg == "-Image" ) || ( arg == "-I" ) ||  ( arg == "--I" ) )
    {
      if (cnt+1 < argc)
      {
        sci_putenv("SCIRUN_SAVE_IMAGE_FILE",argv[cnt+1]);
        cnt++;
      }
    }  
    else if ( ( arg == "--convert" ) || ( arg == "-convert" ) || ( arg == "-c" ) ||  ( arg == "--c" ) )
    {
      sci_putenv("SCIRUN_CONVERT_NET_TO_SRN","1");
      doing_convert_ = true;
    }
    else if ( ( arg == "--eai" ) || ( arg == "-eai" ))
    {
      sci_putenv("SCIRUN_EXTERNAL_APPLICATION_INTERFACE","1");
    }
    else if ( ( arg == "--noeai" ) || ( arg == "-noeai" ))
    {
      sci_putenv("NO_SCIRUN_EXTERNAL_APPLICATION_INTERFACE","1");
    }
    else if ( ( arg == "--regression" ) || ( arg == "-regression" ) || ( arg == "-r" ) ||  ( arg == "--r" ) )
    {
      sci_putenv("SCI_REGRESSION_TESTING","1");
    }
    else if ( ( arg == "--Regressionimagedir" ) || ( arg == "-Regressionimagedir" ) || ( arg == "--regressionimagedir" ) || ( arg == "-regressionimagedir" ) || ( arg == "-R" ) ||  ( arg == "--R" ) )
    {
      if (cnt+1 < argc)
      {
        sci_putenv("SCI_REGRESSION_IMAGE_DIR",argv[cnt+1]);
        cnt++;
      }
    }
    else if ( ( arg == "--datadir" ) || ( arg == "-datadir" ) || ( arg == "-d" ) ||  ( arg == "--d" ) )
    {
      if (cnt+1 < argc)
      {
        sci_putenv("SCIRUN_DATA",argv[cnt+1]);
        cnt++;
      }
    }
    else if ( arg == "--nosplash" )
    {
      sci_putenv("SCIRUN_NOSPLASH", "1");
    }
    else if ( (arg == "--logfile") || (arg == "--log") || (arg == "-l") || (arg == "--l"))
    {
      if (cnt+1 < argc)
      {
        sci_putenv("SCIRUN_LOGFILE",argv[cnt+1]);
        cnt++;
      }
    }
    else if (ends_with(string_tolower(arg),".srn") || ends_with(string_tolower(arg),".srn#"))
    {
      if (!validFile(arg))
      {
        std::cerr << "Couldn't find net file " << arg
            << ".\nNo such file or directory.  Exiting." 
            << std::endl;
        exit(0);
      }
      sci_putenv("SCIRUN_NETFILE", arg);
    }
#ifdef BUILD_BIOMESH3D_REMOTE_SUPPORT
			else if ( ( arg == "--server_port" ) || ( arg == "-server_port" ) || ( arg == "-s_p" ) ||  ( arg == "--s_p" ) )
    {
      int port;
      if ((cnt+1 < argc) && string_to_int(argv[cnt+1], port)) 
      {
        if (port < 1024 || port > 65535) 
        {
          std::cerr << "Server port must be in range 1024-65535\n";
          exit(0);
        }
        cnt++;
      } 
      else 
      {
        port = 0;
      }
      sci_putenv("SCIRUN_SERVER_PORT",to_string(port));
    }    
			else if ( ( arg == "--port_range" ) || ( arg == "-port_range" ) || ( arg == "-p_r" ) ||  ( arg == "--p_r" ) )
    {
				sci_putenv("SCIRUN_PORT_RANGE",to_string(argv[cnt+1]));
				cnt++;
			}
			else if ( ( arg == "--ip_address" ) || ( arg == "-ip_address" ) || ( arg == "-ip" ) ||  ( arg == "--ip" ) )
        {
				sci_putenv("SCIRUN_IP_ADDRESS",to_string(argv[cnt+1]));
				cnt++;
        }
			else if ( ( arg == "--protocol" ) || ( arg == "-protocol" ) )
			{
				sci_putenv("SCIRUN_PROTOCOL",to_string(argv[cnt+1]));
        cnt++;
      } 
			else if (ends_with(string_tolower(arg),".srn") || ends_with(string_tolower(arg),".srn#"))
      {
				if (!validFile(arg))
				{
					std::cerr << "Couldn't find net file " << arg
						<< ".\nNo such file or directory.  Exiting." 
						<< std::endl;
					exit(0);
      }
				sci_putenv("SCIRUN_NETFILE", arg);
    }
#endif
    else if ( ( arg == "--timeout" ) || ( arg == "-timeout" ) || ( arg == "-t" ) ||  ( arg == "--t" ) )
    {
      int timeout;
      if ((cnt+1 < argc) && string_to_int(argv[cnt+1], timeout)) 
      {
        cnt++;
      } 
      else 
      {
        timeout = 250;
      }
      sci_putenv("SCIRUN_TIMEOUT",to_string(timeout));
    }
    else if (arg[0] == '+')
    {
      std::string key = arg.substr(1);
      std::string value;
      std::string::size_type loc = arg.find("=");
      if (loc != std::string::npos)
      { 
        key = arg.substr(1,loc-1);
        value = arg.substr(loc+1);
      }
      sci_putenv(key,value);
    }
    else
    {
      if (!validFile(arg))
      {
        std::cerr << "Couldn't find net file " << arg
            << ".\nNo such file or directory.  Exiting." 
            << std::endl;
        exit(0);
      }
      std::string filename(string_tolower(arg));
      if (!ends_with(filename,".net") && !ends_with(filename,".app") &&
        !ends_with(filename,".srn#"))
      {
        std::cerr << "Valid net files end with .srn, .app, " 
                        << "(or .net prior to v1.25.2) exiting." << std::endl;
        exit(0);
      }

      if (found && !powerapp)
      {
        usage();
      }

      // determine if it is a PowerApp
      if (ends_with(arg,".app")) 
      {
        powerapp = true;
        found = cnt;
      } 
      else if(!powerapp) 
      {
        found = cnt;
      }
    }
    cnt++;
  }
  return found;
}


class SCIRunKiller : public Runnable
{
  public:
    void run()
    {
      // Wait until we did time out
      Time::waitFor(static_cast<double>(timeout_));
      // Notify user that we are exiting
      std::cout << "TIMEOUT: KILLING SCIRUN PROCESS.";
      std::cout.flush();
      std::cerr.flush();
      // Exit with error code
      Thread::exitAll(1);
    }

    SCIRunKiller(int seconds) :
      timeout_(seconds)
    {}

  private:
    int timeout_;
};


// Services start up... 
void
start_eai() {
  // Create a database of all available services. The next piece of code
  // Scans both the SCIRun as well as the Packages directories to find
  // Services that need to be started. Services allow communication with
  // thirdparty software and are Threads that run asychronicly with
  // with the rest of SCIRun. Since the thirdparty software may be running
  // on a different platform it allows for connecting to remote machines
  // and running the service on a different machine 
  ServiceDBHandle servicedb = new ServiceDB;     
  // load all services and find all makers
  servicedb->loadpackages();
  // activate all services
  servicedb->activateall();
  
  // Services are started and created by the ServiceManager, 
  // which will be launched here
  // Two competing managers will be started, 
  // one for purely internal usage and one that
  // communicates over a socket. 
  // The latter will only be created if a port is set.
  // If the current instance of SCIRun should not provide any services 
  // to other instances of SCIRun over the internet, 
  // the second manager will not be launched
  

  
  IComAddress internaladdress("internal","servicemanager");

// Only build log file if needed for debugging  
#ifdef DEBUG  
  const char *chome = sci_getenv("HOME");
  std::string scidir("");
  if (chome)
    scidir = chome+std::string("/SCIRun/");

  // A log file is not necessary but handy for debugging purposes
  ServiceLogHandle internallogfile = 
    new ServiceLog(scidir+"scirun_internal_servicemanager.log");
  
  ServiceManager* internal_service_manager = 
    new ServiceManager(servicedb, internaladdress, internallogfile); 
#else
  ServiceManager* internal_service_manager = 
    new ServiceManager(servicedb, internaladdress); 
#endif

  Thread* t_int = 
    new Thread(internal_service_manager, "internal service manager",
		  0, Thread::NotActivated);
  t_int->setStackSize(1024*20);
  t_int->activate(false);
  t_int->detach();
  
  
  // Use the following environment setting to switch on IPv6 support
  // Most machines should be running a dual-host stack for the internet
  // connections, so it should not hurt to run in IPv6 mode. In most case
  // ipv4 address will work as well.
  // It might be useful
  std::string ipstr(sci_getenv_p("SCIRUN_SERVICE_IPV6")?"ipv6":"");
  
  // Start an external service as well
  const char *serviceport_str = sci_getenv("SCIRUN_SERVICE_PORT");
  // If its not set in the env, we're done
  if (!serviceport_str) return;
  
  // The protocol for conencting has been called "scirun"
  // In the near future this should be replaced with "sciruns" for
  // a secure version which will run over ssl. 
  
  // A log file is not necessary but handy for debugging purposes

  IComAddress externaladdress("scirun","",serviceport_str,ipstr);

#ifdef DEBUG
  ServiceLogHandle externallogfile = 
    new ServiceLog(scidir+"scirun_external_servicemanager.log"); 
  
  ServiceManager* external_service_manager = 
    new ServiceManager(servicedb,externaladdress,externallogfile); 
#else
  ServiceManager* external_service_manager = 
    new ServiceManager(servicedb,externaladdress); 

#endif

  Thread* t_ext = 
    new Thread(external_service_manager,"external service manager",
		  0, Thread::NotActivated);
  t_ext->setStackSize(1024*20);
  t_ext->activate(false);
  t_ext->detach();
}  

int
MainThread::run_tcl()
{
// ----------------------------------------------------
// Generate the TCL interface and Network Editor

  std::cerr << "create network editor" << std::endl;
  NetworkEditor* editor = new NetworkEditor(net_);

  TCLInterface::execute(std::string("global env; set env(TCL_LIBRARY) \"")+sci_getenv("TCL_LIBRARY")+"\"");
  TCLInterface::execute(std::string("global tcl_library; set tcl_library \"")+sci_getenv("TCL_LIBRARY")+"\"");
  TCLInterface::execute(std::string("global env; set env(TK_LIBRARY) \"")+sci_getenv("TK_LIBRARY")+"\"");
  TCLInterface::execute(std::string("global tk_library; set tk_library \"")+sci_getenv("TK_LIBRARY")+"\"");

  std::cerr << "open sciruntcl" << std::endl;
  TCLInterface::execute("package require sciruntcl");    

  std::cerr << "start network editor" << std::endl;
  TCLInterface::execute("makeNetworkEditor");

  // If the user doesnt have a .scirunrc file, or it is out of date,
  // provide them with a default one

  std::cerr << "copy and parse .scirunrc" << std::endl;
  if (!sci_getenv_p("SCIRUN_RC_PARSED")) 
  {
    copy_and_parse_scirunrc();
  }

  // Determine if we are loading an app.
  const bool powerapp_p = (startnetno_ && ends_with(argv_[startnetno_],".app"));
  if (!powerapp_p)
  {
    TCLInterface::eval("set PowerApp 0");
    // Wait for the main window to display before continuing the startup.
    TCLInterface::eval("wm deiconify .");
    TCLInterface::eval("tkwait visibility $minicanvas");
    TCLInterface::eval("showProgress 1 0 1");
  }
  else
  { // If loading an app, don't wait.
    TCLInterface::eval("set PowerApp 1");
    if (argv_[startnetno_+1])
    {
      TCLInterface::eval("set PowerAppSession {"+std::string(argv_[startnetno_+1])+"}");
    }
    // Determine which standalone and set splash.
    if(strstr(argv_[startnetno_], "BioTensor"))
    {
      TCLInterface::eval("set splashImageFile $bioTensorSplashImageFile");
      TCLInterface::eval("showProgress 1 3221 1");
    }
    else if(strstr(argv_[startnetno_], "BioFEM"))
    {
      TCLInterface::eval("set splashImageFile $bioFEMSplashImageFile");
      TCLInterface::eval("showProgress 1 426 1");
    }
    else if(strstr(argv_[startnetno_], "BioImage"))
    {
      // Need to make a BioImage splash screen.
      TCLInterface::eval("set splashImageFile $bioImageSplashImageFile");
      TCLInterface::eval("showProgress 1 759 1");
    }
    else if(strstr(argv_[startnetno_], "FusionViewer"))
    {
      // Need to make a FusionViewer splash screen.
      TCLInterface::eval("set splashImageFile $fusionViewerSplashImageFile");
      TCLInterface::eval("showProgress 1 310 1");
    }
  }
      
  packageDB = new PackageDB();
  // load the packages
  packageDB->loadPackage(!sci_getenv_p("SCIRUN_LOAD_MODULES_ON_STARTUP"));  

  if (!powerapp_p)
  {
    TCLInterface::eval("hideProgress");
  }
  
  // Run the license dialog code
  TCLInterface::eval("licenseDialog 1");
  
  // Activate "File" menu sub-menus once packages are all loaded.
  TCLInterface::eval("activate_file_submenus");

  // Test for correctly initialized OpenGLContext
  bool initialized = false;
#if defined(HAVE_X11)
#  if DEBUG
  std::cerr << "Using X11OpenGLContext." << std::endl;
#  endif
  X11OpenGLContext *context = new X11OpenGLContext(0, 0, 0, 10, 10, 0, false);
  initialized = context->initialized();
  delete context;
  
#elif defined(_WIN32)
#  if DEBUG
  std::cerr << "Using Win32OpenGLContext." << std::endl;
#  endif
  Win32OpenGLContext *context = new Win32OpenGLContext(0, 0, 0, 10, 10, false, false);
  initialized = context->initialized();
  delete context;
#endif

  if (!initialized)
  {
    TCLInterface::eval("UpdateGraphicsDriversDialog");
  }

  TCLInterface::start_processing_events();

  //! Unlock the GUI so we can start processing callbacks from other threads
  TCLInterface::unlock();

#ifdef HAVE_OPENGL
  SLIVR::ShaderProgramARB::init_shaders_supported();
#endif

  bool loaded_network = false;
  // load network from an xml file...
 
  NetworkIO* netio = editor->get_networkio();

#ifdef DEBUG
// All objects that have been created up to now, will not be deleted before
// our final check, hence reset counter to start counting from here
  debug_tag_default_number_of_objects();
#endif
  
  if(sci_getenv("SCIRUN_NETFILE") != 0)
  {
    std::cout<< "loading scirun network file: "<<sci_getenv("SCIRUN_NETFILE")<<std::endl;
    netio->load_net(std::string(sci_getenv("SCIRUN_NETFILE")));
    loaded_network = true;
  }
  // Load the Network file specified from the command line
  std::string fname; 

  if (startnetno_)
  {
    fname = std::string(argv_[startnetno_]);
    std::cout<< "loading scirun network file: "<<fname<<std::endl;
    TCLInterface::eval("loadnet {"+fname+"}");
    loaded_network = true;
  }

  if (loaded_network &&
      (sci_getenv_p("SCIRUN_EXECUTE_ON_STARTUP") || 
       sci_getenv_p("SCI_REGRESSION_TESTING"))) 
  {
    TCLInterface::eval("netedit scheduleall");
  }

  return TCLInterface::OK;
}

}

using namespace SCIRun;

int
main(int argc, char *argv[], char **environment) 
{

#ifdef _WIN32
	_set_invalid_parameter_handler(InvalidParameterHandler);
#else
	XInitThreads();
#endif

  // Setup the SCIRun key/value environment
  create_sci_environment(environment, argv[0]);
  sci_putenv("SCIRUN_VERSION", SCIRUN_VERSION);
  sci_putenv("SCIRUN_RCFILE_SUBVERSION", SCIRUN_RCFILE_SUBVERSION);
  sci_putenv("SCIRUN_SVN_REVISION", SCIRUN_SVN_REVISION);
  sci_putenv("SCIRUN_SVN_URL", SCIRUN_SVN_URL);

  // Parse the command line arguments to find a network to execute
  const int startnetno = parse_args( argc, argv );

  if (sci_getenv_p("SCI_REGRESSION_TESTING"))
  {
    std::cout << "CTEST_FULL_OUTPUT" << std::endl;
  }
  SCIRunInit();
  
  // Substitution often fails when the directory does not have a separator at the
  // end. Adding it here if it is not present.
  if (sci_getenv_p("SCIRUN_DATA"))
  {
    std::string datadir =  sci_getenv("SCIRUN_DATA");
#ifdef _WIN32
    convertToUnixPath(datadir);
#endif
    if (datadir[datadir.size()-1] != '/') datadir += "/";
    sci_putenv("SCIRUN_DATA",datadir);
#if DEBUG
    std::cout << "Using SCIRUN_DATA: " << datadir << std::endl;
#endif
  }

// ----------------------------------------------------
// Provide information on the SCIRun version
#if DEBUG
  std::cout << "DEBUG build: ";
#endif
  std::string revision = sci_getenv("SCIRUN_SVN_REVISION");
  std::string version  = sci_getenv("SCIRUN_VERSION");
  std::cout << "SCIRun version " << version << " revision " << revision << std::endl;


#ifndef _WIN32
  // Now split off a process for running external processes (not on windows)
  if (!sci_getenv_p("NO_SCIRUN_EXTERNAL_APPLICATION_INTERFACE"))
  {
	  if (!sci_getenv_p("DISABLE_MATLAB_PROCESS")) 
	  {
		systemcallmanager_ = new SystemCallManager();
		systemcallmanager_->create();
		start_eai();
	  }
  }
#else
  if (!sci_getenv_p("NO_SCIRUN_EXTERNAL_APPLICATION_INTERFACE"))
  {
    start_eai();
  }
#endif

  //! create the EventManager thread.
  Thread *emt = new Thread(new EventManager(), "Event Manager Thread");
  emt->detach();

  Network* net=new Network();

  // Activate the scheduler.  Arguments and return values are meaningless
  Thread* scheduler=new Thread(new Scheduler(net), "Scheduler");
  scheduler->setDaemon(true);
  scheduler->detach();

  // Start the TCL thread, takes care of loading packages and networks
  Thread* tcl= new Thread(new MainThread(argc, argv, net, startnetno),
                         "TCL main event loop",0,Thread::Activated,1024*1024);
  tcl->detach();
        
  // When doing regressions, make thread to kill ourselves after timeout
  if (sci_getenv_p("SCI_REGRESSION_TESTING")) 
  {
    sci_putenv("SCIRUN_GUI_UseGuiFetch","off");
    sci_putenv("SCIRUN_GUI_MoveGuiToMouse","off");

    int seconds = 250;
    int tmp;
    const char *timeout = sci_getenv("SCIRUN_REGRESSION_TESTING_TIMEOUT");
    if (timeout && string_to_int(timeout, tmp)) seconds = tmp;
    
    SCIRunKiller *kill = new SCIRunKiller(seconds);
    Thread *tkill = new Thread(kill, "Regression Testing: kill SCIRun after timeout");
    tkill->detach();
  }

  // For scripting so we can kill SCIRun if it hangs
  if (sci_getenv_p("SCIRUN_TIMEOUT")) 
  {
    const char *timeout = sci_getenv("SCIRUN_TIMEOUT");
    if (timeout)
    { 
      int seconds;
      string_to_int(timeout, seconds);
      
      SCIRunKiller *kill = new SCIRunKiller(seconds);
      Thread *tkill = new Thread(kill, "Kill SCIRun after timeout");
      tkill->detach();
    }
  }

#ifdef _WIN32
  // windows has a semantic problem with atexit(), so we wait here instead.
  HANDLE forever = CreateSemaphore(0,0,1,"forever");
  WaitForSingleObject(forever,INFINITE);
#endif

  Semaphore wait("main wait", 0);
  wait.down();

  return 0;
}




