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

#include <tcl.h>
#include <tk.h>
#include <itcl.h>
#include <itk.h>

#include <Core/Util/Environment.h>
#include <Core/Util/Debug.h>

#include <Dataflow/GuiInterface/TCLInterface.h>
#include <Dataflow/GuiInterface/TCLThread.h>

#include <Dataflow/GuiInterface/share.h>

#include <iostream>
#include <stdlib.h>

namespace SCIRun {

static
int
x_error_handler(Display* dpy, XErrorEvent* error)
{
#ifndef _WIN32
    char msg[200];
    XGetErrorText(dpy, error->error_code, msg, 200);
    std::cerr << "X Error: " << msg << std::endl;
    abort();
#endif
    return 0; // Never reached...
}

extern "C" SCISHARE Tcl_Interp* the_interp;

extern "C" SCISHARE Tcl_Interp* the_interp;

extern "C" SCISHARE int OpenGLCmd _ANSI_ARGS_((ClientData clientData,
                                      Tcl_Interp *interp, 
                                      int argc, CONST84 char **argv));
extern "C" SCISHARE int BevelCmd _ANSI_ARGS_((ClientData clientData,
                                              Tcl_Interp *interp, 
                                              int argc, CONST84 char **argv));
extern "C" SCISHARE int Tk_RangeObjCmd _ANSI_ARGS_((ClientData clientData, 
                                                    Tcl_Interp *interp, 
                                                    int objc, Tcl_Obj *CONST objv[])); 
extern "C" SCISHARE int Tk_CursorCmd _ANSI_ARGS_((ClientData clientData,
                                                  Tcl_Interp *interp, 
                                                  int argc, CONST84 char **argv));
extern "C" SCISHARE int BLineInit _ANSI_ARGS_((void));
extern "C" int Blt_SafeInit _ANSI_ARGS_((Tcl_Interp *interp));
extern "C" int Blt_Init _ANSI_ARGS_((Tcl_Interp *interp));

static
int
exitproc(ClientData, Tcl_Interp*, int, CONST84 char* [])
{
  Thread::exitAll(0);
  return TCL_OK; // not reached
}

static
int
objectlist(ClientData, Tcl_Interp*, int, CONST84 char* [])
{
  debug_print_objects();
  return TCL_OK;
}

static TCLThread* tcl_thread_ = 0;

int
appInit(Tcl_Interp *interp)
{
  the_interp = interp;

  if (Tcl_Init(the_interp) == TCL_ERROR) 
  {
    std::cerr << "\nTcl_Init() failed with TCL_ERROR="<< TCL_ERROR <<std::endl;
    return TCL_ERROR;
  }

  if(Tk_Init(the_interp) == TCL_ERROR)
  {    
    std::cerr << "Tk_Init failed: " << Tcl_GetStringResult(the_interp) << std::endl;
#ifndef WIN32
    // If the error message doesn't mention DISPLAY, we should.
    if(std::string(Tcl_GetStringResult(the_interp)).find("DISPLAY") == std::string::npos)
    {
      std::cerr << "In most case, this means your DISPLAY environment "
                << "variable is not set correctly." << std::endl;
    }
#endif
    Thread::exitAll(-1);
  }

  Tcl_StaticPackage(the_interp, "Tk", Tk_Init, Tk_SafeInit);

  Tcl_CreateCommand(the_interp, "opengl", OpenGLCmd, 
                    (ClientData) Tk_MainWindow(the_interp), 0);

  Tcl_CreateCommand(the_interp, "bevel", BevelCmd,
                    (ClientData) Tk_MainWindow(the_interp), 0);

  Tcl_CreateObjCommand(the_interp, "range", (Tcl_ObjCmdProc *)Tk_RangeObjCmd,
                       (ClientData) 0, 0);

  Tcl_CreateCommand(the_interp, "cursor", Tk_CursorCmd,
                    (ClientData) Tk_MainWindow(the_interp), 0);

  BLineInit();

  Tcl_CreateCommand(the_interp, "exit", exitproc, 0, 0);
  Tcl_CreateCommand(the_interp, "objectlist",  objectlist, 0, 0);  
  
  return (tcl_thread_->run_tcl());
}

TCLThread::TCLThread(int argc, char** argv) :
  argc_(argc), argv_(argv)
{
  // Setup the error handler to catch errors...
  // The default one exits, and makes it very hard to 
  // track down errors.  We need core dumps!
  XSetErrorHandler(x_error_handler);
}

TCLThread::~TCLThread()
{
}

void
TCLThread::run()
{
  //! Mark that we are inside the thread running TCL
  //! This will disable the GUI locks() called on this thread
  //! Instead the TCLInterface::lock will be treated as a counter to check
  //! availability to execute jobs from other threads
  TCLInterface::mark_as_tcl_thread();

  //! Lock TCL, in case of none experimental thread this blocks threads that 
  //! are trying to send jobs to TCL, in case of the experimental thread, it
  //! will block these jobs when it is receiving them from the mailbox
  TCLInterface::lock();
  
  tcl_thread_ = this;

  if (sci_getenv_p("SCIRUN_NOGUI"))
    Tcl_Main(1, argv_, appInit);
  else
    Tk_Main(1, argv_, appInit);
}

}
