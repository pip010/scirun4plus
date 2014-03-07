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
//    File   : main.cc
//    Author : McKay Davis
//    Date   : Tue May 30 21:38:23 MDT 2006

#include <string.h>
#include <stdlib.h>

#include <sci_debug.h>
#include <sci_defs/x11_defs.h>
#include <wxchar_fix.h>

#include <Core/Util/Environment.h>
#include <Core/Skinner/XMLIO.h>
#include <Core/Skinner/Skinner.h>
#include <Core/Skinner/Window.h>
#include <Applications/Seg3D/Seg3DFrame.h>
#include <Applications/Seg3D/Seg3DVersion.h>
#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

#include <Core/Events/EventManager.h>


#include <iostream>

using namespace SCIRun;

#include <wx/wx.h>
#include <wx/splash.h>
#include <wx/help.h>
#include <wx/stdpaths.h>
#include <wx/xrc/xmlres.h>
#include <wx/fs_zip.h>


// Where to look for the data directory
#define SYSDATAPATH "/usr/share/Seg3D/" SEG3D_VERS_STRING "/data"

// Main entry point for application
class Seg3D : public wxApp 
{
public:
  // Called on application startup
  virtual bool OnInit();
  virtual int OnExit();
  
  Skinner::GLWindow *skinner_window_;
  // Thread that manages all Seg3D events
  Thread *event_manager_thread_;
};


// Needed for WxWidgets
DECLARE_APP(Seg3D)
IMPLEMENT_APP_CONSOLE(Seg3D)

// Function called as the program starts
bool
Seg3D::OnInit()
{
#if DEBUG
    std::cerr << "Starting Seg3D" << std::endl;
#endif

  // Look for the data directory in the current directory and the install path.
  const char *datapath = "./data";
  if (validDir(SYSDATAPATH))
  {
    datapath = SYSDATAPATH;
  }
#ifndef __APPLE__
  else if (!validDir(datapath))
  {
    std::cerr << "Seg3D Error: Could not locate Seg3D data directory.  Exiting.\n";
    exit(0);
  }
#endif
#if DEBUG
    std::cerr << __FILE__ << ", " << __LINE__ << " Seg3D data is located in " << datapath << std::endl;
#endif

// SPASH SCREEN CODE
#ifndef DEBUG 
  // Turn off splash screen on debug - it's too hard to debug with the
  // splash screen that won't go anywhere in the way.
  wxString splashScreenPath;
  splashScreenPath = std2wx(datapath + string("/splash-seg3d.png"));
  
#if defined(__APPLE__)
  wxStandardPathsCF * paths = new wxStandardPathsCF();
  splashScreenPath = paths->GetResourcesDir() + wxT("/splash-seg3d.png");
#endif
    
  wxImage::AddHandler( new wxPNGHandler );
  wxBitmap bitmap;            
  if (bitmap.LoadFile(splashScreenPath, wxBITMAP_TYPE_PNG))
  {
    new wxSplashScreen(bitmap,
                       wxSPLASH_CENTRE_ON_SCREEN|wxSPLASH_TIMEOUT,
                       4000, 0, -1, wxDefaultPosition, wxDefaultSize,
                       wxSIMPLE_BORDER|wxFRAME_NO_TASKBAR);
  }
  wxYield();
#endif
// DONE DOING SPLASH SCREEN

  // SETUP ENVIRONMENT
  // TODO: CLEANUP ENVIRONMENT CONTROL AS SEG3D DOES NOT NEED ALL OF SCIRUN
  create_sci_environment(0, (char *)argv[0]);  
#if DEBUG
  std::cout  << __FILE__ << ", " << __LINE__ << " SCI environment has been created." << std::endl;
#endif

  // Need to make the frame before the Painter (since Painter's tools
  // will refer to it.
  
  wxXmlResource::Get()->InitAllHandlers();
  wxFileSystem::AddHandler(new wxZipFSHandler);
  wxString xrcPath = std2wx(datapath + string("/Seg3D.xrc"));

#if DEBUG
  std::cout << __FILE__ << ", " << __LINE__ << " Handlers have been initialized." << std::endl;
#endif

// Since Seg3D is an App on the mac we need to deal with this seperately
#if defined(__APPLE__)
  wxStandardPathsCF * paths2 = new wxStandardPathsCF();
  xrcPath = paths2->GetResourcesDir() + wxT("/Seg3D.xrc");
#endif
  
  //
  if (!wxXmlResource::Get()->Load(xrcPath)) {
    std::cerr << "Seg3D Error: Could not load Seg3D.xrc from " << xrcPath 
	      << ".  Exiting.\n" << std::endl;
    exit(0);
  }

  // Main WXWidgets frame
  Seg3DFrame* frame  = 
    new Seg3DFrame("Seg3D", NULL, _T("Seg3D"),
                   wxDefaultPosition, wxSize(900+PANEL_WIDTH, 750));
#if DEBUG
  std::cout << __FILE__ << ", " << __LINE__ << " Seg3DFrame has been instantiated." << std::endl;
#endif

  // Store Seg3D in global static variable
  Painter::global_seg3dframe_pointer_ = frame;
  Skinner::XMLIO::register_maker<Painter>();
  Skinner::Skinner *skinner = Skinner::load_default_skin(datapath);
  if (!skinner) {
    std::cerr << "Errors encounted loading default skin.\n";
    return false;
  }
#if DEBUG
  std::cout << __FILE__ << ", " << __LINE__ << " Skinner has loaded the default skin." << std::endl;
#endif

  Skinner::GLWindow *skinner_window = skinner->get_windows()[0];
  
  skinner_window_ = skinner_window;
  // Attach the skinner window to it's opengl context.
  skinner_window->attach_to_context(frame->getContext());
#if DEBUG
  std::cout << __FILE__ << ", " << __LINE__ << " Skinner window has attached it's OpenGL context." << std::endl;
#endif

  EventManager *em = new EventManager();
  event_manager_thread_ = new Thread(em, "Event Manager");
#if DEBUG
  std::cout << __FILE__ << ", " << __LINE__ << " Event manager thread started." << std::endl;
#endif

#if defined(__APPLE__)
  frame->Move(wxDefaultPosition);
  Painter::ThrowSkinnerSignal("Painter::Autoview");
#endif
    
  return true;
}


int
Seg3D::OnExit()
{
  delete skinner_window_;
  EventManager::add_event( new QuitEvent() );
  event_manager_thread_->join();
  return wxApp::OnExit();
}


