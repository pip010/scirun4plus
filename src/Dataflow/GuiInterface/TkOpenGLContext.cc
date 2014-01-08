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
 *  TkOpenGLContext.cc:
 *
 *  Written by:
 *   McKay Davis
 *   December 2004
 *
 */

#include <slivr/ShaderProgramARB.h>
#include <Core/Util/StringUtil.h>
#include <Core/Datatypes/Color.h>
#include <Core/Exceptions/InternalError.h>

#include <Core/Thread/Thread.h>
#include <Core/Thread/Mutex.h>
#include <Core/Util/Assert.h>

#include <Dataflow/GuiInterface/TCLInterface.h>
#include <Dataflow/GuiInterface/TkOpenGLContext.h>

#include <sci_gl.h>
#include <sci_glx.h>

#include <sci_defs/opengl_defs.h>

#include <iostream>
#include <set>

#include <tk.h>

#ifdef _WIN32
#  include <tkWinInt.h>
#  include <tkWinPort.h>
#  include <X11\XUtil.h>
#endif

#ifdef _WIN32
#include <tkWinInt.h>
#include <tkIntPlatDecls.h>
#include <windows.h>
#include <sstream>
#include <process.h>
#endif

using namespace SCIRun;

  
extern "C" Tcl_Interp* the_interp;

// TODO: Make this threadsafe
// Threaded code should **NOT** have static variables

#ifndef _WIN32
  static GLXContext first_context = 0;
#else
  static HGLRC first_context = 0;
#endif

  namespace SCIRun
  {
    class TkOpenGLContextPrivate
    {
    public:
      TkOpenGLContextPrivate() : 
          display_(0),
            tkwin_(0),
            mainwin_(0),
            x11_win_(0),
            context_(0),
            vi_(0) {}
#ifdef _WIN32
          Window tkglCallback(TkOpenGLContext* ctx, Tk_Window tkwin, Window parent);
#endif

    public:
      Display *         display_; /* X's token for the window's display. */
      Tk_Window         tkwin_;
      Tk_Window         mainwin_;
      Window            x11_win_;
#ifndef _WIN32
      GLXContext        context_;
#else
      HGLRC             context_;
      HDC               hDC_;
      HWND              hWND_;
#endif
      XVisualInfo*      vi_;
      Colormap          colormap_;
    };
  }

#ifdef _WIN32
HINSTANCE dllHINSTANCE=0;

void
PrintErr(char* func_name)
{
  LPVOID lpMsgBuf;
  DWORD dw = GetLastError(); 
  
  if (dw) {
    FormatMessage(
		  FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		  FORMAT_MESSAGE_FROM_SYSTEM,
		  0,
		  dw,
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		  (LPTSTR) &lpMsgBuf,
		  0, 0 );
    
    fprintf(stderr, 
	    "%s failed with error %ld: %s", 
	    func_name, dw, (char*)lpMsgBuf); 
    LocalFree(lpMsgBuf);
  }
  SetLastError(0);
}

static Window
TkGLMakeWindow(Tk_Window tkwin, Window parent, ClientData data)
{
  TkOpenGLContext *ctx = reinterpret_cast<TkOpenGLContext*>(data);
  return ctx->impl_->tkglCallback(ctx, tkwin, parent);
}

char *tkGlClassName = "TkGL";
bool tkGlClassInitialized = false;
WNDCLASS tkGlClass;
WNDPROC tkWinChildProc=0;

// win32 event loop
LRESULT CALLBACK 
WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{

  LONG result;
  WNDCLASS childClass;

  switch (msg) 
  {
  case WM_WINDOWPOSCHANGED:
    /* Should be processed by DefWindowProc, otherwise a double buffered
       context is not properly resized when the corresponding window is resized.*/
    result = TRUE;
    break;
    
  case WM_NCCREATE:
    result = TRUE;
    break;

  case WM_DESTROY:
    result = TRUE;
    break;
    
  default:
    {
      if (tkWinChildProc == 0) 
      {
        GetClassInfo(Tk_GetHINSTANCE(),TK_WIN_CHILD_CLASS_NAME,
          &childClass);
        tkWinChildProc = childClass.lpfnWndProc;
      }
      result = tkWinChildProc(hWnd, msg, wParam, lParam);
    }
  }
  Tcl_ServiceAll();
  return result;
}

Window
TkOpenGLContextPrivate::tkglCallback(TkOpenGLContext* ctx, Tk_Window tkwin, Window parent) 
{
  HWND parentWin;
  int style;
  HINSTANCE hInstance;

  vi_ = new XVisualInfo();
  display_ = Tk_Display(tkwin_);

  hInstance = Tk_GetHINSTANCE();

  // next register our own window class.... 
  
  // TODO: THread safety: This code needs a mutex
  
  if (!tkGlClassInitialized) 
  {
    tkGlClassInitialized = true;
    tkGlClass.style = CS_HREDRAW | CS_VREDRAW;// | CS_OWNDC;
    tkGlClass.cbClsExtra = 0;
    tkGlClass.cbWndExtra = 0;
    tkGlClass.hInstance = dllHINSTANCE;
    tkGlClass.hbrBackground = 0;
    tkGlClass.lpszMenuName = 0;
    //tkGlClass.lpszClassName = TK_WIN_CHILD_CLASS_NAME;
    //tkGlClass.lpfnWndProc = TkWinChildProc;
    tkGlClass.lpszClassName = "TkGL";
    tkGlClass.lpfnWndProc = WndProc;
    tkGlClass.hIcon = 0;
    tkGlClass.hCursor = 0;

    if (RegisterClass(&tkGlClass) == 0)
      PrintErr("MakeWindow RegisterClass");
  }

  /*
   * Create the window, then ensure that it is at the top of the
   * stacking order.
   */

  int x = Tk_X(tkwin), y = Tk_Y(tkwin), width = Tk_Width(tkwin), height = Tk_Height(tkwin);
  if (width == 0 || height == 0) 
  {
    style = WS_POPUP;
    parentWin = 0;
  }
  else 
  {
    style = WS_CHILD | WS_CLIPCHILDREN | WS_CLIPSIBLINGS;
    parentWin = Tk_GetHWND(parent);
  }

  hWND_ = CreateWindow(tkGlClassName, "SCIRun GL Viewer Screen",
			      style, x, y, width, height,
			      parentWin, 0, dllHINSTANCE, 0);
  PrintErr("CreateWindow");

  if (width != 0 && height != 0)
    SetWindowPos(hWND_, HWND_TOP, 0, 0, 0, 0, SWP_NOACTIVATE | SWP_NOMOVE | SWP_NOSIZE);

  hDC_ = GetDC(hWND_);
  
  /* Just for portability, define the simplest visinfo */
  vi_->visual = DefaultVisual(display_, DefaultScreen(display_));   
  vi_->depth = vi_->visual->bits_per_rgb;
  /*
  * find a colormap
  */

  ctx->screen_number_ = Tk_ScreenNumber(tkwin_);

  colormap_ = DefaultColormap(display_, ctx->screen_number_);

  int result = Tk_SetWindowVisual(tkwin_, vi_->visual, 
				  vi_->depth,colormap_ );
  if (result != 1) throw new InternalError("Cannot set Tk Window Visual", 
					      __FILE__, __LINE__);

  SelectPalette(hDC_, ((TkWinColormap *)colormap_)->palette, TRUE);
  RealizePalette(hDC_);
  //visualid_ = iPixelFormat;

  Window win = Tk_AttachHWND(tkwin_, hWND_);
  XMapWindow(display_, win);
  return win;

}

#endif // _WIN32

TkOpenGLContext::TkOpenGLContext(const std::string &id, int visualid, 
				 int width, int height)
  : OpenGLContext(),
    impl_(new TkOpenGLContextPrivate),
    visualid_(visualid),
    screen_number_(0),
    id_(id)
{

  TCLInterface::obtain_tcl_pause();

  // Get the Main Window from TK to find out which display screen we are using
  impl_->mainwin_ = Tk_MainWindow(the_interp);
  if (!impl_->mainwin_)
  {
    throw new InternalError("Cannot find main Tk window",__FILE__,__LINE__);
  }
  
  impl_->display_ = Tk_Display(impl_->mainwin_);
  if (!impl_->display_) 
  {
    throw new InternalError("Cannot find X Display", __FILE__, __LINE__);
  }

  screen_number_ = Tk_ScreenNumber(impl_->mainwin_);

  if (valid_visuals_.empty()) 
  {
    listvisuals();
  }
  
  if (visualid < 0 || visualid >= static_cast<int>(valid_visuals_.size()))
  {
    visualid_ = 0;
  } 
  else 
  {
    visualid_ = valid_visuals_[visualid];
  }

  if (visualid_)
  {
    int n;
    XVisualInfo temp_vi;
    temp_vi.visualid = visualid_;
    impl_->vi_ = XGetVisualInfo(impl_->display_, VisualIDMask, &temp_vi, &n);
    if(!impl_->vi_ || n!=1) 
    {
      throw new InternalError("Cannot find Visual ID #"+to_string(visualid_), __FILE__, __LINE__);
    }
  }

  impl_->tkwin_ = Tk_CreateWindowFromPath(the_interp, impl_->mainwin_, id_.c_str(),0);

  if (!impl_->tkwin_) 
  {
    throw new InternalError("Cannot create Tk Window", __FILE__, __LINE__);
  }


#ifdef _WIN32

  Tk_ClassProcs *procsPtr;
  procsPtr = (Tk_ClassProcs*) Tcl_Alloc(sizeof(Tk_ClassProcs));
  procsPtr->size             = sizeof(Tk_ClassProcs);
  procsPtr->createProc       = TkGLMakeWindow;
  procsPtr->worldChangedProc = 0;
  procsPtr->modalProc        = 0;
  Tk_SetClassProcs(impl_->tkwin_,procsPtr,(ClientData)this); 

  Tk_GeometryRequest(impl_->tkwin_, width, height);
  Tk_ResizeWindow(impl_->tkwin_, width, height);

#else

  // We need vi_ for the definition of the colormap
  if (!impl_->vi_) 
  {
    throw new InternalError("Cannot find Visual", __FILE__, __LINE__);
  }
  
  impl_->colormap_ = XCreateColormap(impl_->display_, Tk_WindowId(impl_->mainwin_), impl_->vi_->visual, AllocNone);

  Tk_GeometryRequest(impl_->tkwin_, width, height);
  Tk_ResizeWindow(impl_->tkwin_, width, height);

  int result = Tk_SetWindowVisual(impl_->tkwin_, impl_->vi_->visual, impl_->vi_->depth, impl_->colormap_);
  if (result != 1) 
  {
    throw new InternalError("Tk_SetWindowVisual failed",__FILE__,__LINE__);
  }
  

#endif  

  Tk_MakeWindowExist(impl_->tkwin_);
  if (Tk_WindowId(impl_->tkwin_) == 0) 
  {
    throw InternalError("Tk_MakeWindowExist failed", __FILE__, __LINE__);
  }

  impl_->x11_win_ = Tk_WindowId(impl_->tkwin_);
  if (!impl_->x11_win_) 
  {
    throw new InternalError("Cannot get Tk X11 window ID",__FILE__,__LINE__);
  }
  
  TCLInterface::release_tcl_pause();

  XSync(impl_->display_, False);
}

void TkOpenGLContext::restackWindow()
{
  Tk_RestackWindow(impl_->tkwin_, Above, 0);
}

void TkOpenGLContext::define_cursor(const std::string& name)
{
  Tk_Window &tkwin = impl_->tkwin_;
  Tk_DefineCursor(tkwin, Tk_GetCursor(the_interp, tkwin, ccast_unsafe(name)));
}

void TkOpenGLContext::xsync()
{
  XSync(impl_->display_, 0);
}

// TODO - NEED TO IMPLEMENT THIS PROPERLY
TkOpenGLContext::~TkOpenGLContext()
{
//   lock();
//   release();
// #ifdef _WIN32
//   if (context_ != first_context)
//     wglDeleteContext(context_);

//   // Remove the TkGL Context from the window so it won't delete it in
//   // the callback.
//   SetWindowLong(hWND_, 0, (LONG) 0);
//   DestroyWindow(hWND_);
// #else // X11
//   glXDestroyContext(display_, context_);
//   XSync(display_, False);
// #endif
//   unlock();
}

void
TkOpenGLContext::lock()
{
  TCLInterface::lock();
}

void
TkOpenGLContext::unlock()
{
  TCLInterface::unlock();
}

bool
TkOpenGLContext::make_current()
{
  bool result = true;

  TCLInterface::obtain_tcl_pause();

  XSync(impl_->display_, False);

  try
  {
#ifdef _WIN32
    if (!impl_->context_)
    {
      // Lock the TCLThread, so it will not communicate with the Windows manager
      
      // describe and set the pixel format
      DWORD dwFlags = PFD_DRAW_TO_WINDOW | PFD_STEREO | PFD_SUPPORT_OPENGL | 
                      PFD_DOUBLEBUFFER | PFD_GENERIC_ACCELERATED;

      // WARNING, if any of the following fields are bad, it can cause the opengl driver
      // to (bizarrely) revert back to version (the generic) 1.1 (instead of, eg, 2.0.X).
      PIXELFORMATDESCRIPTOR pfd = 
      { 
        sizeof(PIXELFORMATDESCRIPTOR),  
        1,                     // version number 
        dwFlags,
        PFD_TYPE_RGBA,         // RGBA type 
        24, // color depth
        8, 0, 8, 0, 8, 0,  // color bits  
        0, 0,  // alpha buffer 
        0,// accumulation buffer 
        0, 0, 0, 0,// accum bits 
        32,  // 32-bit z-buffer 
        0,// no stencil buffer 
        0, // no auxiliary buffer 
        PFD_MAIN_PLANE,        // main layer 
        0,                     // reserved 
        0, 0, 0                // layer masks ignored 
      }; 

      int iPixelFormat = ChoosePixelFormat(impl_->hDC_, &pfd);
      if (iPixelFormat == 0)
      { 
        throw new InternalError("Could not choose PixelFormat",__FILE__,__LINE__);
      }
    
      if (!SetPixelFormat(impl_->hDC_, iPixelFormat, &pfd))
      { 
        throw new InternalError("Could not set PixelFormat",__FILE__,__LINE__);
      }

      /* Get the actual pixel format */
      if (!DescribePixelFormat(impl_->hDC_, iPixelFormat, sizeof(pfd), &pfd))
      { 
        throw new InternalError("Could not get PixelFormat",__FILE__,__LINE__);
      }

      if ((impl_->context_ = wglCreateContext(impl_->hDC_)) == 0)
      { 
        throw new InternalError("Cannot create OpenGL Context",__FILE__,__LINE__);
      }

      if (!first_context) 
      {
        if ((first_context = wglCreateContext(impl_->hDC_)) == 0) 
        {
          throw new InternalError("Cannot create default OpenGL Context",__FILE__,__LINE__);
        }
      }

      if (!(wglShareLists(first_context,impl_->context_))) 
      {
        throw new InternalError("OpenGL Context sharing is broken",__FILE__,__LINE__);
      }
    }
    
    HGLRC current = wglGetCurrentContext();
    
    if (current != impl_->context_) 
    {
      result = wglMakeCurrent(impl_->hDC_,impl_->context_);
      if (!result)
      {
        throw new InternalError("wglMakeCurrent failed to retrieve current context",__FILE__,__LINE__);
      }
    }
#else

    // Make current is being called by the thread that runs the OpenGL rendering
    // This thread needs to generate the context to in order to follow proper 
    // OpenGL threading rules
    
    // TODO: Need to make this thread safe                                                 
    if (!impl_->context_)
    {
      if (!first_context) 
      {
        first_context = glXCreateContext(impl_->display_, impl_->vi_, 0, 1);
      }
      impl_->context_ = glXCreateContext(impl_->display_, impl_->vi_, first_context, 1);

      if (!impl_->context_)
      { 
        throw new InternalError("Cannot create GLX Context",__FILE__,__LINE__);
      }
    }
    
    result = glXMakeCurrent(impl_->display_, impl_->x11_win_, impl_->context_);

    if (!result)
    {
      throw new InternalError("Cannot create GLX Context current",__FILE__,__LINE__);
    }
#endif

    // Unlock the TCLThread
    TCLInterface::release_tcl_pause();
  }
  catch(...)
  {
    std::cerr  << "TkOpenGLContext error\n";
    TCLInterface::release_tcl_pause();
  }


  return result;
}


void
TkOpenGLContext::release()
{
  // A TCL pause will force the TCL Thread to be locked so it cannot call the
  // windowing system. This will force the wgl and glX functions to be called
  // thread safe
  TCLInterface::obtain_tcl_pause();

  XSync(impl_->display_, False);

#ifdef _WIN32
  if (wglGetCurrentContext() != 0)
    wglMakeCurrent(0,0);
#else 
  glXMakeCurrent(impl_->display_, None, 0);
#endif

  // Allow TCL Thread to continue
  TCLInterface::release_tcl_pause();
}


int
TkOpenGLContext::width()
{
  return Tk_Width(impl_->tkwin_);
}


int
TkOpenGLContext::height()
{
  return Tk_Height(impl_->tkwin_);
}


void
TkOpenGLContext::swap()
{  
  // TODO:  Thread safe?
#ifdef _WIN32
  SwapBuffers(impl_->hDC_);
#else // X11
  glXSwapBuffers(impl_->display_, impl_->x11_win_);
#endif
}

#define GETCONFIG(attrib) \
if(glXGetConfig(display, &vinfo[i], attrib, &value) != 0){\
  std::cerr << "Error getting attribute: " << #attrib << std::endl; \
  TCLInterface::unlock(); \
  return std::string(""); \
}

bool
TkOpenGLContext::has_shaders()
{
  return(SLIVR::ShaderProgramARB::shaders_supported());
}


bool
TkOpenGLContext::initialized()
{
#ifdef HAVE_OPENGL
  return(SLIVR::ShaderProgramARB::initialized());
#else
  return (true);
#endif
}


std::string
TkOpenGLContext::listvisuals()
{
#ifdef _WIN32
  valid_visuals_.clear();
  return "";
#else // X11
  Tk_Window topwin=Tk_MainWindow(the_interp);
  if(!topwin)
  {
    std::cerr << "Unable to locate main window!\n";
    TCLInterface::unlock();
    return std::string("");
  }

  Display *display = Tk_Display(topwin);
  int screen = Tk_ScreenNumber(topwin);
  valid_visuals_.clear();
  std::vector<std::string> visualtags;
  std::vector<int> scores;
  int nvis;
  XVisualInfo* vinfo=XGetVisualInfo(display, 0, 0, &nvis);
  if(!vinfo)
  {
    std::cerr << "XGetVisualInfo failed";
    TCLInterface::unlock();
    return std::string("");
  }

  for (int i=0;i<nvis;i++)
  {
    int score=0;
    int value;
    GETCONFIG(GLX_USE_GL);
    if(!value)
      continue;
    GETCONFIG(GLX_RGBA);
    if(!value)
      continue;
    GETCONFIG(GLX_LEVEL);
    if(value != 0)
      continue;
    if(vinfo[i].screen != screen)
      continue;

    std::string tag = "id=" + to_string(static_cast<unsigned int>(vinfo[i].visualid));
    valid_visuals_.push_back(vinfo[i].visualid);

    GETCONFIG(GLX_DOUBLEBUFFER);
    if (value)
    {
      score+=200;
      tag += "double, ";
    }
    else
    {
      tag += "single, ";
    }
    GETCONFIG(GLX_STEREO);
    if (value)
    {
      score+=1;
      tag += "stereo, ";
    }
    tag += "rgba=";
    GETCONFIG(GLX_RED_SIZE);
    tag+=to_string(value)+":";
    score+=value;
    GETCONFIG(GLX_GREEN_SIZE);
    tag+=to_string(value)+":";
    score+=value;
    GETCONFIG(GLX_BLUE_SIZE);
    tag+=to_string(value)+":";
    score+=value;
    GETCONFIG(GLX_ALPHA_SIZE);
    tag+=to_string(value);
    score+=value;
    GETCONFIG(GLX_DEPTH_SIZE);
    tag += ", depth=" + to_string(value);
    score+=value*5;
    GETCONFIG(GLX_STENCIL_SIZE);
    score += value * 2;
    tag += ", stencil="+to_string(value);
    tag += ", accum=";
    GETCONFIG(GLX_ACCUM_RED_SIZE);
    tag += to_string(value) + ":";
    GETCONFIG(GLX_ACCUM_GREEN_SIZE);
    tag += to_string(value) + ":";
    GETCONFIG(GLX_ACCUM_BLUE_SIZE);
    tag += to_string(value) + ":";
    GETCONFIG(GLX_ACCUM_ALPHA_SIZE);
    tag += to_string(value);
    tag += ", score=" + to_string(score);
    
    visualtags.push_back(tag);
    scores.push_back(score);
  }
  for(int i=0; i < int(scores.size())-1; i++)
  {
    for(int j=i+1; j < int(scores.size()); j++)
    {
      if(scores[i] < scores[j])
      {
        std::swap(scores[i], scores[j]);
        std::swap(visualtags[i], visualtags[j]);
        std::swap(valid_visuals_[i], valid_visuals_[j]);
      }
    }
  }
  std::string ret_val;
  for (unsigned int k = 0; k < visualtags.size(); ++k) 
  {
    ret_val = ret_val + "{" + visualtags[k] +"} ";
  }
  return ret_val;
#endif
}

#ifdef _WIN32
BOOL WINAPI DllMain(HINSTANCE hinstance, DWORD reason, LPVOID reserved)
{
  switch (reason) {
  case DLL_PROCESS_ATTACH:
    dllHINSTANCE = hinstance; break;
  default: break;
  }
  return TRUE;
}
#endif
