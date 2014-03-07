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
//    File   : WXOpenGLContext.cc
//    Author : David Brayford
//    Date   : July 2007

#include <Applications/Seg3D/WXOpenGLContext.h>

#include <Core/Events/EventManager.h>
#include <Core/Events/BaseEvent.h>
#include <Core/Events/keysyms.h>
#include <slivr/ShaderProgramARB.h>


using namespace SCIRun;

//TODO - add shared context

BEGIN_EVENT_TABLE(WXOpenGLContext, wxGLCanvas)
    //EVT_SIZE(WXOpenGLContext::OnSize) // wxwidgets wiki -> this handler buggy
    EVT_LEFT_DOWN(WXOpenGLContext::OnMousePress)
    EVT_MIDDLE_DOWN(WXOpenGLContext::OnMousePress)
    EVT_RIGHT_DOWN(WXOpenGLContext::OnMousePress)
    EVT_LEFT_UP(WXOpenGLContext::OnMouseRelease)
    EVT_MIDDLE_UP(WXOpenGLContext::OnMouseRelease)
    EVT_RIGHT_UP(WXOpenGLContext::OnMouseRelease)
    EVT_MOTION(WXOpenGLContext::OnMouseMove )
    EVT_MOUSEWHEEL(WXOpenGLContext::OnMouseWheelMove )
    EVT_ENTER_WINDOW(WXOpenGLContext::OnMouseEnter)
    EVT_LEAVE_WINDOW(WXOpenGLContext::OnMouseLeave)
    EVT_PAINT(WXOpenGLContext::OnPaint)
    EVT_ERASE_BACKGROUND(WXOpenGLContext::OnEraseBackground)
    EVT_CHAR(WXOpenGLContext::OnChar)
END_EVENT_TABLE()


static int context_args[] =
{ WX_GL_RGBA, WX_GL_DOUBLEBUFFER, WX_GL_DEPTH_SIZE, 16, 0 };


WXOpenGLContext::WXOpenGLContext(const std::string& target, wxWindow *parent,
                                 wxWindowID id,
                                 const wxPoint& pos, const wxSize& size,
                                 long style, const wxString& name)
  : wxGLCanvas(parent, id, pos, size, style|wxFULL_REPAINT_ON_RESIZE, name,
               context_args),
    target(target)
{
  wxWindow* frame = parent;
  while (dynamic_cast<wxFrame*>(frame) == 0 && frame->GetParent() != 0)
    frame = frame->GetParent();
  // Could be NULL if the parent is a dialog.
  parentWindow = dynamic_cast<wxFrame*>(frame);
}

WXOpenGLContext::~WXOpenGLContext()
{
}



bool
WXOpenGLContext::make_current()
{
  SetCurrent();

  // this has to be here (or like here), only one thread can
  // call make_current on a single opengl window, since WX does not
  // provide a 'release'.
  static bool shaders_initialized = false;

  if (!shaders_initialized) {
    SLIVR::ShaderProgramARB::init_shaders_supported();    
    shaders_initialized = true;
  }

  return true;
}

void
WXOpenGLContext::release()
{
}

int
WXOpenGLContext::width()
{
  int width;
  GetSize(&width, 0);
  return width;
}

int
WXOpenGLContext::height()
{
  int height;
  GetSize(0, &height);
  return height;
}

void
WXOpenGLContext::swap()
{
  SwapBuffers();
}

bool
WXOpenGLContext::has_shaders()
{
  return(SLIVR::ShaderProgramARB::shaders_supported());
}

bool
WXOpenGLContext::initialized()
{
  return(SLIVR::ShaderProgramARB::initialized());
}

void
WXOpenGLContext::lock()
{
#ifndef _WIN32
  wxMutexGuiEnter();
#endif
}

void
WXOpenGLContext::unlock()
{
#ifndef _WIN32
  wxMutexGuiLeave();
#endif
}

void WXOpenGLContext::OnSize(wxSizeEvent& event)
{
}

void WXOpenGLContext::OnEraseBackground(wxEraseEvent& event)
{
}

// David Brayford

unsigned int WXOpenGLContext::ModifierKeyToSkinner(wxMouseEvent& event)
{
  unsigned int m = 0; 
  // Key Press with function button down
  if( event.ShiftDown() )
  {
    m |= EventModifiers::SHIFT_E;
  }
  if( event.ControlDown() )
  {
    m |= EventModifiers::CONTROL_E;
  }
  if( event.AltDown() )
  {
    m |= EventModifiers::ALT_E;
  }
  //  if ( event.CmdDown() )
  //  {
  //    m |= EventModifiers::M1_E; // Apple key
  //  }
  return m;
}

unsigned int
WXOpenGLContext::MouseDownToSkinner(wxMouseEvent& event, unsigned int mod)
{
  unsigned int state = 0;
  if ( event.LeftIsDown() )
  {
#if defined(__APPLE__)
    if (mod & EventModifiers::ALT_E)
    {
      state |= PointerEvent::BUTTON_2_E;
    }
    else if (mod & EventModifiers::CONTROL_E)
    {
      state |= PointerEvent::BUTTON_3_E;
    }
    else
    {
      state |= PointerEvent::BUTTON_1_E;
    }
#else
    state |= PointerEvent::BUTTON_1_E;

	if (event.AltDown() )
	{
		state |= PointerEvent::BUTTON_2_E;
	}
#endif
  }
  else if ( event.MiddleIsDown() )
  {
    state |= PointerEvent::BUTTON_2_E;
  }
  else if ( event.RightIsDown() )
  {
    state |= PointerEvent::BUTTON_3_E;
  }

  return state;
}


unsigned int
WXOpenGLContext::MouseUpToSkinner(wxMouseEvent& event, unsigned int mod)
{
  unsigned int state = 0;
  if ( event.LeftUp() )
  {
#if defined(__APPLE__)
    if (mod & EventModifiers::ALT_E)
    {
      state |= PointerEvent::BUTTON_2_E;
    }
    else if (mod & EventModifiers::CONTROL_E)
    {
      state |= PointerEvent::BUTTON_3_E;
    }
    else
    {
      state |= PointerEvent::BUTTON_1_E;
    }
#else
    state |= PointerEvent::BUTTON_1_E;
#endif
  }
  else if ( event.MiddleUp() )
  {
    state |= PointerEvent::BUTTON_2_E;
  }
  else if ( event.RightUp() )
  {
    state |= PointerEvent::BUTTON_3_E;
  }

  return state;
}


void WXOpenGLContext::OnMousePress(wxMouseEvent& event)
{
  unsigned int state = PointerEvent::BUTTON_PRESS_E;
  unsigned int modifiers = ModifierKeyToSkinner(event);
  unsigned int button = MouseDownToSkinner(event, modifiers);
  
  if (button == 0) return; // Unknown event;
  state |= button;

  const int x = static_cast<int>( event.GetX() );
  const int y = static_cast<int>( height() - 1 - event.GetY() );
  const unsigned long time = (unsigned long)( event.GetTimestamp() );
  PointerEvent *sci_event = new PointerEvent(state, x, y, "", time);
  sci_event->set_modifiers( modifiers );
  EventManager::add_event(sci_event);
}


void WXOpenGLContext::OnMouseRelease(wxMouseEvent& event)
{
  unsigned int state = PointerEvent::BUTTON_RELEASE_E;
  unsigned int modifiers = ModifierKeyToSkinner(event);
  unsigned int button = MouseUpToSkinner(event, modifiers);
  
  if (button == 0) return; // Unknown event;
  state |= button;

  const int x = static_cast<int>( event.GetX() );
  const int y = static_cast<int>( height() - 1 - event.GetY() );
  const unsigned long time = (unsigned long)( event.GetTimestamp() );
  PointerEvent *sci_event = new PointerEvent(state, x, y, "", time);
  sci_event->set_modifiers( modifiers );
  EventManager::add_event(sci_event);
}


void WXOpenGLContext::OnMouseMove(wxMouseEvent& event)
{
  unsigned int state = PointerEvent::MOTION_E;
  unsigned int modifiers = ModifierKeyToSkinner(event);
  unsigned int button = MouseDownToSkinner(event, modifiers);
  
  //if (button == 0) return; // Unknown event;
  state |= button;

  const int x = static_cast<int>( event.GetX() );
  const int y = static_cast<int>( height() - 1 - event.GetY() );	
  const unsigned long time = (unsigned long)( event.GetTimestamp() );
  PointerEvent *sci_event = new PointerEvent(state, x, y, "", time);
  sci_event->set_modifiers( modifiers );
  EventManager::add_event(sci_event);
}

void WXOpenGLContext::OnMouseWheelMove(wxMouseEvent& event)
{
	unsigned int state = PointerEvent::BUTTON_PRESS_E;

	if( event.GetWheelRotation() > 0 )
	{
		state |= PointerEvent::BUTTON_4_E;
	}
	else
	{
		state |= PointerEvent::BUTTON_5_E;
	}

	const int x = static_cast<int>( event.GetX() );
	const int y = static_cast<int>( height() - 1 - event.GetY() );	
	const unsigned long time = (unsigned long)( event.GetTimestamp() );
	PointerEvent *sci_event = new PointerEvent(state, x, y, "", time);
	sci_event->set_modifiers( ModifierKeyToSkinner(event) );
	EventManager::add_event(sci_event);
}

void WXOpenGLContext::OnMouseEnter(wxMouseEvent& event)
{
  unsigned int state = PointerEvent::MOTION_E;

  if( event.Entering() )
  {
    SetFocus();

    const int x = static_cast<int>( event.GetX() );
    const int y = static_cast<int>( height() - 1 - event.GetY() );	
    const unsigned long time = (unsigned long)( event.GetTimestamp() );
    PointerEvent *sci_event = new PointerEvent(state, x, y, "", time);
    sci_event->set_modifiers( ModifierKeyToSkinner(event) );
    EventManager::add_event(sci_event);
  }
}

void WXOpenGLContext::OnMouseLeave(wxMouseEvent& event)
{
  unsigned int state = PointerEvent::MOTION_E;

  if( event.Leaving() )
  {	
    const int x = static_cast<int>( event.GetX() );
    const int y = static_cast<int>( height() - 1 - event.GetY() );	
    const unsigned long time = (unsigned long)( event.GetTimestamp() );
    PointerEvent *sci_event = new PointerEvent(state, x, y, "", time);
    sci_event->set_modifiers( ModifierKeyToSkinner(event) );
    EventManager::add_event(sci_event);
  }
}

void WXOpenGLContext::closeWindow()
{
  if (parentWindow) 
    parentWindow->Close();
}


void WXOpenGLContext::OnPaint( wxPaintEvent& event )
{
  // Must always be here
  // wxClientDC on windows causes bad hangs and crashes
  // Always use wxPaintDC, doesn't hurt other platforms.
  wxPaintDC dc(this);
	
  // Pass to seg3d to do the drawing - it can draw to opengl
  // and swapbuffers on its own 
  WindowEvent *sci_event = new WindowEvent(WindowEvent::REDRAW_E);
  int time = static_cast<int>( event.GetTimestamp() );
  sci_event->set_time(time);
  
  EventManager::add_event(sci_event);
}

unsigned int WXOpenGLContext::KeyToSkinner(int wxkeycode)
{
	unsigned long keyval = wxkeycode;

	switch ( wxkeycode )
	{
	case WXK_CLEAR:
		keyval = SCIRun_Clear;
	break;
	case WXK_LEFT:
		keyval = SCIRun_Left;
	break;
	case WXK_RIGHT:
		keyval = SCIRun_Right;
	break;
	case WXK_UP:
		keyval = SCIRun_Up;
	break;
	case WXK_DOWN:
		keyval = SCIRun_Down;
	break;
	case WXK_HOME:
		keyval = SCIRun_Home;
	break;
	case WXK_END:
		keyval = SCIRun_End;
	break;
	case WXK_INSERT:
		keyval = SCIRun_Insert;
	break;
	case WXK_DELETE:
		keyval = SCIRun_Delete;
	break;
	case WXK_NEXT:
		keyval = SCIRun_Next;
	break;
	case WXK_SHIFT:
		keyval = 0;
	break;
	case WXK_CONTROL:
		keyval = 0;
	break;
	case WXK_MENU:
		keyval = 0;
	break;
	case WXK_ALT:
		keyval = 0;
	break;
	case WXK_CAPITAL:
                keyval = 2; //caps;
	break;
	case WXK_PAUSE:
		keyval = SCIRun_Pause;
	break;
	case WXK_PRINT:
		keyval = SCIRun_Print;
	break;
	case WXK_SELECT:
		keyval = SCIRun_Select;
	break;
	case WXK_EXECUTE:
		keyval = SCIRun_Execute;
	break;
	case WXK_HELP:
		keyval = SCIRun_Help;
	break;
	case WXK_NUMLOCK:
		keyval = SCIRun_Num_Lock;
	break;
	case WXK_SCROLL:
		keyval = SCIRun_Scroll_Lock;
	break;
	case WXK_BACK:
		keyval = SCIRun_BackSpace;
	break;
	case WXK_RETURN:
		keyval = SCIRun_Return;
	break;
	case WXK_ESCAPE:
		keyval = SCIRun_Escape;
	break;
	case WXK_TAB:
		keyval = SCIRun_Tab;
	break;
	default:
	break;
	}

        return keyval;
}


void WXOpenGLContext::OnChar( wxKeyEvent &event )
{
    unsigned int modifiers = 0;
    if (event.ShiftDown())    modifiers |= EventModifiers::SHIFT_E;
    if (event.ControlDown())  modifiers |= EventModifiers::CONTROL_E;
    if (event.AltDown())      modifiers |= EventModifiers::ALT_E;

    const unsigned int keyval = KeyToSkinner( event.GetKeyCode() );

    KeyEvent *sci_event =
      new KeyEvent(KeyEvent::KEY_PRESS_E, modifiers, keyval);

    int time = static_cast<int>( event.GetTimestamp() );
    sci_event->set_time(time);

    EventManager::add_event(sci_event);

    event.Skip();
}
