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
//    File   : OSXEventSpawner.cc
//    Author : McKay Davis
//    Date   : Thu Sep  7 20:34:32 2006

#include <sci_defs/bits_defs.h>

#if defined(__APPLE__) && !defined(SCI_64BITS)

#include <Core/Geom/X11Lock.h>
#include <Core/Events/OSXEventSpawner.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/BaseEvent.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/Timer.h>
#include <Core/Util/Assert.h>
#include <Core/Util/Environment.h>
#include <Core/Thread/Semaphore.h>
#include <iostream>

#include <Carbon/Carbon.h>
#include <Core/Events/keysyms.h>

namespace SCIRun {

static pascal OSStatus
window_event_callback(EventHandlerCallRef nextHandler,
               EventRef theEvent,
               void *userData)
{
  ASSERT (GetEventClass(theEvent) == kEventClassWindow);

  // Being lazy here, redraw handles everything, but this could
  // be optimized for performance.
  WindowEvent *event = new WindowEvent(WindowEvent::REDRAW_E);
  event->set_target(*((std::string *)userData));
  EventManager::add_event(event);

  return (CallNextEventHandler(nextHandler, theEvent));
}


static pascal OSStatus
window_close_event_callback(EventHandlerCallRef nextHandler,
                            EventRef theEvent,
                            void *userData)
{
  ASSERT (GetEventClass(theEvent) == kEventClassWindow);
  event_handle_t event = new WindowEvent(WindowEvent::DESTROY_E,
                                         *((std::string *)userData));
  EventManager::add_event(event);
  return 0;
}


static pascal OSStatus
mouse_event_callback(EventHandlerCallRef nextHandler,
                     EventRef theEvent,
                     void *userData)
{
  // Type of mouse event
  unsigned int state = 0;
  switch (GetEventKind(theEvent)) {
  case kEventMouseDown    : state = PointerEvent::BUTTON_PRESS_E; break;
  case kEventMouseUp      : state = PointerEvent::BUTTON_RELEASE_E; break;
  case kEventMouseDragged : state = PointerEvent::MOTION_E; break;
  case kEventMouseMoved   : state = PointerEvent::MOTION_E; break;
  case kEventMouseWheelMoved: {
    EventMouseWheelAxis axis;
    GetEventParameter(theEvent, kEventParamMouseWheelAxis, typeMouseWheelAxis,
                      NULL, sizeof(axis), NULL, &axis);
    // We only handle the Y for now
    if (axis != 1) {
      return (CallNextEventHandler(nextHandler, theEvent));
    }

    static bool down = true;
    // For some reason, kEventMouseWheelMoved events come in pairs, so we'll
    // just assume one is for button down, and the next is for button up
    state = (down ? PointerEvent::BUTTON_PRESS_E :
             PointerEvent::BUTTON_RELEASE_E);
    down = !down;

    SInt32 delta;
    GetEventParameter(theEvent, kEventParamMouseWheelDelta, typeSInt32,
                      NULL, sizeof(delta), NULL, &delta);

    if (delta > 0) {
      state |= PointerEvent::BUTTON_4_E;
    } else if (delta < 0) {
      state |= PointerEvent::BUTTON_5_E;
    }
  } break;

  default: break;
  }


  UInt32 carbon_modifiers = GetCurrentEventKeyModifiers();
  unsigned int modifiers = 0;
  if (carbon_modifiers & 0x0100) modifiers |= EventModifiers::M1_E;//Apple key
  if (carbon_modifiers & 0x0200) modifiers |= EventModifiers::SHIFT_E;
  if (carbon_modifiers & 0x0800) modifiers |= EventModifiers::ALT_E;
  if (carbon_modifiers & 0x1000) modifiers |= EventModifiers::CONTROL_E;

  // Get Buttons
  // TODO, change to mouse chord for multiple buttons
  EventMouseButton button;
  GetEventParameter(theEvent, kEventParamMouseButton, typeMouseButton,
                    NULL, sizeof(EventMouseButton), NULL, &button);

  switch (button) {
  case 1: {
    if (modifiers & EventModifiers::ALT_E) {
      state |= PointerEvent::BUTTON_2_E;
    } else if (modifiers & EventModifiers::M1_E) {
      state |= PointerEvent::BUTTON_3_E;
    } else if (modifiers & EventModifiers::CONTROL_E) {
      state |= PointerEvent::BUTTON_3_E;
    } else {
      state |= PointerEvent::BUTTON_1_E;
    }
  } break;
    // Swap buttons 2 and 3 for mac
  case 3: state |= PointerEvent::BUTTON_2_E; break;
  case 2: state |= PointerEvent::BUTTON_3_E; break;
  case 4: state |= PointerEvent::BUTTON_4_E; break;
  case 5: state |= PointerEvent::BUTTON_5_E; break;
  default: break;
  }


  // Get mouse cursor location
  Point	point;
  GetEventParameter(theEvent, kEventParamMouseLocation, typeQDPoint,
                    NULL, sizeof(Point), NULL, &point);

  WindowPtr window;
  GetEventParameter(theEvent, kEventParamWindowRef, typeWindowRef,
                    NULL, sizeof(Point), NULL, &window);


  RgnHandle win_reg = NewRgn();
  GetWindowRegion(window, kWindowGlobalPortRgn, win_reg);
  Rect rect;
  GetRegionBounds(win_reg, &rect);
  DisposeRgn(win_reg);

  //Create the event
  PointerEvent *sci_event = new PointerEvent();
  sci_event->set_pointer_state(state);
  sci_event->set_x(point.h - rect.left);
  sci_event->set_y(rect.bottom - point.v);



  sci_event->set_modifiers(modifiers);

  // Send it to our own event manager
  sci_event->set_target(*((std::string *)userData));
  EventManager::add_event(sci_event);

  return (CallNextEventHandler(nextHandler, theEvent));
}


static pascal OSStatus
key_event_callback(EventHandlerCallRef nextHandler,
                     EventRef theEvent,
                     void *userData)
{
  KeyEvent *sci_event = new KeyEvent();

  switch (GetEventKind(theEvent)) {
  case kEventRawKeyDown: sci_event->set_key_state(KeyEvent::KEY_PRESS_E);break;
  case kEventRawKeyUp: sci_event->set_key_state(KeyEvent::KEY_RELEASE_E);break;
  default: break;
  }

  UInt32 keycode;
  GetEventParameter(theEvent, kEventParamKeyCode, typeUInt32,
                    NULL, sizeof(keycode), NULL, &keycode);

  char keychar;
  GetEventParameter(theEvent, kEventParamKeyMacCharCodes, typeChar,
                    NULL, sizeof(keychar), NULL, &keychar);


  int keyval = keychar;
  switch (keychar) {
  case 0x0D: keyval = SCIRun_Return; break;
  case 0x09: keyval = SCIRun_Tab; break;
  case 0x08: keyval = SCIRun_BackSpace; break;
  case 0x1C: keyval = SCIRun_Left; break;
  case 0x1D: keyval = SCIRun_Right; break;
  case 0x1E: keyval = SCIRun_Up; break;
  case 0x1F: keyval = SCIRun_Down; break;

  default: break;
  }

  sci_event->set_keyval(keyval);

  UInt32 carbon_modifiers = GetCurrentEventKeyModifiers();
  unsigned int modifiers = 0;
  if (carbon_modifiers & 0x0100) modifiers |= EventModifiers::M1_E;//Apple key
  if (carbon_modifiers & 0x0200) modifiers |= EventModifiers::SHIFT_E;
  if (carbon_modifiers & 0x0800) modifiers |= EventModifiers::ALT_E;
  if (carbon_modifiers & 0x1000) modifiers |= EventModifiers::CONTROL_E;
  sci_event->set_modifiers(modifiers);

  // Send it to our own event manager
  sci_event->set_target(*((std::string *)userData));
  EventManager::add_event(sci_event);

  return (CallNextEventHandler(nextHandler, theEvent));
}


OSXEventSpawner::OSXEventSpawner(const std::string &target,
                                 WindowPtr window) :
  EventSpawner(target)
{
  X11Lock::lock();

  InstallStandardEventHandler(GetWindowEventTarget(window));

  static EventTypeSpec window_flags[] = {
    { kEventClassWindow, kEventWindowUpdate },
    { kEventClassWindow, kEventWindowBoundsChanged }
    /*
      { kEventClassWindow, kEventWindowDrawContent },
      { kEventClassWindow, kEventWindowShown },
      { kEventClassWindow, kEventWindowHidden },
      { kEventClassWindow, kEventWindowActivated },
      { kEventClassWindow, kEventWindowDeactivated },
      { kEventClassWindow, kEventWindowClose },

    */
  };

  EventHandlerUPP window_callback = NewEventHandlerUPP(window_event_callback);
  InstallWindowEventHandler(window,
                            window_callback,
                            GetEventTypeCount(window_flags),
                            window_flags,
                            &target_, 0L);
  DisposeEventHandlerUPP(window_callback);


  static EventTypeSpec mouse_flags[] = {
    { kEventClassMouse, kEventMouseDown },
    { kEventClassMouse, kEventMouseUp },
    { kEventClassMouse, kEventMouseMoved },
    { kEventClassMouse, kEventMouseDragged },
    { kEventClassMouse, kEventMouseWheelMoved }
  };


  EventHandlerUPP mouse_callback = NewEventHandlerUPP(mouse_event_callback);
  InstallWindowEventHandler(window,
                            mouse_callback,
                            GetEventTypeCount(mouse_flags),
                            mouse_flags,
                            &target_, 0L);
  DisposeEventHandlerUPP(mouse_callback);

  static EventTypeSpec key_flags[] = {
    { kEventClassKeyboard, kEventRawKeyDown },
    { kEventClassKeyboard, kEventRawKeyUp }
  };

  EventHandlerUPP key_callback = NewEventHandlerUPP(key_event_callback);
  InstallWindowEventHandler(window,
                            key_callback,
                            GetEventTypeCount(key_flags),
                            key_flags,
                            &target_, 0L);
  DisposeEventHandlerUPP(key_callback);


  static EventTypeSpec close_flags[] = {
    { kEventClassWindow, kEventWindowClose }
  };

  EventHandlerUPP window_close_callback =
    NewEventHandlerUPP(window_close_event_callback);
  InstallWindowEventHandler(window,
                            window_close_callback,
                            GetEventTypeCount(close_flags),
                            close_flags,
                            &target_, 0L);
  DisposeEventHandlerUPP(window_close_callback);


  X11Lock::unlock();
}


OSXEventSpawner::~OSXEventSpawner()
{
}


bool
OSXEventSpawner::iterate()
{
  EventRef theEvent;
  EventTargetRef theTarget = GetEventDispatcherTarget();
  while (ReceiveNextEvent(0,NULL,kEventDurationForever,
                          true, &theEvent) == noErr) {
    SendEventToEventTarget(theEvent,theTarget);
    ReleaseEvent(theEvent);
  }

  return true;
}

}

#endif
