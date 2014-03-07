#!/usr/bin/env python
#
#  For more information, please see: http://software.sci.utah.edu
# 
# The MIT License
#
# Copyright (c) 2009 Scientific Computing and Imaging Institute,
# University of Utah.
#
#
# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the 
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
#

from wxPython.glcanvas import wxGLCanvas
from wxPython.wx import *
from OpenGL.GLU import *
from OpenGL.GL import *
import sys,math
import threading
import sr_py
import time


class M3DGLCanvas(wxGLCanvas):
    def __init__(self, parent, id):
	wxGLCanvas.__init__(self, parent, id)
	EVT_PAINT(self, self.OnPaint)
	self.init = 0

        EVT_LEFT_DOWN(self, self.pointer_down)
        EVT_LEFT_UP(self, self.pointer_up)
        EVT_MIDDLE_DOWN(self, self.pointer_down)
        EVT_MIDDLE_UP(self, self.pointer_up)
        EVT_RIGHT_DOWN(self, self.pointer_down)
        EVT_RIGHT_UP(self, self.pointer_up)
        EVT_MOTION(self, self.pointer_motion)

        # set up callbacks for the context operations.
        self.sci_context = sr_py.CallbackOpenGLContext()
        self.sci_context.set_pymake_current_func(self.make_current)
        self.sci_context.set_pyswap_func(self.swap)
        self.sci_context.set_pywidth_func(self.width)
        self.sci_context.set_pyheight_func(self.height)
        self.sci_context.set_pyrelease_func(self.release)

        return


    def make_current(self) :
        wxMutexGuiEnter()
        self.SetCurrent()
        wxMutexGuiLeave()
        return 1

    def swap(self) :
        wxMutexGuiEnter()
        self.SwapBuffers()
        wxMutexGuiLeave()

    def width(self) :
        wxMutexGuiEnter()
        w, h = self.GetClientSize()
        wxMutexGuiLeave()
        return w

    def height(self) :
        wxMutexGuiEnter()
        w, h = self.GetClientSize()
        wxMutexGuiLeave()
        return h

    def release(self) :
        pass

    def OnPaint(self,event):
        dc = wxPaintDC(self) # <<-- this is completely retarted
	if not self.init:
            
	    sr_py.run_viewer_thread(self.sci_context)
            self.init = 1
            
            wxYield()
            #wxYield()
            #time.sleep(1)
	    
	#self.OnDraw()
	return

    def OnDraw(self):
        return
	
    def InitGL(self):
        return

    def pointer_event(self, e, event) :
        e.set_x(event.m_x)
        e.set_y(event.m_y)
        
        mask = e.get_modifiers()
        mask = self.get_sr_py_modifier_mask(event, mask) 
        e.set_modifiers(mask)
        
        state = e.get_pointer_state()
        state, n = self.get_sr_py_pointer_modifier_mask(event, state)
        e.set_pointer_state(state)
        sr_py.add_pointer_event(e)
        
    def pointer_down(self, event) :
        e = sr_py.PointerEvent()
        #e.set_time(long(event.time))
        e.set_pointer_state(sr_py.PointerEvent.BUTTON_PRESS_E)
        self.pointer_event(e, event)

    def pointer_up(self, event) :
        e = sr_py.PointerEvent()
        #e.set_time(long(event.time))
        e.set_pointer_state(sr_py.PointerEvent.BUTTON_RELEASE_E)
        self.pointer_event(e, event)
        
    def pointer_motion(self, event) :
        e = sr_py.PointerEvent()
        #e.set_time(long(event.time))
        e.set_pointer_state(sr_py.PointerEvent.MOTION_E)
        self.pointer_event(e, event)
        

    def get_sr_py_modifier_mask(self, event, mask) :
        if event.m_shiftDown :
            mask  |= sr_py.EventModifiers.SHIFT_E
        if event.m_controlDown :
            mask  |= sr_py.EventModifiers.CONTROL_E
        if event.m_altDown :
            mask  |= sr_py.EventModifiers.ALT_E
        if event.m_metaDown :
            mask  |= sr_py.EventModifiers.M1_E
        return mask

    def get_sr_py_pointer_modifier_mask(self, event, mask) :
        n = 0
        if event.m_leftDown :
            mask  |= sr_py.PointerEvent.BUTTON_1_E
            n = 1
        if event.m_middleDown :
            mask  |= sr_py.PointerEvent.BUTTON_2_E
            n = 2
        if event.m_rightDown :
            mask  |= sr_py.PointerEvent.BUTTON_3_E
            n = 3
        return mask, n
