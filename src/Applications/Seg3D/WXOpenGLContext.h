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
//    File   : WXOpenGLContext.h
//    Author : David Brayford
//    Date   : May 2008

#ifndef CORE_EVENTS_WIN32EVENTSPAWNER_H
#define CORE_EVENTS_WIN32EVENTSPAWNER_H

#include <Core/Geom/OpenGLContext.h>

#include <wx/wx.h>
#include <wx/glcanvas.h>

#if !wxUSE_GLCANVAS
#error Need WX GL Canvas!
#endif

#include <Applications/Seg3D/share.h>

namespace SCIRun {

class WXOpenGLContext : public wxGLCanvas, public OpenGLContext
{
public:
  WXOpenGLContext(const string& target,
                  wxWindow *parent, wxWindowID id = wxID_ANY,
                  const wxPoint& pos = wxDefaultPosition,
                  const wxSize& size = wxDefaultSize, long style = 0,
                  const wxString& name = wxT("WXOpenGLContext"));

  ~WXOpenGLContext();

  // To make compliant with the existing OpenGLContext.
  virtual bool			make_current();
  virtual void			release();
  virtual int			width();
  virtual int			height();
  virtual void			swap();
  virtual bool			has_shaders();
  virtual bool			initialized();
  virtual void                  lock();
  virtual void                  unlock();

  // Called when Seg3D wants to shut the window down - it (currently)
  // only has access to the GL portion.
  void closeWindow();

protected:
  void OnPaint(wxPaintEvent& event);
  void OnSize(wxSizeEvent& event);
  void OnEraseBackground(wxEraseEvent& event);
  void OnMousePress(wxMouseEvent& event);
  void OnMouseRelease(wxMouseEvent& event);
  void OnMouseMove(wxMouseEvent& event);
  void OnMouseWheelMove(wxMouseEvent& event);
  void OnMouseEnter(wxMouseEvent& event);
  void OnMouseLeave(wxMouseEvent& event);
  void OnChar(wxKeyEvent &event);

  unsigned int ModifierKeyToSkinner(wxMouseEvent& event);
  unsigned int MouseDownToSkinner(wxMouseEvent& event, unsigned int mods);
  unsigned int MouseUpToSkinner(wxMouseEvent& event, unsigned int mods);
  unsigned int KeyToSkinner(int wxkeycode);

private:
  wxFrame* parentWindow;
  std::string target;
  DECLARE_EVENT_TABLE();
};


}
#endif
