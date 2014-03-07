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


#include <wx/wx.h>
#include <wx/image.h>
#include <wx/colordlg.h>
#include <OGLCanvas.h>

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

// begin wxGlade: ::dependencies
#include <wx/tglbtn.h>
// end wxGlade

using namespace SLIVR_EXAMPLES;
class MainWindow: public wxFrame {
public:
    // begin wxGlade: MainWindow::ids
    // end wxGlade

  MainWindow(wxWindow* parent, int id, const wxString& title, const wxPoint& pos=wxDefaultPosition, const wxSize& size=wxDefaultSize, long style=wxDEFAULT_FRAME_STYLE);

private:
    // begin wxGlade: MainWindow::methods
    void set_properties();
    void do_layout();
    // end wxGlade
  void on_quit(wxCommandEvent& event);
  void on_info(wxCommandEvent& event);
  void on_color_menu(wxCommandEvent& event);
  void on_slider_change(wxScrollEvent& event);
  void on_spin_toggle(wxCommandEvent & event);
  void on_blend_bits_toggle(wxCommandEvent & event);
  void on_shading_toggle(wxCommandEvent & event);
  void on_blend_mode_toggle(wxCommandEvent & event);
  void on_2d_cmap_toggle(wxCommandEvent & event);
  void on_timer(wxTimerEvent &event);

protected:
    // begin wxGlade: MainWindow::attributes
    wxStaticBox* samp_sizer_staticbox;
    wxMenuBar* main_window_menubar;
    wxStatusBar* main_window_statusbar;
    wxToggleButton* spin_tog;
    wxToggleButton* blend32_tog;
    wxToggleButton* shading_tog;
    wxToggleButton* bmode_tog;
    wxToggleButton* use_2d_cmap_tog;
    wxSlider* sampling_rate;
    wxTextCtrl* sample_rate_text;
    OGLCanvas* glcanvas_;
    // end wxGlade
  wxTimer *timer_;
  wxColourDialog *cpick_;
  DECLARE_EVENT_TABLE();
}; // wxGlade: end class


#endif // MAINWINDOW_H
