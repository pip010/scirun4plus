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

#include "MainWindow.h"

#define QUIT_MENU_ID 51
#define INFO_MENU_ID 52
#define COLOR_MENU_ID 53

#define SPIN_TOG_ID 101
#define BLEND_BITS_TOG_ID 102
#define SHADING_TOG_ID 103
#define BLEND_MODE_ID 104
#define SAMPLE_RATE_ID 105
#define SAMPLE_RATE_TEXT_ID 106
#define TIMER_ID 107
#define CMAP_2D_ID 108

BEGIN_EVENT_TABLE(MainWindow, wxFrame)
  EVT_MENU(QUIT_MENU_ID, MainWindow::on_quit)
  EVT_MENU(INFO_MENU_ID, MainWindow::on_info)
  EVT_MENU(COLOR_MENU_ID, MainWindow::on_color_menu)
  EVT_COMMAND_SCROLL(SAMPLE_RATE_ID, MainWindow::on_slider_change)
  EVT_TOGGLEBUTTON(SPIN_TOG_ID, MainWindow::on_spin_toggle)
  EVT_TOGGLEBUTTON(BLEND_BITS_TOG_ID, MainWindow::on_blend_bits_toggle)
  EVT_TOGGLEBUTTON(SHADING_TOG_ID, MainWindow::on_shading_toggle)
  EVT_TOGGLEBUTTON(BLEND_MODE_ID, MainWindow::on_blend_mode_toggle)
  EVT_TOGGLEBUTTON(CMAP_2D_ID, MainWindow::on_2d_cmap_toggle)
  EVT_TIMER(TIMER_ID, MainWindow::on_timer)
END_EVENT_TABLE()

MainWindow::MainWindow(wxWindow* parent, int id, const wxString& title, const wxPoint& pos, const wxSize& size, long style):
    wxFrame(parent, id, title, pos, size, wxDEFAULT_FRAME_STYLE)
{
  int attribs[] = {WX_GL_RGBA, WX_GL_DOUBLEBUFFER, WX_GL_DEPTH_SIZE, 24, 0};

    // begin wxGlade: MainWindow::MainWindow
    samp_sizer_staticbox = new wxStaticBox(this, -1, wxT("Sample Rate"));
    main_window_menubar = new wxMenuBar();
    SetMenuBar(main_window_menubar);
    wxMenu* wxglade_tmp_menu_1 = new wxMenu();
    wxglade_tmp_menu_1->Append(QUIT_MENU_ID, wxT("Quit"), wxEmptyString, wxITEM_NORMAL);
    main_window_menubar->Append(wxglade_tmp_menu_1, wxT("File"));
    wxMenu* wxglade_tmp_menu_2 = new wxMenu();
    wxglade_tmp_menu_2->Append(COLOR_MENU_ID, wxT("Color"), wxEmptyString, wxITEM_NORMAL);
    main_window_menubar->Append(wxglade_tmp_menu_2, wxT("Selection"));
    wxMenu* wxglade_tmp_menu_3 = new wxMenu();
    wxglade_tmp_menu_3->Append(INFO_MENU_ID, wxT("Info"), wxEmptyString, wxITEM_NORMAL);
    main_window_menubar->Append(wxglade_tmp_menu_3, wxT("About"));
    main_window_statusbar = CreateStatusBar(1, 0);
    spin_tog = new wxToggleButton(this, SPIN_TOG_ID, wxT("Spin"));
    blend32_tog = new wxToggleButton(this, BLEND_BITS_TOG_ID, wxT("32 bit Blend"));
    shading_tog = new wxToggleButton(this, SHADING_TOG_ID, wxT("Shading"));
    bmode_tog = new wxToggleButton(this, BLEND_MODE_ID, wxT("MIP Mode"));
    use_2d_cmap_tog = new wxToggleButton(this, CMAP_2D_ID, wxT("2D ColorMap"));
    sampling_rate = new wxSlider(this, SAMPLE_RATE_ID, 225, 25, 1000, wxDefaultPosition, wxDefaultSize, wxSL_VERTICAL);
    sample_rate_text = new wxTextCtrl(this, SAMPLE_RATE_TEXT_ID, wxT("4.5"));
    glcanvas_ = new OGLCanvas(this, attribs);

    set_properties();
    do_layout();
    // end wxGlade

    timer_ = new wxTimer();
    timer_->SetOwner(this, TIMER_ID);
    timer_->Start(300);

    cpick_ = new wxColourDialog(this);
}


void MainWindow::set_properties()
{
    // begin wxGlade: MainWindow::set_properties
    SetTitle(wxT("SLIVR Example"));
    SetSize(wxSize(573, 780));
    int main_window_statusbar_widths[] = { -1 };
    main_window_statusbar->SetStatusWidths(1, main_window_statusbar_widths);
    const wxString main_window_statusbar_fields[] = {
        wxT("main_window_statusbar")
    };
    for(int i = 0; i < main_window_statusbar->GetFieldsCount(); ++i) {
        main_window_statusbar->SetStatusText(main_window_statusbar_fields[i], i);
    }
    // end wxGlade
}


void MainWindow::do_layout()
{
    // begin wxGlade: MainWindow::do_layout
    wxBoxSizer* sizer_1 = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* sizer_2 = new wxBoxSizer(wxVERTICAL);
    wxStaticBoxSizer* samp_sizer = new wxStaticBoxSizer(samp_sizer_staticbox, wxVERTICAL);
    sizer_2->Add(spin_tog, 1, wxEXPAND|wxADJUST_MINSIZE, 0);
    sizer_2->Add(blend32_tog, 1, wxEXPAND|wxADJUST_MINSIZE, 0);
    sizer_2->Add(shading_tog, 1, wxEXPAND|wxADJUST_MINSIZE, 0);
    sizer_2->Add(bmode_tog, 1, wxEXPAND|wxADJUST_MINSIZE, 0);
    sizer_2->Add(use_2d_cmap_tog, 1, wxEXPAND|wxADJUST_MINSIZE, 0);
    samp_sizer->Add(sampling_rate, 1, wxEXPAND, 0);
    samp_sizer->Add(sample_rate_text, 0, wxEXPAND|wxADJUST_MINSIZE, 0);
    sizer_2->Add(samp_sizer, 25, wxEXPAND, 0);
    sizer_1->Add(sizer_2, 2, wxEXPAND, 0);
    sizer_1->Add(glcanvas_, 10, wxEXPAND, 0);
    SetSizer(sizer_1);
    Layout();
    // end wxGlade
}

void 
MainWindow::on_quit(wxCommandEvent& event)
{
  Close(TRUE);
  exit(0);
}

void 
MainWindow::on_info(wxCommandEvent& event)
{
  cerr << "Test app for volume rendering version 0.0.0.0.0.1" << endl;
}

void 
MainWindow::on_color_menu(wxCommandEvent& event)
{
  cerr << "color_menu" << endl;
  int rval = cpick_->ShowModal();
  wxColourData &cd = cpick_->GetColourData();
  wxColour &c = cd.GetColour();
  cerr << "chose color: ( " << (int)c.Red() << ", " << (int)c.Green()
       << ", " << (int)c.Blue() << ")" << endl;
}



void 
MainWindow::on_slider_change(wxScrollEvent & event)
{
  float val = (float)event.GetPosition() / 1000.0 * 20.0;
  glcanvas_->set_sampling_rate(val);
  wxString str;
  str << val;
  sample_rate_text->SetValue(str);
  // update the view
//   wxPaintEvent paint_me;
//   wxPostEvent(glcanvas_, paint_me)
    ;
}

void 
MainWindow::on_spin_toggle(wxCommandEvent & event)
{
  if (spin_tog->GetValue()) {
    glcanvas_->set_spin_p(true);
  } else {
    glcanvas_->set_spin_p(false);
  }
}

void 
MainWindow::on_blend_bits_toggle(wxCommandEvent & event)
{
  if (blend32_tog->GetValue()) {
    glcanvas_->set_blend_num_bits_32();
  } else {
    glcanvas_->set_blend_num_bits_8();
  }
}
void 
MainWindow::on_shading_toggle(wxCommandEvent & event)
{
  if (shading_tog->GetValue()) {
    glcanvas_->set_shading(true);
  } else {
    glcanvas_->set_shading(false);
  }
}

void 
MainWindow::on_blend_mode_toggle(wxCommandEvent & event)
{
  if (bmode_tog->GetValue()) {
    glcanvas_->set_blend_mode_over(false);
  } else {
    glcanvas_->set_blend_mode_over(true);
  }
}

void 
MainWindow::on_2d_cmap_toggle(wxCommandEvent & event)
{
  if (use_2d_cmap_tog->GetValue()) {
    glcanvas_->set_use_2d_cmap(true);
  } else {
    glcanvas_->set_use_2d_cmap(false);
  }
}

void 
MainWindow::on_timer(wxTimerEvent & event)
{
  //cerr << "timer event...." << endl;
  wxPaintEvent paint_me;
  wxPostEvent(glcanvas_, paint_me);
}
