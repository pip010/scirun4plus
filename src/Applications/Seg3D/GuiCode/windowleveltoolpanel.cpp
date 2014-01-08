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

/////////////////////////////////////////////////////////////////////////////
// Name:        windowleveltoolpanel.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     Wed 09 Apr 2008 17:18:05 MDT
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"
#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

////@begin includes
////@end includes

#include "windowleveltoolpanel.h"
#include <stdlib.h>
#include <sstream>

////@begin XPM images
////@end XPM images


/*!
 * WindowLevelToolPanel type definition
 */

IMPLEMENT_DYNAMIC_CLASS( WindowLevelToolPanel, wxPanel )


/*!
 * WindowLevelToolPanel event table definition
 */

BEGIN_EVENT_TABLE( WindowLevelToolPanel, wxPanel )

////@begin WindowLevelToolPanel event table entries
////@end WindowLevelToolPanel event table entries

  EVT_BUTTON( XRCID("CLOSE_BUTTON"), WindowLevelToolPanel::OnCloseButtonClick )
  EVT_COMMAND_SCROLL( XRCID("ID_WINDOW_SLIDER"), WindowLevelToolPanel::OnWindowSlider )
  EVT_COMMAND_SCROLL( XRCID("ID_LEVEL_SLIDER"), WindowLevelToolPanel::OnLevelSlider )
  EVT_TEXT( XRCID("ID_WINDOW_TEXTCTRL"), WindowLevelToolPanel::OnWindowText )
  EVT_TEXT( XRCID("ID_LEVEL_TEXTCTRL"), WindowLevelToolPanel::OnLevelText )
  EVT_COMMAND (wxID_ANY, wxEVT_WINDOWLEVELTOOL_CHANGE, WindowLevelToolPanel::OnWindowLevelToolChange )
#ifndef _WIN32
  EVT_PAINT(WindowLevelToolPanel::OnPaint)
#endif
END_EVENT_TABLE()


/*!
 * WindowLevelToolPanel constructors
 */

WindowLevelToolPanel::WindowLevelToolPanel()
{
  Init();
}

WindowLevelToolPanel::WindowLevelToolPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
  Init();
  Create(parent, id, pos, size, style);
}


/*!
 * WindowLevelToolPanel creator
 */

bool WindowLevelToolPanel::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin WindowLevelToolPanel creation
  SetParent(parent);
  CreateControls();
  if (GetSizer())
  {
    GetSizer()->SetSizeHints(this);
  }
////@end WindowLevelToolPanel creation
  return true;
}


/*!
 * WindowLevelToolPanel destructor
 */

WindowLevelToolPanel::~WindowLevelToolPanel()
{
  if( mHistoWindow )
    delete mHistoWindow;
  
  if( histogram_image_data_ )
    free( histogram_image_data_);

////@begin WindowLevelToolPanel destruction
////@end WindowLevelToolPanel destruction
}


/*!
 * Member initialisation
 */

void WindowLevelToolPanel::Init()
{
////@begin WindowLevelToolPanel member initialisation
    mHistoWindow = NULL;
    mMinValue = NULL;
    mMaxValue = NULL;
    mWindowValue = NULL;
    mWindowSlider = NULL;
    mLevelValue = NULL;
    mLevelSlider = NULL;
////@end WindowLevelToolPanel member initialisation

    minval_ = 0.0;
    maxval_ = 0.0;
    window_ = AIR_NAN;
    level_ = AIR_NAN;

    nbins_ = 200;
    histogram_image_data_ = NULL;
}


/*!
 * Control creation for WindowLevelToolPanel
 */

void WindowLevelToolPanel::CreateControls()
{    
  ////@begin WindowLevelToolPanel content construction
  if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_WINDOWLEVELTOOLPANEL")))
    wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));

  mMinValue = XRCCTRL(*this, "ID_MIN_TEXTCTRL", wxTextCtrl);
  mMaxValue = XRCCTRL(*this, "ID_MAX_TEXTCTRL", wxTextCtrl);
  mWindowValue = XRCCTRL(*this, "ID_WINDOW_TEXTCTRL", wxTextCtrl);
  mWindowSlider = XRCCTRL(*this, "ID_WINDOW_SLIDER", wxSlider);
  mLevelValue = XRCCTRL(*this, "ID_LEVEL_TEXTCTRL", wxTextCtrl);
  mLevelSlider = XRCCTRL(*this, "ID_LEVEL_SLIDER", wxSlider);
  // Set validators
  if (FindWindow(XRCID("ID_MIN_TEXTCTRL")))
    FindWindow(XRCID("ID_MIN_TEXTCTRL"))->
      SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
  if (FindWindow(XRCID("ID_MAX_TEXTCTRL")))
    FindWindow(XRCID("ID_MAX_TEXTCTRL"))->
      SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
  if (FindWindow(XRCID("ID_WINDOW_TEXTCTRL")))
    FindWindow(XRCID("ID_WINDOW_TEXTCTRL"))->
      SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
  if (FindWindow(XRCID("ID_LEVEL_TEXTCTRL")))
    FindWindow(XRCID("ID_LEVEL_TEXTCTRL"))->
      SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
  ////@end WindowLevelToolPanel content construction

#ifndef _WIN32
  // Create custom windows not generated automatically here.
  ////@begin WindowLevelToolPanel content initialisation
  mHistoWindow = new wxWindow(this, (wxWindowID) -1,
			      wxPoint(0, 0), wxSize(nbins_, 125),
			      wxSUNKEN_BORDER, wxString(wxT("Histogram")));
  ////@end WindowLevelToolPanel content initialisation

  // Storage for the histogram image.
  histogram_image_data_ = (unsigned char *)
    malloc( sizeof(unsigned char) * nbins_ * 125 * 3 );
#endif
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void WindowLevelToolPanel::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


void WindowLevelToolPanel::CreateHistogramImage()
{
  // Make sure there is valid histrogram data.
  if (histogram_.size() != nbins_)
    return;

  // Create the image data;
  for( unsigned int i=0; i<nbins_ * 125 * 3; ++i )
    histogram_image_data_[i] = 0;

  float max_values = 1;

  // Ignore the first and last bin as it usually contains lots of
  // nothing. This gives a better histogram.
  for( unsigned int i=1; i<histogram_.size()-1; ++i)
    max_values = max_values > histogram_[i] ? max_values : histogram_[i];

  // Load the histogram data
  for( unsigned int i=0; i<histogram_.size(); ++i)
  {
    unsigned int height =
      (unsigned int) (125.0 * (1.0 - (float) histogram_[i] / max_values));

    if( height > 125 ) 
      height = 125;

    for( unsigned int j=0; j<height; ++j)
    {
      histogram_image_data_[(j*nbins_+i)*3  ] = 150;
      histogram_image_data_[(j*nbins_+i)*3+1] = 150;
      histogram_image_data_[(j*nbins_+i)*3+2] = 150;
    }
  }

  // Now do the window and level bounds.(Divide by 1000 ticks into 200
  // bins == 5)
  unsigned int window_lower = (unsigned int)
    ((nbins_- 1.0) *
     (float) (mLevelSlider->GetValue() - mWindowSlider->GetValue()/2 ) / 1000.0);

  if( window_lower < 0 )
    window_lower = 0;

  if( window_lower > (nbins_-1) )
    window_lower = (nbins_-1);

  unsigned int window_upper = (unsigned int)
    ((nbins_-1.0) *
     (float) (mLevelSlider->GetValue() + mWindowSlider->GetValue()/2 ) / 1000.0);

  if( window_upper < 0 )
    window_upper = 0;

  if( window_upper > (nbins_-1) )
    window_upper = (nbins_-1);

  unsigned int level = (unsigned int)
    ((nbins_-1.0) * (float) mLevelSlider->GetValue() / 1000.0);

  if( level < 0 )
    level = 0;

  if( level > (nbins_-1) )
    level = (nbins_-1);

  for( unsigned int j=0; j<125; ++j)
  {
    histogram_image_data_[(j*nbins_+window_lower)*3  ] = 0;
    histogram_image_data_[(j*nbins_+window_lower)*3+1] = 150;
    histogram_image_data_[(j*nbins_+window_lower)*3+2] = 0;
  }

  for( unsigned int j=0; j<125; ++j)
  {
    histogram_image_data_[(j*nbins_+window_upper)*3  ] = 0;
    histogram_image_data_[(j*nbins_+window_upper)*3+1] = 150;
    histogram_image_data_[(j*nbins_+window_upper)*3+2] = 0;
  }

  for( unsigned int j=0; j<125; ++j)
  {
    histogram_image_data_[(j*nbins_+level)*3  ] = 150;
    histogram_image_data_[(j*nbins_+level)*3+1] = 0;
    histogram_image_data_[(j*nbins_+level)*3+2] = 0;
  }

  // Anytime the histogram gets changed throw a paint event so that it
  // shows up on the inital view.
  wxPaintEvent pe(0);
  OnPaint( pe );
}


void WindowLevelToolPanel::OnWindowLevelToolChange( wxCommandEvent& event )
{
  WindowLevelToolChangeStruct *info =
    (WindowLevelToolChangeStruct *)event.GetClientData();

  minval_ = info->minval;
  maxval_ = info->maxval;
  mMinValue->Clear(); *mMinValue << info->minval;
  mMaxValue->Clear(); *mMaxValue << info->maxval;
  mWindowValue->Clear(); *mWindowValue << info->window;
  mLevelValue->Clear(); *mLevelValue << info->level;

  // Update the slider
  double window = (info->window) / (maxval_ - minval_);
  if (window < 0.0) window = 0.0;
  if (window > 1.0) window = 1.0;
  mWindowSlider->SetValue((int)(1000.0 * window));

  double level = (info->level - minval_) / (maxval_ - minval_);
  if (level < 0.0) level = 0.0;
  if (level > 1.0) level = 1.0;
  mLevelSlider->SetValue((int)(1000.0 * level));

#ifndef _WIN32
  if( info->histogram.size() == nbins_ ) {
    histogram_ = info->histogram;
    CreateHistogramImage();
  }
#endif
  delete info;
}

void WindowLevelToolPanel::OnWindowSlider( wxScrollEvent& event )
{
  double window =
    mWindowSlider->GetValue() / 1000.0 * (maxval_ - minval_);

  if (window != window_)
  {
    window_ = window;
    mWindowValue->Clear();
    *mWindowValue << window_;

    PushToSkinner();
  }
}

void WindowLevelToolPanel::OnLevelSlider( wxScrollEvent& event )
{
  double level =
    mLevelSlider->GetValue() / 1000.0 * (maxval_ - minval_) + minval_;
  if (level != level_)
  {
    level_ = level;
    mLevelValue->Clear();
    *mLevelValue << level_;

    PushToSkinner();
  }
}

void WindowLevelToolPanel::OnWindowText( wxCommandEvent& event )
{
  if (mWindowValue && mWindowValue->IsModified())
  {
    const string str = wx2std(mWindowValue->GetValue());
    std::istringstream sstr(str);
    sstr >> window_;

    // Update the slider
    double window = (window_) / (maxval_ - minval_);
    if (window < 0.0) window = 0.0;
    if (window > 1.0) window = 1.0;
    mWindowSlider->SetValue((int)(1000.0 * window));

    PushToSkinner();
  }
}

void WindowLevelToolPanel::OnLevelText( wxCommandEvent& event )
{
  if (mLevelValue && mLevelValue->IsModified())
  {
    const string str = wx2std(mLevelValue->GetValue());
    std::istringstream sstr(str);
    sstr >> level_;

    // Update the slider
    double level = (level_ - minval_) / (maxval_ - minval_);
    if (level < 0.0) level = 0.0;
    if (level > 1.0) level = 1.0;
    mLevelSlider->SetValue((int)(1000.0 * level));

    PushToSkinner();
  }
}

void WindowLevelToolPanel::OnPaint( wxPaintEvent &WXUNUSED(event) )
{
  wxPaintDC dc( mHistoWindow );
  PrepareDC( dc );

  wxImage image;

  image.Create(nbins_, 125, true);

  unsigned char * imageData = image.GetData();

  if( histogram_image_data_ )
    memcpy(imageData, histogram_image_data_, nbins_ * 125 * 3 *sizeof(unsigned char));

  wxBitmap bitmap = wxBitmap( image );

  if (bitmap.Ok())
  {      
    dc.DrawBitmap( bitmap, 0, 0 );
  }
}

/*!
 * should we show tooltips?
 */

bool WindowLevelToolPanel::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap WindowLevelToolPanel::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin WindowLevelToolPanel bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end WindowLevelToolPanel bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon WindowLevelToolPanel::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin WindowLevelToolPanel icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end WindowLevelToolPanel icon retrieval
}


void WindowLevelToolPanel::PushToSkinner()
{
  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent("Painter::UpdateWindowLevelTool");
  tsse->add_var("Painter::windowleveltool_window",
                wx2std(mWindowValue->GetValue()));
  tsse->add_var("Painter::windowleveltool_level",
                wx2std(mLevelValue->GetValue()));
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);

#ifndef _WIN32
  CreateHistogramImage();
#endif
}
