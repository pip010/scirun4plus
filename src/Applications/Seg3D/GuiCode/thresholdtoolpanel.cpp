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
// Name:        thresholdtoolpanel.cpp
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

#include "thresholdtoolpanel.h"
#include <stdlib.h>
#include <sstream>

////@begin XPM images
////@end XPM images


/*!
 * ThresholdToolPanel type definition
 */

IMPLEMENT_DYNAMIC_CLASS( ThresholdToolPanel, wxPanel )


/*!
 * ThresholdToolPanel event table definition
 */

BEGIN_EVENT_TABLE( ThresholdToolPanel, wxPanel )

////@begin ThresholdToolPanel event table entries
////@end ThresholdToolPanel event table entries

  EVT_BUTTON( XRCID("CLEAR_SEEDS_BUTTON"), ThresholdToolPanel::OnClearSeedsButtonClick )
  EVT_BUTTON( XRCID("START_BUTTON"), ThresholdToolPanel::OnStartButtonClick )
  EVT_BUTTON( XRCID("CLOSE_BUTTON"), ThresholdToolPanel::OnCloseButtonClick )
  EVT_COMMAND_SCROLL( XRCID("ID_LOWER_SLIDER"), ThresholdToolPanel::OnLowerSlider )
  EVT_COMMAND_SCROLL( XRCID("ID_UPPER_SLIDER"), ThresholdToolPanel::OnUpperSlider )
  EVT_TEXT( XRCID("ID_LOWER_TEXTCTRL"), ThresholdToolPanel::OnLowerText )
  EVT_TEXT( XRCID("ID_UPPER_TEXTCTRL"), ThresholdToolPanel::OnUpperText )
  EVT_COMMAND (wxID_ANY, wxEVT_THRESHOLDTOOL_CHANGE, ThresholdToolPanel::OnThresholdToolChange )
#ifndef _WIN32
  EVT_PAINT(ThresholdToolPanel::OnPaint)
#endif
END_EVENT_TABLE()


/*!
 * ThresholdToolPanel constructors
 */

ThresholdToolPanel::ThresholdToolPanel()
{
    Init();
}

ThresholdToolPanel::ThresholdToolPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * ThresholdToolPanel creator
 */

bool ThresholdToolPanel::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin ThresholdToolPanel creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
////@end ThresholdToolPanel creation
    return true;
}


/*!
 * ThresholdToolPanel destructor
 */

ThresholdToolPanel::~ThresholdToolPanel()
{
  if( mHistoWindow )
    delete mHistoWindow;

  if( histogram_image_data_ )
    free( histogram_image_data_);

////@begin ThresholdToolPanel destruction
////@end ThresholdToolPanel destruction
}


/*!
 * Member initialisation
 */

void ThresholdToolPanel::Init()
{
////@begin ThresholdToolPanel member initialisation
    mHistoWindow = NULL;
    mMinValue = NULL;
    mMaxValue = NULL;
    mLowerValue = NULL;
    mLowerSlider = NULL;
    mUpperValue = NULL;
    mUpperSlider = NULL;
////@end ThresholdToolPanel member initialisation

    minval_ = 0.0;
    maxval_ = 0.0;
    lower_ = AIR_NAN;
    upper_ = AIR_NAN;

    nbins_ = 200;
    histogram_image_data_ = NULL;
}


/*!
 * Control creation for ThresholdToolPanel
 */

void ThresholdToolPanel::CreateControls()
{    
  ////@begin ThresholdToolPanel content construction
  if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_THRESHOLDTOOLPANEL")))
    wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));

  mMinValue = XRCCTRL(*this, "ID_MIN_TEXTCTRL", wxTextCtrl);
  mMaxValue = XRCCTRL(*this, "ID_MAX_TEXTCTRL", wxTextCtrl);
  mLowerValue = XRCCTRL(*this, "ID_LOWER_TEXTCTRL", wxTextCtrl);
  mLowerSlider = XRCCTRL(*this, "ID_LOWER_SLIDER", wxSlider);
  mUpperValue = XRCCTRL(*this, "ID_UPPER_TEXTCTRL", wxTextCtrl);
  mUpperSlider = XRCCTRL(*this, "ID_UPPER_SLIDER", wxSlider);
  // Set validators
  if (FindWindow(XRCID("ID_MIN_TEXTCTRL")))
    FindWindow(XRCID("ID_MIN_TEXTCTRL"))->
      SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
  if (FindWindow(XRCID("ID_MAX_TEXTCTRL")))
    FindWindow(XRCID("ID_MAX_TEXTCTRL"))->
      SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
  if (FindWindow(XRCID("ID_LOWER_TEXTCTRL")))
    FindWindow(XRCID("ID_LOWER_TEXTCTRL"))->
      SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
  if (FindWindow(XRCID("ID_UPPER_TEXTCTRL")))
    FindWindow(XRCID("ID_UPPER_TEXTCTRL"))->
      SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
  ////@end ThresholdToolPanel content construction

#ifndef _WIN32
  // Create custom windows not generated automatically here.
  ////@begin ThresholdToolPanel content initialisation
  mHistoWindow = new wxWindow(this, (wxWindowID) -1,
			      wxPoint(0, 0), wxSize(nbins_, 125),
			      wxSUNKEN_BORDER, wxString(wxT("Histogram")));
  ////@end ThresholdToolPanel content initialisation

  // Storage for the histogram image.
  histogram_image_data_ = (unsigned char *)
    malloc( sizeof(unsigned char) * nbins_ * 125 * 3 );
#endif
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLEAR_SEEDS_BUTTON
 */
void ThresholdToolPanel::OnClearSeedsButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::ThrowSkinnerSignal("Painter::SetLayer");
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void ThresholdToolPanel::OnStartButtonClick( wxCommandEvent& event )
{
  wxBeginBusyCursor();
  SCIRun::Painter::ThrowSkinnerSignal("Painter::FinishTool");
  wxEndBusyCursor();
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void ThresholdToolPanel::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


void ThresholdToolPanel::CreateHistogramImage()
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

  // Now do the threshold bounds. (Divide by 1000 ticks into nbins_ bins == 5)
  unsigned int lower = (unsigned int)
    ((nbins_-1.0) * (float) mLowerSlider->GetValue() / 1000.0);
  unsigned int upper = (unsigned int)
    ((nbins_-1.0) * (float) mUpperSlider->GetValue() / 1000.0);

  if( lower < 0 )
    lower = 0;

  if( lower > (nbins_-1) )
    lower = (nbins_-1);

  if( upper < 0 )
    upper = 0;

  if( upper > (nbins_-1) )
    upper = (nbins_-1);

  for( unsigned int j=0; j<125; ++j)
  {
    histogram_image_data_[(j*nbins_+lower)*3  ] = 150;
    histogram_image_data_[(j*nbins_+lower)*3+1] = 0;
    histogram_image_data_[(j*nbins_+lower)*3+2] = 0;
  }

  for( unsigned int j=0; j<125; ++j)
  {
    histogram_image_data_[(j*nbins_+upper)*3  ] = 150;
    histogram_image_data_[(j*nbins_+upper)*3+1] = 0;
    histogram_image_data_[(j*nbins_+upper)*3+2] = 0;
  }

  // Anytime the histogram gets changed throw a paint event so that it
  // shows up on the inital view.
  wxPaintEvent pe(0);
  OnPaint( pe );
}


void ThresholdToolPanel::OnThresholdToolChange( wxCommandEvent& event )
{
  ThresholdToolChangeStruct *info =
    (ThresholdToolChangeStruct *)event.GetClientData();

  minval_ = info->minval;
  maxval_ = info->maxval;
  mMinValue->Clear(); *mMinValue << info->minval;
  mMaxValue->Clear(); *mMaxValue << info->maxval;
  mLowerValue->Clear(); *mLowerValue << info->lower;
  mUpperValue->Clear(); *mUpperValue << info->upper;

  // Update the slider
  double lowern = (info->lower - minval_) / (maxval_ - minval_);
  if (lowern < 0.0) lowern = 0.0;
  if (lowern > 1.0) lowern = 1.0;
  mLowerSlider->SetValue((int)(1000.0 * lowern));

  double uppern = (info->upper - minval_) / (maxval_ - minval_);
  if (uppern < 0.0) uppern = 0.0;
  if (uppern > 1.0) uppern = 1.0;
  mUpperSlider->SetValue((int)(1000.0 * uppern));

#ifndef _WIN32
  if( info->histogram.size() == nbins_ ) {
    histogram_ = info->histogram;
    CreateHistogramImage();
  }
#endif
  delete info;
}

void ThresholdToolPanel::OnLowerSlider( wxScrollEvent& event )
{
  double lower =
    mLowerSlider->GetValue() / 1000.0 * (maxval_ - minval_) + minval_;
  if (lower != lower_)
  {
    lower_ = lower;
    mLowerValue->Clear();
    *mLowerValue << lower_;

    PushToSkinner();
  }
}

void ThresholdToolPanel::OnUpperSlider( wxScrollEvent& event )
{
  double upper =
    mUpperSlider->GetValue() / 1000.0 * (maxval_ - minval_) + minval_;
  if (upper != upper_)
  {
    upper_ = upper;
    mUpperValue->Clear();
    *mUpperValue << upper_;

    PushToSkinner();
  }
}

void ThresholdToolPanel::OnLowerText( wxCommandEvent& event )
{
  if (mLowerValue && mLowerValue->IsModified())
  {
    const string str = wx2std(mLowerValue->GetValue());
    std::istringstream sstr(str);
    sstr >> lower_;

    // Update the slider
    double lowern = (lower_ - minval_) / (maxval_ - minval_);
    if (lowern < 0.0) lowern = 0.0;
    if (lowern > 1.0) lowern = 1.0;
    mLowerSlider->SetValue((int)(1000.0 * lowern));

    PushToSkinner();
  }
}

void ThresholdToolPanel::OnUpperText( wxCommandEvent& event )
{
  if (mUpperValue && mUpperValue->IsModified())
  {
    const string str = wx2std(mUpperValue->GetValue());
    std::istringstream sstr(str);
    sstr >> upper_;

    // Update the slider
    double uppern = (upper_ - minval_) / (maxval_ - minval_);
    if (uppern < 0.0) uppern = 0.0;
    if (uppern > 1.0) uppern = 1.0;
    mUpperSlider->SetValue((int)(1000.0 * uppern));

    PushToSkinner();
  }
}

void ThresholdToolPanel::OnPaint( wxPaintEvent &event )
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

bool ThresholdToolPanel::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap ThresholdToolPanel::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin ThresholdToolPanel bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end ThresholdToolPanel bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon ThresholdToolPanel::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin ThresholdToolPanel icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end ThresholdToolPanel icon retrieval
}


void ThresholdToolPanel::PushToSkinner()
{
  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent("Painter::UpdateThresholdTool");
  tsse->add_var("Painter::thresholdtool_lower",
                wx2std(mLowerValue->GetValue()));
  tsse->add_var("Painter::thresholdtool_upper",
                wx2std(mUpperValue->GetValue()));
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);

#ifndef _WIN32
  CreateHistogramImage();
#endif
}
