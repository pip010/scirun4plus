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
// Name:        optionlessfilter.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     Wed 19 Dec 2007 14:11:29 MST
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"
#include <Applications/Seg3D/Painter.h>

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

////@begin includes
////@end includes

#include "optionlessfilter.h"
#include "seg3devents.h"

////@begin XPM images
////@end XPM images


/*!
 * OptionlessFilter type definition
 */

IMPLEMENT_DYNAMIC_CLASS( OptionlessFilter, wxPanel )


/*!
 * OptionlessFilter event table definition
 */

BEGIN_EVENT_TABLE( OptionlessFilter, wxPanel )

////@begin OptionlessFilter event table entries
////@end OptionlessFilter event table entries

    EVT_BUTTON( XRCID("START_BUTTON"), OptionlessFilter::OnStartButtonClick )
    EVT_BUTTON( XRCID("CLOSE_BUTTON"), OptionlessFilter::OnCloseButtonClick )

    EVT_COMMAND( wxID_ANY, wxEVT_SET_PROGRESS, OptionlessFilter::OnSetProgress)

END_EVENT_TABLE()


/*!
 * OptionlessFilter constructors
 */

OptionlessFilter::OptionlessFilter()
{
    Init();
}

OptionlessFilter::OptionlessFilter( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * OptionLessFilter creator
 */

bool OptionlessFilter::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin OptionlessFilter creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end OptionlessFilter creation
    return true;
}


/*!
 * OptionlessFilter destructor
 */

OptionlessFilter::~OptionlessFilter()
{
////@begin OptionlessFilter destruction
////@end OptionlessFilter destruction
}


/*!
 * Member initialisation
 */

void OptionlessFilter::Init()
{
////@begin OptionlessFilter member initialisation
    mPercentage = NULL;
////@end OptionlessFilter member initialisation

    disabler_ = NULL;
    skinner_callback_ = "Painter::FinishTool";
    show_progress_ = false;
}


/*!
 * Control creation for OptionLessFilter
 */

void OptionlessFilter::CreateControls()
{    
////@begin OptionlessFilter content construction
    if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_OPTIONLESSFILTER")))
        wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
    mPercentage = XRCCTRL(*this, "ID_GAUGE", wxGauge);
////@end OptionlessFilter content construction

    // Create custom windows not generated automatically here.
////@begin OptionlessFilter content initialisation
////@end OptionlessFilter content initialisation
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void OptionlessFilter::OnStartButtonClick( wxCommandEvent& event )
{
  if (show_progress_)
  {
    SCIRun::Painter::ThrowSkinnerSignal(skinner_callback_, false);
  }
  else
  {
    // Call synchronously and wait.
    wxBeginBusyCursor();
    SCIRun::Painter::ThrowSkinnerSignal(skinner_callback_, true);
    wxEndBusyCursor();
  }
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void OptionlessFilter::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


void OptionlessFilter::OnSetProgress( wxCommandEvent& event )
{
  int progress = event.GetInt();
  
  // start_progress() sends -1
  if (progress < 0)
  {
    disabler_ = new wxWindowDisabler();
#ifndef _WIN32
    wxBeginBusyCursor();
#endif
    progress = 0;
  }
  // finish_progress() sends 101
  if (progress > 100)
  {
    if (disabler_) { delete disabler_; disabler_ = 0; }
#ifndef _WIN32
    wxEndBusyCursor();
#endif
    progress = 100;
  }
  if (progress == 0 || progress > mPercentage->GetValue())
  {
    mPercentage->SetValue(progress);
  }
}


void OptionlessFilter::SetSkinnerCallback(const char *callback)
{
  skinner_callback_ = callback;
}

void OptionlessFilter::SetShowProgress(bool show_progress)
{
  show_progress_ = show_progress;
  if (show_progress_)
  {
    mPercentage->Show();
  }
  else
  {
    mPercentage->Hide();
  }
}

/*!
 * Should we show tooltips?
 */

bool OptionlessFilter::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap OptionlessFilter::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin OptionlessFilter bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end OptionlessFilter bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon OptionlessFilter::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin OptionlessFilter icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end OptionlessFilter icon retrieval
}
