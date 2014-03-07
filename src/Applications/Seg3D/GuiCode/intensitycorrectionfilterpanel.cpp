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
// Name:        intensitycorrectionfilterpanel.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     Tue 08 Jul 2008 16:45:24 MDT
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

#include "intensitycorrectionfilterpanel.h"
#include "seg3devents.h"
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

////@begin XPM images
////@end XPM images


/*!
 * IntensityCorrectionFilterPanel type definition
 */

IMPLEMENT_DYNAMIC_CLASS( IntensityCorrectionFilterPanel, wxPanel )


/*!
 * IntensityCorrectionFilterPanel event table definition
 */

BEGIN_EVENT_TABLE( IntensityCorrectionFilterPanel, wxPanel )

////@begin IntensityCorrectionFilterPanel event table entries
////@end IntensityCorrectionFilterPanel event table entries

    EVT_BUTTON( XRCID("START_BUTTON"), IntensityCorrectionFilterPanel::OnStartButtonClick )
    EVT_BUTTON( XRCID("CLOSE_BUTTON"), IntensityCorrectionFilterPanel::OnCloseButtonClick )

    EVT_COMMAND( wxID_ANY, wxEVT_SET_PROGRESS, IntensityCorrectionFilterPanel::OnSetProgress)

END_EVENT_TABLE()


/*!
 * IntensityCorrectionFilterPanel constructors
 */

IntensityCorrectionFilterPanel::IntensityCorrectionFilterPanel()
{
    Init();
}

IntensityCorrectionFilterPanel::IntensityCorrectionFilterPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * HistoEqFilter creator
 */

bool IntensityCorrectionFilterPanel::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin IntensityCorrectionFilterPanel creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
////@end IntensityCorrectionFilterPanel creation
    return true;
}


/*!
 * IntensityCorrectionFilterPanel destructor
 */

IntensityCorrectionFilterPanel::~IntensityCorrectionFilterPanel()
{
////@begin IntensityCorrectionFilterPanel destruction
////@end IntensityCorrectionFilterPanel destruction
}


/*!
 * Member initialisation
 */

void IntensityCorrectionFilterPanel::Init()
{
////@begin IntensityCorrectionFilterPanel member initialisation
    mOrder = NULL;
    mEdge = NULL;
    mPercentage = NULL;
////@end IntensityCorrectionFilterPanel member initialisation

    disabler_ = NULL;
}


/*!
 * Control creation for HistoEqFilter
 */

void IntensityCorrectionFilterPanel::CreateControls()
{    
////@begin IntensityCorrectionFilterPanel content construction
    if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_INTENSITYCORRECTIONFILTER")))
        wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
    mOrder = XRCCTRL(*this, "ID_ORDER", wxSpinCtrl);
    mEdge = XRCCTRL(*this, "ID_EDGE", wxTextCtrl);
    mPercentage = XRCCTRL(*this, "ID_GAUGE1", wxGauge);
    // Set validators
    if (FindWindow(XRCID("ID_EDGE")))
        FindWindow(XRCID("ID_EDGE"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
////@end IntensityCorrectionFilterPanel content construction

    // Create custom windows not generated automatically here.
////@begin IntensityCorrectionFilterPanel content initialisation
////@end IntensityCorrectionFilterPanel content initialisation
}

void IntensityCorrectionFilterPanel::OnStartButtonClick( wxCommandEvent& event )
{
  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent("Painter::FinishTool");
  tsse->add_var("Painter::IntensityCorrection::order",
                SCIRun::to_string(mOrder->GetValue()));
  tsse->add_var("Painter::IntensityCorrection::::edge",
                wx2std(mEdge->GetValue()));
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void IntensityCorrectionFilterPanel::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}

void IntensityCorrectionFilterPanel::OnSetProgress( wxCommandEvent& event )
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


/*!
 * Should we show tooltips?
 */

bool IntensityCorrectionFilterPanel::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap IntensityCorrectionFilterPanel::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin IntensityCorrectionFilterPanel bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end IntensityCorrectionFilterPanel bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon IntensityCorrectionFilterPanel::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin IntensityCorrectionFilterPanel icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end IntensityCorrectionFilterPanel icon retrieval
}
