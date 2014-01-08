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
// Name:        itkspeedtopathiterateneighborhoodfilter.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     12/03/2008 15:18:27
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"
#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/Seg3DFrame.h>
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

////@begin includes
////@end includes

#include "itkspeedtopathiterateneighborhoodfilter.h"
#include "seg3devents.h"

////@begin XPM images
////@end XPM images


/*!
 * ITKSpeedToPathIterateNeighborhoodFilter type definition
 */

IMPLEMENT_DYNAMIC_CLASS( ITKSpeedToPathIterateNeighborhoodFilter, wxPanel )


/*!
 * ITKSpeedToPathIterateNeighborhoodFilter event table definition
 */

BEGIN_EVENT_TABLE( ITKSpeedToPathIterateNeighborhoodFilter, wxPanel )


EVT_BUTTON( XRCID("START_BUTTON"), ITKSpeedToPathIterateNeighborhoodFilter::OnStartButtonClick )
EVT_BUTTON( XRCID("CLOSE_BUTTON"), ITKSpeedToPathIterateNeighborhoodFilter::OnCloseButtonClick )
EVT_BUTTON( XRCID("CLEAR_SEEDS_BUTTON"), ITKSpeedToPathIterateNeighborhoodFilter::OnClearSeedsButtonClick )
EVT_CHECKBOX( XRCID("SEED_CHECKBOX"), ITKSpeedToPathIterateNeighborhoodFilter::OnAddSeedsCheck )
EVT_BUTTON( XRCID("SPEED_BUTTON"), ITKSpeedToPathIterateNeighborhoodFilter::OnCreateSpeedButtonClick )

EVT_COMMAND( wxID_ANY, wxEVT_SET_PROGRESS, ITKSpeedToPathIterateNeighborhoodFilter::OnSetProgress)


////@begin ITKSpeedToPathIterateNeighborhoodFilter event table entries
////@end ITKSpeedToPathIterateNeighborhoodFilter event table entries

END_EVENT_TABLE()


/*!
 * ITKSpeedToPathIterateNeighborhoodFilter constructors
 */

ITKSpeedToPathIterateNeighborhoodFilter::ITKSpeedToPathIterateNeighborhoodFilter()
{
    Init();
}

ITKSpeedToPathIterateNeighborhoodFilter::ITKSpeedToPathIterateNeighborhoodFilter( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * ITKSpeedToPathIterateNeighborhoodFilter creator
 */

bool ITKSpeedToPathIterateNeighborhoodFilter::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin ITKSpeedToPathIterateNeighborhoodFilter creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
////@end ITKSpeedToPathIterateNeighborhoodFilter creation
    return true;
}


/*!
 * ITKSpeedToPathIterateNeighborhoodFilter destructor
 */

ITKSpeedToPathIterateNeighborhoodFilter::~ITKSpeedToPathIterateNeighborhoodFilter()
{
////@begin ITKSpeedToPathIterateNeighborhoodFilter destruction
////@end ITKSpeedToPathIterateNeighborhoodFilter destruction
}


/*!
 * Member initialisation
 */

void ITKSpeedToPathIterateNeighborhoodFilter::Init()
{
////@begin ITKSpeedToPathIterateNeighborhoodFilter member initialisation
    mTermination = NULL;
    mStepLengthFactor = NULL;
    mSeedPoints = NULL;
    mPercentage = NULL;
////@end ITKSpeedToPathIterateNeighborhoodFilter member initialisation
}


/*!
 * Control creation for ITKSpeedToPathIterateNeighborhoodFilter
 */

void ITKSpeedToPathIterateNeighborhoodFilter::CreateControls()
{    
////@begin ITKSpeedToPathIterateNeighborhoodFilter content construction
    if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_ITKSPEEDTOPATHITERATENEIGHBORHOODFILTER")))
        wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
    mTermination = XRCCTRL(*this, "TEXT_TERMINATION", wxTextCtrl);
    mStepLengthFactor = XRCCTRL(*this, "ID_TEXT_STEP_LENGTH_FACTOR", wxTextCtrl);
    mSeedPoints = XRCCTRL(*this, "SEED_CHECKBOX", wxCheckBox);
    mPercentage = XRCCTRL(*this, "ID_GAUGE1", wxGauge);
    // Set validators
    if (FindWindow(XRCID("TEXT_TERMINATION")))
        FindWindow(XRCID("TEXT_TERMINATION"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
    if (FindWindow(XRCID("ID_TEXT_STEP_LENGTH_FACTOR")))
        FindWindow(XRCID("ID_TEXT_STEP_LENGTH_FACTOR"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
////@end ITKSpeedToPathIterateNeighborhoodFilter content construction

    // Create custom windows not generated automatically here.
////@begin ITKSpeedToPathIterateNeighborhoodFilter content initialisation
////@end ITKSpeedToPathIterateNeighborhoodFilter content initialisation
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for SPEED_BUTTON
 */
void ITKSpeedToPathIterateNeighborhoodFilter::OnCreateSpeedButtonClick( wxCommandEvent& event )
{
    SCIRun::ThrowSkinnerSignalEvent *tsse = new SCIRun::ThrowSkinnerSignalEvent("Painter::SetDataLayer");
  tsse->add_var( "ITKSpeedToPathIterateNeighborhoodImageFilterTool::mTermination", SCIRun::to_string( mTermination->GetValue() ) );
  tsse->add_var( "ITKSpeedToPathIterateNeighborhoodImageFilterTool::mStepLengthFactor", SCIRun::to_string( mStepLengthFactor->GetValue() ) );
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLEAR_SEEDS_BUTTON
 */
void ITKSpeedToPathIterateNeighborhoodFilter::OnClearSeedsButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::ThrowSkinnerSignal("Painter::SetLayer");
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for SEED_CHECKBOX
 */
void ITKSpeedToPathIterateNeighborhoodFilter::OnAddSeedsCheck( wxCommandEvent& event )
{
  SCIRun::ThrowSkinnerSignalEvent *tsse = new SCIRun::ThrowSkinnerSignalEvent("Painter::RedrawAll");
  tsse->add_var( "ITKSpeedToPathIterateNeighborhoodImageFilterTool::mSeedPoints",
                 SCIRun::to_string( mSeedPoints->GetValue() ) );
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void ITKSpeedToPathIterateNeighborhoodFilter::OnStartButtonClick( wxCommandEvent& event )
{
  SCIRun::ThrowSkinnerSignalEvent *tsse = new SCIRun::ThrowSkinnerSignalEvent("Painter::FinishTool");
  tsse->add_var( "ITKSpeedToPathIterateNeighborhoodImageFilterTool::mTermination", SCIRun::to_string( mTermination->GetValue() ) );
  tsse->add_var( "ITKSpeedToPathIterateNeighborhoodImageFilterTool::mStepLengthFactor", SCIRun::to_string( mStepLengthFactor->GetValue() ) );
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void ITKSpeedToPathIterateNeighborhoodFilter::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


void ITKSpeedToPathIterateNeighborhoodFilter::OnSetProgress( wxCommandEvent &event)
{
  int progress = event.GetInt();
  
  if (progress < 0)
  {
    disabler_ = new wxWindowDisabler();
#ifndef _WIN32
    wxBeginBusyCursor();
#endif
    progress = 0;
  }
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

bool ITKSpeedToPathIterateNeighborhoodFilter::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap ITKSpeedToPathIterateNeighborhoodFilter::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin ITKSpeedToPathIterateNeighborhoodFilter bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end ITKSpeedToPathIterateNeighborhoodFilter bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon ITKSpeedToPathIterateNeighborhoodFilter::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin ITKSpeedToPathIterateNeighborhoodFilter icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end ITKSpeedToPathIterateNeighborhoodFilter icon retrieval
}
