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
// Name:        itkspeedtopathregularstepgradientdescentfilter.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     12/03/2008 15:08:43
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

#include "itkspeedtopathregularstepgradientdescentfilter.h"
#include "seg3devents.h"

////@begin XPM images
////@end XPM images


/*!
 * ITKSpeedToPathRegularStepGradientDescentFilter type definition
 */

IMPLEMENT_DYNAMIC_CLASS( ITKSpeedToPathRegularStepGradientDescentFilter, wxPanel )


/*!
 * ITKSpeedToPathRegularStepGradientDescentFilter event table definition
 */

BEGIN_EVENT_TABLE( ITKSpeedToPathRegularStepGradientDescentFilter, wxPanel )

////@begin ITKSpeedToPathRegularStepGradientDescentFilter event table entries
////@end ITKSpeedToPathRegularStepGradientDescentFilter event table entries

EVT_BUTTON( XRCID("START_BUTTON"), ITKSpeedToPathRegularStepGradientDescentFilter::OnStartButtonClick )
EVT_BUTTON( XRCID("CLOSE_BUTTON"), ITKSpeedToPathRegularStepGradientDescentFilter::OnCloseButtonClick )
EVT_BUTTON( XRCID("CLEAR_SEEDS_BUTTON"), ITKSpeedToPathRegularStepGradientDescentFilter::OnClearSeedsButtonClick )
EVT_CHECKBOX( XRCID("SEED_CHECKBOX"), ITKSpeedToPathRegularStepGradientDescentFilter::OnAddSeedsCheck )
EVT_BUTTON( XRCID("SPEED_BUTTON"), ITKSpeedToPathRegularStepGradientDescentFilter::OnCreateSpeedButtonClick )

EVT_COMMAND( wxID_ANY, wxEVT_SET_PROGRESS, ITKSpeedToPathRegularStepGradientDescentFilter::OnSetProgress)

END_EVENT_TABLE()


/*!
 * ITKSpeedToPathRegularStepGradientDescentFilter constructors
 */

ITKSpeedToPathRegularStepGradientDescentFilter::ITKSpeedToPathRegularStepGradientDescentFilter()
{
    Init();
}

ITKSpeedToPathRegularStepGradientDescentFilter::ITKSpeedToPathRegularStepGradientDescentFilter( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * ITKSpeedToPathRegularStepGradientDescentFilter creator
 */

bool ITKSpeedToPathRegularStepGradientDescentFilter::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin ITKSpeedToPathRegularStepGradientDescentFilter creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
////@end ITKSpeedToPathRegularStepGradientDescentFilter creation
    return true;
}


/*!
 * ITKSpeedToPathRegularStepGradientDescentFilter destructor
 */

ITKSpeedToPathRegularStepGradientDescentFilter::~ITKSpeedToPathRegularStepGradientDescentFilter()
{
////@begin ITKSpeedToPathRegularStepGradientDescentFilter destruction
////@end ITKSpeedToPathRegularStepGradientDescentFilter destruction
}


/*!
 * Member initialisation
 */

void ITKSpeedToPathRegularStepGradientDescentFilter::Init()
{
////@begin ITKSpeedToPathRegularStepGradientDescentFilter member initialisation
    mIterations = NULL;
    mTermination = NULL;
    mStepLengthFactor = NULL;
    mStepLengthRelax = NULL;
    mSeedPoints = NULL;
    mPercentage = NULL;
////@end ITKSpeedToPathRegularStepGradientDescentFilter member initialisation
	
}


/*!
 * Control creation for ITKSpeedToPathRegularStepGradientDescentFilter
 */

void ITKSpeedToPathRegularStepGradientDescentFilter::CreateControls()
{    
////@begin ITKSpeedToPathRegularStepGradientDescentFilter content construction
    if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_ITKSPEEDTOPATHREGULARSTEPGRADIENTDESCENTFILTER")))
        wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
    mIterations = XRCCTRL(*this, "SPIN_ITERATIONS", wxSpinCtrl);
    mTermination = XRCCTRL(*this, "TEXT_TERMINATION", wxTextCtrl);
    mStepLengthFactor = XRCCTRL(*this, "ID_TEXT_STEP_LENGTH_FACTOR", wxTextCtrl);
    mStepLengthRelax = XRCCTRL(*this, "ID_TEXT_STEP_LENGTH_RELAX", wxTextCtrl);
    mSeedPoints = XRCCTRL(*this, "SEED_CHECKBOX", wxCheckBox);
    mPercentage = XRCCTRL(*this, "ID_GAUGE1", wxGauge);
    // Set validators
    if (FindWindow(XRCID("TEXT_TERMINATION")))
        FindWindow(XRCID("TEXT_TERMINATION"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
    if (FindWindow(XRCID("ID_TEXT_STEP_LENGTH_FACTOR")))
        FindWindow(XRCID("ID_TEXT_STEP_LENGTH_FACTOR"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
    if (FindWindow(XRCID("ID_TEXT_STEP_LENGTH_RELAX")))
        FindWindow(XRCID("ID_TEXT_STEP_LENGTH_RELAX"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
////@end ITKSpeedToPathRegularStepGradientDescentFilter content construction

    // Create custom windows not generated automatically here.
////@begin ITKSpeedToPathRegularStepGradientDescentFilter content initialisation
////@end ITKSpeedToPathRegularStepGradientDescentFilter content initialisation
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for SPEED_BUTTON
 */
void ITKSpeedToPathRegularStepGradientDescentFilter::OnCreateSpeedButtonClick( wxCommandEvent& event )
{
  SCIRun::ThrowSkinnerSignalEvent *tsse = new SCIRun::ThrowSkinnerSignalEvent("Painter::SetDataLayer");
  tsse->add_var( "ITKSpeedToPathRegularStepGradientDescentImageFilterTool::mIterations", SCIRun::to_string( mIterations->GetValue() ) );
  tsse->add_var( "ITKSpeedToPathRegularStepGradientDescentImageFilterTool::mTermination", SCIRun::to_string( mTermination->GetValue() ) );
  tsse->add_var( "ITKSpeedToPathRegularStepGradientDescentImageFilterTool::mStepLengthFactor", SCIRun::to_string( mStepLengthFactor->GetValue() ) );
  tsse->add_var( "ITKSpeedToPathRegularStepGradientDescentImageFilterTool::mStepLengthRelax", SCIRun::to_string( mStepLengthRelax->GetValue() ) );
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLEAR_SEEDS_BUTTON
 */
void ITKSpeedToPathRegularStepGradientDescentFilter::OnClearSeedsButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::ThrowSkinnerSignal("Painter::SetLayer");
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for SEED_CHECKBOX
 */
void ITKSpeedToPathRegularStepGradientDescentFilter::OnAddSeedsCheck( wxCommandEvent& event )
{
  SCIRun::ThrowSkinnerSignalEvent *tsse = new SCIRun::ThrowSkinnerSignalEvent("Painter::RedrawAll");
  tsse->add_var( "ITKSpeedToPathRegularStepGradientDescentImageFilterTool::mSeedPoints",
                 SCIRun::to_string( mSeedPoints->GetValue() ) );
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void ITKSpeedToPathRegularStepGradientDescentFilter::OnStartButtonClick( wxCommandEvent& event )
{
  SCIRun::ThrowSkinnerSignalEvent *tsse = new SCIRun::ThrowSkinnerSignalEvent("Painter::FinishTool");
  tsse->add_var( "ITKSpeedToPathRegularStepGradientDescentImageFilterTool::mIterations", SCIRun::to_string( mIterations->GetValue() ) );
  tsse->add_var( "ITKSpeedToPathRegularStepGradientDescentImageFilterTool::mTermination", SCIRun::to_string( mTermination->GetValue() ) );
  tsse->add_var( "ITKSpeedToPathRegularStepGradientDescentImageFilterTool::mStepLengthFactor", SCIRun::to_string( mStepLengthFactor->GetValue() ) );
  tsse->add_var( "ITKSpeedToPathRegularStepGradientDescentImageFilterTool::mStepLengthRelax", SCIRun::to_string( mStepLengthRelax->GetValue() ) );
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void ITKSpeedToPathRegularStepGradientDescentFilter::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


void ITKSpeedToPathRegularStepGradientDescentFilter::OnSetProgress( wxCommandEvent &event)
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

bool ITKSpeedToPathRegularStepGradientDescentFilter::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap ITKSpeedToPathRegularStepGradientDescentFilter::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin ITKSpeedToPathRegularStepGradientDescentFilter bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end ITKSpeedToPathRegularStepGradientDescentFilter bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon ITKSpeedToPathRegularStepGradientDescentFilter::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin ITKSpeedToPathRegularStepGradientDescentFilter icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end ITKSpeedToPathRegularStepGradientDescentFilter icon retrieval
}
