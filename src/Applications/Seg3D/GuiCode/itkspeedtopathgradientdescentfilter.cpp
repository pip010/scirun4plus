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
// Name:        itkspeedtopathgradientdescentfilter.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     12/03/2008 13:58:25
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

#include "itkspeedtopathgradientdescentfilter.h"
#include "seg3devents.h"

////@begin XPM images
////@end XPM images


/*!
 * ITKSpeedToPathGradientDescentFilter type definition
 */

IMPLEMENT_DYNAMIC_CLASS( ITKSpeedToPathGradientDescentFilter, wxPanel )


/*!
 * ITKSpeedToPathGradientDescentFilter event table definition
 */

BEGIN_EVENT_TABLE( ITKSpeedToPathGradientDescentFilter, wxPanel )

EVT_BUTTON( XRCID("START_BUTTON"), ITKSpeedToPathGradientDescentFilter::OnStartButtonClick )
EVT_BUTTON( XRCID("CLOSE_BUTTON"), ITKSpeedToPathGradientDescentFilter::OnCloseButtonClick )
EVT_BUTTON( XRCID("CLEAR_SEEDS_BUTTON"), ITKSpeedToPathGradientDescentFilter::OnClearSeedsButtonClick )
EVT_CHECKBOX( XRCID("SEED_CHECKBOX"), ITKSpeedToPathGradientDescentFilter::OnAddSeedsCheck )
EVT_BUTTON( XRCID("SPEED_BUTTON"), ITKSpeedToPathGradientDescentFilter::OnCreateSpeedButtonClick )

EVT_COMMAND( wxID_ANY, wxEVT_SET_PROGRESS, ITKSpeedToPathGradientDescentFilter::OnSetProgress)

////@begin ITKSpeedToPathGradientDescentFilter event table entries
////@end ITKSpeedToPathGradientDescentFilter event table entries

END_EVENT_TABLE()


/*!
 * ITKSpeedToPathGradientDescentFilter constructors
 */

ITKSpeedToPathGradientDescentFilter::ITKSpeedToPathGradientDescentFilter()
{
    Init();
}

ITKSpeedToPathGradientDescentFilter::ITKSpeedToPathGradientDescentFilter( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * ITKSpeedToPathGradientDescentFilter creator
 */

bool ITKSpeedToPathGradientDescentFilter::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin ITKSpeedToPathGradientDescentFilter creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
////@end ITKSpeedToPathGradientDescentFilter creation
    return true;
}


/*!
 * ITKSpeedToPathGradientDescentFilter destructor
 */

ITKSpeedToPathGradientDescentFilter::~ITKSpeedToPathGradientDescentFilter()
{
////@begin ITKSpeedToPathGradientDescentFilter destruction
////@end ITKSpeedToPathGradientDescentFilter destruction
}


/*!
 * Member initialisation
 */

void ITKSpeedToPathGradientDescentFilter::Init()
{
////@begin ITKSpeedToPathGradientDescentFilter member initialisation
    mIterations = NULL;
    mTermination = NULL;
    mSeedPoints = NULL;
    mPercentage = NULL;
////@end ITKSpeedToPathGradientDescentFilter member initialisation
}


/*!
 * Control creation for ITKSpeedToPathGradientDescentFilter
 */

void ITKSpeedToPathGradientDescentFilter::CreateControls()
{    
////@begin ITKSpeedToPathGradientDescentFilter content construction
    if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_ITKSPEEDTOPATHGRADIENTDESCENTFILTER")))
        wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
    mIterations = XRCCTRL(*this, "SPIN_ITERATIONS", wxSpinCtrl);
    mTermination = XRCCTRL(*this, "TEXT_TERMINATION", wxTextCtrl);
    mSeedPoints = XRCCTRL(*this, "SEED_CHECKBOX", wxCheckBox);
    mPercentage = XRCCTRL(*this, "ID_GAUGE1", wxGauge);
    // Set validators
    if (FindWindow(XRCID("TEXT_TERMINATION")))
        FindWindow(XRCID("TEXT_TERMINATION"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
////@end ITKSpeedToPathGradientDescentFilter content construction

    // Create custom windows not generated automatically here.
////@begin ITKSpeedToPathGradientDescentFilter content initialisation
////@end ITKSpeedToPathGradientDescentFilter content initialisation
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for SPEED_BUTTON
 */
void ITKSpeedToPathGradientDescentFilter::OnCreateSpeedButtonClick( wxCommandEvent& event )
{
	  SCIRun::ThrowSkinnerSignalEvent *tsse = new SCIRun::ThrowSkinnerSignalEvent("Painter::SetDataLayer");
  tsse->add_var( "ITKSpeedToPathGradientDescentImageFilterTool::mIterations", SCIRun::to_string( mIterations->GetValue() ) );
  tsse->add_var( "ITKSpeedToPathGradientDescentImageFilterTool::mTermination", SCIRun::to_string( mTermination->GetValue() ) );
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLEAR_SEEDS_BUTTON
 */
void ITKSpeedToPathGradientDescentFilter::OnClearSeedsButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::ThrowSkinnerSignal("Painter::SetLayer");
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for SEED_CHECKBOX
 */
void ITKSpeedToPathGradientDescentFilter::OnAddSeedsCheck( wxCommandEvent& event )
{
  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent("Painter::RedrawAll");
  tsse->add_var( "ITKSpeedToPathGradientDescentImageFilterTool::mSeedPoints",
                 SCIRun::to_string( mSeedPoints->GetValue() ) );
  SCIRun::Painter::ThrowSkinnerSignal(tsse);
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void ITKSpeedToPathGradientDescentFilter::OnStartButtonClick( wxCommandEvent& event )
{
  SCIRun::ThrowSkinnerSignalEvent *tsse = new SCIRun::ThrowSkinnerSignalEvent("Painter::FinishTool");
  tsse->add_var( "ITKSpeedToPathGradientDescentImageFilterTool::mIterations", SCIRun::to_string( mIterations->GetValue() ) );
  tsse->add_var( "ITKSpeedToPathGradientDescentImageFilterTool::mTermination", SCIRun::to_string( mTermination->GetValue() ) );
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void ITKSpeedToPathGradientDescentFilter::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


void ITKSpeedToPathGradientDescentFilter::OnSetProgress( wxCommandEvent &event)
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

bool ITKSpeedToPathGradientDescentFilter::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap ITKSpeedToPathGradientDescentFilter::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin ITKSpeedToPathGradientDescentFilter bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end ITKSpeedToPathGradientDescentFilter bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon ITKSpeedToPathGradientDescentFilter::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin ITKSpeedToPathGradientDescentFilter icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end ITKSpeedToPathGradientDescentFilter icon retrieval
}
