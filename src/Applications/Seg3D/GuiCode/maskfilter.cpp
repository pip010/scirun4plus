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
// Name:        maskfilter.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     Fri 28 Dec 2007 11:08:49 MST
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

////@begin includes
////@end includes

#include "maskfilter.h"
#include <Applications/Seg3D/Painter.h>
#include "seg3devents.h"

////@begin XPM images
////@end XPM images


/*!
 * MaskFilter type definition
 */

IMPLEMENT_DYNAMIC_CLASS( MaskFilter, wxPanel )


/*!
 * MaskFilter event table definition
 */

BEGIN_EVENT_TABLE( MaskFilter, wxPanel )

////@begin MaskFilter event table entries
////@end MaskFilter event table entries

EVT_BUTTON( XRCID("START_BUTTON"), MaskFilter::OnStartButtonClick )
EVT_BUTTON( XRCID("CLOSE_BUTTON"), MaskFilter::OnCloseButtonClick )
EVT_BUTTON( XRCID("ID_SET_MASK_BUTTON"), MaskFilter::OnSetMaskLayer )
EVT_RADIOBUTTON(XRCID("CURRENT_MIN_BUTTON"), MaskFilter::OnSetCurrentMin )
EVT_RADIOBUTTON(XRCID("NEW_MIN_BUTTON"), MaskFilter::OnSetNewMin )
EVT_RADIOBUTTON(XRCID("ZERO_BUTTON"), MaskFilter::OnSetZero )
EVT_RADIOBUTTON(XRCID("NEW_MAX_BUTTON"), MaskFilter::OnSetNewMax )
EVT_RADIOBUTTON(XRCID("CURRENT_MAX_BUTTON"), MaskFilter::OnSetCurrentMax )

EVT_COMMAND( wxID_ANY, wxEVT_SET_PROGRESS, MaskFilter::OnSetProgress)

END_EVENT_TABLE()


/*!
 * MaskFilter constructors
 */

MaskFilter::MaskFilter()
{
  Init();
}

MaskFilter::MaskFilter( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
  Init();
  Create(parent, id, pos, size, style);
}


/*!
 * MaskFilter creator
 */

bool MaskFilter::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
  ////@begin MaskFilter creation
  SetParent(parent);
  CreateControls();
  if (GetSizer())
  {
    GetSizer()->SetSizeHints(this);
  }
  Centre();
  ////@end MaskFilter creation
  return true;
}


/*!
 * MaskFilter destructor
 */

MaskFilter::~MaskFilter()
{
  ////@begin MaskFilter destruction
  ////@end MaskFilter destruction
}


/*!
 * Member initialisation
 */

void MaskFilter::Init()
{
  ////@begin MaskFilter member initialisation
  ////@end MaskFilter member initialisation

  skinner_callback_ = "Painter::FinishTool";

  mMaskValue_ = 2;
}


/*!
 * Control creation for MaskFilter
 */

void MaskFilter::CreateControls()
{    
  ////@begin MaskFilter content construction
  if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_MASKFILTER")))
    wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
  ////@end MaskFilter content construction

  // Create custom windows not generated automatically here.
  ////@begin MaskFilter content initialisation
  ////@end MaskFilter content initialisation
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void MaskFilter::OnStartButtonClick( wxCommandEvent& event )
{
  wxBeginBusyCursor();

  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent(skinner_callback_);

  tsse->add_var("MaskFilter::MaskValue", SCIRun::to_string(mMaskValue_));
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);

  wxEndBusyCursor();
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void MaskFilter::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


void MaskFilter::OnSetProgress( wxCommandEvent& event )
{
  //mPercentage->SetValue(event.GetInt());
}


void MaskFilter::OnSetMaskLayer( wxCommandEvent& event )
{
  SCIRun::Painter::ThrowSkinnerSignal("Painter::SetMaskLayer");
}

void MaskFilter::SetSkinnerCallback(const char *callback)
{
  skinner_callback_ = callback;
}

/*!
 * Should we show tooltips?
 */

bool MaskFilter::ShowToolTips()
{
  return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap MaskFilter::GetBitmapResource( const wxString& name )
{
  // Bitmap retrieval
  ////@begin MaskFilter bitmap retrieval
  wxUnusedVar(name);
  return wxNullBitmap;
  ////@end MaskFilter bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon MaskFilter::GetIconResource( const wxString& name )
{
  // Icon retrieval
  ////@begin MaskFilter icon retrieval
  wxUnusedVar(name);
  return wxNullIcon;
  ////@end MaskFilter icon retrieval
}
