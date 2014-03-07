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
// Name:        masklogicalfilter.cpp
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

#include "masklogicalfilter.h"
#include <Applications/Seg3D/Painter.h>
#include "seg3devents.h"

////@begin XPM images
////@end XPM images


/*!
 * MaskLogicalFilter type definition
 */

IMPLEMENT_DYNAMIC_CLASS( MaskLogicalFilter, wxPanel )


/*!
 * MaskLogicalFilter event table definition
 */

BEGIN_EVENT_TABLE( MaskLogicalFilter, wxPanel )

////@begin MaskLogicalFilter event table entries
////@end MaskLogicalFilter event table entries

EVT_BUTTON( XRCID("START_BUTTON"), MaskLogicalFilter::OnStartButtonClick )
EVT_BUTTON( XRCID("CLOSE_BUTTON"), MaskLogicalFilter::OnCloseButtonClick )
EVT_BUTTON( XRCID("ID_SET_MASK_BUTTON"), MaskLogicalFilter::OnSetMaskLayer )

EVT_COMMAND( wxID_ANY, wxEVT_SET_PROGRESS, MaskLogicalFilter::OnSetProgress)

END_EVENT_TABLE()


/*!
 * MaskLogicalFilter constructors
 */

MaskLogicalFilter::MaskLogicalFilter()
{
  Init();
}

MaskLogicalFilter::MaskLogicalFilter( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
  Init();
  Create(parent, id, pos, size, style);
}


/*!
 * MaskLogicalFilter creator
 */

bool MaskLogicalFilter::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
  ////@begin MaskLogicalFilter creation
  SetParent(parent);
  CreateControls();
  if (GetSizer())
  {
    GetSizer()->SetSizeHints(this);
  }
  Centre();
  ////@end MaskLogicalFilter creation
  return true;
}


/*!
 * MaskLogicalFilter destructor
 */

MaskLogicalFilter::~MaskLogicalFilter()
{
  ////@begin MaskLogicalFilter destruction
  ////@end MaskLogicalFilter destruction
}


/*!
 * Member initialisation
 */

void MaskLogicalFilter::Init()
{
  ////@begin MaskLogicalFilter member initialisation
  ////@end MaskLogicalFilter member initialisation

  skinner_callback_ = "Painter::FinishTool";
}


/*!
 * Control creation for MaskLogicalFilter
 */

void MaskLogicalFilter::CreateControls()
{    
  ////@begin MaskLogicalFilter content construction
  if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_MASKLOGICALFILTER")))
    wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
  ////@end MaskLogicalFilter content construction

  // Create custom windows not generated automatically here.
  ////@begin MaskLogicalFilter content initialisation
  ////@end MaskLogicalFilter content initialisation
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void MaskLogicalFilter::OnStartButtonClick( wxCommandEvent& event )
{
  wxBeginBusyCursor();

  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent(skinner_callback_);

  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);

  wxEndBusyCursor();
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void MaskLogicalFilter::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


void MaskLogicalFilter::OnSetProgress( wxCommandEvent& event )
{
  //mPercentage->SetValue(event.GetInt());
}


void MaskLogicalFilter::OnSetMaskLayer( wxCommandEvent& event )
{
  SCIRun::Painter::ThrowSkinnerSignal("Painter::SetMaskLayer");
}

void MaskLogicalFilter::SetSkinnerCallback(const char *callback)
{
  skinner_callback_ = callback;
}

/*!
 * Should we show tooltips?
 */

bool MaskLogicalFilter::ShowToolTips()
{
  return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap MaskLogicalFilter::GetBitmapResource( const wxString& name )
{
  // Bitmap retrieval
  ////@begin MaskLogicalFilter bitmap retrieval
  wxUnusedVar(name);
  return wxNullBitmap;
  ////@end MaskLogicalFilter bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon MaskLogicalFilter::GetIconResource( const wxString& name )
{
  // Icon retrieval
  ////@begin MaskLogicalFilter icon retrieval
  wxUnusedVar(name);
  return wxNullIcon;
  ////@end MaskLogicalFilter icon retrieval
}
