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
// Name:        histoeqfilter.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     Thu 21 Feb 2008 13:36:14 MST
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

#include "histoeqfilter.h"
#include "seg3devents.h"
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

////@begin XPM images
////@end XPM images


/*!
 * HistoEqFilter type definition
 */

IMPLEMENT_DYNAMIC_CLASS( HistoEqFilter, wxPanel )


/*!
 * HistoEqFilter event table definition
 */

BEGIN_EVENT_TABLE( HistoEqFilter, wxPanel )

////@begin HistoEqFilter event table entries
////@end HistoEqFilter event table entries

    EVT_BUTTON( XRCID("START_BUTTON"), HistoEqFilter::OnStartButtonClick )
    EVT_BUTTON( XRCID("CLOSE_BUTTON"), HistoEqFilter::OnCloseButtonClick )

END_EVENT_TABLE()


/*!
 * HistoEqFilter constructors
 */

HistoEqFilter::HistoEqFilter()
{
    Init();
}

HistoEqFilter::HistoEqFilter( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * HistoEqFilter creator
 */

bool HistoEqFilter::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin HistoEqFilter creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
////@end HistoEqFilter creation
    return true;
}


/*!
 * HistoEqFilter destructor
 */

HistoEqFilter::~HistoEqFilter()
{
////@begin HistoEqFilter destruction
////@end HistoEqFilter destruction
}


/*!
 * Member initialisation
 */

void HistoEqFilter::Init()
{
////@begin HistoEqFilter member initialisation
    mIgnore = NULL;
    mAlpha = NULL;
////@end HistoEqFilter member initialisation
}


/*!
 * Control creation for HistoEqFilter
 */

void HistoEqFilter::CreateControls()
{    
////@begin HistoEqFilter content construction
    if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_HISTOEQFILTER")))
        wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
    mIgnore = XRCCTRL(*this, "IGNORE", wxSpinCtrl);
    mAlpha = XRCCTRL(*this, "ID_TEXTCTRL", wxTextCtrl);
    // Set validators
    if (FindWindow(XRCID("ID_TEXTCTRL")))
        FindWindow(XRCID("ID_TEXTCTRL"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
////@end HistoEqFilter content construction

    // Create custom windows not generated automatically here.
////@begin HistoEqFilter content initialisation
////@end HistoEqFilter content initialisation
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void HistoEqFilter::OnStartButtonClick( wxCommandEvent& event )
{
  // Call synchronously and wait.
  wxBeginBusyCursor();
  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent("Painter::HistoEqFilter");
  tsse->add_var("Painter::HistoEq::bins",
                SCIRun::to_string(mIgnore->GetValue()));
  tsse->add_var("Painter::HistoEq::alpha",
                wx2std(mAlpha->GetValue()));
  SCIRun::Painter::ThrowSkinnerSignal(tsse);
  wxEndBusyCursor();
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void HistoEqFilter::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


/*!
 * Should we show tooltips?
 */

bool HistoEqFilter::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap HistoEqFilter::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin HistoEqFilter bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end HistoEqFilter bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon HistoEqFilter::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin HistoEqFilter icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end HistoEqFilter icon retrieval
}
