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
// Name:        polylinetoolpanel.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     Wed 09 Apr 2008 16:29:38 MDT
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

#include "polylinetoolpanel.h"

////@begin XPM images
////@end XPM images


/*!
 * PolylineToolPanel type definition
 */

IMPLEMENT_DYNAMIC_CLASS( PolylineToolPanel, wxPanel )


/*!
 * PolylineToolPanel event table definition
 */

BEGIN_EVENT_TABLE( PolylineToolPanel, wxPanel )

////@begin PolylineToolPanel event table entries
////@end PolylineToolPanel event table entries
    EVT_BUTTON( XRCID("CLEAR_SEEDS_BUTTON"), PolylineToolPanel::OnClearSeedsButtonClick )
    EVT_BUTTON( XRCID("START_BUTTON"), PolylineToolPanel::OnStartButtonClick )
    EVT_BUTTON( XRCID("ERASE_BUTTON"), PolylineToolPanel::OnEraseButtonClick )
    EVT_BUTTON( XRCID("CLOSE_BUTTON"), PolylineToolPanel::OnCloseButtonClick )

END_EVENT_TABLE()


/*!
 * PolylineToolPanel constructors
 */

PolylineToolPanel::PolylineToolPanel()
{
    Init();
}

PolylineToolPanel::PolylineToolPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * PolylineToolPanel creator
 */

bool PolylineToolPanel::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin PolylineToolPanel creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
////@end PolylineToolPanel creation
    return true;
}


/*!
 * PolylineToolPanel destructor
 */

PolylineToolPanel::~PolylineToolPanel()
{
////@begin PolylineToolPanel destruction
////@end PolylineToolPanel destruction
}


/*!
 * Member initialisation
 */

void PolylineToolPanel::Init()
{
////@begin PolylineToolPanel member initialisation
////@end PolylineToolPanel member initialisation
}


/*!
 * Control creation for PolylineToolPanel
 */

void PolylineToolPanel::CreateControls()
{    
////@begin PolylineToolPanel content construction
    if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_POLYLINETOOL")))
        wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
////@end PolylineToolPanel content construction

    // Create custom windows not generated automatically here.
////@begin PolylineToolPanel content initialisation
////@end PolylineToolPanel content initialisation
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLEAR_SEEDS_BUTTON
 */
void PolylineToolPanel::OnClearSeedsButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::ThrowSkinnerSignal("Painter::SetLayer");
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void PolylineToolPanel::OnStartButtonClick( wxCommandEvent& event )
{
  wxBeginBusyCursor();
  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent("Painter::FinishTool");
  tsse->add_var("Painter::polylinetool_erase", "0");
  SCIRun::Painter::ThrowSkinnerSignal(tsse);
  wxEndBusyCursor();
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for ERASE_BUTTON
 */
void PolylineToolPanel::OnEraseButtonClick( wxCommandEvent& event )
{
  wxBeginBusyCursor();
  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent("Painter::FinishTool");
  tsse->add_var("Painter::polylinetool_erase", "1");
  SCIRun::Painter::ThrowSkinnerSignal(tsse);
  wxEndBusyCursor();
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void PolylineToolPanel::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


/*!
 * Should we show tooltips?
 */

bool PolylineToolPanel::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap PolylineToolPanel::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin PolylineToolPanel bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end PolylineToolPanel bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon PolylineToolPanel::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin PolylineToolPanel icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end PolylineToolPanel icon retrieval
}
