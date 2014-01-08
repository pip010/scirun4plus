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
// Name:        medianfiltertool.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     12/12/2007 10:09:58
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

#include "medianfiltertool.h"
#include "seg3devents.h"

////@begin XPM images
////@end XPM images


/*!
 * MedianFilterTool type definition
 */

IMPLEMENT_DYNAMIC_CLASS( MedianFilterTool, wxPanel )


/*!
 * MedianFilterTool event table definition
 */

BEGIN_EVENT_TABLE( MedianFilterTool, wxPanel )

////@begin MedianFilterTool event table entries
////@end MedianFilterTool event table entries

    EVT_BUTTON( XRCID("START_BUTTON"), MedianFilterTool::OnStartButtonClick )
    EVT_BUTTON( XRCID("CLOSE_BUTTON"), MedianFilterTool::OnCloseButtonClick )
    EVT_COMMAND( wxID_ANY, wxEVT_SET_PROGRESS, MedianFilterTool::OnSetProgress)

END_EVENT_TABLE()


/*!
 * MedianFilterTool constructors
 */

MedianFilterTool::MedianFilterTool()
{
    Init();
}

MedianFilterTool::MedianFilterTool( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * MedianFilterTool creator
 */

bool MedianFilterTool::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin MedianFilterTool creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end MedianFilterTool creation
    return true;
}


/*!
 * MedianFilterTool destructor
 */

MedianFilterTool::~MedianFilterTool()
{
////@begin MedianFilterTool destruction
////@end MedianFilterTool destruction
}


/*!
 * Member initialisation
 */

void MedianFilterTool::Init()
{
////@begin MedianFilterTool member initialisation
    mRadius = NULL;
    mPercentage = NULL;
////@end MedianFilterTool member initialisation
}


/*!
 * Control creation for MedianFilterTool
 */

void MedianFilterTool::CreateControls()
{    
////@begin MedianFilterTool content construction
    if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_MEDIANFILTERTOOL")))
        wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
    mRadius = XRCCTRL(*this, "ID_Radius", wxSpinCtrl);
    mPercentage = XRCCTRL(*this, "ID_GAUGE", wxGauge);
////@end MedianFilterTool content construction

    // Disable the progress bar for now.
    mPercentage->Hide();

    // Create custom windows not generated automatically here.
////@begin MedianFilterTool content initialisation
////@end MedianFilterTool content initialisation
}


/*!
 * Should we show tooltips?
 */

bool MedianFilterTool::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap MedianFilterTool::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin MedianFilterTool bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end MedianFilterTool bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon MedianFilterTool::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin MedianFilterTool icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end MedianFilterTool icon retrieval
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void MedianFilterTool::OnStartButtonClick( wxCommandEvent& event )
{
  wxBeginBusyCursor();
  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent("Painter::MedianFilterVolume");
  tsse->add_var("MedianFilterVolume::radius",
                SCIRun::to_string(mRadius->GetValue()));
  SCIRun::Painter::ThrowSkinnerSignal(tsse);
  wxEndBusyCursor();
}


void MedianFilterTool::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}

void
MedianFilterTool::OnSetProgress( wxCommandEvent &event)
{	
  mPercentage->SetValue(event.GetInt());
}
