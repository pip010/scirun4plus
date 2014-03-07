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
// Name:        movescaletoolpanel.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     Mon 21 Apr 2008 14:27:39 MDT
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

#include "movescaletoolpanel.h"
#include "seg3devents.h"
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

////@begin XPM images
////@end XPM images


/*!
 * MoveScaleToolPanel type definition
 */

IMPLEMENT_DYNAMIC_CLASS( MoveScaleToolPanel, wxPanel )


/*!
 * MoveScaleToolPanel event table definition
 */

BEGIN_EVENT_TABLE( MoveScaleToolPanel, wxPanel )

////@begin MoveScaleToolPanel event table entries
////@end MoveScaleToolPanel event table entries

  EVT_BUTTON( XRCID("GET_BUTTON"), MoveScaleToolPanel::OnGetButtonClick )
  EVT_BUTTON( XRCID("SET_BUTTON"), MoveScaleToolPanel::OnSetButtonClick )
  EVT_BUTTON( XRCID("SET_ALL_BUTTON"), MoveScaleToolPanel::OnSetAllButtonClick )
  EVT_BUTTON( XRCID("CLOSE_BUTTON"), MoveScaleToolPanel::OnCloseButtonClick )
  EVT_COMMAND( wxID_ANY, wxEVT_MOVESCALETOOL_CHANGE, MoveScaleToolPanel::OnMoveScaleToolChange )

END_EVENT_TABLE()


/*!
 * MoveScaleToolPanel constructors
 */

MoveScaleToolPanel::MoveScaleToolPanel()
{
    Init();
}

MoveScaleToolPanel::MoveScaleToolPanel( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * MoveScaleToolPanel creator
 */

bool MoveScaleToolPanel::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin MoveScaleToolPanel creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
////@end MoveScaleToolPanel creation
    return true;
}


/*!
 * MoveScaleToolPanel destructor
 */

MoveScaleToolPanel::~MoveScaleToolPanel()
{
////@begin MoveScaleToolPanel destruction
////@end MoveScaleToolPanel destruction
}


/*!
 * Member initialisation
 */

void MoveScaleToolPanel::Init()
{
////@begin MoveScaleToolPanel member initialisation
    mXOrigin = NULL;
    mYOrigin = NULL;
    mZOrigin = NULL;
    mXSpacing = NULL;
    mYSpacing = NULL;
    mZSpacing = NULL;
////@end MoveScaleToolPanel member initialisation
}


/*!
 * Control creation for MoveScaleToolPanel
 */

void MoveScaleToolPanel::CreateControls()
{    
////@begin MoveScaleToolPanel content construction
    if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_MOVESCALETOOLPANEL")))
        wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
    mXOrigin = XRCCTRL(*this, "ID_TEXTCTRL5", wxTextCtrl);
    mYOrigin = XRCCTRL(*this, "ID_TEXTCTRL4", wxTextCtrl);
    mZOrigin = XRCCTRL(*this, "ID_TEXTCTRL3", wxTextCtrl);
    mXSpacing = XRCCTRL(*this, "ID_TEXTCTRL2", wxTextCtrl);
    mYSpacing = XRCCTRL(*this, "ID_TEXTCTRL1", wxTextCtrl);
    mZSpacing = XRCCTRL(*this, "ID_TEXTCTRL", wxTextCtrl);
    // Set validators
    if (FindWindow(XRCID("ID_TEXTCTRL5")))
        FindWindow(XRCID("ID_TEXTCTRL5"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
    if (FindWindow(XRCID("ID_TEXTCTRL4")))
        FindWindow(XRCID("ID_TEXTCTRL4"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
    if (FindWindow(XRCID("ID_TEXTCTRL3")))
        FindWindow(XRCID("ID_TEXTCTRL3"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
    if (FindWindow(XRCID("ID_TEXTCTRL2")))
        FindWindow(XRCID("ID_TEXTCTRL2"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
    if (FindWindow(XRCID("ID_TEXTCTRL1")))
        FindWindow(XRCID("ID_TEXTCTRL1"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
    if (FindWindow(XRCID("ID_TEXTCTRL")))
        FindWindow(XRCID("ID_TEXTCTRL"))->SetValidator( wxTextValidator(wxFILTER_NUMERIC) );
////@end MoveScaleToolPanel content construction

    // Create custom windows not generated automatically here.
////@begin MoveScaleToolPanel content initialisation
////@end MoveScaleToolPanel content initialisation
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for GET_BUTTON
 */
void MoveScaleToolPanel::OnGetButtonClick( wxCommandEvent& event )
{
  wxBeginBusyCursor();
  SCIRun::ThrowSkinnerSignalEvent *tsse = new SCIRun::ThrowSkinnerSignalEvent("Painter::MoveScaleLayerUpdateGUI");
  SCIRun::Painter::ThrowSkinnerSignal(tsse);
  wxEndBusyCursor();
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for SET_BUTTON
 */
void MoveScaleToolPanel::OnSetButtonClick( wxCommandEvent& event )
{
  wxBeginBusyCursor();
  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent("Painter::MoveScaleLayer");
  tsse->add_var("Painter::movescale::spacing::x",
                wx2std(mXSpacing->GetValue()));
  tsse->add_var("Painter::movescale::spacing::y",
                wx2std(mYSpacing->GetValue()));
  tsse->add_var("Painter::movescale::spacing::z",
                wx2std(mZSpacing->GetValue()));
  tsse->add_var("Painter::movescale::origin::x",
                wx2std(mXOrigin->GetValue()));
  tsse->add_var("Painter::movescale::origin::y",
                wx2std(mYOrigin->GetValue()));
  tsse->add_var("Painter::movescale::origin::z",
                wx2std(mZOrigin->GetValue()));
  SCIRun::Painter::ThrowSkinnerSignal(tsse);
  wxEndBusyCursor();
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for SET_ALL_BUTTON
 */
void MoveScaleToolPanel::OnSetAllButtonClick( wxCommandEvent& event )
{
  wxBeginBusyCursor();
  SCIRun::ThrowSkinnerSignalEvent *tsse =
    new SCIRun::ThrowSkinnerSignalEvent("Painter::MoveScaleAllLayers");
  tsse->add_var("Painter::movescale::spacing::x",
                wx2std(mXSpacing->GetValue()));
  tsse->add_var("Painter::movescale::spacing::y",
                wx2std(mYSpacing->GetValue()));
  tsse->add_var("Painter::movescale::spacing::z",
                wx2std(mZSpacing->GetValue()));
  tsse->add_var("Painter::movescale::origin::x",
                wx2std(mXOrigin->GetValue()));
  tsse->add_var("Painter::movescale::origin::y",
                wx2std(mYOrigin->GetValue()));
  tsse->add_var("Painter::movescale::origin::z",
                wx2std(mZOrigin->GetValue()));
  SCIRun::Painter::ThrowSkinnerSignal(tsse);
  wxEndBusyCursor();
}


/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void MoveScaleToolPanel::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}


void MoveScaleToolPanel::OnMoveScaleToolChange( wxCommandEvent& event )
{
  MoveScaleToolChangeStruct *info =
    (MoveScaleToolChangeStruct *)event.GetClientData();

  mXSpacing->Clear(); *mXSpacing << info->spacingx;
  mYSpacing->Clear(); *mYSpacing << info->spacingy;
  mZSpacing->Clear(); *mZSpacing << info->spacingz;

  mXOrigin->Clear(); *mXOrigin << info->originx;
  mYOrigin->Clear(); *mYOrigin << info->originy;
  mZOrigin->Clear(); *mZOrigin << info->originz;

  delete info;
}


/*!
 * Should we show tooltips?
 */

bool MoveScaleToolPanel::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap MoveScaleToolPanel::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin MoveScaleToolPanel bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end MoveScaleToolPanel bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon MoveScaleToolPanel::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin MoveScaleToolPanel icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end MoveScaleToolPanel icon retrieval
}
