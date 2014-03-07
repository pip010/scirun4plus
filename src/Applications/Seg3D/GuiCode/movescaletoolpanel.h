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
// Name:        movescaletoolpanel.h
// Purpose:     
// Author:      
// Modified by: 
// Created:     Mon 21 Apr 2008 14:27:39 MDT
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _MOVESCALETOOLPANEL_H_
#define _MOVESCALETOOLPANEL_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/xrc/xmlres.h"
#include "wx/valtext.h"
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_MOVESCALETOOLPANEL 10012
#define SYMBOL_MOVESCALETOOLPANEL_STYLE 0
#define SYMBOL_MOVESCALETOOLPANEL_TITLE _("MoveScaleToolPanel")
#define SYMBOL_MOVESCALETOOLPANEL_IDNAME ID_MOVESCALETOOLPANEL
#define SYMBOL_MOVESCALETOOLPANEL_SIZE wxSize(200, -1)
#define SYMBOL_MOVESCALETOOLPANEL_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * MoveScaleToolPanel class declaration
 */

class MoveScaleToolPanel: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( MoveScaleToolPanel )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    MoveScaleToolPanel();
    MoveScaleToolPanel( wxWindow* parent, wxWindowID id = SYMBOL_MOVESCALETOOLPANEL_IDNAME, const wxPoint& pos = SYMBOL_MOVESCALETOOLPANEL_POSITION, const wxSize& size = SYMBOL_MOVESCALETOOLPANEL_SIZE, long style = SYMBOL_MOVESCALETOOLPANEL_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_MOVESCALETOOLPANEL_IDNAME, const wxPoint& pos = SYMBOL_MOVESCALETOOLPANEL_POSITION, const wxSize& size = SYMBOL_MOVESCALETOOLPANEL_SIZE, long style = SYMBOL_MOVESCALETOOLPANEL_STYLE );

    /// Destructor
    ~MoveScaleToolPanel();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

    void OnGetButtonClick( wxCommandEvent& event );
    void OnSetButtonClick( wxCommandEvent& event );
    void OnSetAllButtonClick( wxCommandEvent& event );
    void OnCloseButtonClick( wxCommandEvent& event );
    void OnMoveScaleToolChange( wxCommandEvent& event );

////@begin MoveScaleToolPanel event handler declarations

////@end MoveScaleToolPanel event handler declarations

////@begin MoveScaleToolPanel member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end MoveScaleToolPanel member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin MoveScaleToolPanel member variables
    wxTextCtrl* mXOrigin;
    wxTextCtrl* mYOrigin;
    wxTextCtrl* mZOrigin;
    wxTextCtrl* mXSpacing;
    wxTextCtrl* mYSpacing;
    wxTextCtrl* mZSpacing;
////@end MoveScaleToolPanel member variables
};

#endif
    // _MOVESCALETOOLPANEL_H_
