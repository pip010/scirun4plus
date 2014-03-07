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
// Name:        histoeqfilter.h
// Purpose:     
// Author:      
// Modified by: 
// Created:     Thu 21 Feb 2008 13:36:14 MST
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _HISTOEQFILTER_H_
#define _HISTOEQFILTER_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/xrc/xmlres.h"
#include "wx/spinctrl.h"
#include "wx/valtext.h"
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class wxSpinCtrl;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_HISTOEQFILTER 10039
#define SYMBOL_HISTOEQFILTER_STYLE 0
#define SYMBOL_HISTOEQFILTER_TITLE _("HistoEqFilter")
#define SYMBOL_HISTOEQFILTER_IDNAME ID_HISTOEQFILTER
#define SYMBOL_HISTOEQFILTER_SIZE wxSize(200, -1)
#define SYMBOL_HISTOEQFILTER_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * HistoEqFilter class declaration
 */

class HistoEqFilter: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( HistoEqFilter )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    HistoEqFilter();
    HistoEqFilter( wxWindow* parent, wxWindowID id = SYMBOL_HISTOEQFILTER_IDNAME, const wxPoint& pos = SYMBOL_HISTOEQFILTER_POSITION, const wxSize& size = SYMBOL_HISTOEQFILTER_SIZE, long style = SYMBOL_HISTOEQFILTER_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_HISTOEQFILTER_IDNAME, const wxPoint& pos = SYMBOL_HISTOEQFILTER_POSITION, const wxSize& size = SYMBOL_HISTOEQFILTER_SIZE, long style = SYMBOL_HISTOEQFILTER_STYLE );

    /// Destructor
    ~HistoEqFilter();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

    void OnStartButtonClick( wxCommandEvent& event );
    void OnCloseButtonClick( wxCommandEvent& event );

////@begin HistoEqFilter event handler declarations

////@end HistoEqFilter event handler declarations

////@begin HistoEqFilter member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end HistoEqFilter member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin HistoEqFilter member variables
    wxSpinCtrl* mIgnore;
    wxTextCtrl* mAlpha;
////@end HistoEqFilter member variables
};

#endif
    // _HISTOEQFILTER_H_
