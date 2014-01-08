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
// Name:        optionlessfilter.h
// Purpose:     
// Author:      
// Modified by: 
// Created:     Wed 19 Dec 2007 14:11:29 MST
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _OPTIONLESSFILTER_H_
#define _OPTIONLESSFILTER_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/xrc/xmlres.h"
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
#define ID_OPTIONLESSFILTER 10069
#define SYMBOL_OPTIONLESSFILTER_STYLE wxTAB_TRAVERSAL
#define SYMBOL_OPTIONLESSFILTER_TITLE _("OptionlessFilter")
#define SYMBOL_OPTIONLESSFILTER_IDNAME ID_OPTIONLESSFILTER
#define SYMBOL_OPTIONLESSFILTER_SIZE wxSize(200, -1)
#define SYMBOL_OPTIONLESSFILTER_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * OptionlessFilter class declaration
 */

class OptionlessFilter: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( OptionlessFilter )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    OptionlessFilter();
    OptionlessFilter( wxWindow* parent, wxWindowID id = SYMBOL_OPTIONLESSFILTER_IDNAME, const wxPoint& pos = SYMBOL_OPTIONLESSFILTER_POSITION, const wxSize& size = SYMBOL_OPTIONLESSFILTER_SIZE, long style = SYMBOL_OPTIONLESSFILTER_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_OPTIONLESSFILTER_IDNAME, const wxPoint& pos = SYMBOL_OPTIONLESSFILTER_POSITION, const wxSize& size = SYMBOL_OPTIONLESSFILTER_SIZE, long style = SYMBOL_OPTIONLESSFILTER_STYLE );

    /// Destructor
    ~OptionlessFilter();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

    void OnStartButtonClick( wxCommandEvent& event );
    void OnCloseButtonClick( wxCommandEvent& event );
    void OnSetProgress( wxCommandEvent& event );

    void SetSkinnerCallback(const char *callback = "Painter::FinishTool");
    void SetShowProgress(bool show_progress);
  

////@begin OptionlessFilter event handler declarations

////@end OptionlessFilter event handler declarations

////@begin OptionlessFilter member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end OptionlessFilter member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin OptionlessFilter member variables
    wxGauge* mPercentage;
////@end OptionlessFilter member variables

    const char *skinner_callback_;
    bool show_progress_;
    wxWindowDisabler *disabler_;
};

#endif // _OPTIONLESSFILTER_H_
