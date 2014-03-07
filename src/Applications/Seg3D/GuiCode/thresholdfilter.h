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
// Name:        thresholdfilter.h
// Purpose:     
// Author:      
// Modified by: 
// Created:     19/12/2007 11:54:43
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _THRESHOLDFILTER_H_
#define _THRESHOLDFILTER_H_


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
#define ID_THRESHOLDFILTER 10039
#define SYMBOL_THRESHOLDFILTER_STYLE 0
#define SYMBOL_THRESHOLDFILTER_TITLE _("ThresholdFilter")
#define SYMBOL_THRESHOLDFILTER_IDNAME ID_THRESHOLDFILTER
#define SYMBOL_THRESHOLDFILTER_SIZE wxSize(200, -1)
#define SYMBOL_THRESHOLDFILTER_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * ThresholdFilter class declaration
 */

class ThresholdFilter: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( ThresholdFilter )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    ThresholdFilter();
    ThresholdFilter( wxWindow* parent, wxWindowID id = SYMBOL_THRESHOLDFILTER_IDNAME, const wxPoint& pos = SYMBOL_THRESHOLDFILTER_POSITION, const wxSize& size = SYMBOL_THRESHOLDFILTER_SIZE, long style = SYMBOL_THRESHOLDFILTER_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_THRESHOLDFILTER_IDNAME, const wxPoint& pos = SYMBOL_THRESHOLDFILTER_POSITION, const wxSize& size = SYMBOL_THRESHOLDFILTER_SIZE, long style = SYMBOL_THRESHOLDFILTER_STYLE );

    /// Destructor
    ~ThresholdFilter();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin ThresholdFilter event handler declarations

////@end ThresholdFilter event handler declarations

////@begin ThresholdFilter member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end ThresholdFilter member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin ThresholdFilter member variables
    wxTextCtrl* mMinValue;
    wxTextCtrl* mMaxValue;
    wxGauge* mPercentage;
////@end ThresholdFilter member variables
};

#endif
    // _THRESHOLDFILTER_H_
