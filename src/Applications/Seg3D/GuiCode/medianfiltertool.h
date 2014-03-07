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
// Name:        medianfiltertool.h
// Purpose:     
// Author:      
// Modified by: 
// Created:     12/12/2007 10:09:58
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _MEDIANFILTERTOOL_H_
#define _MEDIANFILTERTOOL_H_


/*!
 * Includes
 */

////@begin includes
#include "wx/xrc/xmlres.h"
#include "wx/spinctrl.h"
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
#define ID_MEDIANFILTERTOOL 10062
#define SYMBOL_MEDIANFILTERTOOL_STYLE wxTAB_TRAVERSAL
#define SYMBOL_MEDIANFILTERTOOL_TITLE _("MedianFilterTool")
#define SYMBOL_MEDIANFILTERTOOL_IDNAME ID_MEDIANFILTERTOOL
#define SYMBOL_MEDIANFILTERTOOL_SIZE wxSize(200, -1)
#define SYMBOL_MEDIANFILTERTOOL_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * MedianFilterTool class declaration
 */

class MedianFilterTool: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( MedianFilterTool )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    MedianFilterTool();
    MedianFilterTool( wxWindow* parent, wxWindowID id = SYMBOL_MEDIANFILTERTOOL_IDNAME, const wxPoint& pos = SYMBOL_MEDIANFILTERTOOL_POSITION, const wxSize& size = SYMBOL_MEDIANFILTERTOOL_SIZE, long style = SYMBOL_MEDIANFILTERTOOL_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_MEDIANFILTERTOOL_IDNAME, const wxPoint& pos = SYMBOL_MEDIANFILTERTOOL_POSITION, const wxSize& size = SYMBOL_MEDIANFILTERTOOL_SIZE, long style = SYMBOL_MEDIANFILTERTOOL_STYLE );

    /// Destructor
    ~MedianFilterTool();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin MedianFilterTool event handler declarations

////@end MedianFilterTool event handler declarations

////@begin MedianFilterTool member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end MedianFilterTool member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

    void OnStartButtonClick( wxCommandEvent& event );
    void OnCloseButtonClick( wxCommandEvent& event );
    void OnSetProgress( wxCommandEvent& event );

////@begin MedianFilterTool member variables
    wxSpinCtrl* mRadius;
    wxGauge* mPercentage;
////@end MedianFilterTool member variables
};

#endif
    // _MEDIANFILTERTOOL_H_
