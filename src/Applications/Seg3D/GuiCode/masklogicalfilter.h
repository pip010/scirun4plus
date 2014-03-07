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
// Name:        masklogicalfilter.h
// Purpose:     
// Author:      
// Modified by: 
// Created:     Fri 28 Dec 2007 11:08:49 MST
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _MASKLOGICALFILTER_H_
#define _MASKLOGICALFILTER_H_


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
#define ID_MASKLOGICALFILTER 10070
#define SYMBOL_MASKLOGICALFILTER_STYLE wxTAB_TRAVERSAL
#define SYMBOL_MASKLOGICALFILTER_TITLE _("MaskLogicalFilter")
#define SYMBOL_MASKLOGICALFILTER_IDNAME ID_MASKLOGICALFILTER
#define SYMBOL_MASKLOGICALFILTER_SIZE wxSize(200, -1)
#define SYMBOL_MASKLOGICALFILTER_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * Masklogicalfilter class declaration
 */

class MaskLogicalFilter: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( Masklogicalfilter )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    MaskLogicalFilter();
    MaskLogicalFilter( wxWindow* parent, wxWindowID id = SYMBOL_MASKLOGICALFILTER_IDNAME, const wxPoint& pos = SYMBOL_MASKLOGICALFILTER_POSITION, const wxSize& size = SYMBOL_MASKLOGICALFILTER_SIZE, long style = SYMBOL_MASKLOGICALFILTER_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_MASKLOGICALFILTER_IDNAME, const wxPoint& pos = SYMBOL_MASKLOGICALFILTER_POSITION, const wxSize& size = SYMBOL_MASKLOGICALFILTER_SIZE, long style = SYMBOL_MASKLOGICALFILTER_STYLE );

    /// Destructor
    ~MaskLogicalFilter();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

    void OnStartButtonClick( wxCommandEvent& event );
    void OnCloseButtonClick( wxCommandEvent& event );
    void OnSetProgress( wxCommandEvent& event );
    void OnSetMaskLayer( wxCommandEvent& event );
    void SetSkinnerCallback(const char *callback = "Painter::FinishTool");

////@begin Masklogicalfilter event handler declarations

////@end Masklogicalfilter event handler declarations

////@begin Masklogicalfilter member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end Masklogicalfilter member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin Masklogicalfilter member variables
////@end Masklogicalfilter member variables

    const char *skinner_callback_;

    unsigned int mMaskValue_;
};

#endif
    // _MASKLOGICALFILTER_H_
