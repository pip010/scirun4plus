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
// Name:        maskfilter.h
// Purpose:     
// Author:      
// Modified by: 
// Created:     Fri 28 Dec 2007 11:08:49 MST
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _MASKFILTER_H_
#define _MASKFILTER_H_


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
#define ID_MASKFILTER 10069
#define SYMBOL_MASKFILTER_STYLE wxTAB_TRAVERSAL
#define SYMBOL_MASKFILTER_TITLE _("MaskFilter")
#define SYMBOL_MASKFILTER_IDNAME ID_MASKFILTER
#define SYMBOL_MASKFILTER_SIZE wxSize(200, -1)
#define SYMBOL_MASKFILTER_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * MaskFilter class declaration
 */

class MaskFilter: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( MaskFilter )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    MaskFilter();
    MaskFilter( wxWindow* parent, wxWindowID id = SYMBOL_MASKFILTER_IDNAME, const wxPoint& pos = SYMBOL_MASKFILTER_POSITION, const wxSize& size = SYMBOL_MASKFILTER_SIZE, long style = SYMBOL_MASKFILTER_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_MASKFILTER_IDNAME, const wxPoint& pos = SYMBOL_MASKFILTER_POSITION, const wxSize& size = SYMBOL_MASKFILTER_SIZE, long style = SYMBOL_MASKFILTER_STYLE );

    /// Destructor
    ~MaskFilter();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

    void OnStartButtonClick( wxCommandEvent& event );
    void OnCloseButtonClick( wxCommandEvent& event );
    void OnSetProgress( wxCommandEvent& event );
    void OnSetMaskLayer( wxCommandEvent& event );
    void OnSetCurrentMin( wxCommandEvent& event ) { mMaskValue_ = 0; };
    void OnSetNewMin( wxCommandEvent& event ) { mMaskValue_ = 1; };
    void OnSetZero( wxCommandEvent& event ) { mMaskValue_ = 2; };
    void OnSetNewMax( wxCommandEvent& event ) { mMaskValue_ = 3; };
    void OnSetCurrentMax( wxCommandEvent& event ) { mMaskValue_ = 4; };

    void SetSkinnerCallback(const char *callback = "Painter::FinishTool");

////@begin MaskFilter event handler declarations

////@end MaskFilter event handler declarations

////@begin MaskFilter member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end MaskFilter member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin MaskFilter member variables
////@end MaskFilter member variables

    const char *skinner_callback_;

    unsigned int mMaskValue_;
};

#endif
    // _MASKFILTER_H_
