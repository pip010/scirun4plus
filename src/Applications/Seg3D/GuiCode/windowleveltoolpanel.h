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
// Name:        windowleveltoolpanel.h
// Purpose:     
// Author:      
// Modified by: 
// Created:     Wed 09 Apr 2008 17:18:05 MDT
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _WINDOWLEVELTOOLPANEL_H_
#define _WINDOWLEVELTOOLPANEL_H_


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
#define ID_WINDOWLEVELTOOLPANEL 10013
#define SYMBOL_WINDOWLEVELTOOLPANEL_STYLE 0
#define SYMBOL_WINDOWLEVELTOOLPANEL_TITLE _("WindowlevelToolPanel")
#define SYMBOL_WINDOWLEVELTOOLPANEL_IDNAME ID_WINDOWLEVELTOOLPANEL
#define SYMBOL_WINDOWLEVELTOOLPANEL_SIZE wxSize(200, -1)
#define SYMBOL_WINDOWLEVELTOOLPANEL_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * WindowlevelToolPanel class declaration
 */

class WindowLevelToolPanel: public wxPanel
{    
  DECLARE_DYNAMIC_CLASS( WindowLevelToolPanel )
    DECLARE_EVENT_TABLE()

    public:
  /// Constructors
  WindowLevelToolPanel();
  WindowLevelToolPanel( wxWindow* parent, wxWindowID id = SYMBOL_WINDOWLEVELTOOLPANEL_IDNAME, const wxPoint& pos = SYMBOL_WINDOWLEVELTOOLPANEL_POSITION, const wxSize& size = SYMBOL_WINDOWLEVELTOOLPANEL_SIZE, long style = SYMBOL_WINDOWLEVELTOOLPANEL_STYLE );

  /// Creation
  bool Create( wxWindow* parent, wxWindowID id = SYMBOL_WINDOWLEVELTOOLPANEL_IDNAME, const wxPoint& pos = SYMBOL_WINDOWLEVELTOOLPANEL_POSITION, const wxSize& size = SYMBOL_WINDOWLEVELTOOLPANEL_SIZE, long style = SYMBOL_WINDOWLEVELTOOLPANEL_STYLE );

  /// Destructor
  ~WindowLevelToolPanel();

  /// Initialises member variables
  void Init();

  /// Creates the controls and sizers
  void CreateControls();
  void CreateHistogramImage();
  void OnCloseButtonClick( wxCommandEvent& event );
  void OnWindowLevelToolChange( wxCommandEvent& event );
  void OnWindowSlider( wxScrollEvent& event );
  void OnLevelSlider( wxScrollEvent& event );
  void OnWindowText ( wxCommandEvent& event );
  void OnLevelText ( wxCommandEvent& event );
  void OnPaint( wxPaintEvent &event );

  ////@begin WindowLevelToolPanel event handler declarations

  ////@end WindowLevelToolPanel event handler declarations

  ////@begin WindowLevelToolPanel member function declarations

  /// Retrieves bitmap resources
  wxBitmap GetBitmapResource( const wxString& name );

  /// Retrieves icon resources
  wxIcon GetIconResource( const wxString& name );
  ////@end WindowLevelToolPanel member function declarations

  /// Should we show tooltips?
  static bool ShowToolTips();

  void PushToSkinner();

  ////@begin WindowLevelToolPanel member variables
  wxWindow* mHistoWindow;
  wxTextCtrl* mMinValue;
  wxTextCtrl* mMaxValue;
  wxTextCtrl* mWindowValue;
  wxSlider* mWindowSlider;
  wxTextCtrl* mLevelValue;
  wxSlider* mLevelSlider;
  ////@end WindowLevelToolPanel member variables

  double minval_;
  double maxval_;
  double level_;
  double window_;
  vector<int> histogram_;

  unsigned int nbins_;

  unsigned char * histogram_image_data_;
};

#endif
    // _WINDOWLEVELTOOLPANEL_H_
