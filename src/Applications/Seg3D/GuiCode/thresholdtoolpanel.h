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
// Name:        thresholdtoolpanel.h
// Purpose:     
// Author:      
// Modified by: 
// Created:     Wed 09 Apr 2008 17:18:05 MDT
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _THRESHOLDTOOLPANEL_H_
#define _THRESHOLDTOOLPANEL_H_


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
#define ID_THRESHOLDTOOLPANEL 10012
#define SYMBOL_THRESHOLDTOOLPANEL_STYLE 0
#define SYMBOL_THRESHOLDTOOLPANEL_TITLE _("ThresholdToolPanel")
#define SYMBOL_THRESHOLDTOOLPANEL_IDNAME ID_THRESHOLDTOOLPANEL
#define SYMBOL_THRESHOLDTOOLPANEL_SIZE wxSize(200, -1)
#define SYMBOL_THRESHOLDTOOLPANEL_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * ThresholdToolPanel class declaration
 */

class ThresholdToolPanel: public wxPanel
{    
  DECLARE_DYNAMIC_CLASS( ThresholdToolPanel )
    DECLARE_EVENT_TABLE()

    public:
  /// Constructors
  ThresholdToolPanel();
  ThresholdToolPanel( wxWindow* parent, wxWindowID id = SYMBOL_THRESHOLDTOOLPANEL_IDNAME, const wxPoint& pos = SYMBOL_THRESHOLDTOOLPANEL_POSITION, const wxSize& size = SYMBOL_THRESHOLDTOOLPANEL_SIZE, long style = SYMBOL_THRESHOLDTOOLPANEL_STYLE );

  /// Creation
  bool Create( wxWindow* parent, wxWindowID id = SYMBOL_THRESHOLDTOOLPANEL_IDNAME, const wxPoint& pos = SYMBOL_THRESHOLDTOOLPANEL_POSITION, const wxSize& size = SYMBOL_THRESHOLDTOOLPANEL_SIZE, long style = SYMBOL_THRESHOLDTOOLPANEL_STYLE );

  /// Destructor
  ~ThresholdToolPanel();

  /// Initialises member variables
  void Init();

  /// Creates the controls and sizers
  void CreateControls();
  void CreateHistogramImage();
  void OnClearSeedsButtonClick( wxCommandEvent& event );
  void OnStartButtonClick( wxCommandEvent& event );
  void OnCloseButtonClick( wxCommandEvent& event );
  void OnThresholdToolChange( wxCommandEvent& event );
  void OnLowerSlider( wxScrollEvent& event );
  void OnUpperSlider( wxScrollEvent& event );
  void OnLowerText ( wxCommandEvent& event );
  void OnUpperText ( wxCommandEvent& event );
  void OnPaint( wxPaintEvent &event );

  ////@begin ThresholdToolPanel event handler declarations

  ////@end ThresholdToolPanel event handler declarations

  ////@begin ThresholdToolPanel member function declarations

  /// Retrieves bitmap resources
  wxBitmap GetBitmapResource( const wxString& name );

  /// Retrieves icon resources
  wxIcon GetIconResource( const wxString& name );
  ////@end ThresholdToolPanel member function declarations

  /// Should we show tooltips?
  static bool ShowToolTips();

  void PushToSkinner();

  ////@begin ThresholdToolPanel member variables
  wxWindow* mHistoWindow;
  wxTextCtrl* mMinValue;
  wxTextCtrl* mMaxValue;
  wxTextCtrl* mLowerValue;
  wxSlider* mLowerSlider;
  wxTextCtrl* mUpperValue;
  wxSlider* mUpperSlider;
  ////@end ThresholdToolPanel member variables

  double minval_;
  double maxval_;
  double upper_;
  double lower_;
  vector<int> histogram_;

  unsigned int nbins_;
  unsigned char * histogram_image_data_;
};

#endif
    // _THRESHOLDTOOLPANEL_H_
