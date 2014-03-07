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
// Name:        itkspeedtopathgradientdescentfilter.h
// Purpose:     
// Author:      
// Modified by: 
// Created:     12/03/2008 13:58:25
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITKSPEEDTOPATHGRADIENTDESCENTFILTER_H_
#define _ITKSPEEDTOPATHGRADIENTDESCENTFILTER_H_


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
#define ID_ITKSPEEDTOPATHGRADIENTDESCENTFILTER 10013
#define SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_STYLE 0
#define SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_TITLE _("ITKSpeedToPathGradientDescentFilter")
#define SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_IDNAME ID_ITKSPEEDTOPATHGRADIENTDESCENTFILTER
#define SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_SIZE wxSize(200, -1)
#define SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * ITKSpeedToPathGradientDescentFilter class declaration
 */

class ITKSpeedToPathGradientDescentFilter: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( ITKSpeedToPathGradientDescentFilter )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    ITKSpeedToPathGradientDescentFilter();
    ITKSpeedToPathGradientDescentFilter( wxWindow* parent, wxWindowID id = SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_IDNAME, const wxPoint& pos = SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_POSITION, const wxSize& size = SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_SIZE, long style = SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_IDNAME, const wxPoint& pos = SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_POSITION, const wxSize& size = SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_SIZE, long style = SYMBOL_ITKSPEEDTOPATHGRADIENTDESCENTFILTER_STYLE );

    /// Destructor
    ~ITKSpeedToPathGradientDescentFilter();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();
	void OnStartButtonClick( wxCommandEvent& event );
	void OnCloseButtonClick( wxCommandEvent& event );
	void OnClearSeedsButtonClick( wxCommandEvent& event );
	void OnAddSeedsCheck( wxCommandEvent& event );
	void OnCreateSpeedButtonClick( wxCommandEvent& event );
	void OnSetProgress( wxCommandEvent& event );

////@begin ITKSpeedToPathGradientDescentFilter event handler declarations

////@end ITKSpeedToPathGradientDescentFilter event handler declarations

////@begin ITKSpeedToPathGradientDescentFilter member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end ITKSpeedToPathGradientDescentFilter member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin ITKSpeedToPathGradientDescentFilter member variables
    wxSpinCtrl* mIterations;
    wxTextCtrl* mTermination;
    wxCheckBox* mSeedPoints;
    wxGauge* mPercentage;
////@end ITKSpeedToPathGradientDescentFilter member variables
    
	wxWindowDisabler *disabler_;
};

#endif
    // _ITKSPEEDTOPATHGRADIENTDESCENTFILTER_H_
