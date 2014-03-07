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

#ifndef _ITKDISCRETEGAUSSIANIMAGEFILTER_H_
#define _ITKDISCRETEGAUSSIANIMAGEFILTER_H_


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
#define ID_ITKDISCRETEGAUSSIANIMAGEFILTER 10039
#define SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_STYLE 0
#define SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_TITLE _("itkDiscreteGaussianImageFilter")
#define SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_IDNAME ID_ITKDISCRETEGAUSSIANIMAGEFILTER
#define SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_SIZE wxSize(200, -1)
#define SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_POSITION wxDefaultPosition
////@end control identifiers



/*!
 * itkDiscreteGaussianImageFilter class declaration
 */

class itkDiscreteGaussianImageFilter: public wxPanel
{    
    DECLARE_DYNAMIC_CLASS( itkDiscreteGaussianImageFilter )
    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    itkDiscreteGaussianImageFilter();
    itkDiscreteGaussianImageFilter( wxWindow* parent, wxWindowID id = SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_IDNAME, const wxPoint& pos = SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_POSITION, const wxSize& size = SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_SIZE, long style = SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_IDNAME, const wxPoint& pos = SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_POSITION, const wxSize& size = SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_SIZE, long style = SYMBOL_ITKDISCRETEGAUSSIANIMAGEFILTER_STYLE );

    /// Destructor
    ~itkDiscreteGaussianImageFilter();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin itkDiscreteGaussianImageFilter event handler declarations

////@end itkDiscreteGaussianImageFilter event handler declarations

////@begin itkDiscreteGaussianImageFilter member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
////@end itkDiscreteGaussianImageFilter member function declarations

	void OnStartButtonClick( wxCommandEvent& event );
	void OnCloseButtonClick( wxCommandEvent& event );
	void OnSetProgress( wxCommandEvent &event);

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin itkDiscreteGaussianImageFilter member variables
    wxSpinCtrl* mVariance;
    wxSpinCtrl* mMaxError;
    wxGauge* mPercentage;
////@end itkDiscreteGaussianImageFilter member variables

    wxWindowDisabler *disabler_;
};

#endif
