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

#ifndef _ITKCANNYEDGEDETECTIONIMAGEFILTER_H_
#define _ITKCANNYEDGEDETECTIONIMAGEFILTER_H_


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
#define ID_ITKCANNYEDGEDETECTIONIMAGEFILTER 10040
#define SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_STYLE 0
#define SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_TITLE _("itkCannyEdgeDetectionImageFilter")
#define SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_IDNAME ID_ITKCANNYEDGEDETECTIONIMAGEFILTER
#define SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_SIZE wxSize(200, -1)
#define SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_POSITION wxDefaultPosition
////@end control identifiers



/*!
 * itkCannyEdgeDetectionImageFilter class declaration
 */

class itkCannyEdgeDetectionImageFilter: public wxPanel
{    
  DECLARE_DYNAMIC_CLASS( itkCannyEdgeDetectionImageFilter )
    DECLARE_EVENT_TABLE()

    public:
  /// Constructors
  itkCannyEdgeDetectionImageFilter();
  itkCannyEdgeDetectionImageFilter( wxWindow* parent, wxWindowID id = SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_IDNAME, const wxPoint& pos = SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_POSITION, const wxSize& size = SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_SIZE, long style = SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_STYLE );

  /// Creation
  bool Create( wxWindow* parent, wxWindowID id = SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_IDNAME, const wxPoint& pos = SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_POSITION, const wxSize& size = SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_SIZE, long style = SYMBOL_ITKCANNYEDGEDETECTIONIMAGEFILTER_STYLE );

  /// Destructor
  ~itkCannyEdgeDetectionImageFilter();

  /// Initialises member variables
  void Init();

  /// Creates the controls and sizers
  void CreateControls();

  ////@begin itkCannyEdgeDetectionImageFilter event handler declarations

  ////@end itkCannyEdgeDetectionImageFilter event handler declarations

  ////@begin itkCannyEdgeDetectionImageFilter member function declarations

  /// Retrieves bitmap resources
  wxBitmap GetBitmapResource( const wxString& name );

  /// Retrieves icon resources
  wxIcon GetIconResource( const wxString& name );
  ////@end itkCannyEdgeDetectionImageFilter member function declarations

  void OnStartButtonClick( wxCommandEvent& event );
  void OnCloseButtonClick( wxCommandEvent& event );
  void OnSetProgress( wxCommandEvent &event);

  /// Should we show tooltips?
  static bool ShowToolTips();

  ////@begin itkCannyEdgeDetectionImageFilter member variables
  wxSpinCtrl* mVariance;
  wxSpinCtrl* mMaxError;
  wxSpinCtrl* mThreshold;
  wxGauge* mPercentage;
  ////@end itkCannyEdgeDetectionImageFilter member variables

  wxWindowDisabler *disabler_;
};

#endif
