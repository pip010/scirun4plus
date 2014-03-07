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

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"
#include <Applications/Seg3D/Painter.h>

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

////@begin includes
////@end includes

#include "itkDiscreteGaussianImageFilter.h"
#include "seg3devents.h"
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

////@begin XPM images
////@end XPM images

/*!
 * itkDiscreteGaussianImageFilter type definition
 */

IMPLEMENT_DYNAMIC_CLASS( itkDiscreteGaussianImageFilter, wxPanel )


/*!
 * itkDiscreteGaussianImageFilter event table definition
 */

BEGIN_EVENT_TABLE( itkDiscreteGaussianImageFilter, wxPanel )

////@begin itkDiscreteGaussianImageFilter event table entries
////@end itkDiscreteGaussianImageFilter event table entries
  

  EVT_BUTTON( XRCID("START_BUTTON"), itkDiscreteGaussianImageFilter::OnStartButtonClick )
  EVT_BUTTON( XRCID("CLOSE_BUTTON"), itkDiscreteGaussianImageFilter::OnCloseButtonClick )

  EVT_COMMAND( wxID_ANY, wxEVT_SET_PROGRESS, itkDiscreteGaussianImageFilter::OnSetProgress)

END_EVENT_TABLE()


/*!
 * itkDiscreteGaussianImageFilter constructors
 */

itkDiscreteGaussianImageFilter::itkDiscreteGaussianImageFilter()
{
    Init();
}

itkDiscreteGaussianImageFilter::itkDiscreteGaussianImageFilter( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, pos, size, style);
}


/*!
 * itkDiscreteGaussianImageFilter creator
 */

bool itkDiscreteGaussianImageFilter::Create( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style )
{
////@begin itkDiscreteGaussianImageFilter creation
    SetParent(parent);
    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
////@end itkDiscreteGaussianImageFilter creation
    return true;
}


/*!
 * itkDiscreteGaussianImageFilter destructor
 */

itkDiscreteGaussianImageFilter::~itkDiscreteGaussianImageFilter()
{
////@begin itkDiscreteGaussianImageFilter destruction
////@end itkDiscreteGaussianImageFilter destruction
}


/*!
 * Member initialisation
 */

void itkDiscreteGaussianImageFilter::Init()
{
////@begin itkDiscreteGaussianImageFilter member initialisation
    mVariance = NULL;
    mMaxError = NULL;
    mPercentage = NULL;
////@end itkDiscreteGaussianImageFilter member initialisation
}


/*!
 * Control creation for itkDiscreteGaussianImageFilter
 */

void itkDiscreteGaussianImageFilter::CreateControls()
{    
////@begin itkDiscreteGaussianImageFilter content construction
    if (!wxXmlResource::Get()->LoadPanel(this, GetParent(), _T("ID_ITKDISCRETEGAUSSIANIMAGEFILTER")))
        wxLogError(wxT("Missing wxXmlResource::Get()->Load() in OnInit()?"));
    mVariance = XRCCTRL(*this, "VARIANCE", wxSpinCtrl);
    mMaxError = XRCCTRL(*this, "MAXERROR", wxSpinCtrl);
    mPercentage = XRCCTRL(*this, "ID_GAUGE", wxGauge);
////@end itkDiscreteGaussianImageFilter content construction

    // Create custom windows not generated automatically here.
////@begin itkDiscreteGaussianImageFilter content initialisation
////@end itkDiscreteGaussianImageFilter content initialisation
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for START_BUTTON
 */
void itkDiscreteGaussianImageFilter::OnStartButtonClick( wxCommandEvent& event )
{
  SCIRun::ThrowSkinnerSignalEvent *tsse = new SCIRun::ThrowSkinnerSignalEvent("Painter::FinishTool");
  tsse->add_var("itkDiscreteGaussianImageFilter::mVariance",SCIRun::to_string(mVariance->GetValue()));
  tsse->add_var("itkDiscreteGaussianImageFilter::mMaxError",SCIRun::to_string(mMaxError->GetValue()));
  SCIRun::Painter::ThrowSkinnerSignal(tsse, false);
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for CLOSE_BUTTON
 */
void itkDiscreteGaussianImageFilter::OnCloseButtonClick( wxCommandEvent& event )
{
  SCIRun::Painter::global_seg3dframe_pointer_->HideTool();
}

void itkDiscreteGaussianImageFilter::OnSetProgress( wxCommandEvent &event)
{	
  int progress = event.GetInt();
  
  if (progress < 0)
  {
    disabler_ = new wxWindowDisabler();
    wxBeginBusyCursor();
    progress = 0;
  }
  if (progress > 100)
  {
    if (disabler_) { delete disabler_; disabler_ = 0; }
    wxEndBusyCursor();
    progress = 100;
  }
  if (progress == 0 || progress > mPercentage->GetValue())
  {
    mPercentage->SetValue(progress);
  }
}

/*!
 * Should we show tooltips?
 */

bool itkDiscreteGaussianImageFilter::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap itkDiscreteGaussianImageFilter::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin itkDiscreteGaussianImageFilter bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end itkDiscreteGaussianImageFilter bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon itkDiscreteGaussianImageFilter::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin itkDiscreteGaussianImageFilter icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end itkDiscreteGaussianImageFilter icon retrieval
}
