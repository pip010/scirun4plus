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
// Name:        seg3devents.cpp
// Purpose:     
// Author:      
// Modified by: 
// Created:     10/08/2007 11:53:59
// RCS-ID:      
// Copyright:   
// Licence:     
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif


#include "seg3devents.h"


DEFINE_EVENT_TYPE(wxEVT_CURSOR_INFORMATION_CHANGE)
DEFINE_EVENT_TYPE(wxEVT_BRUSH_RADIUS_CHANGE)
DEFINE_EVENT_TYPE(wxEVT_CROP_SET_RANGE)
DEFINE_EVENT_TYPE(wxEVT_CROP_SET_MINMAX)
DEFINE_EVENT_TYPE(wxEVT_SET_PROGRESS)
DEFINE_EVENT_TYPE(wxEVT_THRESHOLDTOOL_CHANGE)
DEFINE_EVENT_TYPE(wxEVT_WINDOWLEVELTOOL_CHANGE)
DEFINE_EVENT_TYPE(wxEVT_MOVESCALETOOL_CHANGE)
DEFINE_EVENT_TYPE(wxEVT_VOLUME_INFO_PANEL)
DEFINE_EVENT_TYPE(wxEVT_MEASUREMENT_UPDATE)
DEFINE_EVENT_TYPE(wxEVT_ARITHMETIC_SET_INPUTS)
DEFINE_EVENT_TYPE(wxEVT_OPACITY_PANEL)
DEFINE_EVENT_TYPE(wxEVT_OPACITY_SET_VALUE)

