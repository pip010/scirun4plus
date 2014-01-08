//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : VolumeFilter.cc
//    Author : McKay Davis
//    Date   : Feburary 2005

#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/VolumeFilter.h>

#include <sci_comp_warn_fixes.h>


namespace SCIRun {


VolumeFilterBase::VolumeFilterBase(NrrdVolumeHandle volume,
                                   Skinner::Variables *vars)
  : Skinner::Parent(vars ? vars : new Skinner::Variables("VolumeFilter", 0)),
    volume_(volume),
    abort_(false),
    progress_(0.0),
    progress_offset_(0.0),
    progress_show_(true),
    passes_(1)
{
}
  
    
VolumeFilterBase::~VolumeFilterBase()
{
}


void
VolumeFilterBase::operator()(NrrdDataHandle nrrdh)
{
}

void
VolumeFilterBase::start(NrrdDataHandle nrrdh)
{
}


void
VolumeFilterBase::update()
{
}


void
VolumeFilterBase::update_progress()
{
  const int progress = (int)(progress_ * 100 + progress_offset_);

  if( progress_show_ )
    volume_->painter_->update_progress(progress);
}


}
