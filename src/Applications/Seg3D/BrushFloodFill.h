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
//    File   : BrushFloodFill.h
//    Author : McKay Davis
//    Date   : Sat Oct 14 15:51:56 2006

#ifndef SEG3D_BrushFloodFill_h
#define SEG3D_BrushFloodFill_h

#include <sci_defs/insight_defs.h>

#include <string>
#include <Core/Events/BaseTool.h>
#include <Core/Datatypes/ITKDatatype.h>
#include <Core/Thread/ThreadLock.h>
#include <Applications/Seg3D/VolumeFilter.h>
#include <Applications/Seg3D/SeedTool.h>

namespace SCIRun {

class Painter;
class NrrdVolume;
  
class BrushFloodFill : public SeedTool
{
public:
  BrushFloodFill(Painter *painter);
  unsigned int          value_;
  bool                  erase_;
  virtual void          run_filter();
  propagation_state_e   process_event(event_handle_t);

  BaseTool::propagation_state_e flood_fill_slice(bool erase);
  //BaseTool::propagation_state_e flood_fill_volume(bool erase);

  static void flood_fill(Painter *painter,
                         NrrdDataHandle dnrrd, unsigned int dlabel,
                         NrrdDataHandle snrrd, unsigned int slabel,
                         NrrdDataHandle mnrrd, unsigned int mlabel,
                         bool erase, const vector<vector<int> > &seeds);
};
  
  
}

#endif
