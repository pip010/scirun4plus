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
//    File   : ITKSpeedToPathIterateNeighborhoodImageFilterTool.h
//    Author : David Brayford
//    Date   : May 2008

#ifndef SEG3D_SpeedToPath_h
#define SEG3D_SpeedToPath_h

#define NO_DUMMY_DECL

#include <sci_defs/insight_defs.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <itksys/SystemTools.hxx>

#include <Core/Events/BaseTool.h>
#include <Core/Datatypes/ITKDatatype.h>
#include <Core/Thread/ThreadLock.h>
#include <Applications/Seg3D/VolumeFilter.h>
#include <Applications/Seg3D/PolylineTool.h>

namespace SCIRun {

class Painter;
class NrrdVolume;

class SpeedToPathTool : public PolylineTool
{
public:
  SpeedToPathTool(const string &name, Painter *painter);
  virtual ~SpeedToPathTool();

  propagation_state_e   process_event(event_handle_t);

protected:
  virtual propagation_state_e   run_filter(bool erase);
  virtual void			draw_gl(SliceWindow &window);
  virtual void			seed_change_callback(int index);
  virtual int   compute_nearest_segment_index(const Point &curpos, int axis);

  // Compute the path from seed[i] to seed[i+1] and put it in paths[i].
  // Wraps to zero if i+1 overflows (it's a loop).
  virtual void                  compute_path(size_t seed_index);

  void generate_speed_function();

  NrrdVolumeHandle		source_volume_;
  vector<vector<Point> >        paths_;

  struct ltptr
  {
    bool operator()(const void* a, const void *b) const
    {
      return a < b;
    }
  };

  typedef map<void *, int, ltptr> window_slice_map_t;
  window_slice_map_t window_slice_map_;
};
  
  
}

#endif
