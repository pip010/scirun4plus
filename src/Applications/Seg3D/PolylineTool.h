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
//    File   : PolylineTool.h
//    Author : Michael Callahan
//    Date   : March 2008

#ifndef SEG3D_PolylineTool_h
#define SEG3D_PolylineTool_h

#include <sci_defs/insight_defs.h>

#include <string>
#include <Core/Events/BaseTool.h>

namespace SCIRun {

class Painter;
class NrrdVolume;


class PolylineTool : public virtual BaseTool, public PointerTool
{
public:
  PolylineTool(const string &name, Painter *painter);

  propagation_state_e   pointer_down(int, int, int, unsigned int, int);
  propagation_state_e   pointer_up(int, int, int, unsigned int, int);
  propagation_state_e   pointer_motion(int, int, int, unsigned int, int);

  propagation_state_e   process_event(event_handle_t);

  virtual propagation_state_e   run_filter(bool erase);

protected:
  typedef vector<Point> seeds_t;

  virtual void  draw_gl(SliceWindow &window);
  virtual void	seed_change_callback(int index);
  virtual int   compute_nearest_segment_index(const Point &curpos, int axis);

  void                  compute_bbox(const vector<vector<float> > &seeds,
                                     float &x0, float &y0,
                                     float &x1, float &y1);

  bool                  check_on_boundary(const vector<vector<float> > &seeds,
                                          float x, float y);
  bool                  check_crossings(const vector<vector<float> > &seeds,
                                        float x, float y);
  void                  rasterize(NrrdDataHandle dnrrd, unsigned int dlabel,
                                  NrrdDataHandle mnrrd, unsigned int mlabel,
                                  const vector<vector<float> > &seeds,
                                  bool erase);

  Painter *             painter_;

  seeds_t               seeds_;
  ThreadLock            seed_lock_;

  // For 2D seed operations, do the operation where the seed was put.
  SliceWindow          *last_seed_window_;

private:
  int                   last_seed_index_;
  seeds_t               drag_seeds_;
  Point                 drag_point_;
};
  
  
}

#endif // PolylineTool_h
