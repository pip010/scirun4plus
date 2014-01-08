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


/*
 *  ColorMappedNrrdTextureObj.h
 *
 *  Written by:
 *   McKay Davis
 *   School of Computing
 *   University of Utah
 *   November, 2005
 *
 */

#ifndef SCIRun_Dataflow_Modules_Render_ColorMappedNrrdTextureObj_h
#define SCIRun_Dataflow_Modules_Render_ColorMappedNrrdTextureObj_h

#include <string>
#include <vector>

#include <Core/Datatypes/NrrdData.h>
#include <Core/Datatypes/ColorMap.h>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/BBox.h>
#include <Core/Containers/LockingHandle.h>
#include <Core/Thread/Mutex.h>
#include <Core/Geom/share.h>
namespace SCIRun {

  class SCISHARE ColorMappedNrrdTextureObj : public UsedWithLockingHandle<Mutex>
{
public:
  ColorMappedNrrdTextureObj(NrrdDataHandle &nrrd_handle,
                            ColorMapHandle &cmap, unsigned int label = 0);
  ~ColorMappedNrrdTextureObj();

  void                  set_colormap(ColorMapHandle &cmap);
  void                  set_clut_minmax(float min, float max);
  void                  set_dirty() { dirty_ = true; }
  void                  set_coords(const Point &, 
                                   const Vector &,
                                   const Vector &);  
  void			draw_quad(bool noalpha = false);
  void                  get_bounds(BBox&);
  void                  apply_colormap(int, int, int, int, int border=0);
  void                  set_label_color(const Color &c) { label_color_ = c; }
  void                  set_opacity(float op) { opacity_ = op; } 
  void                  set_threshold(double lower, double upper);

  const Point &         get_min() { return min_; }

private:
  friend class GeomColorMappedNrrdTextureObj;
  bool			bind(int x, int y);

  template <class T> 
  void			apply_colormap_to_data(unsigned char *dst,
                                               T *src,
                                               int row_width,
                                               int region_start,
                                               int region_width,
                                               int region_height,
                                               const unsigned char *rgba,
                                               int ncolors,
                                               float scale, float bias);
  template <class T> 
  void			apply_label_to_data(unsigned char *dst,
                                            T *src,
                                            int row_width,
                                            int region_start,
                                            int region_width,
                                            int region_height,
                                            unsigned int mask);

  template <class T> 
  void			apply_threshold_to_data(unsigned char *dst,
                                                T *src,
                                                int row_width,
                                                int region_start,
                                                int region_width,
                                                int region_height);

  typedef std::vector<std::pair<unsigned int, unsigned int> > divisions_t;

  ColorMapHandle        colormap_;
  unsigned char *       colormap_rgba_;
  Color                 label_color_;
  NrrdDataHandle	nrrd_handle_;
  bool  		dirty_;
  std::vector<bool>		texture_id_dirty_;
  std::vector<unsigned int>  texture_id_;
  divisions_t           xdiv_;
  divisions_t           ydiv_;
  float                 clut_min_;
  float                 clut_max_;
  unsigned char *       data_;
  unsigned int          label_;
  float                 opacity_;
  Point                 min_;
  Vector                xdir_;
  Vector                ydir_;
  double                lower_threshold_;
  double                upper_threshold_;
  bool                  thresholding_;
};


typedef LockingHandle<ColorMappedNrrdTextureObj> ColorMappedNrrdTextureObjHandle;


}

  
#endif
