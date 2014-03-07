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
//    File   : VolumeSlice.h
//    Author : McKay Davis
//    Date   : Fri Oct 13 15:50:22 2006

#ifndef Seg3D_VolumeSlice_h
#define Seg3D_VolumeSlice_h

#include <slivr/Plane.h>
#include <Core/Geom/GeomObj.h>
#include <Core/Geom/ColorMappedNrrdTextureObj.h>
#include <Core/Geom/NrrdBitmaskOutline.h>
#include <map>
#include <vector>


namespace SCIRun {

class Painter;
class NrrdVolume;
class SliceWindow;

class VolumeSlice {
public:

  VolumeSlice(NrrdVolume *, 
              const SLIVR::Plane &,
              NrrdDataHandle nrrd ,
              unsigned int label);
  
  // For LockingHandle
  Mutex                                 lock;
  int                                   ref_cnt;

  void                                  draw();
  void                                  set_tex_dirty();
  void                                  set_region_dirty(int x1, int y1,
                                                         int x2, int y2,
                                                         int border);

  unsigned int                          get_axis();
  const SLIVR::Plane &                  get_plane() { return plane_; }
  GeomHandle                      get_geom_texture() { return geom_texture_; }

  NrrdVolume *                          volume_;
  NrrdDataHandle                        nrrd_handle_;
  
private:
  void                                  bind();
  void                                  extract_nrrd_slice_from_volume();

  NrrdBitmaskOutlineHandle              outline_;
  ColorMappedNrrdTextureObjHandle       texture_;
  GeomHandle                            geom_texture_;

  SLIVR::Plane                          plane_;
  unsigned int                          label_;
  bool                                  tex_dirty_;

public:  
  Point                                 pos_;
  Vector                                xdir_;
  Vector                                ydir_;
};


typedef LockingHandle<VolumeSlice> VolumeSliceHandle;
typedef std::vector<VolumeSliceHandle> VolumeSlices_t;
typedef map<string, VolumeSlices_t> VolumeSliceGroups_t;    
typedef map<int, VolumeSlices_t> VolumeSlicesAlongAxis_t;
typedef map<int, VolumeSlicesAlongAxis_t> VolumeSliceCache_t;


}

#endif
  
