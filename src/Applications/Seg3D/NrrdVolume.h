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
//    File   : NrrdVolume.h
//    Author : McKay Davis
//    Date   : Fri Oct 13 16:03:50 2006


#ifndef SEG3D_NrrdVolume
#define SEG3D_NrrdVolume

#include <vector>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Util/StringUtil.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/NrrdData.h>
#include <Core/Skinner/Variables.h>
#include <Core/Datatypes/NrrdToITK.h>
#include <Applications/Seg3D/VolumeSlice.h>
#include <Applications/Seg3D/VolumeOps.h>
#include <sci_defs/insight_defs.h>
#include <Core/Geom/IndexedGroup.h>
#include <Core/Thread/ThreadLock.h>

using std::vector;


namespace SCIRun {

class Painter;
class VolumeSlice;
class NrrdVolume;
class LayerButton;
class GeomSkinnerVarSwitch;
typedef LockingHandle<NrrdVolume> NrrdVolumeHandle;
typedef vector<NrrdVolumeHandle>	NrrdVolumes;

class NrrdVolume { 
public:
  // For LockingHandle<NrrdVolume> NrrdVolumeHandle typedef after this class
  ThreadLock          lock;
  unsigned int        ref_cnt;

  // Constructor
  NrrdVolume	      (Painter *painter, 
                       const string &name,
                       NrrdDataHandle &nrrdh,
                       const unsigned int label = 0);

  ~NrrdVolume();

  ITKDatatypeHandle   get_itk_image() 
  { return nrrd_to_itk_image(nrrd_handle_); }
    
  bool                write(const string &filename,
			    bool quantize = true, bool nrrd_on_error = true);
  bool                writeWithHeader(const string &filename, 
				      const string &hdr_filename,
				      const string &series_desc,
				      bool quantize = true, 
				      bool nrrd_on_error = true);
  static string       create_nrrd_header(const Nrrd *nrrd);

  void                set_label(unsigned int label);
  void                set_label_color(const Color &c);
  // Old way of doing this, retain compatability with 1.8.1 and prior.
  void                set_label_color_legacy(unsigned int index);
  static Color        get_next_label_color();

  void                set_opacity(const double opacity);
  double              get_opacity();
  
  void                set_nrrd(NrrdDataHandle &);
  NrrdDataHandle      get_nrrd_strip_subaxis();  // Makes a copy of the nrrd.
  void                set_dirty() { dirty_ = true; }

//   NrrdVolumeHandle    create_label_volume(unsigned int label=1, 
//                                           NrrdDataHandle nrrdh = 0);
  NrrdVolumeHandle    create_child_label_volume(unsigned int label=0);

  unsigned int        compute_label_mask(unsigned int label = 0);

  // Generates a VolumeSlice class if Plane intersects the volume,
  // Returns 0 if the Plane does not intersect the volume
  VolumeSliceHandle   get_volume_slice(const SLIVR::Plane &);

  void                reset_data_range();

  // Methods to transform between index and world space
private:
  void                build_index_to_world_matrix();
public:
  void                rebuild_transform();
  Point               index_to_world(const vector<int> &index);
  Point               index_to_point(const vector<double> &index);
  vector<int>         world_to_index(const Point &p);
  vector<int>         world_to_index(const SLIVR::Point &p);
  vector<double>      point_to_index(const Point &p);
  vector<double>      point_to_index(const SLIVR::Point &p);
  vector<double>      vector_to_index(const Vector &v);
  vector<double>      vector_to_index(const SLIVR::Vector &v);
  Vector              index_to_vector(const vector<double> &);
  bool                index_valid(const vector<int> &index);
  Point               center(int axis = -1, int slice = -1);
  Point               min(int axis = -1, int slice = -1);
  Point               max(int axis = -1, int slice = -1);

  // Voxel getter/setters
  template<class T>
  void                get_value(const vector<int> &index, T &value);
  template<class T>
  void                set_value(const vector<int> &index, T value);
  template<class T>
  void                and_value(const vector<int> &index, T value);
  template<class T>
  void                or_value(const vector<int> &index, T value);


  Vector              scale();
  double              scale(unsigned int axis);

  vector<int>         max_index();
  int                 max_index(unsigned int axis);

  bool                inside_p(const Point &p);
  ColorMapHandle      get_colormap();
  Color               get_label_color();
  MaterialHandle      get_label_material();
  GeomIndexedGroup*   get_geom_group();
  void                set_geom_switch(Skinner::Var<bool> v);

  NrrdDataHandle      extract_label_as_bit();

  GeomHandle          isosurface_label();

  void                change_type_from_float_to_bit(float val = 0);
  void                change_type_from_bit_to_float(float val = 1);
  void                clear();
  size_t              numbytes();
  int                 bit();

  void                set_slices_dirty();
  void                reset_clut();

  bool                visible();
  
  // This next function must be called inside an opengl context, because
  // it may call a geometry class desctructor that is deletes GL elements
  // that are bound to a gl context, and thus the context must be current.
  void                purge_unused_slices();

private:
  NrrdVolumeHandle    root();
public:
  int                 depth();

  // Do this when deleting a label volume so that it's internal
  // allocation bitplane is defeferenced.
  void                unparent();

  // True if this is the only label volume that is using this
  // particular nrrd.
  bool                only_one_label();

  Painter *           painter_;
private:
  NrrdVolumeHandle    parent_;
  NrrdVolumes         children_;
public:
  NrrdDataHandle      nrrd_handle_;
  string              name_;
  string              filename_;
  string              full_path_;  // Used to guess which cmap2 to load.

protected:
  double	      opacity_;
public:
  double              clut_min_;
  double              clut_max_;
  double              data_min_;
  double              data_max_;

  unsigned int        label_;

  DenseMatrix         transform_; 
  LayerButton *       button_;
  bool                tmp_visible_; // Don't use outside rebuild_layer_buttons
  bool                dirty_;

private:
  ColorMapHandle      colormap_;
  Color               label_color_;
  GeomHandle          geom_switch_;
  VolumeSlices_t      all_slices_;
  vector<int>         stub_axes_;
};


template<class T>
void
NrrdVolume::get_value(const vector<int> &index, T &value) 
{
  ASSERT(index_valid(index));
  nrrd_handle_->lock.lock();
  VolumeOps::nrrd_get_value(nrrd_handle_->nrrd_, index, value);
  nrrd_handle_->lock.unlock();
}


template <class T>
void
NrrdVolume::set_value(const vector<int> &index, T value) {
  ASSERT(index_valid(index));
  VolumeOps::nrrd_set_value_unmasked(nrrd_handle_->nrrd_, index, value);
}


template <class T>
void
NrrdVolume::and_value(const vector<int> &index, T value) {
  ASSERT(index_valid(index));
  VolumeOps::nrrd_and_value_unmasked(nrrd_handle_->nrrd_, index, value);
}


template <class T>
void
NrrdVolume::or_value(const vector<int> &index, T value) {
  ASSERT(index_valid(index));
  VolumeOps::nrrd_or_value_unmasked(nrrd_handle_->nrrd_, index, value);
}


}


#endif
