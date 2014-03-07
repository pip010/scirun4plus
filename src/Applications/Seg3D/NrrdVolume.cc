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
//    File   : NrrdVolume.cc
//    Author : McKay Davis
//    Date   : Fri Oct 13 15:06:57 2006

#include <Applications/Seg3D/NrrdVolume.h>
#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/VolumeOps.h>
#include <sci_comp_warn_fixes.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <math.h>
#include <map>
#include <typeinfo>
#include <iostream>
#include <sci_gl.h>
#include <algorithm>
#include <vector>
#include <sstream>

#include <Core/Datatypes/Field.h>
#include <Core/Exceptions/GuiException.h>
#include <Core/Geom/OpenGLViewport.h>
#include <Core/Geom/GeomQuads.h>
#include <Core/Math/MiscMath.h>
#include <Core/Math/MinMax.h>
#include <Core/Util/Environment.h>
#include <Core/Geom/FontManager.h>
#include <Core/Skinner/Variables.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/SceneGraphEvent.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Skinner/GeomSkinnerVarSwitch.h>
#include <itkImageSeriesWriter.h>
#include <itkImageSeriesReader.h>
#include <itkNumericSeriesFileNames.h>
#include <itkGDCMImageIO.h>
#include <itkMetaDataObject.h>
#include <gdcmUtil.h>
#include <Core/Skinner/Histogram.h>
#include <Applications/Seg3D/Seg3DwxGuiUtils.h>

namespace SCIRun {



static SLIVR::ColorMap *
create_greyscale_fixed()
{
  float rgba[256*4];
  for (int c = 0; c < 256*4; ++c) {
    rgba[c] = (c % 4 == 3) ? 1.0f : (c/4) / 255.0f;
  }

  return new SLIVR::ColorMap(rgba);
}


static SLIVR::ColorMap *
create_label_colormap(const Color &c, const double opacity, int position)
{
  float rgba[256*4];

  // Clear it all
  for (int i = 0; i < 256*4; ++i)
  {
    rgba[i] = 0.0;
  }

  rgba[position*4+0] = c.r();
  rgba[position*4+1] = c.g();
  rgba[position*4+2] = c.b();
  rgba[position*4+3] = opacity;
  
  return new SLIVR::ColorMap(rgba);
}

static double static_hue = 280.0;
static double static_saturation = 1.05;


Color
NrrdVolume::get_next_label_color()
{
  double hue = static_hue + 360.0 / 6.4;
  if (hue > 360.0) hue -= 360.0;
  double sat = static_saturation - 0.05;
  if (sat < 0.50) sat = 1.0;
  return HSVColor(hue, sat, 1.0);
}


NrrdVolume::NrrdVolume(Painter *painter,
                       const string &name,
                       NrrdDataHandle &nrrd,
                       unsigned int label) :
  lock(("NrrdVolume "+name).c_str()),
  ref_cnt(0),
  painter_(painter),
  parent_(0),
  children_(0),
  nrrd_handle_(0),
  name_(name),
  filename_(name),
  full_path_(""),
  opacity_(1.0),
  clut_min_(0.0),
  clut_max_(1.0),
  data_min_(0),
  data_max_(1.0),
  label_(label),
  transform_(),
  tmp_visible_(true),
  label_color_(0.0, 0.0, 0.0),
  stub_axes_()
{
  if (label_)
  {
    static_hue += 360.0 / 6.4;
    if (static_hue > 360.0) static_hue -= 360.0;
    static_saturation -= 0.05;
    if (static_saturation < 0.50) static_saturation = 1.0;

    set_label_color( HSVColor(static_hue, static_saturation, 1.0) );
  }
  else
  {
    SLIVR::ColorMap *cmtmp = create_greyscale_fixed();
    colormap_ = new SCIRun::ColorMap(*cmtmp);
    delete cmtmp;
  }

  set_nrrd(nrrd);
}


NrrdVolume::~NrrdVolume()
{
  nrrd_handle_ = 0;
}

void
NrrdVolume::set_label(unsigned int label)
{
  label_ = label;
}


void
NrrdVolume::set_label_color(const Color &c)
{
  label_color_ = c;
  if (label_) {
    SLIVR::ColorMap *cmtmp =
      create_label_colormap(label_color_, opacity_, bit()+1);
    colormap_ = new SCIRun::ColorMap(*cmtmp);
    delete cmtmp;
  }
}


void
NrrdVolume::set_label_color_legacy(unsigned int index)
{
  // This ugly formula converts index into a hue.
  // It isn't particularly meaningful except that it's what was in
  // the code originally to get a reasonable spread of hue values.
  const double hue = 90.0 + (360.0 / 255.0 * index + 180.0) * (bit()+1);
  const Color c = HSVColor(hue, 1.0, 1.0);
  set_label_color(c);
}


double
NrrdVolume::get_opacity()
{
  return opacity_;
}


void
NrrdVolume::set_opacity(double opacity)
{
  opacity_ = opacity;
}


void
NrrdVolume::set_nrrd(NrrdDataHandle &nrrd_handle)
{
  nrrd_handle_ = nrrd_handle;
  nrrd_handle_->lock.lock();
  Nrrd *n = nrrd_handle_->nrrd_;

  stub_axes_.clear();
  if (n->axis[0].size > 4) {
    nrrdAxesInsert(n, n, 0);
    n->axis[0].min = 0.0;
    n->axis[0].max = 1.0;
    n->axis[0].spacing = 1.0;
    stub_axes_.push_back(0);
  }

  if (n->dim == 3) {
    nrrdAxesInsert(n, n, 3);
    n->axis[3].min = 0.0;
    n->axis[3].max = 1.0;
    n->axis[3].spacing = 1.0;
    stub_axes_.push_back(3);
  }


  for (unsigned int a = 0; a < n->dim; ++a) {

    if (n->axis[a].center == nrrdCenterUnknown)
      n->axis[a].center = nrrdCenterNode;

    if (n->axis[a].kind == nrrdKindUnknown)
      n->axis[a].kind = nrrdKindDomain;

    if (n->axis[a].min > n->axis[a].max)
      std::swap(n->axis[a].min,n->axis[a].max);
    if (n->axis[a].spacing < 0.0)
      n->axis[a].spacing *= -1.0;
  }

  nrrd_handle_->lock.unlock();
  build_index_to_world_matrix();
  dirty_ = true;
  reset_data_range();
  // Reset the 3D geometry bounding box.
  geom_switch_ = 0;
}


void
NrrdVolume::reset_data_range()
{
  NrrdRange range;
  nrrdRangeSet(&range, nrrd_handle_->nrrd_, 0);
  if (data_min_ != range.min || data_max_ != range.max) {
    data_min_ = range.min;
    data_max_ = range.max;
  }
  reset_clut();
  opacity_ = 1.0;
}


NrrdDataHandle
NrrdVolume::get_nrrd_strip_subaxis()
{
  NrrdDataHandle nrrd_handle = nrrd_handle_;
  nrrd_handle.detach();
  NrrdDataHandle nrrd2_handle = new NrrdData();

  //   nrrdBasicInfoCopy(nrrd->nrrd, nrrd_handle_->nrrd_,0);
  //   nrrdAxisInfoCopy(nrrd->nrrd, nrrd_handle_->nrrd_, 0,0);
  //   nrrd->nrrd->data = nrrd_handle_->nrrd_->data;

  for (int s = stub_axes_.size()-1; s >= 0 ; --s) {
    nrrdAxesDelete(nrrd2_handle->nrrd_, nrrd_handle->nrrd_, stub_axes_[s]);
    nrrd_handle = nrrd2_handle;
  }

  while (nrrd_handle->nrrd_->dim && nrrd_handle->nrrd_->axis[0].size == 1) {
    nrrdAxesDelete(nrrd2_handle->nrrd_, nrrd_handle->nrrd_, 0);
    nrrd_handle = nrrd2_handle;
  }

  nrrdKeyValueCopy(nrrd_handle->nrrd_, nrrd_handle_->nrrd_);

  //  unsigned long ptr = (unsigned long)(&painter_);
  //  nrrdKeyValueAdd(nrrd_handle->nrrd_,
  //                  "progress_ptr", to_string(ptr).c_str());

  return nrrd_handle;
}



Point
NrrdVolume::center(int axis, int slice)
{
  nrrd_handle_->lock.lock();
  Nrrd *n = nrrd_handle_->nrrd_;
  vector<int> index(n->dim,0);
  for (unsigned int a = 0; a < index.size(); ++a)
    index[a] = n->axis[a].size/2;
  if (axis >= 0 && axis < int(index.size()))
    index[axis] = Clamp(slice, 0, n->axis[axis].size-1);
  ASSERT(index_valid(index));
  nrrd_handle_->lock.unlock();
  return index_to_world(index);
}


Point
NrrdVolume::min(int axis, int slice)
{
  nrrd_handle_->lock.lock();
  Nrrd *n = nrrd_handle_->nrrd_;
  vector<int> index(n->dim,0);
  if (axis >= 0 && axis < int(index.size()))
    index[axis] = Clamp(slice, 0, n->axis[axis].size-1);
  ASSERT(index_valid(index));
  nrrd_handle_->lock.unlock();
  return index_to_world(index);
}


Point
NrrdVolume::max(int axis, int slice)
{
  vector<int> index = max_index();
  nrrd_handle_->lock.lock();
  Nrrd *n = nrrd_handle_->nrrd_;
  if (axis >= 0 && axis < int(index.size()))
    index[axis] = Clamp(slice, 0, n->axis[axis].size-1);
  ASSERT(index_valid(index));
  nrrd_handle_->lock.unlock();
  return index_to_world(index);
}


Vector
NrrdVolume::scale()
{
  vector<int> index_zero(nrrd_handle_->nrrd_->dim,0);
  vector<int> index_one(nrrd_handle_->nrrd_->dim,1);
  return index_to_world(index_one) - index_to_world(index_zero);
}


double
NrrdVolume::scale(unsigned int axis)
{
  ASSERT(axis >= 0 && (unsigned int) axis < nrrd_handle_->nrrd_->dim);
  return scale()[axis];
}



vector<int>
NrrdVolume::max_index()
{
  vector<int> max_index(nrrd_handle_->nrrd_->dim,0);
  for (unsigned int a = 0; a < nrrd_handle_->nrrd_->dim; ++a)
    max_index[a] = nrrd_handle_->nrrd_->axis[a].size;
  return max_index;
}


int
NrrdVolume::max_index(unsigned int axis)
{
  ASSERT(axis >= 0 && (unsigned int) axis < nrrd_handle_->nrrd_->dim);
  return max_index()[axis];
}


bool
NrrdVolume::inside_p(const Point &p)
{
  return index_valid(world_to_index(p));
}


Point
NrrdVolume::index_to_world(const vector<int> &index)
{
  unsigned int dim = index.size()+1;
  ColumnMatrix index_matrix(dim);
  ColumnMatrix world_coords(dim);
  for (unsigned int i = 0; i < dim-1; ++i)
    index_matrix[i] = index[i];
  index_matrix[dim-1] = 1.0;
  DenseMatrix transform = transform_;
  transform.mult(index_matrix, world_coords);
  Point return_val;
  for (int i = 1; i < 4; ++i)
    return_val(i-1) = world_coords[i];
  return return_val;
}


Point
NrrdVolume::index_to_point(const vector<double> &index)
{
  unsigned int dim = index.size()+1;
  ColumnMatrix index_matrix(dim);
  ColumnMatrix world_coords(dim);
  for (unsigned int i = 0; i < dim-1; ++i)
    index_matrix[i] = index[i];
  index_matrix[dim-1] = 1.0;
  DenseMatrix transform = transform_;
  transform.mult(index_matrix, world_coords);
  Point return_val;
  for (int i = 1; i < 4; ++i)
    return_val(i-1) = world_coords[i];
  return return_val;
}


vector<int>
NrrdVolume::world_to_index(const Point &p)
{
  DenseMatrix transform = transform_;
  ColumnMatrix index_matrix(transform.ncols());
  ColumnMatrix world_coords(transform.nrows());
  for (int i = 0; i < transform.nrows(); ++i)
    if (i > 0 && i < 4)
      world_coords[i] = p(i-1)-transform.get(i,transform.ncols()-1);
    else
      world_coords[i] = 0.0;;
  transform.solve(world_coords, index_matrix, 1);
  vector<int> return_val(index_matrix.nrows()-1);
  for (unsigned int i = 0; i < return_val.size(); ++i) {
    return_val[i] = Floor(index_matrix[i]);
  }
  return return_val;
}


vector<int>
NrrdVolume::world_to_index(const SLIVR::Point &p)
{
  DenseMatrix transform = transform_;
  ColumnMatrix index_matrix(transform.ncols());
  ColumnMatrix world_coords(transform.nrows());
  for (int i = 0; i < transform.nrows(); ++i)
    if (i > 0 && i < 4)
      world_coords[i] = p(i-1)-transform.get(i,transform.ncols()-1);
    else
      world_coords[i] = 0.0;;
  transform.solve(world_coords, index_matrix, 1);
  vector<int> return_val(index_matrix.nrows()-1);
  for (unsigned int i = 0; i < return_val.size(); ++i) {
    return_val[i] = Floor(index_matrix[i]);
  }
  return return_val;
}


vector<double>
NrrdVolume::point_to_index(const Point &p)
{
  DenseMatrix transform = transform_;
  ColumnMatrix index_matrix(transform.ncols());
  ColumnMatrix world_coords(transform.nrows());
  for (int i = 0; i < transform.nrows(); ++i)
    if (i > 0 && i < 4)
      world_coords[i] = p(i-1)-transform.get(i,transform.ncols()-1);
    else
      world_coords[i] = 0.0;;
  transform.solve(world_coords, index_matrix, 1);
  vector<double> return_val(index_matrix.nrows()-1);
  for (unsigned int i = 0; i < return_val.size(); ++i) {
    return_val[i] = index_matrix[i];
  }
  return return_val;
}


vector<double>
NrrdVolume::point_to_index(const SLIVR::Point &p)
{
  DenseMatrix transform = transform_;
  ColumnMatrix index_matrix(transform.ncols());
  ColumnMatrix world_coords(transform.nrows());
  for (int i = 0; i < transform.nrows(); ++i)
    if (i > 0 && i < 4)
      world_coords[i] = p(i-1)-transform.get(i,transform.ncols()-1);
    else
      world_coords[i] = 0.0;;
  transform.solve(world_coords, index_matrix, 1);
  vector<double> return_val(index_matrix.nrows()-1);
  for (unsigned int i = 0; i < return_val.size(); ++i) {
    return_val[i] = index_matrix[i];
  }
  return return_val;
}


vector<double>
NrrdVolume::vector_to_index(const Vector &v)
{
  Point zero(0,0,0);
  vector<double> zero_idx = point_to_index(zero);
  vector<double> idx = point_to_index(v.asPoint());
  for (unsigned int i = 0; i < zero_idx.size(); ++i)
    idx[i] = idx[i] - zero_idx[i];
  return idx;
}


vector<double>
NrrdVolume::vector_to_index(const SLIVR::Vector &v)
{
  Point zero(0,0,0);
  vector<double> zero_idx = point_to_index(zero);
  vector<double> idx = point_to_index(v.asPoint());
  for (unsigned int i = 0; i < zero_idx.size(); ++i)
    idx[i] = idx[i] - zero_idx[i];
  return idx;
}

Vector
NrrdVolume::index_to_vector(const vector<double> &index)
{
  vector<double> zero_index(index.size(),0.0);
  return index_to_point(index) - index_to_point(zero_index);
}


void
NrrdVolume::build_index_to_world_matrix()
{
  nrrd_handle_->lock.lock();
  Nrrd *nrrd = nrrd_handle_->nrrd_;
  int dim = nrrd->dim+1;
  DenseMatrix matrix(dim, dim);
  matrix.zero();
  for (int i = 0; i < dim-1; ++i) {
    if (airExists(nrrd->axis[i].spacing)) {
      if (nrrd->axis[i].spacing == 0.0) {
        nrrd->axis[i].spacing = 1.0;
      }
      matrix.put(i,i,nrrd->axis[i].spacing);
    } else if (airExists(nrrd->axis[i].min) && airExists(nrrd->axis[i].max)) {
      if (nrrd->axis[i].min == nrrd->axis[i].max) {
        nrrd->axis[i].spacing = 1.0;
        matrix.put(i,i,1.0);
      } else {
        matrix.put(i,i,((nrrd->axis[i].max-nrrd->axis[i].min)/
                        nrrd->axis[i].size));
      }
    } else {
      matrix.put(i,i, 1.0);
    }

    if (airExists(nrrd->axis[i].min))
      matrix.put(i, nrrd->dim, nrrd->axis[i].min);
  }

  if (nrrd->axis[0].size != 1) {
    matrix.put(2,nrrd->dim, nrrd->axis[2].min+nrrd->axis[2].size*matrix.get(2,2));
    matrix.put(2,2,-matrix.get(2,2));
  }

  matrix.put(dim-1, dim-1, 1.0);

  transform_ = matrix;
  nrrd_handle_->lock.unlock();
}


void
NrrdVolume::rebuild_transform()
{
  build_index_to_world_matrix();
  geom_switch_ = 0;  // Reset the geometry position as well.
}


bool
NrrdVolume::index_valid(const vector<int> &index)
{
  unsigned int dim = nrrd_handle_->nrrd_->dim;
  if (index.size() != dim) return false;
  for (unsigned int a = 0; a < dim; ++a)
    if (index[a] < 0 ||
	(unsigned int) index[a] >= nrrd_handle_->nrrd_->axis[a].size) {
      return false;
    }
  return true;
}


unsigned int
NrrdVolume::compute_label_mask(unsigned int label)
{
  ASSERT((label & label_) == 0);
  label |= label_;
  for (unsigned int i = 0; i < children_.size(); ++i) {
    label |= children_[i]->compute_label_mask(label);
  }
  return label;
}


NrrdVolumeHandle
NrrdVolume::create_child_label_volume(unsigned int label)
{
  if (!label) {
    const unsigned char max_bit = sizeof(label_type)*8;
    unsigned char bit = 0;
    unsigned int used_labels = root()->compute_label_mask();
    while (bit < max_bit && (used_labels & (1 << bit))) ++bit;
    if (bit == max_bit) {
      // Cannot create child label volume!
      return 0;
    }
    label = 1 << bit;
  }

  NrrdVolumeHandle vol =
    new NrrdVolume(painter_,
                   root()->name_ + " " + to_string(label),
                   nrrd_handle_, label);
  vol->parent_ = this;
  children_.push_back(vol);
  return vol;
}


template <class T>
static bool
write_itk_image(ITKDatatypeHandle &itk_image_h, const string &filename)
{
  // create a new writer
  typedef itk::Image < T, 3 > ImageType;
  typedef itk::ImageFileWriter< ImageType > FileWriterType;
  typename FileWriterType::Pointer writer = FileWriterType::New();
  ImageType *img =
    dynamic_cast<ImageType *>(itk_image_h->data_.GetPointer());
  ASSERT(img);

  // set writer
  writer->SetFileName( filename.c_str() );
  writer->SetInput(img);
  writer->SetUseCompression(true);

  writer->Update();

  return true;
}


static pair<string, string>
split_extension(const string &filename)
{
  string base;
  vector<string> fileext = split_string(filename, '.');
  for (int i = 0; i < (int)fileext.size()-1; i++)
  {
    if (i > 0) base = base + ".";
    base = base + fileext[i];
  }
  string ext;
  if (fileext.size() > 1)
  {
    ext = fileext.back();
  }
  return pair<string, string>(base, ext);
}

// img_center, img_width are used for writing DICOM headers, 
// img_center becomes the window center
// img_width becomes the window width
template <class T>
static bool
write_itk_series(ITKDatatypeHandle &itk_image_h, const string &filename,
		 const vector<string> &hdr_series, const string &series_desc,
		 bool dicom, float img_center, float img_width)
{
  // Create a new writer
  typedef itk::Image < T, 3 > ImageType;
  typedef itk::Image < T, 2 > Image2DType;
  typedef itk::ImageSeriesWriter< ImageType, Image2DType > FileWriterType;
  typedef itk::ImageSeriesReader< itk::Image<T, 3> > FileReaderType;

  typename FileWriterType::Pointer writer = FileWriterType::New();
  typename FileReaderType::Pointer reader;
  itk::GDCMImageIO::Pointer gdcmio;
  typename FileReaderType::DictionaryArrayType dict_array;
  ImageType *img =
    dynamic_cast<ImageType *>(itk_image_h->data_.GetPointer());
  ASSERT(img);

  if (dicom)
  { 
    gdcmio = itk::GDCMImageIO::New();
    if (!hdr_series.empty())
    {
      reader = FileReaderType::New();
      reader->SetFileNames(hdr_series); 
      reader->SetImageIO( gdcmio );
      try {
	reader->Update();
      } catch (itk::ExceptionObject & err ) {
	cerr << "NrrdVolume::write_itk_series - ITK ExceptionObject caught!" << std::endl;
	cerr << err.GetDescription() << std::endl;
	return false;
      }
      dict_array = *(reader->GetMetaDataDictionaryArray());
      typename FileReaderType::DictionaryArrayType::iterator dict_iter;
      // now keep the study UID, but generate new series and SOP UIDs.  ITK only lets us
      // keep all 4 UIDs or discard all 4 UIDs, hence this bit of hackery.
      // series UID should be the same across the series, SOP UID should be different for
      // each slice
      std::string uid_prefix = gdcmio->GetUIDPrefix();
      std::string seriesUID = gdcm::Util::CreateUniqueUID(uid_prefix);
      // DICOM standard only allows 64 characters for series description
      std::string short_desc(series_desc, 0, 64); 
      std::ostringstream ss;
      ss << img_center;
      std::string img_center_str(ss.str());
      ss.str("");
      ss << img_width;
      std::string img_width_str(ss.str());
      ss.str("");
      for (dict_iter = dict_array.begin(); dict_iter != dict_array.end(); dict_iter++)
      {
	// this works for GDCM_MAJOR_VERSION < 2, see itkGDCMIO.cxx for how to 
	// generate UID for GDCM_MAJOR_VERSION >= 2
	std::string sopUID = gdcm::Util::CreateUniqueUID(uid_prefix);
	itk::EncapsulateMetaData<std::string>(**dict_iter, "0020|000e", seriesUID); // [Series Instance UID]
	itk::EncapsulateMetaData<std::string>(**dict_iter, "0008|0018", sopUID); // [SOP Instance UID]
	itk::EncapsulateMetaData<std::string>(**dict_iter, "0002|0003", sopUID); // [Media Stored SOP Instance UID]
	itk::EncapsulateMetaData<std::string>(**dict_iter, "0028|1050", img_center_str); // [ Window Center ]
	itk::EncapsulateMetaData<std::string>(**dict_iter, "0028|1051", img_width_str); // [ Window Width ]
	if (short_desc != "")
	{
	  itk::EncapsulateMetaData<std::string>(**dict_iter, "0008|103e", short_desc); // [Series Description]
	}
	// unset the following values so that itkGDCMImageIO will fill them in correctly
	itk::EncapsulateMetaData<std::string>(**dict_iter, "0028|0100", ""); // [Bits Allocated]
	itk::EncapsulateMetaData<std::string>(**dict_iter, "0028|0101", ""); // [Bits Stored]
	itk::EncapsulateMetaData<std::string>(**dict_iter, "0028|0102", ""); // [High Bit]
	itk::EncapsulateMetaData<std::string>(**dict_iter, "0028|0103", ""); // [Pixel Representation]
      }
      // now tell the Dicom writer to use the provided UIDs
      gdcmio->KeepOriginalUIDOn();
      writer->SetMetaDataDictionaryArray( &dict_array );
    }
    writer->SetImageIO( gdcmio ); 
  }

  // Compute the slice range.
  typename ImageType::RegionType region = img->GetLargestPossibleRegion();
  typename ImageType::IndexType  start = region.GetIndex();
  typename ImageType::SizeType   size  = region.GetSize();
  const size_t firstSlice = start[2];
  const size_t lastSlice = start[2] + size[2] - 1;

  // Create the name generator
  typedef itk::NumericSeriesFileNames NameGeneratorType;
  typename NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();

  // Split up the filename, add in the numbers.
  pair<string, string> file_base_ext = split_extension(filename);
  // TODO: Need to fix this for VC++ lnk 2019 error.
  int numerals = (int)Ceil(log10((float)(lastSlice - firstSlice + 1)));
  const string format = file_base_ext.first +
    "%0" + to_string(numerals) + "d" +
    "." + file_base_ext.second;
  nameGenerator->SetSeriesFormat(format.c_str());

  nameGenerator->SetStartIndex(firstSlice);
  nameGenerator->SetEndIndex(lastSlice);
  nameGenerator->SetIncrementIndex(1);

  // Set writer
  writer->SetFileNames( nameGenerator->GetFileNames() );
  writer->SetInput(img);

  writer->Update();

  return true;
}


static bool
vff_writer(NrrdDataHandle nin, const string &filename)
{
  NrrdDataHandle nout;

  if (nin->nrrd_->type == nrrdTypeUChar)
  {
    nout = nin;
  }
  else
  {
    nout = new NrrdData();
    if (nrrdConvert(nout->nrrd_, nin->nrrd_, nrrdTypeShort))
    {
      char *err = biffGetDone(NRRD);
      string errstr(err);
      free(err);
      throw "VFF Writer converting to short: " + errstr;
    }
  }

  ofstream vffFileStream(filename.c_str(), ios::binary);

  if (! vffFileStream)
  {
    throw "Could not open file " + filename;
    return false;
  }

  vffFileStream << "ncaa\n";
  vffFileStream << "rank=3;\n";
  vffFileStream << "type=raster;\n";
  vffFileStream << "size=" <<
    nin->nrrd_->axis[1].size << " " <<
    nin->nrrd_->axis[2].size << " " <<
    nin->nrrd_->axis[3].size << ";\n";
  vffFileStream << "bits=" <<
    (nin->nrrd_->type == nrrdTypeUChar?"8":"16") << ";\n";
  vffFileStream << "format=slice;\n";
  vffFileStream << "spacing=" <<
    nin->nrrd_->axis[1].spacing << " " <<
    nin->nrrd_->axis[2].spacing << " " <<
    nin->nrrd_->axis[3].spacing << ";\n";
  vffFileStream << "origin=" <<
    nin->nrrd_->axis[1].min << " " <<
    nin->nrrd_->axis[2].min << " " <<
    nin->nrrd_->axis[3].min << ";\n"; // TODO: check - failed
  vffFileStream << "elementsize=" <<
    nin->nrrd_->axis[1].size << ";\n"; // TODO: check - failed
  vffFileStream << (char)12 << '\n';

  if (nout->nrrd_->type != nrrdTypeUChar && AIR_ENDIAN != airEndianBig)
  {
    nrrdSwapEndian(nout->nrrd_);
  }

  vffFileStream.write((const char *)(nout->nrrd_->data),
                      VolumeOps::nrrd_data_size(nout));

  vffFileStream.close();
  
  return true;
}


bool
NrrdVolume::write(const string &ofname,
		  bool quantize,
		  bool nrrd_on_error)
{
  const string fname = substituteTilde(ofname);

  const string lfname = string_tolower(fname);

  if (ends_with(lfname, ".vff"))
  {
    return vff_writer(nrrd_handle_, fname);
  }

  else if (ends_with(lfname, ".nrrd"))
  {
    // Restrict TEEM access: it is not thread safe
    NrrdGuard guard;

    NrrdIoState *nio = nrrdIoStateNew();
    nio->encoding = nrrdEncodingArray[1]; // raw
    nio->format = nrrdFormatArray[1];     // nrrd
    nio->endian = AIR_ENDIAN;             // endian 

    if (nrrd_handle_->nrrd_->dim == 4 || nrrd_handle_->nrrd_->axis[0].size == 1)
    {
      Nrrd* tmp = nrrdNew();
      nrrdBasicInfoCopy(tmp,nrrd_handle_->nrrd_,0);
      tmp->data = nrrd_handle_->nrrd_->data;
      
      for (int i=0; i< 3; i++)
      {
        tmp->axis[i].size = nrrd_handle_->nrrd_->axis[i+1].size;
        tmp->axis[i].center = nrrd_handle_->nrrd_->axis[i+1].center;
        tmp->axis[i].kind = nrrd_handle_->nrrd_->axis[i+1].kind;
        tmp->axis[i].min = nrrd_handle_->nrrd_->axis[i+1].min;
        tmp->axis[i].max = nrrd_handle_->nrrd_->axis[i+1].max;
        tmp->axis[i].spacing = nrrd_handle_->nrrd_->axis[i+1].spacing;
        tmp->axis[i].thickness = nrrd_handle_->nrrd_->axis[i+1].thickness;
        for (int d=0; d<tmp->spaceDim; d++) tmp->axis[i].spaceDirection[d] = nrrd_handle_->nrrd_->axis[i+1].spaceDirection[d];
      }      
      tmp->dim = 3;
      
      if (nrrdSave(fname.c_str(), tmp, nio)) 
      {
        char *err = biffGet(NRRD);      
        cerr << "Error writing nrrd " << fname << ": "<< err << endl;
        free(err);
        biffDone(NRRD);
        nrrdNix(tmp);
        return false;
      }      
      nrrdNix(tmp);
    }
    else
    {
      if (nrrdSave(fname.c_str(), nrrd_handle_->nrrd_, nio)) 
      {
        char *err = biffGet(NRRD);      
        cerr << "Error writing nrrd " << fname << ": "<< err << endl;
        free(err);
        biffDone(NRRD);
        return false;
      }
    }
    
    return true;
  }
  // First try to save as a Matlab matrix
  else if (ends_with(lfname, ".mat") &&
      // TODO: This appears to make a copy of the nrrd and then remove
      // an axis from it.  Double check to see if the copy is really
      // needed for the matlab writer.
      painter_->MatlabNrrd_writer(nrrd_handle_, fname.c_str())) {
    return true;
  }

  else if (ends_with(lfname, ".tiff") ||
	   ends_with(lfname, ".tif") ||
	   ends_with(lfname, ".bmp") ||
	   ends_with(lfname, ".dcm") ||
	   ends_with(lfname, ".dicom") ||
	   ends_with(lfname, ".ima") ||
	   ends_with(lfname, ".png") ||
	   ends_with(lfname, ".jpeg") ||
	   ends_with(lfname, ".jpg"))
  {
    try {
      const bool is_dicom =
        ends_with(lfname, ".dcm") || ends_with(lfname, ".dicom") || ends_with(lfname, ".ima");
      const vector<string> header;

      // TODO: Need a general output degrader (see vff_writer) and UI.
      if (nrrd_handle_->nrrd_->type == nrrdTypeUChar)
      {
        ITKDatatypeHandle img = nrrd_to_itk_image(nrrd_handle_);
        return write_itk_series<unsigned char>(img, fname, header, "", is_dicom, 
					       (data_min_ + data_max_) / 2.0, data_max_ - data_min_);
      }
      else if (ends_with(lfname, ".tiff") ||
               ends_with(lfname, ".tif") ||
               ends_with(lfname, ".dcm") ||
               ends_with(lfname, ".dicom") ||
	       ends_with(lfname, ".ima") || 
               ends_with(lfname, ".png"))
      {
        // These formats support 16 bit output.
        // TODO: dicom should, but the results were broken.
        // NOTE: Maybe dicom requires signed shorts like vff?
	NrrdDataHandle tmp;

	if( quantize )
	  tmp = Skinner::Histogram::UnuQuantize(nrrd_handle_, 0.0, 0.0, 16);
	else
	  tmp = Skinner::Histogram::UnuQuantize(nrrd_handle_, 0.0, 65535.0, 16);

	NrrdRange range;
	nrrdRangeSet(&range, nrrd_handle_->nrrd_, 0);
        ITKDatatypeHandle img = nrrd_to_itk_image(tmp);
        return write_itk_series<unsigned short>(img, fname, header, "", is_dicom,
						(range.min + range.max) / 2.0, range.max - range.min);
      }
      else
      {
        NrrdDataHandle tmp;

	if( quantize )
	  tmp = Skinner::Histogram::UnuQuantize(nrrd_handle_, 0.0, 0.0, 8);
	else
	  tmp = Skinner::Histogram::UnuQuantize(nrrd_handle_, 0.0, 255, 8);

	NrrdRange range;
	nrrdRangeSet(&range, nrrd_handle_->nrrd_, 0);
        ITKDatatypeHandle img = nrrd_to_itk_image(tmp);
        return write_itk_series<unsigned char>(img, fname, header, "", is_dicom, 
					       (range.min + range.max) / 2.0, range.max - range.min);	
      }
    }
    catch  ( itk::ExceptionObject & err )
    {
      cerr << "NrrdVolume::write ITK::ExceptionObject caught" << endl;
      cerr << err.GetDescription() << endl;
      painter_->set_status_safe(err.GetDescription());
      return false;
    }
    return true;
  }

  else
  {
    try {
      ITKDatatypeHandle img = get_itk_image();
      switch (nrrd_handle_->nrrd_->type) {
	
      case LabelNrrdType:
	return write_itk_image<label_type>(img, fname);
	break;

      case nrrdTypeFloat:
	return write_itk_image<float>(img, fname);
	break;

      default:
	throw "Unsupported nrrd format write.";
      }
    }
    catch  ( itk::ExceptionObject & err )
    {
      const string msg = err.GetDescription();
      if (!strncmp(msg.c_str(), " Could not create IO object for file", 36))
      {
	const string nfname = changeExtension(ofname, "nrrd");
	
	// Ask the user if they want to use .nrrd here.
	wxMessageDialog dialog(Painter::global_seg3dframe_pointer_,
			       std2wx("Unsupported output format.\nTry " 
				      + nfname + "?"),
			       _T("Unsupported format"),
			       wxYES_NO | wxICON_ERROR);
	if (dialog.ShowModal() == wxID_YES)
	{
	  return write(nfname, false);
	}
	return false;
      }
      cerr << "NrrdVolume::write ITK::ExceptionObject caught" << endl;
      cerr << msg << endl;
      painter_->set_status_safe(err.GetDescription());
      return false;
    }
  }

  return false;
}
bool
NrrdVolume::writeWithHeader(const string &ofname,
			    const string &hdr_filename,
			    const string &series_desc,
			    bool quantize,
			    bool nrrd_on_error)
{
  const string fname = substituteTilde(ofname);

  const string lfname = string_tolower(fname);

  const vector<string> hdr_series = painter_->get_filename_series(hdr_filename);

  if (ends_with(lfname, ".dcm") ||
      ends_with(lfname, ".dicom") ||
      ends_with(lfname, ".ima"))
  {
    try 
    {
      const bool is_dicom = true;

      // TODO: Need a general output degrader (see vff_writer) and UI.
      if (nrrd_handle_->nrrd_->type == nrrdTypeUChar)
      {
	NrrdRange range;
	nrrdRangeSet(&range, nrrd_handle_->nrrd_, 0);
	ITKDatatypeHandle img = nrrd_to_itk_image(nrrd_handle_);
        return write_itk_series<unsigned char>(img, fname, hdr_series, series_desc, is_dicom,
					       (range.min + range.max) / 2.0, range.max - range.min);	
      }
      else
      {
        // These formats support 16 bit output.
        // TODO: dicom should, but the results were broken.
        // NOTE: Maybe dicom requires signed shorts like vff?
	NrrdDataHandle tmp;

	if( quantize )
	  tmp = Skinner::Histogram::UnuQuantize(nrrd_handle_, 0.0, 0.0, 16);
	else
	  tmp = Skinner::Histogram::UnuQuantize(nrrd_handle_, 0.0, 65535.0, 16);

	NrrdRange range;
	nrrdRangeSet(&range, tmp->nrrd_, 0);
        ITKDatatypeHandle img = nrrd_to_itk_image(tmp);
        return write_itk_series<unsigned short>(img, fname, hdr_series, series_desc, is_dicom, 
					       (range.min + range.max) / 2.0, range.max - range.min);	
      }
    }
    catch  ( itk::ExceptionObject & err )
    {
      cerr << "NrrdVolume::write ITK::ExceptionObject caught" << endl;
      cerr << err.GetDescription() << endl;
      painter_->set_status_safe(err.GetDescription());
      return false;
    }
    return true;
  }

  else
  {
    // Ask the user if they want to use .nrrd here.
    wxMessageDialog dialog(Painter::global_seg3dframe_pointer_,
			   std2wx("File " + fname + " unsupported output format.  " +
				  " Can only export dicom headers for files with " +
				  ".dcm, .dicom, or .ima extensions"),
			   _T("Unsupported format"),
			   wxOK | wxICON_ERROR);
    dialog.ShowModal();
    cerr << "NrrdVolume::writeWithHeader expected dicom format, but received " << lfname << endl;
    painter_->set_status_safe("NrrdVolume::writeWithHeader received unexpected file format");
    return false;
  }
  return false;
}


string
NrrdVolume::create_nrrd_header(const Nrrd *nrrd)
{
  ostringstream ostrm;

  ostrm << "NRRD0004\n";
  ostrm << "# Complete NRRD file format specification at:\n";
  ostrm << "# http://teem.sourceforge.net/nrrd/format.html\n";
  ostrm << "# Generated by Seg3D nrrd session writer.\n";
  ostrm << "type: " << nrrd_type_to_string(nrrd->type) << "\n";
  ostrm << "dimension: 3\n";
  ostrm << "sizes: "
	<< nrrd->axis[1].size << " "
	<< nrrd->axis[2].size << " "
	<< nrrd->axis[3].size << "\n";
  ostrm << "spacings: "
	<< nrrd->axis[1].spacing << " "
	<< nrrd->axis[2].spacing << " "
	<< nrrd->axis[3].spacing << "\n";
  ostrm << "axis mins: "
	<< nrrd->axis[1].min << " "
	<< nrrd->axis[2].min << " "
	<< nrrd->axis[3].min << "\n";
  ostrm << "axis maxs: "
	<< nrrd->axis[1].max << " "
	<< nrrd->axis[2].max << " "
	<< nrrd->axis[3].max << "\n";
  ostrm << "centerings: node node node\n";
  ostrm << "kinds: domain domain domain\n";

  if( AIR_ENDIAN == airEndianBig )
    ostrm << "endian: big\n";
  else
    ostrm << "endian: little\n";

  ostrm << "encoding: raw\n";
  ostrm << "\n";

  return ostrm.str();
}


VolumeSliceHandle
NrrdVolume::get_volume_slice(const SLIVR::Plane &plane)
{
  NrrdVolumeHandle parent = root();
  VolumeSlices_t::iterator siter = parent->all_slices_.begin();
  VolumeSlices_t::iterator send = parent->all_slices_.end();
  for (; siter != send; ++siter) {
    if ((*siter)->get_plane() == plane &&
        (*siter)->volume_->label_ == label_) {
      return (*siter);
    }
  }

  if (this != parent.get_rep()) {
    VolumeSlices_t::iterator siter = parent->all_slices_.begin();
    VolumeSlices_t::iterator send = parent->all_slices_.end();
    for (; siter != send; ++siter) {
      if ((*siter)->get_plane() == plane) {
        parent->all_slices_.push_back
          (new VolumeSlice(this, plane, (*siter)->nrrd_handle_, label_));
        return parent->all_slices_.back();
      }
    }
  }

  parent->all_slices_.push_back(new VolumeSlice(this, plane,0,label_));
  return parent->all_slices_.back();
}


void
NrrdVolume::purge_unused_slices()
{
  NrrdVolumeHandle parent = root();

  VolumeSlices_t new_all_slices;
  if (!dirty_) {
    for (unsigned int j = 0; j < parent->all_slices_.size(); ++j) {
      if (parent->all_slices_[j]->ref_cnt > 1) {
        new_all_slices.push_back(parent->all_slices_[j]);
      }
    }
  }
  dirty_ = false;
  parent->all_slices_ = new_all_slices;
}


ColorMapHandle
NrrdVolume::get_colormap()
{
  return colormap_;
}


Color
NrrdVolume::get_label_color()
{
  return label_color_;
}


MaterialHandle
NrrdVolume::get_label_material()
{
  MaterialHandle mat = new Material(label_color_);
  mat->transparency = opacity_;
  return mat;
}


GeomIndexedGroup *
NrrdVolume::get_geom_group()
{
  if (!geom_switch_.get_rep()) {
    GeomIndexedGroup *grp = new GeomIndexedGroup();
    geom_switch_ = new GeomSkinnerVarSwitch(grp, button_->layer_visible_);
    event_handle_t add_geom_switch_event =
      new SceneGraphEvent(geom_switch_, name_);
    EventManager::add_event(add_geom_switch_event);
  }
  GeomSkinnerVarSwitch *sw =
    dynamic_cast<GeomSkinnerVarSwitch *>(geom_switch_.get_rep());
  GeomIndexedGroup *ret_val =
    dynamic_cast<GeomIndexedGroup *>(sw->get_child().get_rep());

  return ret_val;
}


void
NrrdVolume::set_geom_switch(Skinner::Var<bool> v)
{
  get_geom_group(); // Force geometry creation first.
  GeomSkinnerVarSwitch *gs =
    dynamic_cast<GeomSkinnerVarSwitch *>(geom_switch_.get_rep());
  if (gs) gs->set_var(v);
}


NrrdDataHandle
NrrdVolume::extract_label_as_bit()
{
  Nrrd *n = nrrd_handle_->nrrd_;
  if (!label_ || !n || n->type != LabelNrrdType) return 0;

  NrrdDataHandle newnrrdh = new NrrdData();
  nrrdCopy(newnrrdh->nrrd_, n);
  label_type *data = (label_type *)newnrrdh->nrrd_->data;
  ASSERT(data);

  const size_t count = VolumeOps::nrrd_elem_count(newnrrdh);
  unsigned int bit = this->bit();

  for (size_t i = 0; i < count; ++i) {
    data[i] = (data[i] >> bit) & 1;
  }

  return newnrrdh;
}


int
NrrdVolume::bit()
{
  if (!label_) return -1;
  int lbit = 0;
  while (!(label_ & (1 << lbit))) ++lbit;
  return lbit;
}


void
NrrdVolume::change_type_from_float_to_bit(float val)
{
  NrrdDataHandle nrrdh = VolumeOps::float_to_bit(nrrd_handle_, val, label_);
  set_nrrd(nrrdh);
}


void
NrrdVolume::change_type_from_bit_to_float(float val)
{
  nrrd_handle_ = VolumeOps::bit_to_float(nrrd_handle_, label_, val);
  dirty_ = true;
}


void
NrrdVolume::clear()
{
  VolumeOps::clear_nrrd(nrrd_handle_);
}


size_t
NrrdVolume::numbytes()
{
  return VolumeOps::nrrd_data_size(nrrd_handle_);
}


void
NrrdVolume::set_slices_dirty()
{
  for (VolumeSlices_t::iterator iter = all_slices_.begin();
       iter != all_slices_.end(); ++iter) {
    (*iter)->set_tex_dirty();
  }
}


void
NrrdVolume::reset_clut()
{
  clut_min_ = data_min_; // - (data_max_ - data_min_)/254.0;
  clut_max_ = data_max_;
  set_slices_dirty();
}


bool
NrrdVolume::visible()
{
  return button_ && button_->layer_visible_();
}


NrrdVolumeHandle
NrrdVolume::root()
{
  NrrdVolumeHandle parent = this;
  while (parent->parent_.get_rep()) parent = parent->parent_;
  return parent;
}


int
NrrdVolume::depth()
{
  int depth = label_?1:0;
  return depth;
}


bool
NrrdVolume::only_one_label()
{
  return !(parent_.get_rep() || children_.size());
}


void
NrrdVolume::unparent()
{
  // Remove this reference from the parent's children list.
  if (parent_.get_rep())
  {
    NrrdVolumes pchildren;
    for (size_t i = 0; i < parent_->children_.size(); i++)
    {
      if (parent_->children_[i].get_rep() != this)
      {
        pchildren.push_back(parent_->children_[i]);
      }
    }
    parent_->children_ = pchildren;
  }

  if (children_.size())
  {
    // Set a new parent.
    children_[0]->parent_ = parent_;
    for (size_t i = 1; i < children_.size(); i++)
    {
      // Move the children to the new parent.
      children_[0]->children_.push_back(children_[i]);
      // Set the children's parent to the new parent.
      children_[i]->parent_ = children_[0];
    }
  }
}


GeomHandle
NrrdVolume::isosurface_label()
{
  ASSERT(label_ > 0);

  Nrrd *nrrd = nrrd_handle_->nrrd_;
  label_type *data = (label_type *)nrrd->data;

  const size_t isize = nrrd->axis[1].size;
  const size_t jsize = nrrd->axis[2].size;
  const size_t ksize = nrrd->axis[3].size;

  const size_t ijsize = isize * jsize;

  GeomFastQuads *quads = new GeomTranspQuads;

  for (unsigned int k = 0; k <= ksize; k++)
  {
    for (unsigned int j = 0; j <= jsize; j++)
    {
      for (unsigned int i = 0; i <= isize; i++)
      {
        label_type kval = 0;
        if (k > 0 && i < isize && j < jsize)
        {
          kval = data[(k-1) * ijsize + j * isize + i];
        }
        
        label_type jval = 0;
        if (j > 0 && k < ksize && i < isize)
        {
          jval = data[k * ijsize + (j-1) * isize + i];
        }

        label_type ival = 0;
        if (i > 0 && j < jsize && k < ksize)
        {
          ival = data[k * ijsize + j * isize + (i-1)];
        }
        
        label_type val = 0;
        if (i < isize && j < jsize && k < ksize)
        {
          val = data[k * ijsize + j * isize + i];
        }
         
        if (val != kval && (val & label_ || kval & label_))
        {
          quads->add(Point(i+0, j+0, k+0),
                     Point(i+0, j+1, k+0),
                     Point(i+1, j+1, k+0),
                     Point(i+1, j+0, k+0));
        }

        if (val != jval && (val & label_ || jval & label_))
        {
          quads->add(Point(i+0, j+0, k+0),
                     Point(i+1, j+0, k+0),
                     Point(i+1, j+0, k+1),
                     Point(i+0, j+0, k+1));
        }
        if (val != ival && (val & label_ || ival & label_))
        {
          quads->add(Point(i+0, j+0, k+0),
                     Point(i+0, j+0, k+1),
                     Point(i+0, j+1, k+1),
                     Point(i+0, j+1, k+0));
        }
      }
    }
  }

  return quads;
}


}
