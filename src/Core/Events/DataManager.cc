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
//    File   : DataManager.cc
//    Author : Martin Cole
//    Date   : Tue Sep 12 10:35:09 2006


#include <Core/Events/DataManager.h>
#include <Core/Events/OpenGLBase.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

namespace SCIRun {

  DataManager* DataManager::dm_ = 0;

} // end namespace SCIRun

using namespace SCIRun;


DataManager::DataManager() :
  lock_("DataManager lock"),
  next_id_(1) // init to 1, 0 is invalid.
{
  bundles_.clear();
  colormaps_.clear();
  colormap2s_.clear();
  fields_.clear();
  geoms_.clear();
  mats_.clear();
  nrrds_.clear();
  paths_.clear();
  texs_.clear();
  views_.clear();
}

DataManager::~DataManager()
{
  bundles_.clear();
  colormaps_.clear();
  colormap2s_.clear();
  fields_.clear();
  geoms_.clear();
  mats_.clear();
  nrrds_.clear();
  paths_.clear();
  texs_.clear();
  views_.clear();
}

DataManager*
DataManager::get_dm() 
{
  if (dm_ == 0) {
    dm_ = new DataManager();
  }
  return dm_;
}

BundleHandle
DataManager::get_bundle(size_t id)
{
  BundleHandle bundle;
  lock_.readLock();
  bundle = bundles_[id];
  lock_.readUnlock();
  return bundle;
}

ColorMapHandle
DataManager::get_colormap(size_t id)
{
  ColorMapHandle cm;
  lock_.readLock();
  cm = colormaps_[id];
  lock_.readUnlock();
  return cm;
}

ColorMap2Handle
DataManager::get_colormap2(size_t id)
{
  ColorMap2Handle cm;
  lock_.readLock();
  cm = colormap2s_[id];
  lock_.readUnlock();
  return cm;
}

FieldHandle
DataManager::get_field(size_t id)
{
  FieldHandle fld;
  lock_.readLock();
  fld = fields_[id];
  lock_.readUnlock();
  return fld;
}

DataManager::geom_info_t
DataManager::get_geom(size_t id)
{
  geom_info_t gi;
  lock_.readLock();
  gi = geoms_[id];
  lock_.readUnlock();
  return gi;
}

MatrixHandle
DataManager::get_matrix(size_t id)
{
  MatrixHandle m;
  lock_.readLock();
  m = mats_[id];
  lock_.readUnlock();
  return m;
}


NrrdDataHandle
DataManager::get_nrrd(size_t id)
{
  NrrdDataHandle nd;
  lock_.readLock();
  nd = nrrds_[id];
  lock_.readUnlock();
  return nd;
}

PathHandle
DataManager::get_path(size_t id)
{
  PathHandle handle;
  lock_.readLock();
  handle = paths_[id];
  lock_.readUnlock();
  return handle;
}

TextureHandle
DataManager::get_texture(size_t id)
{
  TextureHandle th;
  lock_.readLock();
  th = texs_[id];
  lock_.readUnlock();
  return th;
}

OpenGLBase*
DataManager::get_view(size_t id)
{
  OpenGLBase *v;
  lock_.readLock();
  v = views_[id];
  lock_.readUnlock();
  return v;
}

size_t
DataManager::add_bundle(BundleHandle bh)
{
  size_t rval = 0;
  lock_.writeLock();
  rval = next_id_;
  bundles_[next_id_++] = bh;
  lock_.writeUnlock();
  return rval;
}

size_t
DataManager::add_colormap(ColorMapHandle cmh)
{
  size_t rval = 0;
  lock_.writeLock();
  rval = next_id_;
  colormaps_[next_id_++] = cmh;
  lock_.writeUnlock();
  return rval;
}

size_t
DataManager::add_colormap2(ColorMap2Handle cmh)
{
  size_t rval = 0;
  lock_.writeLock();
  rval = next_id_;
  colormap2s_[next_id_++] = cmh;
  lock_.writeUnlock();
  return rval;
}

size_t
DataManager::add_field(FieldHandle fh)
{
  size_t rval = 0;
  lock_.writeLock();
  rval = next_id_;
  fields_[next_id_++] = fh;
  lock_.writeUnlock();
  return rval;
}

size_t
DataManager::add_geom(GeomHandle gh, std::string name)
{
  size_t rval = 0;
  geom_info_t p;
  p.first = gh;
  p.second = name;
  lock_.writeLock();
  rval = next_id_;
  geoms_[next_id_++] = p;
  lock_.writeUnlock();
  return rval;
}

size_t
DataManager::add_matrix(MatrixHandle mh)
{
  size_t rval = 0;
  lock_.writeLock();
  rval = next_id_;
  mats_[next_id_++] = mh;
  lock_.writeUnlock();
  return rval;
}

size_t
DataManager::add_nrrd(NrrdDataHandle nh)
{
  size_t rval = 0;
  lock_.writeLock();
  rval = next_id_;
  nrrds_[next_id_++] = nh;
  lock_.writeUnlock();
  return next_id_ - 1;
}

size_t
DataManager::add_path(PathHandle handle)
{
  size_t rval = 0;
  lock_.writeLock();
  rval = next_id_;
  paths_[next_id_++] = handle;
  lock_.writeUnlock();
  return next_id_ - 1;
}

size_t
DataManager::add_texture(TextureHandle th)
{
  size_t rval = 0;
  lock_.writeLock();
  rval = next_id_;
  texs_[next_id_++] = th;
  lock_.writeUnlock();
  return rval;
}


//! \todo {add interface to delete the viewer objects}
size_t
DataManager::add_view(OpenGLBase *v)
{
  size_t rval = 0;
  lock_.writeLock();
  rval = next_id_;
  views_[next_id_++] = v;
  lock_.writeUnlock();
  return rval;
}

size_t       
DataManager::get_selection_target()
{ 
  lock_.readLock();
  size_t rval = sel_fid_;
  lock_.readUnlock();

  return rval; 
}

void 
DataManager::set_selection_target(size_t st)
{ 
  lock_.writeLock();
  sel_fid_ = st; 
  lock_.writeUnlock();
}

size_t
DataManager::load_field(std::string fname)
{
  FieldHandle fld;
  PiostreamPtr stream = auto_istream(fname);
  if (!stream) {
    std::cerr << "Couldn't open file: "<< fname << std::endl;
    return 0;
  }
  Pio(*stream, fld);
  return add_field(fld);
}


size_t
DataManager::load_matrix(std::string fname)
{
  MatrixHandle m;
  PiostreamPtr stream = auto_istream(fname);
  if (!stream) {
    std::cerr << "Couldn't open file: "<< fname << std::endl;
    return 0;
  }
  Pio(*stream, m);
  return add_matrix(m);
}


size_t
DataManager::load_nrrd(std::string fname)
{
  NrrdDataHandle n;
  PiostreamPtr stream = auto_istream(fname);
  if (!stream) {
    std::cerr << "Couldn't open file: "<< fname << std::endl;
    return 0;
  }
  Pio(*stream, n);
  return add_nrrd(n);
}
