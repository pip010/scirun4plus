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
//    File   : DataManager.h
//    Author : Martin Cole
//    Date   : Tue Sep 12 09:55:13 2006


#if !defined(DataManager_h)
#define DataManager_h

#include <Core/Thread/CrowdMonitor.h>
#include <Core/Datatypes/Bundle.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/NrrdData.h>
#include <Core/Events/ToolManager.h>
#include <Core/Events/EventManager.h>
#include <Core/Geom/GeomObj.h>
#include <Core/Datatypes/ColorMap.h>
#include <Core/Geom/Path.h>
#include <Core/Volume/Texture.h>
#include <Core/Volume/ColorMap2.h>
#include <Core/Events/share.h>

namespace SCIRun {
class OpenGLBase;

class SCISHARE DataManager 
{
 public:
  //! name is important, it triggers transparency among other things.
  //! geom_info_t == (handle, name)
  typedef std::pair<GeomHandle, std::string> geom_info_t;

  DataManager();
  virtual ~DataManager();

  // return the valid dm pointer.
  static DataManager* get_dm(); 

  // Keep these in alphebetical order based on get_XYZ()
  BundleHandle         get_bundle(size_t id);
  ColorMapHandle       get_colormap(size_t id);
  ColorMap2Handle      get_colormap2(size_t id);
  FieldHandle          get_field(size_t id);
  geom_info_t          get_geom(size_t id);
  MatrixHandle         get_matrix(size_t id);
  NrrdDataHandle       get_nrrd(size_t id);
  PathHandle           get_path(size_t id);
  TextureHandle        get_texture(size_t id);
  OpenGLBase*          get_view(size_t id);

  //! Store a new handle.
  size_t add_bundle(BundleHandle);
  size_t add_colormap(ColorMapHandle cmh);
  size_t add_colormap2(ColorMap2Handle cmh);
  size_t add_field(FieldHandle);
  size_t add_geom(GeomHandle gh, std::string name);
  size_t add_matrix(MatrixHandle);
  size_t add_nrrd(NrrdDataHandle);
  size_t add_path(PathHandle);
  size_t add_texture(TextureHandle cmh);
  size_t add_view(OpenGLBase *v);

  //! return the id of the 'selected' field.
  size_t get_selection_target();
  void   set_selection_target(size_t st);

  //! load from file.
  size_t load_field(std::string fname);
  size_t load_matrix(std::string fname);
  size_t load_nrrd(std::string fname);

 private:
  static DataManager*                         dm_;

  //! Regulate access to the handle maps.
  CrowdMonitor                                lock_;

  std::map<size_t, BundleHandle>                   bundles_;
  std::map<size_t, ColorMapHandle>                 colormaps_;
  std::map<size_t, ColorMap2Handle>                colormap2s_;
  std::map<size_t, FieldHandle>                    fields_;
  //! stored scene graph events, with associated name.
  std::map<size_t, geom_info_t>                    geoms_;
  std::map<size_t, MatrixHandle>                   mats_;
  std::map<size_t, NrrdDataHandle>                 nrrds_;
  std::map<size_t, PathHandle>                     paths_;
  std::map<size_t, TextureHandle>                  texs_;
  std::map<size_t, OpenGLBase*>                    views_;

  size_t                                      sel_fid_;
  			
  //! next_id_ is shared by all the maps,
  //! separate into additional it becomes an issue.
  size_t                                      next_id_;

};

}

#endif //DataManager_h
