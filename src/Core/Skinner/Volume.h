//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//  License for the specific language governing rights and limitations under
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
//    File   : Volume.h
//    Author : McKay Davis
//    Date   : Thu Feb  8 14:54:30 2007

#ifndef SKINNER_VOLUME_H
#define SKINNER_VOLUME_H

#include <Core/Skinner/Parent.h>
#include <Core/Volume/Texture.h>
#include <Core/Volume/ColorMap2.h>
#include <Core/Geom/GeomObj.h>

namespace SCIRun {
class VolumeRenderer;

namespace Skinner {


class Volume : public Parent {
public:
  Volume (Variables *variables);
  virtual ~Volume();
  CatcherFunction_t           process_event;

private:
  CatcherFunction_t           SetParameters;

  void                        set_nrrd(NrrdDataHandle nrrd);
  bool                        load_colormap2(const string &filename);
  void                        create_default_colormap2(ColorMap2Handle h);
  void                        create_renderer();
  void                        clear_renderer();
  NrrdDataHandle              resize_nrrd(NrrdDataHandle nrrd);

  Var<string>                 filename_;
  Var<string>                 cmap2_filename_;
  Var<double>                 sample_rate_;
  Var<bool>                   shading_;
  Var<double>                 ambient_;
  Var<double>                 diffuse_;
  Var<double>                 specular_;
  Var<double>                 shine_;
  Var<bool>                   mip_;
  Var<bool>                   interpolate_;
  Var<double>                 slice_alpha_;

  TextureHandle               volume_texture_;
  NrrdDataHandle              nrrd_;
  VolumeRenderer *            volume_renderer_;
  vector<ColorMap2Handle>     cm2v_;
  GeomHandle                  volume_renderer_handle_;
  ColorMap2Handle             cm2_;
};


}
}

#endif
