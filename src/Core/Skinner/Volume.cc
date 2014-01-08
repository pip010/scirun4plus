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
//    File   : Volume.cc
//    Author : McKay Davis
//    Date   : Tue Feb  6 17:22:56 2007

#include <Core/Skinner/Volume.h>
#include <Core/Skinner/Histogram.h>
#include <Core/Skinner/GeomSkinnerVarSwitch.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/SceneGraphEvent.h>
#include <Core/Events/ColorMap2Event.h>
#include <Core/Events/LoadVolumeEvent.h>
#include <Core/Util/Environment.h>
#include <Core/Util/FileUtils.h>
#include <Core/Volume/ColorMap2.h>
#include <Core/Volume/CM2Widget.h>
#include <Core/Volume/VolumeRenderer.h>
#include <Core/Datatypes/ColorMap.h>
#include <Core/Datatypes/NrrdData.h>
#include <Core/Math/MiscMath.h>

#include <slivr/TextureRenderer.h>
#include <slivr/Plane.h>
#include <string.h>

#include <boost/filesystem.hpp>

namespace bfs=boost::filesystem;

namespace SCIRun {
namespace Skinner {


Volume::Volume(Variables *vars) :
  Parent(vars),
  filename_(vars, "Volume::filename"),
  sample_rate_(vars, "Volume::sample_rate",1.5),
  shading_(vars, "Volume::lighting::on", true),
  ambient_(vars, "Volume::lighting::ambient",0.5),
  diffuse_(vars, "Volume::lighting::diffuse",0.5),
  specular_(vars, "Volume::lighting::specular",1.0),
  shine_(vars, "Volume::lighting::shine",16.0),
  mip_(vars, "Volume::mip", 0),
  interpolate_(vars, "Volume::interpolate",1),
  slice_alpha_(vars, "Volume::slice_alpha",0.0),
  volume_texture_(0),
  nrrd_(0),
  volume_renderer_(0)
{
  REGISTER_CATCHER_TARGET(Volume::SetParameters);

  cm2_ = new ColorMap2();
  create_default_colormap2(cm2_);
  cm2v_.push_back(cm2_.get_rep());
}


Volume::~Volume()
{
}


BaseTool::propagation_state_e
Volume::process_event(event_handle_t &event)
{
  ColorMap2Event *cm2e = dynamic_cast<ColorMap2Event *>(event.get_rep());
  if (cm2e)
  {
    if (volume_renderer_) { volume_renderer_->lock.lock(); }
    ColorMap2Event::ColorMap2Vector_t cm2s = cm2e->get_data();
    if (cm2s.size())
    {
      cm2_ = cm2s[0];
      cm2v_[0] = cm2_.get_rep();
      if (volume_renderer_) { volume_renderer_->set_colormap2(cm2v_); }
    }
    if (volume_renderer_) { volume_renderer_->lock.unlock(); }
  }

  LoadVolumeEvent *lve = dynamic_cast<LoadVolumeEvent *>(event.get_rep());
  if (lve) {
    NrrdDataHandle nrrd = lve->get_nrrd();
    if (nrrd.get_rep() != nrrd_.get_rep()) {
      if (nrrd.get_rep() == 0)
      {
        clear_renderer();
      }
      else
      {
        set_nrrd(nrrd);
        if (! lve->get_filename().empty() )
        {
          bfs::path cmap2_filename_path = bfs::change_extension(lve->get_filename().c_str(), ".cmap2");
          load_colormap2(cmap2_filename_path.string());
        }
        create_renderer();

        if (lve->get_change_histogram())
        {
          event_handle_t uhe = new UpdateHistogramEvent(nrrd_);
          EventManager::add_event(uhe);
        }
      }
    }
  }

  return Parent::process_event(event);
}


void
Volume::set_nrrd(NrrdDataHandle nin)
{
  nin = resize_nrrd(nin);

  if (nin->nrrd_->type != nrrdTypeUChar) {
    nin = Histogram::UnuQuantize(nin);
  }

  for (unsigned int a = 0; a < nin->nrrd_->dim; ++a) {
    if (nin->nrrd_->axis[a].center == nrrdCenterUnknown) {
      nin->nrrd_->axis[a].center = nrrdCenterNode;
    }

    if (airIsNaN(nin->nrrd_->axis[a].min) ||
        airIsNaN(nin->nrrd_->axis[a].max) ||
        nin->nrrd_->axis[a].min == nin->nrrd_->axis[a].max) {
      if (airIsNaN(nin->nrrd_->axis[a].spacing)) {
        nin->nrrd_->axis[a].spacing = 1.0;
      }
      nin->nrrd_->axis[a].min = 0.0;
      nin->nrrd_->axis[a].max =
        nin->nrrd_->axis[a].spacing * (nin->nrrd_->axis[a].size - 1);

      // Images result in degenerate spacing for third axis, pad to 1.0.
      if (nin->nrrd_->axis[a].size <= 1)
      {
        nin->nrrd_->axis[a].max = 1.0;
      }
    }

    if (nin->nrrd_->axis[a].min > nin->nrrd_->axis[a].max) {
      std::swap(nin->nrrd_->axis[a].min, nin->nrrd_->axis[a].max);
    }

    if (nin->nrrd_->axis[a].spacing < 0.0) {
      nin->nrrd_->axis[a].spacing *= -1.0;
    }
  }
  nrrd_ = nin;
}


bool
Volume::load_colormap2(const string &fn)
{
  // Read the file.
  ColorMap2 *cmap2 = new ColorMap2();
  ColorMap2Handle icmap = cmap2;

  if (!validFile(fn))
  {
    create_default_colormap2(icmap);
  }
  else
  {
    Piostream *stream = auto_istream(fn, 0);
    if (!stream) {
      create_default_colormap2(icmap);
    }
    else
    {
      try {
        Pio(*stream, icmap);
      } catch (...) {
        cerr << "Error loading "+fn+", using default instead." << std::endl;
        create_default_colormap2(icmap);
      }
      delete stream;
    }
  }

  if (volume_renderer_) { volume_renderer_->lock.lock(); }
  cm2_ = icmap;
  cm2v_[0] = cm2_.get_rep();
  if (volume_renderer_)
  {
    volume_renderer_->set_colormap2(cm2v_);
    volume_renderer_->lock.unlock();
  }

  ColorMap2Event *cm2e = new ColorMap2Event(cm2_);
  EventManager::add_event(cm2e);

  return true;
}


void
Volume::create_default_colormap2(ColorMap2Handle icmap2)
{
  icmap2->widgets().clear();
  SLIVR::RectangleCM2Widget *rect =
    new SLIVR::RectangleCM2Widget(SLIVR::CM2_RECTANGLE_1D,
                                  0.25, 0.0, 0.75, 1.0, 1.0);
  rect->set_color(SCIRun::Color(1.0, 1.0, 1.0));
  rect->set_alpha(0.25);
  rect->set_faux(false);
  icmap2->widgets().push_back(rect);
}


void
Volume::create_renderer()
{
  // Volume creates and owns the volume rendering related objects so
  // that it can keep the handles and locking on the slivr side
  // correct.

  const int card_mem = 64;

  // Create the texture as needed.
  if (volume_texture_.get_rep() == 0)
  {
    volume_texture_ = new Texture;
  }

  // Create the volume renderer as needed.
  if (volume_renderer_ == NULL)
  {
    vector<Plane> *planes = new vector<Plane>;
    volume_renderer_ = new VolumeRenderer(NULL, NULL, cm2v_, *planes,
                                          Round(card_mem*1024*1024*0.8));
    Skinner::Var<bool>
      visible(get_vars(),"Painter::volume_visible", true);
    volume_renderer_handle_ =
      new GeomSkinnerVarSwitch(volume_renderer_, visible);
  }

  // Rebuild the contents of the texture and push it to the volume renderer.
  volume_renderer_->lock.lock();
  volume_texture_->build(nrrd_->nrrd_, 0, 0, 256, 0, 0, card_mem);
  volume_renderer_->set_texture(volume_texture_);
  volume_renderer_->lock.unlock();

  // Add/resend this scene graph event.  Even if it is already there
  // this causes it to recompute the bounds properly.
  event_handle_t scene_event =
    new SceneGraphEvent(volume_renderer_handle_, "FOO Transparent");
  SetParameters(scene_event);
  EventManager::add_event(scene_event);
}


void
Volume::clear_renderer()
{
  if (volume_renderer_)
  {
    volume_renderer_->lock.lock();
    volume_renderer_->set_texture(0);
    volume_renderer_->lock.unlock();
  }
  volume_texture_ = 0;
  nrrd_ = 0;

  // Go ahead and send the clear scene event for completeness.
  SceneGraphEvent *raw_scene_event =
    new SceneGraphEvent(0, "FOO Transparent");
  event_handle_t scene_event = raw_scene_event;
  EventManager::add_event(scene_event);
}


BaseTool::propagation_state_e
Volume::SetParameters(event_handle_t& /*event*/)
{
  if (!volume_renderer_) return STOP_E;

  volume_renderer_->lock.lock();

  volume_renderer_->set_shading(shading_);
  volume_renderer_->set_material(ambient_, diffuse_, specular_, shine_);
  volume_renderer_->set_sampling_rate(sample_rate_);
  volume_renderer_->set_slice_alpha(slice_alpha_);
  volume_renderer_->set_interp(interpolate_);
  volume_renderer_->set_mode(mip_ ? SLIVR::TextureRenderer::MODE_MIP :
                          SLIVR::TextureRenderer::MODE_OVER);

  volume_renderer_->lock.unlock();

  return CONTINUE_E;
}


#define BIG 1<<25


static size_t
compute_size(size_t samples[4])
{
  size_t size = 1;
  for (int i = 0; i < 4; i++)
  {
    size *= samples[i];
  }
  return size;
}


static size_t
find_next_size(size_t a)
{
  static size_t chart[12] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 510, 1020, 2042 };
  for (int i = 1; i < 12; i++)
  {
    if (chart[i] >= a) { return chart[i-1]; }
  }
  return chart[11];
}


static void
compute_new_size(size_t samples[4], double spacings[4], size_t big)
{
  size_t size = compute_size(samples);
  while (size > big)
  {
    int mi = 1;
    for (int i = 2; i < 4; i++)
    {
      if (samples[i] != 1 &&
          spacings[i] < spacings[mi] ||
          spacings[i] == spacings[mi] && samples[i] > samples[mi])
      {
        mi = i;
      }
    }
    size_t nsize = find_next_size(samples[mi]); 
    spacings[mi] *= samples[mi] / (double)nsize;
    samples[mi] = nsize;
    size = compute_size(samples);
  }    
}


NrrdDataHandle
Volume::resize_nrrd(NrrdDataHandle nrrd)
{
  size_t i;

  Nrrd *nin = nrrd->nrrd_;

  size_t size = 1;
  for (i = 0; i < nin->dim; i++)
  {
    size *= nin->axis[i].size;
  }

  if (size < BIG)
  {
    return nrrd;
  }

  NrrdResampleInfo *info = nrrdResampleInfoNew();
  NrrdKernel *kern = nrrdKernelAQuartic;
  double p[NRRD_KERNEL_PARMS_NUM];
  memset(p, 0, NRRD_KERNEL_PARMS_NUM * sizeof(double));
  p[0] = 1.0;
  p[1] = 0.0834;

  size_t samples[4];
  double spacings[4];
  for (i = 0; i < 4; i++)
  {
    samples[i] = nin->axis[i].size;
    spacings[i] = nin->axis[i].spacing;
  }

  compute_new_size(samples, spacings, BIG);

  for (int a = 0; a < 4; a++)
  {
    info->kernel[a] = a==0?0:kern;

    info->samples[a] = samples[a];

    memcpy(info->parm[a], p, NRRD_KERNEL_PARMS_NUM * sizeof(double));

    if (info->kernel[a] && 
    	(!(airExists(nin->axis[a].min) && airExists(nin->axis[a].max))))
    {
      nrrdAxisInfoMinMaxSet(nin, a, nin->axis[a].center ? 
                            nin->axis[a].center : nrrdDefaultCenter);
    }

    info->min[a] = nin->axis[a].min;
    info->max[a] = nin->axis[a].max;
  }    
  info->boundary = nrrdBoundaryBleed;
  info->type = nin->type;
  info->renormalize = AIR_TRUE;

  NrrdDataHandle nout_handle = new NrrdData;
  Nrrd *nout = nout_handle->nrrd_;
  if (nrrdSpatialResample(nout, nin, info)) {
    char *err = biffGetDone(NRRD);
    string errstr(err);
    free(err);
    throw "Trouble resampling: " + errstr;
  }
  nrrdResampleInfoNix(info); 

  nrrdKeyValueCopy(nout, nin);

  return nout_handle;
}


} // namespace Skinnner
} // namespace SCIRun
