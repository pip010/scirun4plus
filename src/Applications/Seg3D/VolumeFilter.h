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
//    File   : VolumeFilter.h
//    Author : McKay Davis
//    Date   : Sun Nov  5 21:37:09 2006

#ifndef SEG3D_VolumeFilter_H
#define SEG3D_VolumeFilter_H

#include <sci_comp_warn_fixes.h>
#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/NrrdVolume.h>
#include <Core/Datatypes/NrrdToITK.h>
#include <Core/Skinner/Parent.h>
#include <Core/Skinner/Variables.h>
#include <Core/Util/Assert.h>

#include <sci_defs/insight_defs.h>

#  include <Core/Datatypes/ITKDatatype.h>
#  include <itkCommand.h>

namespace SCIRun {

using SCIRun::ITKDatatypeHandle;


class VolumeFilterBase : public Skinner::Parent
{
public:
  VolumeFilterBase      (NrrdVolumeHandle volume  = 0,
                         Skinner::Variables *vars = 0);
  
  virtual               ~VolumeFilterBase();

  void                  set_volume(NrrdVolumeHandle vol) { volume_ = vol; }
  void                  set_offset(double offset) { progress_offset_ = offset; }
  void                  stop() { abort_ = true; }
  bool                  stopped() { return abort_; }

  virtual void          operator()(NrrdDataHandle nrrdh);
  virtual void          start(NrrdDataHandle nrrdh);
  virtual void          update();
  virtual void          update_progress();
  virtual void          show_progress( bool val = true )
    { progress_show_ = val; };

protected:
  NrrdVolumeHandle      volume_;
  bool                  abort_;
  double                progress_;
  double                progress_offset_;
  int                   passes_;  

  bool                  progress_show_;
};


// todo, integrate scirun progressreporter class to this callback
template <class FilterType>
class VolumeFilter : public VolumeFilterBase
{
  typedef typename      FilterType::Pointer FilterPointer;
public:
  VolumeFilter          (NrrdVolumeHandle volume  = 0,
                         Skinner::Variables *vars = 0);
  virtual               ~VolumeFilter() {};
  FilterPointer         operator->() { return filter_; };
  void                  operator()(NrrdDataHandle nrrdh = 0);
  void                  update();
  void                  start(NrrdDataHandle nrrdh = 0);
private:
  void                  filter_callback(itk::Object *, 
                                        const itk::EventObject &);
  
  void                  filter_callback_const (const itk::Object *, 
                                               const itk::EventObject &);

  FilterPointer         filter_;
  
};


template <class FilterType>  
VolumeFilter<FilterType>::VolumeFilter(NrrdVolumeHandle volume,
                                       Skinner::Variables *vars) 
  : VolumeFilterBase(volume, vars),
    filter_(FilterType::New())
{
}


template <class FilterType>  
void
VolumeFilter<FilterType>::update() 
{
  if (volume_->painter_->filters_.find(this) == 
      volume_->painter_->filters_.end()) {
    volume_->painter_->filters_.insert(this);
  } else {
    // This filter is already running
    return;
  }

  abort_ = false;

  bool error = false;

  try {
    filter_->Update();
  } catch (itk::ExceptionObject &err) {
    if (!abort_) {
      cerr << "ITK Exception: \n";
      err.Print(cerr);
      error = true;
    } 
  } catch (...) {
    cerr << "ITK Filter error!\n";
    error = true;
    return;
  }

  if (volume_->painter_->filters_.find(this) != 
      volume_->painter_->filters_.end()) {
    volume_->painter_->filters_.erase(this);
  }

  typedef typename FilterType::OutputImageType ImageOutT;
  if (!error) {
    SCIRun::ITKDatatypeHandle output_img = new SCIRun::ITKDatatype();
    output_img->data_ = filter_->GetOutput();
    volume_->nrrd_handle_ = 
      itk_image_to_nrrd<typename ImageOutT::PixelType, 3>(output_img);
  }
}

  


template <class FilterType>  
void
VolumeFilter<FilterType>::operator()(NrrdDataHandle nrrd_handle) 
{
  this->start(nrrd_handle);
}

template <class FilterType>  
void
VolumeFilter<FilterType>::start(NrrdDataHandle nrrd_handle) 
{
  if (!nrrd_handle.get_rep()) {
    nrrd_handle = volume_->nrrd_handle_;
  }
  typedef typename FilterType::InputImageType ImageInT;
  typedef typename FilterType::OutputImageType ImageOutT;
  typedef typename 
    itk::MemberCommand< VolumeFilter<FilterType> > RedrawCommandType;

  typename RedrawCommandType::Pointer callback = RedrawCommandType::New();

  callback->SetCallbackFunction
    (this, &VolumeFilter<FilterType>::filter_callback);

  callback->SetCallbackFunction
    (this, &VolumeFilter<FilterType>::filter_callback_const);

  filter_->AddObserver(itk::ProgressEvent(), callback);
  filter_->AddObserver(itk::IterationEvent(), callback);
  
  if (nrrd_handle.get_rep()) {
    ITKDatatypeHandle img_handle = nrrd_to_itk_image(nrrd_handle);  
    ImageInT *imgp = dynamic_cast<ImageInT*>(img_handle->data_.GetPointer());
    if (!imgp) {
      return;
    }
    filter_->SetInput(imgp);
  }

  update();
}


template <class FilterType>
void
VolumeFilter<FilterType>::filter_callback(itk::Object *object,
                                               const itk::EventObject &event)
{
  typedef typename FilterType::InputImageType ImageInT;
  typedef typename FilterType::OutputImageType ImageOutT;
  typedef itk::ProcessObject::Pointer ProcessPointer;

  ProcessPointer process = dynamic_cast<itk::ProcessObject *>(object);
  ASSERT(process);
  if (!process) return;

  progress_ = process->GetProgress();
  if (typeid(itk::ProgressEvent) == typeid(event))
  {
    update_progress();
  }

  if (typeid(itk::IterationEvent) == typeid(event))
  {
    update_progress();

    if (volume_.get_rep()) {
      if (volume_->lock.try_lock()) {
        ITKDatatypeHandle imgh = new ITKDatatype();
        imgh->data_ = filter_->GetOutput();
        
        if (!dynamic_cast<ImageOutT *>(imgh->data_.GetPointer())) return;
        
        typedef typename FilterType::OutputImageType::PixelType OutT;
        volume_->nrrd_handle_ = itk_image_to_nrrd<OutT, 3>(imgh);
        volume_->set_dirty();
        volume_->painter_->extract_all_window_slices();
        volume_->painter_->redraw_all();
        volume_->lock.unlock();        
      }
    }
  }

  if (abort_) {
    process->AbortGenerateDataOn();
  }
}


template <class FilterType>
void
VolumeFilter<FilterType>::
filter_callback_const(const itk::Object *, const itk::EventObject &)
{
  cerr << "filter_callback_const not implemented\n;";
}


}
#endif
