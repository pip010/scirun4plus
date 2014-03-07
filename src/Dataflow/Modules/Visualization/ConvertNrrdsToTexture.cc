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
//    File   : ConvertNrrdsToTexture.cc
//    Author : Milan Ikits
//    Date   : Fri Jul 16 03:28:21 2004


#include <slivr/ShaderProgramARB.h>
#include <slivr/VideoCardInfo.h>

#include <Core/Util/StringUtil.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/NrrdPort.h>
#include <Dataflow/Network/Ports/TexturePort.h>

#include <Dataflow/Modules/Visualization/share.h>

namespace SCIRun {

using namespace SLIVR;

class SCISHARE ConvertNrrdsToTexture : public Module
{
  public:
    ConvertNrrdsToTexture(GuiContext* ctx, const std::string& name = "ConvertNrrdsToTexture",
                          SchedClass sc = Source, const std::string& cat = "Visualization", 
                          const std::string& pack = "SCIRun");
    virtual ~ConvertNrrdsToTexture() {}
    virtual void execute();

  protected:
    TextureHandle tHandle_;

    GuiDouble gui_vminval_;
    GuiDouble gui_vmaxval_;
    GuiDouble gui_gminval_;
    GuiDouble gui_gmaxval_;
    GuiDouble gui_mminval_;
    GuiDouble gui_mmaxval_;

    GuiInt gui_fixed_;
    GuiInt gui_card_mem_;
    GuiInt gui_card_mem_auto_;
    GuiInt gui_uchar_;

    GuiInt gui_histogram_;    
    GuiDouble gui_gamma_;

    int card_mem_;
    int is_uchar_;
    int vnrrd_last_generation_;
    int gnrrd_last_generation_;
    int mnrrd_last_generation_;
    double vminval_, vmaxval_;
    double gminval_, gmaxval_;
    double mminval_, mmaxval_;
};

DECLARE_MAKER(ConvertNrrdsToTexture)
  
ConvertNrrdsToTexture::ConvertNrrdsToTexture(GuiContext* ctx, const std::string& name,
                                             SchedClass sc, const std::string& cat, const std::string& pack) :
  Module(name, ctx, sc, cat, pack),
  gui_vminval_(get_ctx()->subVar("vmin"), 0),
  gui_vmaxval_(get_ctx()->subVar("vmax"), 1),
  gui_gminval_(get_ctx()->subVar("gmin"), 0),
  gui_gmaxval_(get_ctx()->subVar("gmax"), 1),
  gui_mminval_(get_ctx()->subVar("mmin"), 0),
  gui_mmaxval_(get_ctx()->subVar("mmax"), 1),
  gui_fixed_(get_ctx()->subVar("is_fixed"), 0),
  gui_card_mem_(get_ctx()->subVar("card_mem"), 16),
  gui_card_mem_auto_(get_ctx()->subVar("card_mem_auto"), 1),
  gui_uchar_(get_ctx()->subVar("is_uchar"), 1),
  gui_histogram_(get_ctx()->subVar("histogram"),1),
  gui_gamma_(get_ctx()->subVar("gamma"),0.5),
  card_mem_(video_card_memory_size()),
  is_uchar_(1),
  vnrrd_last_generation_(-1),
  gnrrd_last_generation_(-1),
  mnrrd_last_generation_(-1)
{
}


class SortFloatVector : public std::binary_function<size_t,size_t,bool>
{
  public:
    SortFloatVector(const std::vector<float>& hist) : data_(hist)
    {
    }
    
    bool operator()(size_t i1, size_t i2) const
    {
      return (data_[i1] < data_[i2]);
    }

  private:
    const std::vector<float>&      data_;
};   

void
ConvertNrrdsToTexture::execute()
{

  TextureHandle tHandle_;
  bool update = false;

  if (card_mem_ != 0 && gui_card_mem_auto_.get())
    gui_card_mem_.set(card_mem_);
    
  else if (card_mem_ == 0)
    gui_card_mem_auto_.set(0);

  NrrdDataHandle vHandle;
  NrrdDataHandle gHandle;
  NrrdDataHandle mHandle;
  
  get_input_handle("Value Nrrd", vHandle,true);

  Nrrd* nv_nrrd = vHandle->nrrd_;

  if (nv_nrrd->dim != 3 && nv_nrrd->dim != 4) 
  {
    error("Invalid dimension for input value nrrd.");
    return;
  }

  size_t axis_size[4];
  nrrdAxisInfoGet_nva(nv_nrrd, nrrdAxisInfoSize, axis_size);
  if (nv_nrrd->dim == 4 && axis_size[0] != 1 && axis_size[0] != 4) 
  {
    error("Invalid axis size for Normal/Value nrrd.");
    return;
  }

  // The input nrrd type must be unsigned char.
  if (gui_uchar_.get() && vHandle->nrrd_->type != nrrdTypeUChar) 
  {
    error("Normal/Value input nrrd type must be unsigned char.");
    return;
  }

  if( !gui_fixed_.get() )
  {
    // set vmin/vmax
    NrrdRange *range = nrrdRangeNewSet(vHandle->nrrd_, nrrdBlind8BitRangeFalse);

    gui_vminval_.set(range->min);
    gui_vmaxval_.set(range->max);
    nrrdRangeNix(range);
  }

  // Check to see if the input nrrd has changed.
  if( vnrrd_last_generation_ != vHandle->generation  ||
      (gui_vminval_.get() != vminval_) || 
      (gui_vmaxval_.get() != vmaxval_) ||
      tHandle_.get_rep() == 0)
  {
    vnrrd_last_generation_ = vHandle->generation;

    vminval_ = gui_vminval_.get();
    vmaxval_ = gui_vmaxval_.get();

    update = true;
    if (tHandle_.get_rep() == 0) { tHandle_ = new Texture(); }
  }

  // The gradient nrrd input is optional.
  if (get_input_handle("Gradient Magnitude Nrrd", gHandle, false))
  {
    update_state(Executing);
    if (!ShaderProgramARB::shaders_supported()) 
    {
      // TODO: Runtime check, change message to reflect that.
      warning("This machine does not support advanced volume rendering. The gradient nrrd will be ignored.");
      if( gnrrd_last_generation_ != -1 )
      update = true;
      gHandle = 0;
      gnrrd_last_generation_ = -1;

    } 
    else 
    {
      Nrrd* gm_nrrd = gHandle->nrrd_;

      if (gm_nrrd->dim != 3 && gm_nrrd->dim != 4) 
      {
        error("Invalid dimension for input gradient magnitude nrrd.");
        return;
      }
      
      if( gm_nrrd->dim == 4 ) 
      {
        nrrdAxisInfoGet_nva(gm_nrrd, nrrdAxisInfoSize, axis_size);
        if (axis_size[0] != 1) 
        {
          error("Invalid axis size for gradient magnitude nrrd.");
          return;
        }
      }

      // The input nrrd type must be unsigned char.
      if (gui_uchar_.get() && gHandle->nrrd_->type != nrrdTypeUChar) 
      {
        error("Gradient magnitude input nrrd type must be unsigned char.");
        return;
      }

      if( !gui_fixed_.get())
      {
        // set gmin/gmax
        NrrdRange *range =
          nrrdRangeNewSet(gHandle->nrrd_, nrrdBlind8BitRangeFalse);

        gui_gminval_.set(range->min);
        gui_gmaxval_.set(range->max);
      }
	
      // Check to see if the input gradient nrrd has changed.
      if( gnrrd_last_generation_ != gHandle->generation  ||
          (gui_gminval_.get() != gminval_) || 
          (gui_gmaxval_.get() != gmaxval_) ) 
      {
        gnrrd_last_generation_ = gHandle->generation;
          
        gminval_ = gui_gminval_.get();
        gmaxval_ = gui_gmaxval_.get();
          
        update = true;
      }
    }
  } 
  else 
  {
    if( gnrrd_last_generation_ != -1 )
      update = true;

    gnrrd_last_generation_ = -1;
  }


  if( gui_uchar_.get() != is_uchar_ ) 
  {
    is_uchar_ = gui_uchar_.get();

    update = true;
  }

  update_state(Executing);

  if (update) 
  {   
    Nrrd *v = 0;
    Nrrd *g = 0;
    tHandle_ = new Texture;
    if (vHandle.get_rep()) v = vHandle->nrrd_;
    if (gHandle.get_rep()) g = gHandle->nrrd_;
    tHandle_->build(v, g, vminval_, vmaxval_, 
		    gminval_, gmaxval_, gui_card_mem_.get());
        
    tHandle_->value_nrrd_ = vHandle;
    tHandle_->gradient_magnitude_nrrd_ = gHandle;
    tHandle_->mask_nrrd_ = mHandle;
    send_output_handle("Texture", tHandle_, true);
  }

  // Generate a reasonable histogram
  if (gHandle.get_rep() && gui_histogram_.get())
  {
    NrrdDataHandle HistoGramNrrd = new NrrdData;

    // build joint histogram
    size_t sx = 256;
    size_t sy = 256;

    size_t sdim[2]; sdim[0] = sx; sdim[1] = sy;
    nrrdAlloc_nva(HistoGramNrrd->nrrd_ ,nrrdTypeUChar ,2 ,sdim );
    
    int centerdata[2]; centerdata[0] = nrrdCenterNode; centerdata[1] = nrrdCenterNode;
    nrrdAxisInfoSet_nva(HistoGramNrrd->nrrd_, nrrdAxisInfoCenter, centerdata);
  

    size_t sz = sx*sy;
    std::vector<float> hist(sz);
    std::vector<size_t> index(sz);
    
    for (size_t j=0; j<sz; j++) 
    {
      hist[j] = 0.0;
    }
    
    unsigned char* value = static_cast<unsigned char*>(vHandle->nrrd_->data);
    unsigned char* gradient = static_cast<unsigned char*>(gHandle->nrrd_->data);
    unsigned char* data = static_cast<unsigned char*>(HistoGramNrrd->nrrd_->data);

    size_t mul, offset;
    if (vHandle->nrrd_->dim == 4) 
    {
      offset = vHandle->nrrd_->axis[0].size-1;
      mul = vHandle->nrrd_->axis[0].size;
      sz = vHandle->nrrd_->axis[1].size*vHandle->nrrd_->axis[2].size*vHandle->nrrd_->axis[3].size; 
    }
    else
    {
      offset = 0;
      mul = 1;
      sz = vHandle->nrrd_->axis[0].size*vHandle->nrrd_->axis[1].size*vHandle->nrrd_->axis[2].size;       
    }
    
    for (size_t j=0; j<sz; j++)
    {
      unsigned char p = value[mul*j+offset];
      unsigned char q = gradient[j];
      hist[p+q*sx]++; 
    }

    double gamma = gui_gamma_.get();
    sz = sx*sy;
    std::vector<size_t> temp(sz);

    for (size_t j=0; j<sz; j++) 
    {
      temp[j] = j;
    }
    
    std::sort(temp.begin(),temp.end(),SortFloatVector(hist));

    float threshold = hist[temp[static_cast<size_t>(0.99*sz)]]*0.001;
    
    float oldv = 0.0;
    unsigned char old = 0;
    size_t j=0;
    for (; j<sz; j++) 
    {    
      if (hist[temp[j]] > threshold) break;
      data[temp[j]] = 0;
    }
    
    size_t jc = j;
    
    for (; j<sz; j++) 
    {
      if (hist[temp[j]] == oldv) 
      {
        data[temp[j]] = old;
      }
      else
      {
        oldv = hist[temp[j]];
        data[temp[j]] = static_cast<unsigned char>(255.0*pow((double)(j-jc)/(double)(sz-jc-1),gamma)); 
        old = data[temp[j]];
      }
    }

    send_output_handle("JointHistoGram",HistoGramNrrd,true);
  }
}

} // namespace SCIRun
