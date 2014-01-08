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


#include <slivr/VideoCardInfo.h>
#include <slivr/ShaderProgramARB.h>

#include <Core/Util/StringUtil.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/TexturePort.h>

#include <Dataflow/Modules/Visualization/share.h>

namespace SCIRun {

using namespace SLIVR;

class SCISHARE ConvertFieldsToTexture : public Module
{
  public:
    ConvertFieldsToTexture(GuiContext* ctx, const std::string& name = "ConvertFieldsToTexture",
                           SchedClass sc = Source, const std::string& cat = "Visualization", 
                           const std::string& pack = "SCIRun");
    virtual ~ConvertFieldsToTexture() {}
    virtual void execute();

    template<class T>
    void create_scaled_value_nrrd(T* value, unsigned char* nrrd, size_t sz, double fmin, double fmax );

    template<class T>
    void create_scaled_value_normal_nrrd(T* value, Vector* normal, unsigned char* nrrd, size_t sz, double fmin, double fmax );
    void create_scaled_gradient_magnitude_nrrd(Vector* gradient, unsigned char* nrrd, size_t sz, double fmin, double fmax );

  protected:

    GuiDouble gui_vminval_;
    GuiDouble gui_vmaxval_;
    GuiDouble gui_gminval_;
    GuiDouble gui_gmaxval_;

    GuiInt gui_fixed_;
    GuiInt gui_card_mem_;
    GuiInt gui_card_mem_auto_;

    GuiInt gui_histogram_;    
    GuiDouble gui_gamma_;
    
    int card_mem_;
};

DECLARE_MAKER(ConvertFieldsToTexture)
  
ConvertFieldsToTexture::ConvertFieldsToTexture(GuiContext* ctx, const std::string& name,
                                               SchedClass sc, const std::string& cat, const std::string& pack) :
  Module(name, ctx, sc, cat, pack),
  gui_vminval_(get_ctx()->subVar("vmin"), 0),
  gui_vmaxval_(get_ctx()->subVar("vmax"), 1),
  gui_gminval_(get_ctx()->subVar("gmin"), 0),
  gui_gmaxval_(get_ctx()->subVar("gmax"), 1),
  gui_fixed_(get_ctx()->subVar("is_fixed"), 0),
  gui_card_mem_(get_ctx()->subVar("card_mem"), 64),
  gui_card_mem_auto_(get_ctx()->subVar("card_mem_auto"), 1),
  gui_histogram_(get_ctx()->subVar("histogram"),1),
  gui_gamma_(get_ctx()->subVar("gamma"),0.5),
  card_mem_(video_card_memory_size())
{}

template<class T>
void
ConvertFieldsToTexture::create_scaled_value_nrrd(T* value, unsigned char* nrrd, 
            size_t sz, double fmin, double fmax )
{
  T minval = static_cast<T>(fmin); 
  T maxval = static_cast<T>(fmax);

  double factor =  1.0;
  if (maxval > minval) factor = 1.0/(maxval-minval);

  for(size_t i=0; i<sz; i++)
  {
    T val = value[i];
    if (val < minval) val = minval;
    else if (val > maxval) val = maxval;

    nrrd[i] = static_cast<unsigned char>((255*(value[i]-minval))*factor);
    
  }
}

template<class T>
void
ConvertFieldsToTexture::create_scaled_value_normal_nrrd(T* value, Vector* normal, unsigned char* nrrd, 
            size_t sz, double fmin, double fmax )
{
  T minval = static_cast<T>(fmin); 
  T maxval = static_cast<T>(fmax);

  double factor =  1.0;
  if (maxval > minval) factor = 1.0/(maxval-minval);
    
  for(size_t i=0; i<sz; i++)
  {
    T val = value[i];
    if (val < minval) val = minval;
    else if (val > maxval) val = maxval;
    Vector norm = normal[i].normal();

    nrrd[4*i+3] = static_cast<unsigned char>((255*(value[i]-minval))*factor);
    nrrd[i*4] = static_cast<unsigned char>(norm.x()*127.5 + 127.5);
    nrrd[i*4+1] = static_cast<unsigned char>(norm.y()*127.5 + 127.5);
    nrrd[i*4+2] = static_cast<unsigned char>(norm.z()*127.5 + 127.5);
  }
}

void
ConvertFieldsToTexture::create_scaled_gradient_magnitude_nrrd(Vector* gradient, unsigned char* nrrd, 
            size_t sz, double fmin, double fmax )
{
  double factor =  1.0;
  if (fmax > fmin) factor = 1.0/(fmax-fmin);

  for(size_t i=0; i<sz; i++)
  {
    double len = gradient[i].length();
    if (len < fmin) len = fmin;
    else if (len > fmax) len = fmax;

    nrrd[i] = static_cast<unsigned char>((255*(len-fmin))*factor);
  }
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
ConvertFieldsToTexture::execute()
{
  // Fields        
  FieldHandle ValueHandle;
  FieldHandle GradientHandle;
  FieldHandle MaskHandle;
    
  get_input_handle("Value Field",ValueHandle,true);
  get_input_handle("Gradient Magnitude Field",GradientHandle,false);
  
  if (inputs_changed_ ||gui_card_mem_.changed() || gui_card_mem_.changed()  ||
      gui_vminval_.changed() || gui_vmaxval_.changed() || gui_gminval_.changed() ||
      gui_gmaxval_.changed() || gui_fixed_.changed() || 
      gui_histogram_.changed() || !oport_cached("Texture"))
  {
    // Check the memory of the graphics card
    update_state(Executing);
    
    if (card_mem_ != 0 && gui_card_mem_auto_.get()) gui_card_mem_.set(card_mem_);
    else if (card_mem_ == 0) gui_card_mem_auto_.set(0);

    bool has_gradient = false;
    bool has_gradient_magnitude = false;
    
    VField* vfield = ValueHandle->vfield();
    VField* gfield = 0;
    
    if (GradientHandle.get_rep())
    {
      gfield = GradientHandle->vfield();
      if (gfield->is_vector()) has_gradient = true;
      if (gfield->is_scalar()) has_gradient_magnitude = true;
      
      if ((has_gradient == false)&&(has_gradient_magnitude == false)) 
      {
        error("Gradient does not contain scalar or vector data");
        return;
      }
    }

    VMesh* vmesh = ValueHandle->vmesh();
    VMesh::dimension_type vdim;
    vmesh->get_dimensions(vdim);
    
    if (!vmesh->is_latvolmesh())
    { 
      error("Value needs to be a LatVolField");
      return;
    }

    if (!vfield->is_scalar())
    {
      error("Value needs to be a scalar value");
      return;
    }

    VMesh* gmesh = 0;
    VMesh::dimension_type gdim;
   
    double vmin, vmax;  
    double gmin, gmax;  
     
       
    if (has_gradient || has_gradient_magnitude)
    {
      gmesh = GradientHandle->vmesh();
      
      if (!gmesh->is_latvolmesh())
      { 
        error("Gradient needs to be a LatVolField");
        return;
      }    
      
      gmesh->get_dimensions(gdim);
      
      bool samedim = true;
      for (size_t i = 0; i < gdim.size(); i++) if (gdim[i] != vdim[i]) samedim = false;
      
      if (!samedim)
      {
        error("Value and Gradient do not have the same dimension");
        return;
      }
      
      if (gfield->basis_order() == -1)
      {
        error("Gradient does not contain data");
        return;        
      }
      
      if (vfield->basis_order() != gfield->basis_order())
      {
        error("Value and Gradient do not have the same basis order");
        return;    
      }
      
      if( gui_fixed_.get() )
      {
        gmin = gui_gminval_.get();
        gmax = gui_gmaxval_.get();
      }
      else
      {
        gfield->min(gmin);
        gfield->max(gmax);
        gui_gminval_.set(gmin);
        gui_gmaxval_.set(gmax);
      }
    }

    if (vfield->basis_order() == -1)
    {
      error("No data contained in the Value field");
      return;        
    }

    NrrdDataHandle ValueNrrd;
    NrrdDataHandle GradientNrrd;
    
    VMesh::dimension_type dim;
    if (vfield->basis_order() == 1)
    {
      vmesh->get_dimensions(dim);
    }
    else
    {
      vmesh->get_elem_dimensions(dim);
    }

    ValueNrrd = new NrrdData;
    

    if( gui_fixed_.get() )
    {
      vmin = gui_vminval_.get();
      vmax = gui_vmaxval_.get();
    }
    else
    {
      vfield->min(vmin);
      vfield->max(vmax);
      gui_vminval_.set(vmin);
      gui_vmaxval_.set(vmax);      
    }

    if (has_gradient || has_gradient_magnitude)
    {
      GradientNrrd = new NrrdData;
    }
    
    if (!has_gradient)
    {
      size_t sdim[3]; sdim[0] = dim[0]; sdim[1] = dim[1]; sdim[2] = dim[2];
      nrrdAlloc_nva(ValueNrrd->nrrd_ ,nrrdTypeUChar ,3 ,sdim );
      
      int nrrdcenter = nrrdCenterCell;
      if(vfield->basis_order() == 1) nrrdcenter = nrrdCenterNode;
      
      int centerdata[3];
      for (int j=0; j<3; j++) centerdata[j] = nrrdcenter;
      nrrdAxisInfoSet_nva(ValueNrrd->nrrd_, nrrdAxisInfoCenter, centerdata);
      
      ValueNrrd->nrrd_->spaceDim = 3;
      ValueNrrd->nrrd_->space = nrrdSpace3DRightHanded;

      Transform tf; 
      tf = vmesh->get_transform();  
      double Trans[16];
      tf.get(Trans);
      
      for (int p=0;p<3;p++)
      {
        ValueNrrd->nrrd_->spaceOrigin[p] = Trans[3+4*p];
        for (int q=0;q<3;q++)
          ValueNrrd->nrrd_->axis[q].spaceDirection[p] = Trans[q+4*p];
      }

      unsigned char* data = reinterpret_cast<unsigned char *>(ValueNrrd->nrrd_->data);
      
      if (vfield->is_char()) 
      {
        char* ptr = static_cast<char*>(vfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_unsigned_char()) 
      {
        unsigned char* ptr = static_cast<unsigned char*>(vfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_short()) 
      {
        short* ptr = static_cast<short*>(vfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_unsigned_short()) 
      {
        unsigned short* ptr = static_cast<unsigned short*>(vfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_int()) 
      {
        int* ptr = static_cast<int*>(vfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_unsigned_int()) 
      {
        unsigned int* ptr = static_cast<unsigned int*>(vfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_longlong()) 
      {
        long long* ptr = static_cast<long long*>(vfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_unsigned_longlong()) 
      {
        unsigned long long* ptr = static_cast<unsigned long long*>(vfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_float()) 
      {
        float* ptr = static_cast<float*>(vfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_double()) 
      {
        double* ptr = static_cast<double*>(vfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,vfield->num_values(),vmin,vmax);       
      }
    }
    else
    {
      size_t sdim[4]; sdim[0] = 4; sdim[1] = dim[0]; sdim[2] = dim[1]; sdim[3] = dim[2];
      nrrdAlloc_nva(ValueNrrd->nrrd_ ,nrrdTypeUChar, 4 ,sdim);
      
      int nrrdcenter = nrrdCenterCell;
      if(vfield->basis_order() == 1) nrrdcenter = nrrdCenterNode;
      
      int centerdata[4];
      for (int j=0; j<4; j++) centerdata[j] = nrrdcenter;
      nrrdAxisInfoSet_nva(ValueNrrd->nrrd_, nrrdAxisInfoCenter, centerdata);
      
      ValueNrrd->nrrd_->spaceDim = 3;
      ValueNrrd->nrrd_->space = nrrdSpace3DRightHanded;

      Transform tf; 
      tf = vmesh->get_transform();  
      double Trans[16];
      tf.get(Trans);

      for (int p=0;p<3;p++)
      {
        ValueNrrd->nrrd_->spaceOrigin[p] = Trans[3+4*p];
        for (int q=0;q<3;q++)
          ValueNrrd->nrrd_->axis[q+1].spaceDirection[p] = Trans[q+4*p];
      }
    
      unsigned char* data = reinterpret_cast<unsigned char *>(ValueNrrd->nrrd_->data);
      Vector* norm = static_cast<Vector *>(gfield->fdata_pointer());
      
      if (vfield->is_char()) 
      {
        char* ptr = static_cast<char*>(vfield->fdata_pointer());
        create_scaled_value_normal_nrrd(ptr,norm,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_unsigned_char()) 
      {
        unsigned char* ptr = static_cast<unsigned char*>(vfield->fdata_pointer());
        create_scaled_value_normal_nrrd(ptr,norm,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_short()) 
      {
        short* ptr = static_cast<short*>(vfield->fdata_pointer());
        create_scaled_value_normal_nrrd(ptr,norm,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_unsigned_short()) 
      {
        unsigned short* ptr = static_cast<unsigned short*>(vfield->fdata_pointer());
        create_scaled_value_normal_nrrd(ptr,norm,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_int()) 
      {
        int* ptr = static_cast<int*>(vfield->fdata_pointer());
        create_scaled_value_normal_nrrd(ptr,norm,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_unsigned_int()) 
      {
        unsigned int* ptr = static_cast<unsigned int*>(vfield->fdata_pointer());
        create_scaled_value_normal_nrrd(ptr,norm,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_longlong()) 
      {
        long long* ptr = static_cast<long long*>(vfield->fdata_pointer());
        create_scaled_value_normal_nrrd(ptr,norm,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_unsigned_longlong()) 
      {
        unsigned long long* ptr = static_cast<unsigned long long*>(vfield->fdata_pointer());
        create_scaled_value_normal_nrrd(ptr,norm,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_float()) 
      {
        float* ptr = static_cast<float*>(vfield->fdata_pointer());
        create_scaled_value_normal_nrrd(ptr,norm,data,vfield->num_values(),vmin,vmax);       
      }
      if (vfield->is_double()) 
      {
        double* ptr = static_cast<double*>(vfield->fdata_pointer());
        create_scaled_value_normal_nrrd(ptr,norm,data,vfield->num_values(),vmin,vmax);       
      }
    }
    
    
    if (has_gradient)
    {
      size_t sdim[3]; sdim[0] = dim[0]; sdim[1] = dim[1]; sdim[2] = dim[2];
      nrrdAlloc_nva(GradientNrrd->nrrd_ ,nrrdTypeUChar  ,3,sdim );
      
      int nrrdcenter = nrrdCenterCell;
      if(vfield->basis_order() == 1) nrrdcenter = nrrdCenterNode;
      
      int centerdata[3];
      for (int j=0; j<3; j++) centerdata[j] = nrrdcenter;
      nrrdAxisInfoSet_nva(GradientNrrd->nrrd_, nrrdAxisInfoCenter, centerdata);
      
      ValueNrrd->nrrd_->spaceDim = 3;
      ValueNrrd->nrrd_->space = nrrdSpace3DRightHanded;

      Transform tf; 
      tf = vmesh->get_transform();  
      double Trans[16];
      tf.get(Trans);
      
      for (int p=0;p<3;p++)
      {
        GradientNrrd->nrrd_->spaceOrigin[p] = Trans[3+4*p];
        for (int q=0;q<3;q++)
          GradientNrrd->nrrd_->axis[q].spaceDirection[p] = Trans[q+4*p];
      }

      unsigned char* data = reinterpret_cast<unsigned char *>(GradientNrrd->nrrd_->data);

      Vector* ptr = static_cast<Vector*>(gfield->fdata_pointer());
      create_scaled_gradient_magnitude_nrrd(ptr,data,gfield->num_values(),gmin,gmax);       
    }
    else if (has_gradient_magnitude)
    {
      size_t sdim[3]; sdim[0] = dim[0]; sdim[1] = dim[1]; sdim[2] = dim[2];
      nrrdAlloc_nva(GradientNrrd->nrrd_ ,nrrdTypeUChar ,3 ,sdim );
      
      int nrrdcenter = nrrdCenterCell;
      if(vfield->basis_order() == 1) nrrdcenter = nrrdCenterNode;
      
      int centerdata[3];
      for (int j=0; j<3; j++) centerdata[j] = nrrdcenter;
      nrrdAxisInfoSet_nva(GradientNrrd->nrrd_, nrrdAxisInfoCenter, centerdata);
      
      ValueNrrd->nrrd_->spaceDim = 3;
      ValueNrrd->nrrd_->space = nrrdSpace3DRightHanded;

      Transform tf; 
      tf = vmesh->get_transform();  
      double Trans[16];
      tf.get(Trans);
      
      for (int p=0;p<3;p++)
      {
        GradientNrrd->nrrd_->spaceOrigin[p] = Trans[3+4*p];
        for (int q=0;q<3;q++)
          GradientNrrd->nrrd_->axis[q].spaceDirection[p] = Trans[q+4*p];
      }

      unsigned char* data = reinterpret_cast<unsigned char *>(GradientNrrd->nrrd_->data);

      if (gfield->is_char()) 
      {
        char* ptr = static_cast<char*>(gfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,gfield->num_values(),vmin,vmax);       
      }
      if (gfield->is_unsigned_char()) 
      {
        unsigned char* ptr = static_cast<unsigned char*>(gfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,gfield->num_values(),vmin,vmax);       
      }
      if (gfield->is_short()) 
      {
        short* ptr = static_cast<short*>(gfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,gfield->num_values(),vmin,vmax);       
      }
      if (gfield->is_unsigned_short()) 
      {
        unsigned short* ptr = static_cast<unsigned short*>(gfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,gfield->num_values(),vmin,vmax);       
      }
      if (gfield->is_int()) 
      {
        int* ptr = static_cast<int*>(gfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,gfield->num_values(),vmin,vmax);       
      }
      if (gfield->is_unsigned_int()) 
      {
        unsigned int* ptr = static_cast<unsigned int*>(gfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,gfield->num_values(),vmin,vmax);       
      }
      if (gfield->is_longlong()) 
      {
        long long* ptr = static_cast<long long*>(gfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,gfield->num_values(),vmin,vmax);       
      }
      if (gfield->is_unsigned_longlong()) 
      {
        unsigned long long* ptr = static_cast<unsigned long long*>(gfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,gfield->num_values(),vmin,vmax);       
      }
      if (gfield->is_float()) 
      {
        float* ptr = static_cast<float*>(gfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,gfield->num_values(),vmin,vmax);       
      }
      if (gfield->is_double()) 
      {
        double* ptr = static_cast<double*>(gfield->fdata_pointer());
        create_scaled_value_nrrd(ptr,data,gfield->num_values(),vmin,vmax);       
      }   
    }
    
    TextureHandle TextureH = new Texture();
    Nrrd* v = 0;
    Nrrd* g = 0;
    if (ValueNrrd.get_rep()) v = ValueNrrd->nrrd_;
    if (GradientNrrd.get_rep()) g = GradientNrrd->nrrd_;
    
    TextureH->build(v, g, vmin, vmax, 
          gmin, gmax, gui_card_mem_.get());

    TextureH->value_nrrd_ = ValueNrrd;
    TextureH->gradient_magnitude_nrrd_ = GradientNrrd;

    send_output_handle("Texture", TextureH, true);
    
    if (GradientNrrd.get_rep() && gui_histogram_.get())
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
      
      unsigned char* value = static_cast<unsigned char*>(ValueNrrd->nrrd_->data);
      unsigned char* gradient = static_cast<unsigned char*>(GradientNrrd->nrrd_->data);
      unsigned char* data = static_cast<unsigned char*>(HistoGramNrrd->nrrd_->data);

      size_t mul, offset;
      if (ValueNrrd->nrrd_->dim == 4) 
      {
        offset = ValueNrrd->nrrd_->axis[0].size-1;
        mul = ValueNrrd->nrrd_->axis[0].size;
        sz = ValueNrrd->nrrd_->axis[1].size*ValueNrrd->nrrd_->axis[2].size*ValueNrrd->nrrd_->axis[3].size; 
      }
      else
      {
        offset = 0;
        mul = 1;
        sz = ValueNrrd->nrrd_->axis[0].size*ValueNrrd->nrrd_->axis[1].size*ValueNrrd->nrrd_->axis[2].size;       
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
}

} // namespace SCIRun
