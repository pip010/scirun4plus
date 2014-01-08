/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/


#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/ITKDatatypePort.h>
#include <Dataflow/Network/Ports/FieldPort.h>

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkRegion.h>
#include <itkVector.h>

namespace Insight {

using namespace SCIRun;

enum FieldType {LATVOLFIELD, IMAGEFIELD };

class FieldToImage : public Module {

public:
  FieldHandle infield_handle_;

  ITKDatatypeHandle outimage_handle_;
  ITKDatatype* img_;

  FieldToImage(GuiContext*);

  virtual ~FieldToImage() {}

  virtual void execute();

  // Run function will dynamically cast data to determine which
  // instantiation we are working with. The last template type
  // refers to the last template type of the filter intstantiation.
  template< class data > 
  bool run( const FieldHandle &fh );

  bool run_vector( const FieldHandle &fh );
};


DECLARE_MAKER(FieldToImage)
FieldToImage::FieldToImage(GuiContext* ctx)
  : Module("FieldToImage", ctx, Source, "Converters", "Insight")
{
  img_ = new ITKDatatype;
}


template<class data >
bool
FieldToImage::run( const FieldHandle& fh) 
{
  FieldInformation fi(fh);

  if(fi.is_latvolmesh()) 
  {

    typedef itk::Image<data, 3> ImageType;
    VField* vfield = fh->vfield();
    VMesh*  vmesh  = fh->vmesh();

    // create a new itk image
    typename ImageType::Pointer img = ImageType::New(); 

    VField::dimension_type dims;
    vfield->get_values_dimensions(dims);
    
    // set size
    typename ImageType::SizeType fixedSize = {{dims[0], dims[1], dims[2]}};
    img->SetRegions( fixedSize );

    // set origin and spacing
    const BBox bbox = vmesh->get_bounding_box();
    Point mesh_center;
    Vector mesh_size;
    if(bbox.valid()) 
    {
      mesh_center = bbox.center();
      mesh_size = bbox.diagonal();
    }
    else 
    {
      error("No bounding box to get center");
      return false;
    }
    
    double origin[ ImageType::ImageDimension ];
    Point bbox_min = bbox.min();
    origin[0] = bbox_min.x();
    origin[1] = bbox_min.y();
    origin[2] = bbox_min.z();
    
    img->SetOrigin( origin );
    
    double spacing[ ImageType::ImageDimension ];
    spacing[0] = mesh_size.x()/dims[0];
    spacing[1] = mesh_size.y()/dims[1];
    spacing[2] = mesh_size.z()/dims[2];
    
    img->SetSpacing( spacing );

    // set new data container
    data *imageData = vfield->get_values_pointer();
    unsigned long size = (unsigned long)vfield->num_values();

    img->GetPixelContainer()->SetImportPointer(imageData, size, false);

    // send the data downstream
    img_->data_ = img;
    outimage_handle_ = img_;

  }
  else if(fi.is_imagemesh()) 
  {
    typedef itk::Image<data, 2> ImageType;
    VField* vfield = fh->vfield();
    VMesh*  vmesh  = fh->vmesh();

    // create a new itk image
    typename ImageType::Pointer img = ImageType::New(); 

    VField::dimension_type dims;
    vfield->get_values_dimensions(dims);

    // set size
    typename ImageType::SizeType fixedSize = {{dims[0], dims[1]}};
    img->SetRegions( fixedSize );

    // set origin and spacing
    const BBox bbox = vmesh->get_bounding_box();
    Point mesh_center;
    Vector mesh_size;
    if(bbox.valid()) 
    {
      mesh_center = bbox.center();
      mesh_size = bbox.diagonal();
    }
    else 
    {
      error("No bounding box to get center");
      return false;
    }
    
    double origin[ ImageType::ImageDimension ];
    Point bbox_min = bbox.min();
    origin[0] = bbox_min.x();
    origin[1] = bbox_min.y();
    
    img->SetOrigin( origin );
    
    double spacing[ ImageType::ImageDimension ];
    spacing[0] = mesh_size.x()/dims[0];
    spacing[1] = mesh_size.y()/dims[1];
    
    img->SetSpacing( spacing );

    // set new data container
    data* imageData = vfield->get_values_pointer();
    unsigned long size = (unsigned long)vfield->num_values();

    img->GetPixelContainer()->SetImportPointer(imageData, size, false);

    // send the data downstream
    img_->data_ = img;
    outimage_handle_ = img_;
  }

  return true;
}

bool FieldToImage::run_vector( const FieldHandle &fh) 
{
  FieldInformation fi(fh);

  if(fi.is_latvolmesh()) 
  {
    typedef itk::Image<itk::Vector<double>, 3> ImageType;
    typedef itk::ImageRegionIterator< ImageType > IteratorType;

    VField* f = fh->vfield();
    // create a new itk image
    ImageType::Pointer img = ImageType::New(); 

    // set size
    VField:dimension_type dims;
    f->get_values_dimension(dims);
    
    ImageType::SizeType fixedSize = {{dims[0], dims[1], dims[2] }};
    img->SetRegions( fixedSize );

    ImageType::RegionType region;
    
    ImageType::IndexType start;
    
    for(int i=0; i<3; i++) start[i] = 0;

    region.SetSize( fixedSize );
    region.SetIndex( start );
    
    img->SetRegions( region );
    img->Allocate();
    
    // set origin and spacing
    const BBox bbox = f->vmesh()->get_bounding_box();
    Point mesh_center;
    Vector mesh_size;
    if(bbox.valid()) 
    {
      mesh_center = bbox.center();
      mesh_size = bbox.diagonal();
    }
    else 
    {
      error("No bounding box to get center");
      return false;
    }
    
    double origin[ ImageType::ImageDimension ];
    Point bbox_min = bbox.min();
    origin[0] = bbox_min.x();
    origin[1] = bbox_min.y();
    origin[2] = bbox_min.z();
    
    img->SetOrigin( origin );
    
    double spacing[ ImageType::ImageDimension ];
    spacing[0] = mesh_size.x()/f->fdata().dim3();
    spacing[1] = mesh_size.y()/f->fdata().dim2();
    spacing[2] = mesh_size.z()/f->fdata().dim1();
    
    img->SetSpacing( spacing );

    // copy the data
    VMesh::Node::iterator iter, end;
    VMesh* mh = f->vmesh(); 
    
    //LVF::handle_type mh((LVMesh*)(f->mesh().get_rep()));
    mh->begin(iter);
    mh->end(end);

    typedef ImageType::PixelType PixelType;
    PixelType pixel;

    IteratorType img_iter(img, img->GetRequestedRegion());
    img_iter.GoToBegin();

    while(iter != end ) 
    {
      Vector val;
      if (f->get_value(val, *iter)) 
      {	  
        itk::Vector<double> new_val;
        new_val[0] = val[0];
        new_val[1] = val[1];
        new_val[2] = val[2];
        
        img_iter.Set(new_val);
      } 
      ++iter;
      img_iter.operator++();
    }

    // send the data downstream
    img_->data_ = img;
    outimage_handle_ = img_;

  }
  else if(fi.is_imagemesh()) 
  {
    typedef itk::Image< itk::Vector<double>, 2> ImageType;
    typedef itk::ImageRegionIterator< ImageType > IteratorType;
    VField* f = fh->vfield();

    // create a new itk image
    ImageType::Pointer img = ImageType::New(); 

    VField:dimension_type dims;
    f->get_values_dimension(dims);
    // set size
    ImageType::SizeType fixedSize = {{dims[0], dims[1]}};
    img->SetRegions( fixedSize );

    ImageType::RegionType region;
    
    ImageType::IndexType start;
    
    for(int i=0; i<3; i++)
      start[i] = 0;

    region.SetSize( fixedSize );
    region.SetIndex( start );
    
    img->SetRegions( region );
    img->Allocate();

    // set origin and spacing
    const BBox bbox = f->vmesh()->get_bounding_box();
    Point mesh_center;
    Vector mesh_size;
    if(bbox.valid()) 
    {
      mesh_center = bbox.center();
      mesh_size = bbox.diagonal();
    }
    else 
    {
      error("No bounding box to get center");
      return false;
    }
    
    double origin[ ImageType::ImageDimension ];
    Point bbox_min = bbox.min();
    origin[0] = bbox_min.x();
    origin[1] = bbox_min.y();
    
    img->SetOrigin( origin );
    
    double spacing[ ImageType::ImageDimension ];
    spacing[0] = mesh_size.x()/f->fdata().dim2();
    spacing[1] = mesh_size.y()/f->fdata().dim1();
    
    img->SetSpacing( spacing );

     // copy the data
    VMesh::Node::iterator iter, end;
    VMesh* mh= f->vmesh();
    mh->begin(iter);
    mh->end(end);

    typedef ImageType::PixelType PixelType;
    PixelType pixel;

    IteratorType img_iter(img, img->GetRequestedRegion());
    img_iter.GoToBegin();

    while(iter != end ) 
    {
      Vector val;
      if (f->get_value(val, *iter)) 
      {	  
        itk::Vector<double> new_val;
        new_val[0] = val[0];
        new_val[1] = val[1];
        new_val[2] = val[2];
        
        img_iter.Set(new_val);
      } 
      ++iter;
      img_iter.operator++();
    }   
    // send the data downstream
    img_->data_ = img;
    outimage_handle_ = img_;
  }

  return true;
}


void
FieldToImage::execute()
{
  if (!get_input_handle("InputField", infield_handle_)) return;

  // Determine which type of field we are convertion from.
  // Our only options are ImageField, ITKImageField,
  // LatVolField and ITKLatVolField
  VField* vfield = infield_handle_->vfield();
  
  if(vfield->is_double()) { run<double>(infield_handle_); }
  else if(vfield->is_float())) { run<float>(infield_handle_); }
  else if(vfield->is_unsigned_char())) { run<unsigned char>(infield_handle_); }
  else if(vfield->is_char())) { run<char>(infield_handle_); }
  else if(vfield->is_unsigned_short())) { run<unsigned short>(infield_handle_); }
  else if(vfield->is_short())) { run<short>(infield_handle_); }
  else if(vfield->is_unsigned_int())) { run<unsigned int>(infield_handle_); }
  else if(vfield->is_int())) { run<int>(infield_handle_); }  
  else if(run_vector(infield_handle_)) {}
  else {
    error("Unknown type");
    return;
  }

  send_output_handle("OutputImage", outimage_handle_, true);
}


} // End namespace Insight


