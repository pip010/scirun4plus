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
//    File   : ITKSpeedToPathGradientDescentImageFilterTool.cc
//    Author : David Brayford
//    Date   : May 2008

#include <Applications/Seg3D/ITKSpeedToPathGradientDescentImageFilterTool.h>
#include <Applications/Seg3D/GuiCode/itkspeedtopathgradientdescentfilter.h>
#include <Applications/Seg3D/Seg3DFrame.h>
#include <Applications/Seg3D/BrushTool.h>
#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/VolumeOps.h>
#include <Core/Util/Assert.h>

#include <Applications/Seg3D/Seg3DwxGuiUtils.h>
#include <wx/string.h>
#include <wx/variant.h>

#include <itkImageToImageFilter.h>
#include <itkCommand.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkGradientImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkNormalizeImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkNumericTraits.h>
#include <itkTimeProbe.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkPolyLineParametricPath.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkArrivalFunctionToPathFilter.h>
#ifndef _WIN32
#include <itkSpeedFunctionToPathFilter.h>
#endif
#include <itkPathIterator.h>
#include <itkGradientDescentOptimizer.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkTubeSpatialObject.h>
#include <itkTubeSpatialObjectPoint.h>
#include <itkSpatialObjectPoint.h>
#include <itkSpatialObjectWriter.h>
#include <itkBloxBoundaryPointImage.h>

#ifdef ITK_MANUAL_INSTANTIATION
#ifndef _WIN32
#include <itkSpeedFunctionToPathFilter.txx>
#endif
#include <itkArrivalFunctionToPathFilter.txx>
#include <itkImageToPathFilter.txx>
#include <itkIterateNeighborhoodOptimizer.txx>
#include <itkPhysicalCentralDifferenceImageFunction.txx>
#include <itkSingleImageCostFunction.txx>
#include <itkSpeedFunctionToPathFilter.txx>
#endif // ITK_MANUAL_INSTANTIATION

#ifdef _WIN32
#include <itkSpeedFunctionToPathFilter.txx>
#endif
namespace SCIRun {

const static unsigned int Dimension = 2;

typedef unsigned char OutputPixelType;
typedef itk::Image< float, Dimension> ImageType;
typedef itk::Image< OutputPixelType, Dimension> OutputImageType;
typedef itk::PolyLineParametricPath< Dimension > PathType;
typedef itk::SpeedFunctionToPathFilter< ITKImageFloat2D, PathType > PathFilterType;
typedef PathFilterType::CostFunctionType::CoordRepType CoordRepType;
typedef itk::LinearInterpolateImageFunction<ITKImageFloat2D, CoordRepType> InterpolatorType;
typedef itk::GradientDescentOptimizer OptimizerType;
typedef itk::PathIterator< OutputImageType, PathType > PathIteratorType;


ITKSpeedToPathGradientDescentImageFilterTool::ITKSpeedToPathGradientDescentImageFilterTool(Painter *painter) :
  SpeedToPathTool("ITKSpeedToPathGradientDescentImageFilterTool", painter),
  NumberIterations_(painter->get_vars(), "ITKSpeedToPathGradientDescentImageFilterTool::mIterations"),
  TerminationValue_(painter->get_vars(), "ITKSpeedToPathGradientDescentImageFilterTool::mTermination")
{
}


ITKSpeedToPathGradientDescentImageFilterTool::~ITKSpeedToPathGradientDescentImageFilterTool()
{
}


BaseTool::propagation_state_e
ITKSpeedToPathGradientDescentImageFilterTool::pointer_down(int b, int x, int y,
                                                           unsigned int m,
                                                           int t)
{
  if (!painter_->get_vars()->get_bool("ITKSpeedToPathGradientDescentImageFilterTool::mSeedPoints"))
  {
    return CONTINUE_E;
  }

  return PolylineTool::pointer_down(b, x, y, m, t);
}


void 
ITKSpeedToPathGradientDescentImageFilterTool::compute_path(size_t sindex)
{
  const size_t sindex0 = sindex;
  const size_t sindex1 = (sindex+1)%seeds_.size();

  paths_[sindex].clear();

  vector<int> iseed0 = source_volume_->world_to_index(seeds_[sindex0]);
  vector<int> iseed1 = source_volume_->world_to_index(seeds_[sindex1]);
  
  // Out of bounds, push straight line.
  if (!(source_volume_->index_valid(iseed0) &&
        source_volume_->index_valid(iseed1)))
  {
    paths_[sindex].push_back(seeds_[sindex]);
    return;
  }

  // Get speed function slice.
  SliceWindow *window = last_seed_window_;
  SLIVR::Plane plane(SLIVR::Point(window->center_.x(), window->center_.y(),
    window->center_.z()), SLIVR::Vector(window->normal_.x(),
    window->normal_.y(),window->normal_.z()));
  VolumeSliceHandle image_slice;
  image_slice = source_volume_->get_volume_slice(plane);

  ITKDatatypeHandle image_handle =
    nrrd_to_itk_image( image_slice->nrrd_handle_ ); 

  // Replace InputImageType with actual image type non typedef
  ImageType::Pointer pImage =
    dynamic_cast<ImageType*>(image_handle->data_.GetPointer());
    
  if (!pImage) 
  {
    return;
  } 

  InterpolatorType::Pointer interp = InterpolatorType::New();
	
  PathFilterType::CostFunctionType::Pointer cost =
    PathFilterType::CostFunctionType::New();
  cost->SetInterpolator( interp );

  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetNumberOfIterations( NumberIterations_ );
	
  PathFilterType::Pointer filter = PathFilterType::New();
  filter->SetInput( pImage );
  filter->SetCostFunction( cost );
  filter->SetOptimizer( optimizer );
  filter->SetTerminationValue( TerminationValue_ );

  // Add the seed points.
  PathFilterType::PathInfo info;
  PathFilterType::PointType pnt;
  const int axis = image_slice->get_axis();

  // Add starting point.
  if (axis == 1)
  {
    pnt[0] = seeds_[sindex0].y();
    pnt[1] = seeds_[sindex0].z();
  }
  else if (axis == 2)
  {
    pnt[0] = seeds_[sindex0].x();
    pnt[1] = seeds_[sindex0].z();
  }
  else
  {
    pnt[0] = seeds_[sindex0].x();
    pnt[1] = seeds_[sindex0].y();
  }
  info.SetEndPoint( pnt );

  // Add ending point.
  if (axis == 1)
  {
    pnt[0] = seeds_[sindex1].y();
    pnt[1] = seeds_[sindex1].z();
  }
  else if (axis == 2)
  {
    pnt[0] = seeds_[sindex1].x();
    pnt[1] = seeds_[sindex1].z();
  }
  else
  {
    pnt[0] = seeds_[sindex1].x();
    pnt[1] = seeds_[sindex1].y();
  }
  info.SetStartPoint( pnt );
  
  filter->AddPathInfo(info);

  filter->Update();

  // Use straight line if there is no computable path.
  paths_[sindex].push_back(seeds_[sindex]);

  for (size_t i = 0; i < filter->GetNumberOfOutputs(); i++)
  {
    PathType::Pointer path = filter->GetOutput(i);

    for (unsigned int j = 0; j < path->GetVertexList()->Size(); j++)
    {
      const float x = path->GetVertexList()->ElementAt(j)[0];
      const float y = path->GetVertexList()->ElementAt(j)[1];
      const float slice_index = (float)iseed0[axis];

      vector<double> ipnt(4);
      ipnt[0] = 1;
      ipnt[axis] = slice_index;
      if (axis == 1)
      {
        ipnt[2] = x;
        ipnt[3] = y;
      }
      else if (axis == 2)
      {
        ipnt[1] = x;
        ipnt[3] = y;
      }
      else
      {
        ipnt[1] = x;
        ipnt[2] = y;
      }
      paths_[sindex].push_back(source_volume_->index_to_point(ipnt));
    }
  }
}


}
