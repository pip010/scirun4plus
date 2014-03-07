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

#include "itkWin32Header.h"
#include <iostream>
#include <fstream>
#include "itkNumericTraits.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDifferenceImageFilter.h"



#define TEST_DIMENSION_MAX 2

int RegressionTestImage (const char *, const char *, double, unsigned int, unsigned int, int, bool);

int main(int argc, char **argv)
{
 if(argc < 5)
   {
   cerr << "Usage:" << endl;
   cerr << "testRegressionImage, intensityTolerance (good default 2.0), radiusTolerance, numberOfPixelsTolerance, testImage, baselineImage1, [baselineImage2, baselineImage3, ...]" << endl;
   cerr << "Note that if you supply more than one baselineImage, this test will pass if any" << endl;
   cerr << "of them match the testImage" << endl;
   return -1;
   }
 int bestBaselineStatus = 2001;
 int bestBaseline = 2;
 double intensityTolerance = atof( argv[1] );
 unsigned int radiusTolerance = atoi( argv[2] );
 unsigned int numberOfPixelsTolerance = atoi( argv[3] );

 cerr << "intensity tolerance: " << intensityTolerance << " radiusTolerance: " << radiusTolerance << " numberOfPixelsTolerance: " << numberOfPixelsTolerance << "\n";
 try
   {
   if(argc == 5)
     {
     bestBaselineStatus = RegressionTestImage(argv[4], argv[5], intensityTolerance, radiusTolerance, numberOfPixelsTolerance, 0, false);
     }
   else
     {
     for(int i=5;i<argc;i++)
       {
       const int currentStatus = RegressionTestImage(argv[4], argv[i], intensityTolerance,
						     radiusTolerance, numberOfPixelsTolerance,
						     0, false);
       if(currentStatus < bestBaselineStatus)
         {
         bestBaselineStatus = currentStatus;
         bestBaseline = i;
         }
       if(bestBaselineStatus == 0)
         {
         break;
         }
       }
     }
   // generate images of our closest match
   if(bestBaselineStatus == 0)
     {
     RegressionTestImage(argv[4], argv[bestBaseline], intensityTolerance,
			 radiusTolerance, numberOfPixelsTolerance, 1, false);
     }
   else
     {
     RegressionTestImage(argv[4], argv[bestBaseline], intensityTolerance,
			 radiusTolerance, numberOfPixelsTolerance, 1, true);
     }

   }
 catch(const itk::ExceptionObject& e)
   {
   std::cerr << "TestRegressionImages caught an ITK exception:\n";
   std::cerr << e.GetFile() << ":" << e.GetLine() << ":\n"
             << e.GetDescription() << "\n";
   bestBaselineStatus = -1;
   }
 catch(const std::exception& e)
   {
   std::cerr << "TestRegressionImages caught an exception:\n";
   std::cerr << e.what() << "\n";
   bestBaselineStatus = -1;
   }
 catch(...)
   {
   std::cerr << "TestRegressionImages caught an unknown exception!!!\n";
   bestBaselineStatus = -1;
   }
 cout << bestBaselineStatus << endl;
 return bestBaselineStatus;
}

// Regression Testing Code
int RegressionTestImage (const char *testImageFilename, const char *baselineImageFilename,
			 double intensityTolerance, unsigned int radiusTolerance,
			 unsigned int numberOfPixelsTolerance, 
			 int reportErrors, bool differences)
{
 // Use the factory mechanism to read the test and baseline files and convert them to double
 typedef itk::Image<double,TEST_DIMENSION_MAX> ImageType;
 typedef itk::Image<unsigned char,TEST_DIMENSION_MAX> OutputType;
 typedef itk::Image<unsigned char,2> DiffOutputType;
 typedef itk::ImageFileReader<ImageType> ReaderType;

 // Read the baseline file
 ReaderType::Pointer baselineReader = ReaderType::New();
   baselineReader->SetFileName(baselineImageFilename);
 try
   {
   baselineReader->UpdateLargestPossibleRegion();
   }
 catch (itk::ExceptionObject& e)
   {
   std::cerr << "Exception detected while reading " << baselineImageFilename << " : "  << e.GetDescription();
   return 1000;
   }

 // Read the file generated by the test
 ReaderType::Pointer testReader = ReaderType::New();
   testReader->SetFileName(testImageFilename);
 try
   {
   testReader->UpdateLargestPossibleRegion();
   }
 catch (itk::ExceptionObject& e)
   {
   std::cerr << "Exception detected while reading " << testImageFilename << " : "  << e.GetDescription() << std::endl;
   return 1000;
   }

 // The sizes of the baseline and test image must match
 ImageType::SizeType baselineSize;
   baselineSize = baselineReader->GetOutput()->GetLargestPossibleRegion().GetSize();
 ImageType::SizeType testSize;
   testSize = testReader->GetOutput()->GetLargestPossibleRegion().GetSize();

 if (baselineSize != testSize)
   {
   std::cerr << "The size of the Baseline image and Test image do not match!" << std::endl;
   std::cerr << "Baseline image: " << baselineImageFilename
             << " has size " << baselineSize << std::endl;
   std::cerr << "Test image:     " << testImageFilename
             << " has size " << testSize << std::endl;
   return 1;
   }

 // Now compare the two images
 typedef itk::DifferenceImageFilter<ImageType,ImageType> DiffType;
 DiffType::Pointer diff = DiffType::New();
   diff->SetValidInput( baselineReader->GetOutput() );
   diff->SetTestInput( testReader->GetOutput() );
   diff->SetDifferenceThreshold( intensityTolerance );
   diff->SetToleranceRadius( radiusTolerance );
   diff->UpdateLargestPossibleRegion();

 double status = diff->GetTotalDifference();
 const unsigned int numberOfFailedPixels = 
     diff->GetNumberOfPixelsWithDifferences();



 if (reportErrors &&
     ( numberOfFailedPixels > numberOfPixelsTolerance ) )
   {
   typedef itk::RescaleIntensityImageFilter<ImageType,OutputType> RescaleType;
   typedef itk::ImageFileWriter<DiffOutputType> WriterType;
   OutputType::IndexType index; index.Fill(0);
   OutputType::SizeType size; size.Fill(0);

   RescaleType::Pointer rescale = RescaleType::New();
     rescale->SetOutputMinimum(itk::NumericTraits<unsigned char>::NonpositiveMin());
     rescale->SetOutputMaximum(itk::NumericTraits<unsigned char>::max());
     rescale->SetInput(diff->GetOutput());
     rescale->UpdateLargestPossibleRegion();

   WriterType::Pointer writer = WriterType::New();
     writer->SetInput(rescale->GetOutput());
   if(differences)
     {
     // if there are discrepencies, create an diff image
     std::cout << "<DartMeasurement name=\"ImageError\" type=\"numeric/double\">";
     std::cout << status;
     std::cout <<  "</DartMeasurement>" << std::endl;

     std::cout << "<DartMeasurement name=\"NumberOfPixelsWithDifferences\" type=\"numeric/double\">";
     std::cout << numberOfFailedPixels;
     std::cout <<  "</DartMeasurement>" << std::endl;

     ::itk::OStringStream diffName;
       diffName << testImageFilename << ".diff.png";
     try
       {
       rescale->SetInput(diff->GetOutput());
       rescale->Update();
       }
     catch (...)
       {
       std::cerr << "Error during rescale of " << diffName.str() << std::endl;
       }
     writer->SetFileName(diffName.str().c_str());
     try
       {
       writer->Update();
       }
     catch (...)
       {
       std::cerr << "Error during write of " << diffName.str() << std::endl;
       }

     std::cout << "<DartMeasurementFile name=\"DifferenceImage\" type=\"image/png\">";
     std::cout << diffName.str();
     std::cout << "</DartMeasurementFile>" << std::endl;
     }
   ::itk::OStringStream baseName;
   baseName << testImageFilename << ".base.png";
   try
     {
     rescale->SetInput(baselineReader->GetOutput());
     rescale->Update();
     }
   catch (...)
     {
     std::cerr << "Error during rescale of " << baseName.str() << std::endl;
     }
   try
     {
     writer->SetFileName(baseName.str().c_str());
     writer->Update();
     }
   catch (...)
     {
     std::cerr << "Error during write of " << baseName.str() << std::endl;
     }

   std::cout << "<DartMeasurementFile name=\"BaselineImage\" type=\"image/png\">";
   std::cout << baseName.str();
   std::cout << "</DartMeasurementFile>" << std::endl;

   ::itk::OStringStream testName;
   testName << testImageFilename << ".test.png";
   try
     {
     rescale->SetInput(testReader->GetOutput());
     rescale->Update();
     }
   catch (...)
     {
     std::cerr << "Error during rescale of " << testName.str()
               << std::endl;
     }
   try
     {
     writer->SetFileName(testName.str().c_str());
     writer->Update();
     }
   catch (...)
     {
     std::cerr << "Error during write of " << testName.str() << std::endl;
     }

   std::cout << "<DartMeasurementFile name=\"TestImage\" type=\"image/png\">";
   std::cout << testName.str();
   std::cout << "</DartMeasurementFile>" << std::endl;


   }
 return (status != 0) ? 1 : 0;
}
