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

#include <teem/nrrd.h>

#include <iostream>
#include <string>



int main(int argc, char **argv)
{
  if(argc != 4) {
    std::cerr << "Usage:" << std::endl;
    std::cerr << "testRegressionImageTeem testImage baselineImage diffImage" << std::endl;
    return -1;
  }
  const std::string testfilename = argv[1];
  const std::string basefilename = argv[2];
  const std::string difffilename = argv[3];

  // read test image
  Nrrd *ntest = nrrdNew();
  std::string test_regversion;
  if (nrrdLoad(ntest, airStrdup(testfilename.c_str()), 0)) 
  {
    char *err = biffGetDone(NRRD);
    std::cerr << "Read error on '" << testfilename << "': " << err << std::endl;
    free(err);
    return -1;
  }
  
  char* value;
  
  value = nrrdKeyValueGet(ntest,"netversion");
  if (value)
  {
    test_regversion = std::string(value);
    free(value);
  }

  // read base image
  Nrrd *nbase = nrrdNew();
  std::string base_regversion;
  if (nrrdLoad(nbase, airStrdup(basefilename.c_str()), 0))
  {
    NrrdIoState *nio = nrrdIoStateNew();
    // set encoding to be raw
    nio->encoding = nrrdEncodingArray[1];
    // set format to be png
    nio->format = nrrdFormatArray[3];
    // set endian to be endian of machine
    nio->endian = airMyEndian;  
  
    nrrdSave(airStrdup(basefilename.c_str()),ntest,nio);
  
    std::cerr << "No baseline image: copying result as new baseline image" << std::endl;
    return 0;
  }

  value = nrrdKeyValueGet(nbase,"netversion");
  if (value)
  {
    base_regversion = std::string(value);
    free(value);
  }
  

  if (base_regversion != test_regversion)
  {
    NrrdIoState *nio = nrrdIoStateNew();
    // set encoding to be raw
    nio->encoding = nrrdEncodingArray[1];
    // set format to be png
    nio->format = nrrdFormatArray[3];
    // set endian to be endian of machine
    nio->endian = airMyEndian;  
  
    nrrdSave(airStrdup(basefilename.c_str()),ntest,nio);
  
    std::cerr << "Baseline image is out of sync: copying new result as new baseline image" << std::endl;
    return 0; 
  }


  Nrrd *nout1 = nrrdNew();
  NrrdIter *in1 = nrrdIterNew();
  NrrdIter *in2 = nrrdIterNew();
  Nrrd *ntemp;

  // convert nrrds to floats first
  nrrdConvert(ntemp=nrrdNew(), ntest, nrrdTypeFloat);
  nrrdIterSetOwnNrrd(in1, ntemp);  
  nrrdConvert(ntemp=nrrdNew(), nbase, nrrdTypeFloat);
  nrrdIterSetOwnNrrd(in2, ntemp);  

  // unu 2op - testimage.png baseimage.png
  if (nrrdArithIterBinaryOp(nout1, nrrdBinaryOpSubtract, in1, in2)) {
    char *err = biffGetDone(NRRD);
    std::cerr << "Error performing 2op to nrrd: " << err << std::endl;
    free(err);
    return -1;
  }

  // unu 1op abs -i - -o -
  Nrrd *nout2 = nrrdNew();
  if (nrrdArithUnaryOp(nout2, nrrdUnaryOpAbs, nout1))
    {
      char *err = biffGetDone(NRRD);
      std::cerr << "Error performing 1op to nrrd: " << err << std::endl;
      free(err);
      return -1;
    }

  // unu quantize -b 8
  Nrrd *nout3 = nrrdNew();
  NrrdRange *range1 = nrrdRangeNewSet(nout2, nrrdBlind8BitRangeState);
  if (nrrdQuantize(nout3, nout2, range1, 8))
    {
      char *err = biffGetDone(NRRD);
      std::cerr << "Trouble quantizing: " << err << std::endl;
      free(err);
      return -1;
    }

  // if min/max are not both 0 then write out diff image
  NrrdRange *range2 = nrrdRangeNewSet(nout3, nrrdBlind8BitRangeFalse);
  if (range2->min == 0.0 && range2->max == 0.0) {
    std::cout << "No differences\n";
  } else {
    std::cout << "Differences in images (see " << difffilename << " for diff image)\n";
    NrrdIoState *nio = nrrdIoStateNew();
    // set encoding to be raw
    nio->encoding = nrrdEncodingArray[1];
    // set format to be png
    nio->format = nrrdFormatArray[3];
    // set endian to be endian of machine
    nio->endian = airMyEndian;
    
    if (AIR_ENDIAN != nio->endian) {
      nrrdSwapEndian(nout3);
    }

    if (nrrdSave(difffilename.c_str(), nout3, nio)) {
      char *err = biffGet(NRRD);      
      std::cerr << "Error writing nrrd " << difffilename << ": "<< err << std::endl;
      free(err);
      biffDone(NRRD);
      return -1;
    }

    return -1;
/*
     // Disable sending diff images to dashboard; only report that there is a difference in images.
     cout << "<DartMeasurementFile name=\"DifferenceImage\" type=\"image/png\">";
     cout << difffilename;
     cout << "</DartMeasurementFile>" << endl;

     cout << "<DartMeasurementFile name=\"TestImage\" type=\"image/png\">";
     cout << testfilename;
     cout << "</DartMeasurementFile>" << endl;

     cout << "<DartMeasurementFile name=\"BaselineImage\" type=\"image/png\">";
     cout << basefilename;
     cout << "</DartMeasurementFile>" << endl;
*/
  }
  return 0;
}

int RegressionTestImage (const char *testImageFilename, const char *baselineImageFilename)
{

  return 0;
}
