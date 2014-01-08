// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: imagefile.h,v 1.1.1.1 2003/02/12 16:51:52 whitaker Exp $

// File:           imagefile.h
// Author:         Ross T. Whitaker
// Institution:    The University of Tennessee, Knoxville
// Contents:       The VISImageFile class.  This class is used for file i/o
//                 of the VISIm class, which can then be cast as VISImage<T>
//                 class.
// Log of changes: June 28, 1999 -- Added comments

// Three file formats are supported in this library:  fits, tiff, and
// iv.  The fits format is a floating point file format, thus,
// all the information in a VISImage<float> is preserved when using
// this file format.  The tiff format does not preserve floating point
// information.  The iv format is a format used by the 3D rendering application
// ivview, which is available on SGI systems.

// Here is a short example of how to read and write a fits file,
// which contains floats.
//
// VISImage<float> input,output;
// VISImageFile imfile;
// input = VISImage<float>(imfile.read("inputfilename.fit");
// output = input+5.0f;
// imfile.write_fits(output,"outputfilename.fit");
//


#ifndef iris_image_file_h
#define iris_image_file_h

#define USE_FITS
#define USE_TIFF
#define USE_JPEG

#include "image/classtypes.h"
#include "util/defs.h"



class VISImageFile 
{
  protected:
    const char* _filename;
    
  public:
    VISImageFile() {}
    ~VISImageFile() {}
#ifdef USE_TIFF
    int write(const VISIm&, const char*);
    int write_tiff(const VISImageRGBA&, const char*);
    int write_tiff(const VISImage<float>&, const char*);
    int write_tiff(const VISImage<int>&, const char*);
    int write_tiff(const VISImage<byte>&, const char*);
    int write_tiff(const VISImage<short>&, const char*);
#endif
#ifdef USE_FITS
    int write_fits(const VISImage<byte>& image, const char* fname);
    int write_fits(const VISImage<int>& image, const char* fname);
    int write_fits(const VISImage<float>& image, const char* fname);
    int write_fits(const VISImage<short>&, const char*);
// writes a 3D inventor file with points (x, y, f(x, y))
    int write_iv(const VISImage<float>&, const char*);
    int write_iv(const VISImage<float>&,
		 const VISImage<float>&,
		 const VISImage<float>&,
		 const char*);
    int write_iv(const VISImage<float>&,
		 const VISImage<float>&,
		 const VISImage<float>&,
		 const VISImageRGBA&,
		 const char*);
    int write_iv_trans(const VISImage<float>&,
		       const VISImage<float>&,
		       const VISImage<float>&,
		       const VISImage<float>&,
		       const char*);
    int write_iv(const VISImage<float>&, 
		 const VISImageRGBA&, const char*);

    int printFitsError(int fits_status);
#endif

#ifdef USE_JPEG
    int write_jpeg(const VISImageRGBA&, const char*, int quality);
    int write_jpeg(const VISImage<float>, const char*, int quality);
#endif
    VISIm read(const char* fname);
    virtual void update( const char* fname ) { _filename = fname; }

};

#endif







