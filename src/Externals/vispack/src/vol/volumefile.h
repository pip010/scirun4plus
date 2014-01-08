// File:           volumefile.h
// Author:         Ross T. Whitaker
// Institution:    The University of Tennessee, Knoxville
// Contents:       The VISVolumeFile class.  This class is used for file i/o
//                 of the VISVol class, which can then be cast as VISVolume<T>
//                 class.
// Log of changes: July 19, 1999 -- Added comments

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


#ifndef vis_volume_file_h
#define vis_volume_file_h

#include "volume.h"

class VISVolumeFile 
{
  protected:
    const char* _filename;
    //explorer functions no longer in use
    //int write_explorer_lattice(const VISVol vol, 
    //		     FILE* file) const;
    //int write_explorer_data(int index, int& current_index, const VISVol vol, 
    //		     FILE* file) const;
    //int write_explorer_coord(int index, int& current_index, const VISVol vol, 
    //		     FILE* file) const;
    //int write_explorer_values(int index, int& current_index, 
    //		    const VISVolume<byte> vol, FILE* file) const;
    //int write_explorer_values(int index, int& current_index, 
    //		    const VISVolume<float> vol, FILE* file) const;
    //int write_explorer_values(int index, int& current_index, 
    //		    const VISVolume<int> vol, FILE* file) const;
    //int write_explorer_values(int index, int& current_index, 
    //		    const VISVolume<unsigned> vol, FILE* file) const;
    //int write_explorer_values(int index, int& current_index, 
    //		    const VISVolume<rgba> vol, FILE* file) const;
    //int write_explorer_dims(int index, int& current_index, 
    //		    const VISVol vol, FILE* file) const;
    //int write_explorer_bbox(int index, int& current_index, 
    //		    const VISVol vol, FILE* file) const;
  public:
    VISVolumeFile() {}
    ~VISVolumeFile() {}
    //int write(const VISVol&, const char*) const;
    //#ifdef USE_HIPS
    // int write_hips(const VISVol&, const char*) const;
    //#endif
    //int write_avs(const VISVol&, const char*) const;
    //int write_avs(const VISVolume<byte>& volume, const char*) const;
    //int write_explorer(const VISVol&, const char*) const;
    
    //write and read "raw" float file with an ascii header
    //a "raw" format is used because no good widely-used volume file
    //format exists.  The header format is the following
    //volume.width() volume.height() volume.depth()
    int write_float(const VISVolume<float>& volume, const char*) const;
    VISVol read_float(const char*) const;

    int write_byte(const VISVolume<byte>& volume, const char*) const;
  //  int write_ushort(const VISVolume<byte>& volume, const char*) const;
    int write_rgba(const VISVolume<rgba>& volume, const char*) const;

    //write and read a volume of type byte.  The volume is assumed to be
    //binary and the values are stored in the file as bits.  Header is binary
    //The header format is the folowing (where the values are stored in the
    //file as unsigned short's):
    //volume.width() volume.height() volume.depth()
    int write_binary(const VISVolume<byte>& volume, const char*) const;
    VISVol read_binary(const char*) const;
    //VISVol read_avs(const char*) const;
    //uses avs, which is no longer used
    //VISVol read(const char* fname) const;

// figure out what type of volume it is based on the first field in the header
    VISVol read(const char*) const;

    // applies marchingCubes algorithm to a volume.
    // iso_value indicates the value that the algorithm uses to determine
    // the surface location in the volume.
    // fname is a filename of an ivview file created by the function
    // that contains the result.
    int write_marchingCubes(const VISVolume<float>& volume, 
			float iso_value, char* fname) const;

    void writeVTK(const VISVolume<float> &vol_float, const char* filename, boolean swap_bytes=false) const;
    void writeText(const VISVolume<float> &vol_float, const char* filename, const char* header, boolean swap_bytes=false) const;
};

VISVolume<float> volConvertBL(const VISVolume<float> &vol);

template<class T>
int read_raw(char* fname, VISVolume<T> &ret)
{
  FILE* in_file;
   if ((in_file = fopen(fname, "r")) == NULL)
	{
	    printf("ERROR: read raw file open failed: %s\n", fname);
	    return(-1);
	}
   int w = ret.width(), h = ret.height(), d = ret.depth();
   fread((void*)(ret.repRef()->bufferRef()), sizeof(T),
	 w*h*d, 
	 in_file);
   fclose(in_file);
   return(0);
}


// returns -1 is something fails
template<class T>
int write_raw(char* fname, const VISVolume<T> &ret)
{
  FILE* file;
   if ((file = fopen(fname, "w")) == NULL)
	{
	    printf("ERROR: read raw file open failed: %s\n", fname);
	    return(-1);
	}
   int w = ret.width(), h = ret.height(), d = ret.depth();
   if (fwrite((void*)(ret.rep()->buffer()), sizeof(T),
	 w*h*d, 
	      file) != (w*h*d))
     {
       printf("ERROR: write_raw failed to write data: %s\n", fname);
       fclose(file);
       return(-1);
     }
   fclose(file);
   return(0);
}


#endif






