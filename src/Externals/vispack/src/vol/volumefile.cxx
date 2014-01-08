
// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: volumefile.cxx,v 1.3 2003/02/26 22:54:24 whitaker Exp $
/*
 *  Someone needs to comment this
 *
 *  And put an SCCS-ID here
 */


#include "image/imagefile.h"
#include "vol/volumefile.h"
#include "image/image.h"
#include "image/imageRGBA.h"
#include "image/vispacktiffioCC.h"
#include "vol/march.h"
#include <stdio.h>
#ifdef USE_HIPS
#include "image/hips.h"
#endif

// #define DEBUG_PRINT(a) {printf("%s \n", a);}
#define DEBUG_PRINT(a) {}

int VISVolumeFile::write_float(const VISVolume<float>& volume, 
			      const char* fname) const
{
//
// no longer needed since we have moved to ascii in header
// Ross 6-30-98    
//
//    if ((volume.width() > 255)||(volume.height() > 255)
//	||(volume.depth() > 255))
//    {
//	printf("ERROR: 3d file too big for avs format (limit 0-255)\n");
//	exit(-1);
//    }

    FILE* out_file;

    if ((out_file = fopen(fname, "w")) == NULL)
	{
	    printf("ERROR: output file open failed\n");
	    return(-1);
	}

    printf("write float vol width %d, height %d, depth %d\n", 
	   volume.width(), 
	   volume.height(), 
	   volume.depth());
// doing it in ascii now
//    
//    unsigned char uchar_tmp = (unsigned char)volume.width();
//    fwrite((const void*)&uchar_tmp, sizeof(unsigned char), 1, out_file);
//    uchar_tmp = (unsigned char)volume.height();
//    fwrite((const void*)&uchar_tmp, sizeof(unsigned char), 1, out_file);
//    uchar_tmp = (unsigned char)volume.depth();
//    fwrite((const void*)&uchar_tmp, sizeof(unsigned char), 1, out_file);
//

// this stuff now includes the number of bytes in the header
    fprintf(out_file, "4\n %d %d %d \n", volume.width(), volume.height(), 
	    volume.depth());

    fwrite((const void*)(volume.rep()->buffer()), sizeof(float),
	   volume.height()*volume.width()*volume.depth(), 
	   out_file);

    fclose(out_file);

    return(0);

}




int VISVolumeFile::write_byte(const VISVolume<byte>& volume, 
			      const char* fname) const
{
//
// no longer needed since we have moved to ascii in header
// Ross 6-30-98    
//
//    if ((volume.width() > 255)||(volume.height() > 255)
//	||(volume.depth() > 255))
//    {
//	printf("ERROR: 3d file too big for avs format (limit 0-255)\n");
//	exit(-1);
//    }

    FILE* out_file;

    if ((out_file = fopen(fname, "w")) == NULL)
	{
	    printf("ERROR: output file open failed\n");
	    return(-1);
	}

    printf("write byte vol width %d, height %d, depth %d\n", 
	   volume.width(), 
	   volume.height(), 
	   volume.depth());
// doing it in ascii now
//    
//    unsigned char uchar_tmp = (unsigned char)volume.width();
//    fwrite((const void*)&uchar_tmp, sizeof(unsigned char), 1, out_file);
//    uchar_tmp = (unsigned char)volume.height();
//    fwrite((const void*)&uchar_tmp, sizeof(unsigned char), 1, out_file);
//    uchar_tmp = (unsigned char)volume.depth();
//    fwrite((const void*)&uchar_tmp, sizeof(unsigned char), 1, out_file);
//

// this stuff now includes the number of bytes in the header
    fprintf(out_file, "1\n %d %d %d \n", volume.width(), volume.height(), 
	    volume.depth());

    fwrite((const void*)(volume.rep()->buffer()), sizeof(byte),
	   volume.height()*volume.width()*volume.depth(), 
	   out_file);

    fclose(out_file);

    return(0);

}



VISVol VISVolumeFile::read_float(const char* fname) const
{
    VISVolume<float> r;
    FILE* in_file;
    int h, w, d;

    if ((in_file = fopen(fname, "r")) != NULL)
	{

//
// doing it in ascii now
//
//	    fread((void*)&w, 
//		  sizeof(unsigned char), 1, in_file);
//	    fread((void*)&h, 
//		  sizeof(unsigned char), 1, in_file);
//	    fread((void*)&d, 
//		  sizeof(unsigned char), 1, in_file);
//
//

// this is done so that we can now have the # bytes in the header but 
// also backward compatibility
	    fscanf(in_file, "%d \n", &w);
	    if (w == 4)
		fscanf(in_file, "%d %d %d \n", &w, &h, &d);
	    else
		fscanf(in_file, "%d %d \n", &h, &d);

	    r = VISVolume<float>(w, h, d);
	    printf("read float vol width %d, height %d, depth %d\n", 
		   r.width(), 
		   r.height(), 
		   r.depth());
	    fread((void*)(r.repRef()->bufferRef()), sizeof(float),
		  r.height()*r.width()*r.depth(), 
		  in_file);
	    fclose(in_file);
	}
    return(r);
}

VISVol VISVolumeFile::read(const char* fname) const
{
    VISVolume<float> r_float;
    VISVolume<byte> r_byte;
    VISVolume<rgba> r_rgba;
    VISVolume<rgba> r_short;
    byte r, g, b;
    FILE* in_file;
    int h, w, d, bytes;

    if ((in_file = fopen(fname, "r")) != NULL)
	{

//
// doing it in ascii now
//
//	    fread((void*)&w, 
//		  sizeof(unsigned char), 1, in_file);
//	    fread((void*)&h, 
//		  sizeof(unsigned char), 1, in_file);
//	    fread((void*)&d, 
//		  sizeof(unsigned char), 1, in_file);
//
//

// this is done so that we can now have the # bytes in the header but 
// also backward compatibility
	    fscanf(in_file, "%d \n", &bytes);
	    fscanf(in_file, "%d %d %d \n", &w, &h, &d);

	    switch(bytes)
		{
		  case 4:
		    r_float = VISVolume<float>(w, h, d);
		    printf("read float vol width %d, height %d, depth %d\n", 
			   r_float.width(), 
			   r_float.height(), 
			   r_float.depth());
		    fread((void*)(r_float.repRef()->bufferRef()), sizeof(float),
			  r_float.height()*r_float.width()*r_float.depth(), 
			  in_file);
		    fclose(in_file);
		    return(r_float);
		    break;
		  case 2:
		    r_short = VISVolume<float>(w, h, d);
		    printf("read float vol width %d, height %d, depth %d\n", 
			   r_short.width(), r_short.height(), r_short.depth());
		    fread((void*)(r_short.repRef()->bufferRef()), sizeof(short),
			  r_short.height()*r_short.width()*r_short.depth(), 
			  in_file);
		    fclose(in_file);
		    return(r_short);
		    break;
		  case 1:
		    r_byte = VISVolume<byte>(w, h, d);
		    printf("read float vol width %d, height %d, depth %d\n", 
			   r_byte.width(), 
			   r_byte.height(), 
			   r_byte.depth());
		    fread((void*)(r_byte.repRef()->bufferRef()), sizeof(byte),
			  r_byte.height()*r_byte.width()*r_byte.depth(), 
			  in_file);
		    fclose(in_file);
		    return(r_byte);
		    break;
		  case 3:
		    byte  r,g,b;
		    r_rgba = VISVolume<rgba>(w, h, d);
		    printf("read rgba vol width %d, height %d, depth %d\n", 
			   r_rgba.width(), 
			   r_rgba.height(), 
			   r_rgba.depth());
		    for(int i=0;i<d;i++)
			for(int j=0; j<h;j++)
			    for(int k=0; k<w; k++) 
				{	
				    fread((void*)&r,sizeof(byte),1,in_file);
				    fread((void*)&g,sizeof(byte),1,in_file);
				    fread((void*)&b,sizeof(byte),1,in_file);
				    (r_rgba.poke(k,j,i)).r(r);
				    (r_rgba.poke(k,j,i)).g(g);
				    (r_rgba.poke(k,j,i)).b(b);
				}
		    fclose(in_file);
		    return(r_rgba);
		    break;
		}
	}
}


int VISVolumeFile::write_rgba(const VISVolume<rgba>& volume, 
			      const char* fname) const
{
//

    FILE* out_file;

    if ((out_file = fopen(fname, "w")) == NULL)
	{
	    printf("ERROR: output file open failed\n");
	    return(-1);
	}

    printf("write color vol width %d, height %d, depth %d\n", 
	   volume.width(), 
	   volume.height(), 
	   volume.depth());

// this stuff now includes the number of bytes in the header
    fprintf(out_file, "3\n %d %d %d \n", volume.width(), volume.height(), 
	    volume.depth());


    byte  r,g,b;
    for(int i=0;i<volume.depth();i++)
	for(int j=0; j<volume.height();j++)
	    for(int k=0; k<volume.width(); k++) 
		{	
		    r=(volume.peek(k,j,i)).r();
		    g=(volume.peek(k,j,i)).g();
		    b=(volume.peek(k,j,i)).b();
		
		    fwrite((const void *)(&r), sizeof(byte),1, out_file);     
		    fwrite((const void *)(&g), sizeof(byte),1, out_file);
		    fwrite((const void *)(&b), sizeof(byte),1, out_file);
		}

 
    fclose(out_file);

    return(0);

}


int VISVolumeFile::write_binary(const VISVolume<byte>& volume, 
				 const char* fname) const
{
//
// no longer needed since we have moved to ascii in header
// Ross 6-30-98    
//
//    if ((volume.width() > 255)||(volume.height() > 255)
//	||(volume.depth() > 255))
//    {
//	printf("ERROR: 3d file too big for avs format (limit 0-255)\n");
//	exit(-1);
//    }

    FILE* out_file;

    if ((out_file = fopen(fname, "w")) == NULL)
	{
	    printf("ERROR: output file open failed\n");
	    return(-1);
	}

    printf("write float vol width %d, height %d, depth %d\n", 
	   volume.width(), 
	   volume.height(), 
	   volume.depth());
// doing it in ascii now
//    
//    unsigned char uchar_tmp = (unsigned char)volume.width();
//    fwrite((const void*)&uchar_tmp, sizeof(unsigned char), 1, out_file);
//    uchar_tmp = (unsigned char)volume.height();
//    fwrite((const void*)&uchar_tmp, sizeof(unsigned char), 1, out_file);
//    uchar_tmp = (unsigned char)volume.depth();
//    fwrite((const void*)&uchar_tmp, sizeof(unsigned char), 1, out_file);
//
    unsigned short w = volume.width(), h = volume.height(), 
	d = volume.depth();

//    fwrite(out_file, "%d %d %d \n", volume.width(), volume.height(), 
//	   volume.depth());

//    fwrite((const void*)&w, sizeof(unsigned short), 1, out_file);
//    fwrite((const void*)&h, sizeof(unsigned short), 1, out_file);
//    fwrite((const void*)&d, sizeof(unsigned short), 1, out_file);

// zero bytes indicated on the first line is "binary"
    fprintf(out_file, "0\n %d %d %d \n", volume.width(), volume.height(), 
	    volume.depth());
    
    unsigned num_bytes = w*h*d/8 + 1*((w*h*d%8) > 0);

    byte *data = new byte[num_bytes];
    const byte* buffer = volume.rep()->buffer();

    int count = 0;
    int data_local = 0;
    byte byte_tmp = 0;
    for (int i = 0; i < w*h*d; i++)
	{
	    if (buffer[i] > 0)
		byte_tmp++;
	    if (count++ == 7)
		{
		    count = 0;
		    data[data_local++] = byte_tmp;
		    byte_tmp = 0;
		}
	    else 
		byte_tmp <<= 1;
	}
    
    fwrite((const void*)(data), sizeof(byte),
	   num_bytes, 
	   out_file);

    fclose(out_file);

    return(0);

}


VISVol VISVolumeFile::read_binary(const char* fname) const
{
//
// Use if you need smaller file size. Otherwise use read_float
// Ross 6-30-98    
//

  VISVolume<byte> im_r = VISVolume<byte>();
    FILE* in_file;

    if ((in_file = fopen(fname, "r")) == NULL)
	{
	    printf("ERROR: output file open failed\n");
	    return(im_r);
	}

    unsigned short w, h, d;

//    fwrite(in_file, "%d %d %d \n", volume.width(), volume.height(), 
//	   volume.depth());

//    fread((void*)&w, sizeof(unsigned short), 1, in_file);
//    fread((void*)&h, sizeof(unsigned short), 1, in_file);
//    fread((void*)&d, sizeof(unsigned short), 1, in_file);

// this is done so that we can now have the # bytes in the header but 
// also backward compatibility
    fscanf(in_file, "%d \n", &w);
// zero bytes will be considered binary
    if (w == 0)
	fscanf(in_file, "%d %d %d \n", &w, &h, &d);
    else
	fscanf(in_file, "%d %d \n", &h, &d);


    unsigned num_bytes = w*h*d/8 + 1*((w*h*d%8) > 0);

    byte *data = new byte[num_bytes];

    fread((void*)(data), sizeof(byte),
	   num_bytes, 
	   in_file);

    VISVolume<byte> r(w, h, d);

    byte* buffer = r.repRef()->bufferRef();

    int count = 0;
    int data_local = 0;
    byte byte_tmp;
    byte mask  = 1;
    mask <<= 7;

    byte_tmp = data[0];

    for (int i = 0; i < w*h*d; i++)
	{
	    if (byte_tmp&mask)
		buffer[i] = 1;
	    else
		buffer[i] = 0;
	    if (count++ == 7)
		{
		    count = 0;
		    byte_tmp = data[data_local++];
		}
	    else 
		byte_tmp <<= 1;
	}
    
    fclose(in_file);

    return(r);
}


int VISVolumeFile::write_marchingCubes(const VISVolume<float>& volume, 
    float iso_value, char* fname) const{
    
    size s;
    s.x = volume.width();
    s.y = volume.height();
    s.z = volume.depth();
    Marching_Cubes((float*)(volume.rep()->buffer()), s, iso_value, fname);
    return(0);
}


VISVolume<float> volConvertBL(const VISVolume<float> &vol)
{

  int w = vol.width(), h = vol.height(), d = vol.depth();
  int i, j, k, l;
  VISVolume<float> r(w, h, d);
  const char *buf_in = (const char*)(vol.rep())->buffer();
  char *buf_out  = (char*)(r.repRef())->bufferRef();

  for (k = 0; k < w*h*d; k++)
    {
      for (j = 0; j < 4; j++)
	buf_out[3 - j] = buf_in[j];
      buf_out += 4;
      buf_in += 4;
    }

  return(r);
}

void VISVolumeFile::writeVTK(const VISVolume<float> &vol_float, const char* filename, boolean swap_bytes) const
{
  //   float min = (vol_float).min(), max =  (vol_float).max();
  //   float the_max = VISmax(max, (float)fabs(min));
  //   VISVolume<byte> vol_byte((vol_float)*(127.0f/(the_max)) + 128.0f);

  FILE* out_file;
  
  VISVolume<float> vol;
  if (swap_bytes)
    vol = volConvertBL(vol_float);
  else
    vol = vol_float;

  if ((out_file = fopen(filename, "w")) == NULL)
    printf("ERROR WriteVTK: output file open failed\n");

  fprintf(out_file, "# vtk DataFile Version 2.0\n");
  fprintf(out_file, "test volume\n");
  fprintf(out_file, "BINARY\n");
  fprintf(out_file, "DATASET STRUCTURED_POINTS\n");
  fprintf(out_file, "DIMENSIONS %d %d %d\n", vol.width(), vol.height(), vol.depth());
  fprintf(out_file, "ASPECT_RATIO 1.0 1.0 1.0\n");
  fprintf(out_file, "ORIGIN 0.0 0.0 0.0\n");
  fprintf(out_file, "POINT_DATA %d \n", vol.width()*vol.height()*vol.depth());
  fprintf(out_file, "SCALARS scalars float\n");
  fprintf(out_file, "LOOKUP_TABLE default\n");

  fwrite((const void*)(vol.rep()->buffer()), sizeof(float),
	 vol.height()*vol.width()*vol.depth(), 
	 out_file);
  fclose(out_file);
}

void VISVolumeFile::writeText(const VISVolume<float> &vol_float, const char* filename, const char* header, boolean swap_bytes) const
{
  //   float min = (vol_float).min(), max =  (vol_float).max();
  //   float the_max = VISmax(max, (float)fabs(min));
  //   VISVolume<byte> vol_byte((vol_float)*(127.0f/(the_max)) + 128.0f);

  FILE* out_file;
  
  VISVolume<float> vol;
  if (swap_bytes)
    vol = volConvertBL(vol_float);
  else
    vol = vol_float;

  if ((out_file = fopen(filename, "w")) == NULL)
    printf("ERROR WriteVTK: output file open failed\n");

  fprintf(out_file, header);

  fwrite((const void*)(vol.rep()->buffer()), sizeof(float),
	 vol.height()*vol.width()*vol.depth(), 
	 out_file);
  fclose(out_file);
}


    
