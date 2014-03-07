/* *	sccsid "@(#) plyimagefile.C     2.0 2/16/94" */

//
// plyimagefile.h - implementation of the reader for .ply format
//

#ifndef gfx_ply_image_file_h
#define gfx_ply_image_file_h


#include <stdio.h>
#include "image/image.h"


/************************************************************/
class GfxPlyImageFile
{
  protected:
    float _opt_center_x, _opt_center_y, _opt_center_z;
    int _height, _width;

  public:
    GfxPlyImageFile() 
    {
	_opt_center_x = 0.0f;
	_opt_center_y = 0.16;
	_opt_center_z = -1.62;
	_height = 400;
	_width = 512;
    };

    GfxImage<float> read(char* fname);
    // this one also give a confidence image
    GfxImage<float> read(char* fname, GfxImage<float> &confidence);

};

typedef struct A_Vertex {
  float x,y,z;
  float confidence;
  float intensity;
  unsigned char red,grn,blu;
} A_Vertex;


typedef struct RangePnt {
  unsigned char num_pts;
  int *pts;
} RangePnt;

#endif






