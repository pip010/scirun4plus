#include "rangemodel.h"
#include <image/imagefile.h>
#include <string.h>
#include "util/mathutil.h"


// ********************************************************
// ********************************************************
// ******** STAND ALONE FUNCTIONS *****************
// ********************************************************
// ********************************************************


VISVolume<float> volumeResampleSpecial(const VISVolume<float> &vol_in, 
					float scale, 
					float x, float y, float z)
{
    unsigned new_w, new_h, new_d;
     float new_x, new_y, new_z;    
    float tmp_x, tmp_y, tmp_z;    
    unsigned w = vol_in.width(), 
	h = vol_in.height(), d = vol_in.depth();
// for stoning the values at the 8 corners

    int i, j, k;
    
    new_w = (unsigned)w*scale; //why float
    new_h = (unsigned)h*scale;
    new_d = (unsigned)d*scale;


//    new_x = x*scale_x + (float_w - (float)new_w)/2.0;
//    new_y = y*scale_y + (float_h - (float)new_h)/2.0;
//    new_z = z*scale_z + (float_d - (float)new_d)/2.0;

    new_x = (x + 0.5)*scale - 0.5;
    new_y = (y + 0.5)*scale - 0.5;
    new_z = (z + 0.5)*scale - 0.5;

    VISVolume<float> volume_return(new_w, new_h, new_d);
    unsigned int x_lo, x_hi, y_lo, y_hi, z_lo, z_hi;
    float total_neighbors;

    for (k = 0; k < new_d; k++)
	for (j = 0; j < new_h; j++)
	    for (i = 0; i < new_w; i++)
		{
		    tmp_x = ((i - new_x)/scale + x); 
		    tmp_y = ((j - new_y)/scale + y); 
		    tmp_z = ((k - new_z)/scale + z); 		    
		    
		    tmp_x = MAX(tmp_x, 0);
		    tmp_x = MIN(tmp_x, w - 1);

		    tmp_y = MAX(tmp_y, 0);
		    tmp_y = MIN(tmp_y, h - 1);

		    tmp_z = MAX(tmp_z, 0);
		    tmp_z = MIN(tmp_z, d - 1);

		    x_lo = (unsigned int)floor(tmp_x);
		    y_lo = (unsigned int)floor(tmp_y);
		    z_lo = (unsigned int)floor(tmp_z);
		    x_hi = (unsigned int)ceil(tmp_x);
		    y_hi = (unsigned int)ceil(tmp_y); 
		    z_hi = (unsigned int)ceil(tmp_z); 

		    total_neighbors = vol_in.itemAt(x_lo, y_lo, z_lo)
			+ vol_in.itemAt(x_hi, y_lo, z_lo)
			+ vol_in.itemAt(x_lo, y_hi, z_lo)
			+ vol_in.itemAt(x_hi, y_hi, z_lo)
			+ vol_in.itemAt(x_lo, y_lo, z_hi)
			+ vol_in.itemAt(x_hi, y_lo, z_hi)
			+ vol_in.itemAt(x_lo, y_hi, z_hi)
			+ vol_in.itemAt(x_hi, y_hi, z_hi);

		    if (total_neighbors == -8.0f)
			volume_return.at(i, j, k) = -1.0f;
		    else if (total_neighbors == 8.0f)
			volume_return.at(i, j, k) = 1.0f;
		    else
			volume_return.at(i, j, k) = 0.0f;
		}
    return(volume_return);
}



//*****************************************************************************
//*****************************************************************************
//***********  Partion Sh@#$ ************************************************
//*****************************************************************************
//*****************************************************************************

#define TOTAL_PARTITION (M_PI/2.0f)
//#define TOTAL_PARTITION (3.0f*M_PI/5.0f)
//#define TOTAL_PARTITION (M_PI/2.0f)

// !!!! caution, these have changed since the experiments
#define FLAT_REGION (9.0f/10.0f)
//#define FLAT_REGION (0.0f)
#define TRANSITION_REGION (1.0f/10.0f)
//#define TRANSITION_REGION (1.0f)
//#define FLAT_REGION (1.0f/2.0f)
//#define TRANSITION_REGION (1.0f/2.0f)
// used in experiments
// #define FLAT_REGION (1.0f/3.0f)
// #define TRANSITION_REGION (1.0f/3.0f)

// this is a partition of unity used in some of the numerical schemes
float partition(float angle)
{
    float a; 
    a = 1.0f - angle/TOTAL_PARTITION;
    if (a >= (FLAT_REGION + TRANSITION_REGION))
	return(0.0f);
    if (a <= (FLAT_REGION))
	return(1.0f);
// else calculate the transition
    return(0.5f*(cosFast(M_PI*(a - FLAT_REGION)/TRANSITION_REGION) + 1.0f));
}


