// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: volume.cxx,v 1.1.1.1 2003/02/12 16:51:54 whitaker Exp $

/* *	sccsid "@(#) volume.C     2.0 2/16/94" */

#include "vol/volume.h"
#include <math.h>
#include "util/mathutil.h"

#define GAUSS_KERNEL_FACTOR (4.0)


// these are a hack to get sgi templates to work
// otherwise your risk turning up unresolved symbols (for some reason).
VISVolume<unsigned> __vol_unsigned;
VISVolume<int> __vol_int;
VISVolume<rgba> __vol_rgba;
VISVolume<float> __vol_float;
VISVolume<byte> __vol_byte;


VISVolume<float> dz_kernel(int order)
{
    int i, j;
    int h, w, d;
    float previous, next;
    d = 2*((order + 1)/2) + 1;
//    printf("height is %d\n", h);
    h = w = 1;
    VISVolume<float> kernel(w, h, d);
    kernel = 0.0;
    kernel.at(0, 0, (d/2)) = 1.0;
    for (i = 0; i < order/2; i++)
	{
	  previous = kernel.at(0, 0, 1) - 2*kernel.at(0, 0, 0);
	    for (j = 1; j < d - 1; j++)
		{
		    next = kernel.at(0, 0, j - 1)
			+ kernel.at(0, 0, j + 1) - 2*kernel.at(0, 0, j);
		    kernel.at(0, 0, j-1) = previous;
		    previous = next;
		}
	    next = kernel.at(0, 0, j - 1) - 2*kernel.at(0, 0, j);
	    kernel.at(0, 0, j-1) = previous;
	    kernel.at(0, 0, j) = next;	    
	}
    for (i = 0; i < order%2; i++)    
	{
	    previous = -0.5f*kernel.at(0, 0, 1);
	    for (j = 1; j < d - 1; j++)
		{
		    next = 0.5f*kernel.at(0, 0, j - 1)
		      - 0.5f*kernel.at(0, 0, j + 1);
		    kernel.at(0, 0, j-1) = previous;
		    previous = next;
		}
	    next = 0.5f*kernel.at(0, 0, j - 1);
	    kernel.at(0, 0, j-1) = previous;
	    kernel.at(0, 0, j) = next;	    
	}
    return(kernel);
}

VISVolume<float> gauss_slice_kernel(float sigma)
{
    float variance = sigma*sigma;
    int total_kernels =(int) VISmax((float)ceil(2.0f*variance), 2.0f);
    float alpha = 1 - variance/total_kernels;

    VISVolume<float> kernel(1, 1, 3), r(1, 1, 1 + 2*total_kernels);
    kernel.at(0, 0, 0) = (1 - alpha)/2;
    kernel.at(0, 0, 2) = (1 - alpha)/2;
    kernel.at(0, 0, 1) = alpha;
    r = 0;
    r.at(0, 0, total_kernels) = 1.0f;
    
    for (int i = 0; i < total_kernels; i++)
	r = r.convolve(kernel);

    return(r);
}

VISVolume<float> gauss_dz_kernel(int order, float sigma)
{
    VISVolume<float> the_dz, the_gauss, r;
    the_dz = dz_kernel(order);
    the_gauss = gauss_slice_kernel(sigma);
    r = the_gauss.convolve(the_dz);
    for (int i = 0; i < the_dz.width()/2; i++)
	{
	    r.at(0, 0, i) = 0.0;
	    r.at(0, 0, (r.depth() - 1) - i) = 0.0;
	}
    return(r);
}


VISArray< VISVolume<float> >* derivativeVolMasks(unsigned degree)
{
    VISArray< VISVolume<float> > *r;
    r = new VISArray< VISVolume<float> >();
    int h, i, j;
    VISVolume<float> dx_tmp, dy_tmp, dz_tmp, total_tmp;
    
    for (i = 0; i <= degree; i++)
	for (j = 0; j <= (degree - i); j++)
	    for (h = 0; h <= degree - (i + j); h++)
		{
// these checks just optimze things if the kernel is not mixed
		    if ((j == 0)&&(i == 0))
			r->poke(dVolIndex(j, i, h)) = dz_kernel(h);
		    else if ((j == 0)&&(h == 0))
			r->poke(dVolIndex(j, i, h)) = dy_kernel(i);
		    else if ((i == 0)&&(h == 0))
			r->poke(dVolIndex(j, i, h)) = dx_kernel(j);
		    else
			{
			    dx_tmp = dx_kernel(j);
			    dy_tmp = dy_kernel(i);
			    dz_tmp = dz_kernel(h);
			    total_tmp = VISVolume<float>(dx_tmp.width(), 
							 dy_tmp.height(), 
							 dz_tmp.depth()
							 );
//
// this doesn't produce the right kernel, you want to convolve it
//			    for (m = 0; m < total_tmp.depth(); m++)
//				for (k = 0; k < total_tmp.height(); k++)
//				    for (l = 0; l < total_tmp.width(); l++)
//					total_tmp.at(l, k, m) 
//					    = dx_tmp.itemAt(l, 0, 0)*
//					    dy_tmp.itemAt(0, k, 0)*
//					    dz_tmp.itemAt(0, 0, m);
//
			    total_tmp = 0.0;
			    total_tmp.at(total_tmp.width()/2, 
					 total_tmp.height()/2, 
					 total_tmp.depth()/2) = 1.0;
			    r->poke(dVolIndex(j, i, h)) 
				= (((total_tmp.convolve(dx_tmp))
				    .convolve(dy_tmp)).convolve(dz_tmp));
			}
		}
    return(r);
}

unsigned dVolIndex(unsigned c, unsigned b, unsigned a)
{
    return((((a + b + c)*(2*(a + b + c)*(a + b + c) + 3*(a + b + c) + 1))/6
	    + (a + b + c)*(a + b + c + 1)/2)/2
	   + (a + b)*(a + b + 1)/2 + a);
}
