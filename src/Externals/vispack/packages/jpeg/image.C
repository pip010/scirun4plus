
/* *	sccsid "@(#) image.C     2.0 2/16/94" */

#include "image/image.h"
#include <math.h>
#include "util/array.h"
#include "util/mathutil.h"

#ifdef __memory_count
unsigned __MEMORY_SIZE = 0;
#endif

#define GAUSS_KERNEL_FACTOR (4.0)


VISImage<float> dx_kernel(int order)
{
    int i, j;
    int w, h;
    float previous, next;
    w = 2*((order + 1)/2) + 1;
//    printf("widht is %d\n", w);
    h = 1;
    VISImage<float> kernel(w, h, 1);
    kernel = 0.0;
    kernel.at((w/2), 0) = 1.0;
    for (i = 0; i < order/2; i++)
	{
	    previous = kernel.at(1, 0) - 2*kernel.at(0, 0);
	    for (j = 1; j < w - 1; j++)
		{
		    
		    next = kernel.at(j - 1, 0)
			+ kernel.at(j + 1, 0) - 2*kernel.at(j, 0);
		    kernel.at(j-1, 0) = previous;
		    previous = next;
		}
	    next = kernel.at(j - 1, 0) - 2*kernel.at(j, 0);
	    kernel.at(j-1, 0) = previous;
	    kernel.at(j, 0) = next;	    
	}
    for (i = 0; i < order%2; i++)    
	{
	    previous =  0.5*kernel.at(1, 0);
	    for (j = 1; j < w - 1; j++)
		{
		    next = -0.5*kernel.at(j - 1, 0)
			+ 0.5*kernel.at(j + 1, 0);
		    kernel.at(j-1, 0) = previous;
		    previous = next;
		}
	    next = -0.5*kernel.at(j - 1, 0);
	    kernel.at(j-1, 0) = previous;
	    kernel.at(j, 0) = next;	    
	}
    return(kernel);
}



VISImage<float> dy_kernel(int order)
{
    int i, j;
    int h, w;
    float previous, next;
    h = 2*((order + 1)/2) + 1;
//    printf("height is %d\n", h);
    w = 1;
    VISImage<float> kernel(w, h, 1);
    kernel = 0.0;
    kernel.at(0, (h/2)) = 1.0;
    for (i = 0; i < order/2; i++)
	{
	    previous = kernel.at(0, 1) - 2*kernel.at(0, 0);
	    for (j = 1; j < h - 1; j++)
		{
		    next = kernel.at(0, j - 1)
			+ kernel.at(0, j + 1) - 2*kernel.at(0, j);
		    kernel.at(0, j-1) = previous;
		    previous = next;
		}
	    next = kernel.at(0, j - 1) - 2*kernel.at(0, j);
	    kernel.at(0, j-1) = previous;
	    kernel.at(0, j) = next;	    
	}
    for (i = 0; i < order%2; i++)    
	{
	    previous =  0.5*kernel.at(0, 1);
	    for (j = 1; j < h - 1; j++)
		{
		    next = -0.5*kernel.at(0, j - 1)
			+ 0.5*kernel.at(0, j + 1);
		    kernel.at(0, j-1) = previous;
		    previous = next;
		}
	    next = -0.5*kernel.at(0, j - 1);
	    kernel.at(0, j-1) = previous;
	    kernel.at(0, j) = next;	    
	}
    return(kernel);
}




VISImage<float> gauss_row_kernel(float sigma, float window_size)
{
    int i;
    int w, h;
    if (sigma > 0.0)
	{
    w = (2*(int)(window_size*sigma)) + 1;
    h = 1;
    VISImage<float> kernel(w, h, 1);
    int pos,center = (w/2);
    float total = 0.0, value;
    for (i = 0; i < w; i++)
	{
	    pos = i - center;
	    value = exp((float)(-pos*pos)/(2.0*sigma*sigma));
	    kernel.at(i, 0) = value;
	    total += value;
	}
    return(kernel.divAssign(total));
	}
    else 
	{
	    VISImage<float> kernel(1, 1);
	    kernel.at(0, 0) = 1.0f;
	    return(kernel);
	}
}




#define SMALL_SIGMA_THRESH (1.0f)

#define k0 (3.0f/4.0f)
#define k1 (1.0f/8.0f)
#define K_VAR (1.0f/4.0f)

VISImage<float> gauss_row_kernel_small(float sigma)
{
    int i;
    int w;
    int num_iterations;
    float extra;
    float var = sigma*sigma;
    VISImage<float> k(3, 1);
    if (sigma > 0.0)
	{
	    num_iterations = (int)(var/K_VAR);
	    extra = var - K_VAR*(float)num_iterations;
	    if (extra > 0.0f)
		w = 2*(num_iterations + 1) + 1;
	    else
		w = 2*(num_iterations) + 1;

	    VISImage<float> kernel(w, 1, 1);
//	    int pos;
	    int center = (w/2);

	    kernel = 0.0f;
	    kernel.at(center, 0) = 1.0;

	    k.at(1, 0) = k0;
	    k.at(0, 0) = k1;
	    k.at(2, 0) = k1;

	    for (i = 0; i < num_iterations; i++)
		{
		    kernel = kernel.convolve(k);
		}
	    if (extra > 0)
		{
		    k.at(1, 0) = 1 - extra;
		    k.at(0, 0) = k.at(2, 0) = extra/2.0f;
		    kernel = kernel.convolve(k);
		}
	    return(kernel);
	}
    else 
	{
	    VISImage<float> kernel(1, 1);
	    kernel.at(0, 0) = 1.0f;
	    return(kernel);
	}

}

VISImage<float> gauss_col_kernel_small(float sigma)
{
    int i;
    int h;
    int num_iterations;
    float extra;
    float var = sigma*sigma;
    VISImage<float> k(1, 3);
    if (sigma > 0.0)
	{
	    num_iterations = (int)(var/K_VAR);
	    extra = var - K_VAR*(float)num_iterations;
	    if (extra > 0.0f)
		h = 2*(num_iterations + 1) + 1;
	    else
		h = 2*(num_iterations) + 1;

	    VISImage<float> kernel(1, h, 1);
//	    int pos;
	    int center = (h/2);

	    kernel = 0.0f;
	    kernel.at(center, 0) = 1.0;

	    k.at(0, 1) = k0;
	    k.at(0, 0) = k1;
	    k.at(0, 2) = k1;

	    for (i = 0; i < num_iterations; i++)
		{
		    kernel = kernel.convolve(k);
		}
	    if (extra > 0)
		{
		    k.at(0, 1) = 1 - extra;
		    k.at(0, 0) = k.at(0, 2) = extra/2.0f;
		    kernel = kernel.convolve(k);
		}
	    return(kernel);
	}
    else 
	{
	    VISImage<float> kernel(1, 1);
	    kernel.at(0, 0) = 1.0f;
	    return(kernel);
	}

}

#define SMALL_KERNEL_FACTOR (2.0f)

VISImage<float> gauss_col_kernel(float sigma)
{
    float variance = sigma*sigma;
    int total_kernels = max(fceil(SMALL_KERNEL_FACTOR*variance), 2.0f);
    float alpha = 1 - variance/total_kernels;

    VISImage<float> kernel(1, 3), r(1, 1 + 2*total_kernels);
    kernel.at(0, 0) = (1 - alpha)/2;
    kernel.at(0, 2) = (1 - alpha)/2;
    kernel.at(0, 1) = alpha;
    r = 0;
    r.at(0, total_kernels) = 1.0f;
    
    for (int i = 0; i < total_kernels; i++)
	r = r.convolve(kernel);

    return(r);
}

VISImage<float> gauss_row_kernel(float sigma)
{
    float variance = sigma*sigma;
    int total_kernels = max(fceil(SMALL_KERNEL_FACTOR*variance), 2.0f);
    float alpha = 1 - variance/total_kernels;

    VISImage<float> kernel(3, 1), r(1 + 2*total_kernels, 1);
    kernel.at(0, 0) = (1 - alpha)/2;
    kernel.at(2, 0) = (1 - alpha)/2;
    kernel.at(1, 0) = alpha;
    r = 0;
    r.at(total_kernels, 0) = 1.0f;
    
    for (int i = 0; i < total_kernels; i++)
	r = r.convolve(kernel);
    return(r);
}


VISImage<float> gauss_col_kernel(float sigma, float window_size)
{
    VISImage<float> kernel;
    if (sigma > 0.0f)
	{
	    int i;
	    int w, h;
	    h = (2*(int)(window_size*sigma)) + 1;
	    w = 1;
	    kernel = VISImage<float>(w, h, 1);
	    int pos,center = (h/2);
	    float total = 0.0, value;
	    for (i = 0; i < h; i++)
		{
		    pos = i - center;
		    value = exp((float)(-pos*pos)/(2.0*sigma*sigma));
		    kernel.at(0, i) = value;
		    total += value;
		}
	    return(kernel.divAssign(total));
	}
    else 
	{
	    kernel = VISImage<float>(1, 1);
	    kernel.at(0, 0) = 1.0f;
	    return(kernel);
	}

}

VISImage<float> gauss_kernel(float sigma)
{
    VISImage<float> x_kernel = gauss_row_kernel(sigma);
    VISImage<float> y_kernel = gauss_col_kernel(sigma);

    VISImage<float> kernel(x_kernel.width(), y_kernel.height());

    kernel = 0.0f;
    kernel.at(kernel.width()/2, kernel.height()/2) = 1.0f;
    
    return((kernel.convolve(x_kernel)).convolve(y_kernel));
}

VISImage<float> gauss_kernel(float sigma, float window_size)
{
    VISImage<float> x_kernel = gauss_row_kernel(sigma, window_size);
    VISImage<float> y_kernel = gauss_col_kernel(sigma, window_size);

    VISImage<float> kernel(x_kernel.width(), y_kernel.height());

    kernel = 0.0f;
    kernel.at(kernel.width()/2, kernel.height()/2) = 1.0f;
    
    return((kernel.convolve(x_kernel)).convolve(y_kernel));
}


VISImage<float> gauss_dx_kernel(int order, float sigma, float window_size)
{
    VISImage<float> the_dx, the_gauss, r;
    the_dx = dx_kernel(order);
    if (sigma > 0.0)
	{
	    the_gauss = gauss_row_kernel(sigma, window_size);
	    r = the_gauss.convolve(the_dx);
	    for (int i = 0; i < the_dx.width()/2; i++)
		{
		    r.at(i, 0) = 0.0;
		    r.at((r.width() - 1) - i, 0) = 0.0;
		}
	    return(r);
	}
    else 
	return(the_dx);
}

VISImage<float> gauss_dy_kernel(int order, float sigma, float window_size)
{
    VISImage<float> the_dy, the_gauss, r;
    the_dy = dy_kernel(order);
    if (sigma > 0.0)
	{
	    the_gauss = gauss_col_kernel(sigma, window_size);
	    r = the_gauss.convolve(the_dy);
	    for (int i = 0; i < the_dy.height()/2; i++)
		{
		    r.at(0, i) = 0.0;
		    r.at(0, (r.height() - 1) - i) = 0.0;
		}
	    return(r);
	}
    else 
	return(the_dy);
}


VISImage<float> gauss_dx_kernel(int order, float sigma)
{
    VISImage<float> the_dx, the_gauss, tmp, r;
    the_dx = dx_kernel(order);
    int offset = (the_dx.width()/2);
    if (sigma > 0.0)
	{
//	    the_gauss = gauss_row_kernel(sigma, 
//	    GAUSS_KERNEL_FACTOR*((float)order/4.0f + 1.0));
	    the_gauss = gauss_row_kernel(sigma); 
	    tmp = VISImage<float>(the_gauss.width() + 2*offset, 
				   1);
	    tmp = 0.0f;
	    tmp.putROI(the_gauss, offset, 0);
	    r = tmp.convolve(the_dx);

// this gets rid of the part near the edges that might ring
//	    for (int i = 0; i < r.width(); i++)
//		{
//		    r.at(i, 0) = tmp.itemAt(i + the_dx.width()/2, 0);
//		}
	    return(r);
	}
    else 
	return(the_dx);
}

VISImage<float> gauss_dy_kernel(int order, float sigma)
{
    VISImage<float> the_dy, the_gauss, r, tmp;
    the_dy = dy_kernel(order);
    int offset = (the_dy.height()/2);
    if (sigma > 0.0f)
	{
//	    the_gauss = gauss_col_kernel
//		(sigma, GAUSS_KERNEL_FACTOR*((float)order/4.0f + 1.0f));
	    the_gauss = gauss_col_kernel(sigma);

	    tmp = VISImage<float>(1, the_gauss.height() + 2*offset);
	    tmp = 0.0f;
	    tmp.putROI(the_gauss, 0, offset);
	    r = tmp.convolve(the_dy);
	    
// this gets rid of the part near the edges that might ring
//	    for (int i = 0; i < r.height(); i++)
//		r.at(0, i) = tmp.itemAt(0, i + the_dy.height()/2);
	    return(r);
	}
    else 
	return(the_dy);
}

// these next 2 methods assume that that derivatives which cross out
// of image are zero 


VISArray< VISImage<float> >* derivativeMasks(unsigned degree)
{
    VISArray< VISImage<float> > *r;
    r = new VISArray< VISImage<float> >();
    int i, j, k, l;
    VISImage<float> dx_tmp, dy_tmp, total_tmp;
    
    for (i = 0; i <= degree; i++)
	for (j = 0; j <= (degree - i); j++)
	    {
		if (j == 0)
		    r->at(dIndex(j, i)) = dy_kernel(i);
		else if (i == 0)
		    r->at(dIndex(j, i)) = dx_kernel(j);
		else
		    {
			dx_tmp = dx_kernel(i);
			dy_tmp = dy_kernel(j);
			total_tmp = VISImage<float>(dx_tmp.width(), 
						    dy_tmp.height());
			for (k = 0; k < total_tmp.height(); k++)
			    for (l = 0; l < total_tmp.width(); l++)
				total_tmp.at(l, k) = dx_tmp.itemAt(l, 0)*
				    dy_tmp.itemAt(0, k);
			r->at(dIndex(j, i)) = total_tmp;
		    }
	    }
    return(r);
}


VISArray< VISImage<float> >* derivativeMasks(unsigned degree, float scale)
{
    VISArray< VISImage<float> > *r;
    r = new VISArray< VISImage<float> >();
    int i, j, k, l;
    VISImage<float> dx_tmp, dy_tmp, total_tmp;
    
    for (i = 0; i <= degree; i++)
	for (j = 0; j <= (degree - i); j++)
	    {
		if (j == 0)
		    r->at(dIndex(j, i)) = gauss_dy_kernel(i, scale);
		else if (i == 0)
		    r->at(dIndex(j, i)) = gauss_dx_kernel(j, scale);
		else
		    {
			dx_tmp = gauss_dx_kernel(i, scale);
			dy_tmp = gauss_dy_kernel(j, scale);
			total_tmp = VISImage<float>(dx_tmp.width(), 
						    dy_tmp.height());
			for (k = 0; k < total_tmp.height(); k++)
			    for (l = 0; l < total_tmp.width(); l++)
				total_tmp.at(l, k) = dx_tmp.itemAt(l, 0)*
				    dy_tmp.itemAt(0, k);
			r->at(dIndex(j, i)) = total_tmp;
		    }
	    }
    return(r);
}
