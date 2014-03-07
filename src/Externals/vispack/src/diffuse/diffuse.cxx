// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************


#include "diffuse/diffuse.h"
#include "image/imagefile.h"
#include "util/mathutil.h"

void Diffuse::calculateConductance()
{

    VISImage<float> f_x, f_y, 
	dx_half, dy_half, 
	dx, dy, 
	x_kernel, y_kernel, im, cond_mag;
    
    float scale;

    if ((_weights.n() != _num_intermediate_images)||
	(_intermediate_images.n() != _num_intermediate_images))
	{
	    ERROR("diffuse::calculateConductance - array size mismatch\n");
	    return;
	}
    
    im = _intermediate_images.at(0);
    if ((!im.compareSize(_conductance_dx))||
	(!im.compareSize(_conductance_dy)))
	_conductance_dx = _conductance_dy = im.createToSize();
    
    _conductance_dx = 0.0f;
    _conductance_dy = 0.0f;
    int w, h, j, i;

    for (i = 0; i < _num_intermediate_images; i++)
	{
//	    scale = -0.5f/_weights.itemAt(i);
	    scale = 1.0f/_weights.itemAt(i);
	    
	    im = _intermediate_images.at(i);

	    f_x = im.dx();
	    f_y = im.dy();
	    
	    w = im.width();
	    h = im.height();

// need to set the boundaries of the derivatives
	for (j = 0; j < h; j++)
	    {
		f_x.at(0, j) = im.itemAt(1, j) 
		    - im.itemAt(0, j);
		f_x.at(w - 1, j) = im.itemAt(w - 1, j) 
		    - im.itemAt(w - 2, j);
	    }

	for (j = 0; j < w; j++)
	    {
		f_y.at(j, 0) = im.itemAt(j, 1) 
		    - im.itemAt(j, 0);
		f_y.at(j, h - 1) = im.itemAt(j, h - 1) 
		    - im.itemAt(j, h - 2);
	    }
	    
	    dx_half = im.dxHalfForward();
	    dy_half = im.dyHalfForward();

//
// not for now
//	    float grad_mag_average = (f_x.power(2) + f_y.power(2)).average();
//

	    x_kernel = VISImage<float>(3, 1);
	    x_kernel = 0.5f;
	    x_kernel.at(0, 0) = 0.0f;

	    y_kernel = VISImage<float>(1, 3);
	    y_kernel = 0.5f;
	    y_kernel.at(0, 0) = 0.0f;

	    f_x = f_x.convolve(x_kernel);
	    f_y = f_y.convolve(y_kernel);

	    cond_mag = (dx_half.power(2) + f_y.power(2));
	    _conductance_dx += (scale/cond_mag.average())*(cond_mag);

	    cond_mag = dy_half.power(2) + f_x.power(2);
	    _conductance_dy += (scale/cond_mag.average())*(cond_mag);
	}

    float* buffer;

    buffer = (_conductance_dx.repRef())->bufferRef();
    for ( i = 0; i < w*h; i++)
	*(buffer++) = gaussFast(-(*buffer));

    buffer = (_conductance_dy.repRef())->bufferRef();
    for ( i = 0; i < w*h; i++)
	*(buffer++) = gaussFast(-(*buffer));
    
//    _conductance_dx = 1.0f/(1.0f + _conductance_dx);
//    _conductance_dy = 1.0f/(1.0f + _conductance_dy);

}

void Diffuse::update()
{

    if (_change.n() != _images.n())
 	{
	    ERROR("diffuse::update - array size mismatch\n");
	    return;
	}
  
    for (int i = 0; i < _change.n(); i++)
	{
	    _images.at(i) += _change.itemAt(i);
	}
    
}

void Diffuse::calculateChange()
{
    VISImage<float> im;
    for (int i = 0; i < _images.n(); i++)
	{
	    //
	    // this could be made faster by putting the scalar into the 
	    // actual routine
	    //
	    _change.at(i) = 0.125f*
		(_images.itemAt(i)).anisoDiffuse(_conductance_dx, 
						 _conductance_dy);
	}

}


void SphericalRange::calculateIntermediateFeatures()
{
    _intermediate_images.at(0) = _images.itemAt(0);
    _intermediate_images.at(1) = _images.itemAt(1);

    
    VISImage<float> norm, dx_tmp, dy_tmp;

// this is if you don't want 4 features
// range times delta angle
//    dx_tmp = _delta_x*((_images.itemAt(0)).dx()).div(_images.itemAt(0), 0.0f);
//    dy_tmp = _delta_y*((_images.itemAt(0)).dy()).div(_images.itemAt(0), 0.0f);

// this is if you have all four features
// range times delta angle
    dx_tmp = _delta_x*((_images.itemAt(2)).div(_images.itemAt(0), 0.0f));
    dy_tmp = _delta_y*((_images.itemAt(3)).div(_images.itemAt(0), 0.0f));

    norm = (dx_tmp.power(2)
	    + dy_tmp.power(2)
	    + 1.0f).sqrt();
    
    _intermediate_images.at(2) = dx_tmp/norm;
    _intermediate_images.at(3) = dy_tmp/norm;
    _intermediate_images.at(4) = 1.0f/norm;

}

void SphericalRange2::calculateIntermediateFeatures()
{
    _intermediate_images.at(0) = _images.itemAt(0);
    _intermediate_images.at(1) = _images.itemAt(1);

    
    VISImage<float> norm, dx_tmp, dy_tmp;

// this is if you don't want 4 features
// range times delta angle
    dx_tmp = _delta_x*(normal_dx().div(_images.itemAt(0), 0.0f));
    dy_tmp = _delta_y*(normal_dy().div(_images.itemAt(0), 0.0f));
//_delta_x*((_images.itemAt(0)).dx()).div(_images.itemAt(0), 0.0f);
//_delta_y*((_images.itemAt(0)).dy()).div(_images.itemAt(0), 0.0f);

// this is if you have all four features
// range times delta angle
//    dx_tmp = _delta_x*((_images.itemAt(2)).div(_images.itemAt(0), 0.0f));
//    dy_tmp = _delta_y*((_images.itemAt(3)).div(_images.itemAt(0), 0.0f));

    norm = (dx_tmp.power(2)
	    + dy_tmp.power(2)
	    + 1.0f).sqrt();
    
    _intermediate_images.at(2) = dx_tmp/norm;
    _intermediate_images.at(3) = dy_tmp/norm;
    _intermediate_images.at(4) = 1.0f/norm;

}

VISImage<float> SphericalRange2::normal_dx()
{
    int w, h;
    VISImage<float> r(w = _images.itemAt(0).width(), 
		       h = _images.itemAt(0).height());

    VISImage<float> range;

    range = (_images.itemAt(0));

    float forward, back;
    
    int i, j;
    for (j = 0; j < (h); j++)
	for (i = 1; i < (w - 1); i++)
	    {
		forward = 
		    range.itemAt(i + 1, j) - range.itemAt(i, j);

		back = 
		    range.itemAt(i, j) - range.itemAt(i - 1, j);

//		if ((forward*back) < 0)
//		    r.at(i, j) = (forward + back)/2.0f;
//		else
		    if (fabs(forward) > fabs(back))
			r.at(i, j) = forward;
		    else 
			r.at(i, j) = back;
		    
	    }

    for (j = 0; j < (h); j++)    
	{
	    r.at(0, j) = range.itemAt(1, j) 
		- range.itemAt(0, j);

	    r.at(w - 1, j) = range.itemAt(w - 1, j) 
		- range.itemAt(w - 2, j);
	}

    return(r);
}

VISImage<float> SphericalRange2::normal_dy()
{
    int w, h;
    VISImage<float> r(w = _images.itemAt(0).width(), 
		       h = _images.itemAt(0).height());

    VISImage<float> range;

    range = _images.itemAt(0);

    float forward, back;
    
    int i, j;
    for (j = 1; j < (h - 1); j++)
	for (i = 0; i < (w); i++)
	    {
		forward = 
		    range.itemAt(i, j + 1) - range.itemAt(i, j);

		back = 
		    range.itemAt(i, j) - range.itemAt(i, j - 1);

//		if ((forward*back) < 0)
//		    r.at(i, j) = (forward + back)/2.0f;
//		else
		if (fabs(forward) > fabs(back))
		    r.at(i, j) = forward;
		else 
		    r.at(i, j) = back;
	    }

    for (j = 0; j < (w); j++)    
	{
	    r.at(j, 0) = range.itemAt(j, 1) 
		- range.itemAt(j, 0);

	    r.at(j, h - 1) = range.itemAt(j, h - 1) 
		- range.itemAt(j, h - 2);
	}

    return(r);
}
