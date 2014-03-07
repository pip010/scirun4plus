//
// modelnode.h 
//


#ifndef VIS_diffuse_h
#define VIS_diffuse_h

#include "image/image.h"
#include "util/array.h"

class Diffuse
{
//  protected:
  public:
    VISArray< VISImage<float> > _images;
    VISArray< VISImage<float> > _change;
    VISImage<float> _conductance_dx, _conductance_dy;
    int _num_intermediate_images;
    VISArray< VISImage<float> > _intermediate_images;
    VISArray< float > _weights;
    void calculateIntermediateImages();

  public:
    
    const VISImage<float>& image(unsigned i) { return(_images.itemAt(i));}
    void image(const VISImage<float>& im, unsigned i) 
    { _images.at(i) = im;}
    void parameters(VISArray<float> a) {_weights = a;}
    const VISArray<float>& parameters() {return(_weights);}

    virtual void calculateIntermediateFeatures()
    {
	for (int i = 0; i < _images.n(); i++)
	    _intermediate_images.at(i) = _images.itemAt(i);
    }
    
    void calculateConductance();
    void iterate()
	{
	    calculateIntermediateFeatures();
	    calculateConductance();
	    calculateChange();
	    update();
	}

    void iterate(int num_iterations)
	{
	    calculateIntermediateFeatures();
	    calculateConductance();
	    for (int i = 0; i < num_iterations; i++)
		{
		    calculateChange();
		    update();
		}
	}

    void calculateChange();
    void update();
    
};


class Aniso: public Diffuse
{
  protected:

  public:
    virtual void calculateIntermediateFeatures()
    {
	_intermediate_images.at(0) = _images.itemAt(0);
    }
    
    Aniso(VISImage<float> im)
    {
	image(im, 0);
	_num_intermediate_images = 1;
    }

    void conductance(float k)
    {
	_weights.at(0) = k*k;
    }

};

class SphericalRange: public Diffuse
{
  protected:
    
    float _delta_x, _delta_y;

  public:
    virtual void calculateIntermediateFeatures();
    
    SphericalRange(VISImage<float> range, VISImage<float> intensity, 
		   float dx, float dy)
    {

	int w, h, i;
	VISImage<float> im_dx, im_dy;
	
	image(range, 0);
	image(intensity, 1);


	im_dx = range.dx();
	im_dy = range.dy();
	w = range.width();
	h = range.height();


// need to set the boundaries of the derivatives
	for (i = 0; i < h; i++)
	    {
		im_dx.at(0, i) = range.itemAt(1, i) 
		    - range.itemAt(0, i);
		im_dx.at(w - 1, i) = range.itemAt(w - 1, i) 
		    - range.itemAt(w - 2, i);
	    }

	for (i = 0; i < w; i++)
	    {
		im_dy.at(i, 0) = range.itemAt(i, 1) 
		    - range.itemAt(i, 0);
		im_dy.at(i, h - 1) = range.itemAt(i, h - 1) 
		    - range.itemAt(i, h - 2);
	    }

	image(im_dx, 2);
	image(im_dy, 3);

	_delta_x = dx;
	_delta_y = dy;

	_num_intermediate_images = 5;

    }

    void conductance(float k_range, float k_intensity, float k_normal)
    {
	_weights.at(0) = k_range*k_range;
	_weights.at(1) = k_intensity*k_intensity;
	_weights.at(2) = _weights.at(3) = _weights.at(4) = k_normal
	    *k_normal;
    }

};


class SphericalRange2: public Diffuse
{
  protected:
    
    float _delta_x, _delta_y;
    VISImage<float> _edges;


  public:
    virtual void calculateIntermediateFeatures();
    
    SphericalRange2(VISImage<float> range, VISImage<float> intensity, 
		   float dx, float dy)
    {

	int w, h, i;
	VISImage<float> im_dx, im_dy;
	
	image(range, 0);
	image(intensity, 1);

	w = range.width();
	h = range.height();

	_delta_x = dx;
	_delta_y = dy;

	_num_intermediate_images = 5;
    }

    void conductance(float k_range, float k_intensity, float k_normal)
    {
	_weights.at(0) = k_range*k_range;
	_weights.at(1) = k_intensity*k_intensity;
	_weights.at(2) = _weights.at(3) = _weights.at(4) = k_normal
	    *k_normal;
    }

    VISImage<float> normal_dx();
    VISImage<float> normal_dy();

};




#endif








