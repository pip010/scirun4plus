#ifndef __RECONSTRUCTION_FIELD_H__
#define __RECONSTRUCTION_FIELD_H__

#include <cstdlib>
#include "ScalarField.h"
#include "mtxlib.h"
#include "multiDarrays.h"
#include "Kernel.h"
#include <iostream>

typedef float reconfield_datatype;
/************************************************************************/
// NOTE: this is written for 3D only right now!
/************************************************************************/

template <int dim> class ReconstructionField : public ScalarField<dim>
{
public:
  ReconstructionField(const char *df_filename,
                      int kernel_type=APPROXIMATING,
                      float scale_x=1.0, 
                      float scale_y=1.0, 
                      float scale_z=1.0 );
  ~ReconstructionField();

  void computeScalarFieldParams(const vec<dim> &pos, 
                                ScalarFieldParams<dim>&params,
                                bool computeHessian=false) const;

  virtual bool checkBounds(const vec<dim> &pos) const;
  virtual vec<dim> boundingBoxLow() const
  { return(_d_start); };
  virtual vec<dim> boundingBoxHigh() const
  { return(_d_end); };

  enum {INTERPOLATING, APPROXIMATING};

  int xdim () const {return(_xdim);}
  int ydim () const {return(_ydim);}
  int zdim () const {return(_zdim);}

  // protected:
protected:
  int _xdim, _ydim, _zdim;
  vec<dim> _d_start, _d_end;
  float _scale[3];
  array3D<reconfield_datatype> _field;
  Kernel *_kernel, *_kernelD, *_kernelDD;

  void readFile(const char *df_filename);

};

//------------------------------------------------------------------------
// Function    : Constructor and Destructor
// Description : initializing constructor
//------------------------------------------------------------------------
template <int dim>
ReconstructionField<dim>::ReconstructionField( const char *df_filename,
                                               int kernel_type,
                                               float scale_x, 
                                               float scale_y, 
                                               float scale_z )
{
  _scale[0] = scale_x;
  _scale[1] = scale_y;
  _scale[2] = scale_z;
  
  readFile( df_filename );

  // create the reconstruction kernels
  switch ( kernel_type )
  {
    case APPROXIMATING:
    _kernel  = new CubicBSpline();
    _kernelD  = new CubicBSplineD();
    _kernelDD = new CubicBSplineDD();
    break;

    case INTERPOLATING:
    _kernel = new CatmullRom();
    _kernelD  = new CatmullRomD();
    _kernelDD = new CatmullRomDD();
    break;

    default:
    std::cout << "No such kernel type" << std::endl;
    exit( 1 );
    break;
  }
}

template <int dim>
ReconstructionField<dim>::~ReconstructionField()
{
  delete _kernel;
  delete _kernelD;
  delete _kernelDD;
}


#endif // __RECONSTRUCTION_FIELD_H__
