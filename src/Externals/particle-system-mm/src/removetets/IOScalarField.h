//----------------------------------------------------------------------//
// FILE        : Intersection.h                                                
// DESCRIPTION : This surface reads in two distance field volumes and the 
//               first derivatives, along with a sizing field.
//----------------------------------------------------------------------//

#ifndef __IO_SCALAR_FIELD_H__
#define __IO_SCALAR_FIELD_H__

#include "Surface.h"
#include "mtxlib.h"
#include "defines.h"
//#include "svector.h"

#define EPSILON_SQRD 0.00000001

class IOScalarField : public Surface
{
public:
  IOScalarField(const char *filename,
                int kernel_type=APPROXIMATING,
                int main_indicator=0);
  ~IOScalarField();

  enum {INTERPOLATING, APPROXIMATING};

  //--------------------------------------------------------------------
  // evaluate the surface parameters at this position
  bool computeSurfacePointParams(const vector_type &pos, 
                                 SurfacePointParams &params,
                                 bool computeGradient=true,
                                 bool computeHessian=false) const;

  inline void setMainIndicatorFunction(int i)
  { if ( i >=  _num_indicators )
  { std::cout << "PROBLEM: main indicator number too big!" << std::endl;
  exit( 1 ); }
  _main_indicator = i; };

  inline int numIndicators() const 
  { return _num_indicators; };

  void printFunctionValues(const vector_type &pos);

private:
  Surface **_indicators;
  int _num_indicators, _main_indicator;

  float interp_max(const float* const f, int num_values) const;
  void  interp_max_D(const float* const f, const vec<3>* const fx, 
                     int num_values,  vec<3> &D) const;

  float sign_function(float a) const;
  void sign_function_D(float a, const vec<3> &aD, vec<3> &D) const;

  float power_function(float a, int b) const;
};


#endif // __IO_SCALAR_FIELD_H__
