//----------------------------------------------------------------------//
// FILE        : IOScalarField.h                                                
// DESCRIPTION : This function takes a set of indicator functions
//               (ReconstructionField) and computes the
//               inside/outside scalar field
//
//               Assumes that the indicator functions are [0.5,1.0] for
//               material and [0.0, 0.5) for not material. 
//----------------------------------------------------------------------//

#ifndef __IO_SCALAR_FIELD_H__
#define __IO_SCALAR_FIELD_H__

#include <utilExports.h>

#include <Surface.h>
#include <ScalarField.h>
#include <mtxlib.h>

#define EPSILON_SQRD 1.0e-8

class MM_util_SHARE IOScalarField : public Surface
{
public:
  IOScalarField(const char *filename,
                int kernel_type=APPROXIMATING,
                int main_indicator=0);
  ~IOScalarField();

  enum {INTERPOLATING, APPROXIMATING};

  inline void setMainIndicatorFunction(int i)
  { if ( i >=  _num_indicators )
    { std::cout << "PROBLEM: main indicator number too big!" << std::endl;
      exit( 1 ); }
    _main_indicator = i; };
    
  bool computeSurfacePointParams(const vec<3> &pos, 
                                 SurfacePointParams &params,
                                 bool computeGradient=false,
                                 bool computeHessian=false);

  int xdim() const { return _indicators[_main_indicator]->xdim(); }
  int ydim() const { return _indicators[_main_indicator]->ydim(); }
  int zdim() const { return _indicators[_main_indicator]->zdim(); }

  DataCenter xcenter() const { return _indicators[_main_indicator]->xcenter(); }
  DataCenter ycenter() const { return _indicators[_main_indicator]->ycenter(); }
  DataCenter zcenter() const { return _indicators[_main_indicator]->zcenter(); }
  
protected:
  Surface **_indicators;
  int _num_indicators, _main_indicator;
  float _epsilon_sqrd;

  float interp_max(const float* const f, int num_values) const;
  void  interp_max_D(const float* const f, const vec<3>* const fx, 
                     int num_values,  vec<3> &D) const;
  void interp_max_DD(const float* const f, const vec<3>* const fx,
                     const matrix<3,3>* const fxx, int num_values,
                     matrix<3,3> &DD) const;
  void interp_max_D_DD(const float* const f, const vec<3>* const fx,
                       const matrix<3,3>* const fxx, int num_values,
                       vec<3> &D, matrix<3,3> &DD) const;

  float sign_function(float a) const;
  void sign_function_D(float a, const vec<3> &aD, vec<3> &D) const;
  void sign_function_DD(float a, const vec<3> &aD,
                        const matrix<3,3> &aDD, matrix<3,3> &DD) const;
  void sign_function_D_DD(float a, const vec<3> &aD,
                          const matrix<3,3> &aDD, 
                          vec<3> &D, matrix<3,3> &DD) const;

  float power_function(float a, int b) const;
};


#endif // __IO_SCALAR_FIELD_H__
