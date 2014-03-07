#ifndef __SIMPLE_SURFACES_H__
#define __SIMPLE_SURFACES_H__

#include <system/systemExports.h>

#include <system/domain/Surface.h>
#include <features/mtxlib.h>
#include <system/defines.h>

namespace particle_sys 
{

/************************************************************************/
//   These simple surfaces are always centered at the origin. The       //
//              domain is defined by the radius!                        //
/************************************************************************/


/************************************************************************/
//                            SPHERE                                    //
/************************************************************************/
class System_SHARE SimpleSphere : public Surface
{
public:
  SimpleSphere(float radius);
  ~SimpleSphere(){};

  bool computeSurfacePointParams(const vector_type &pos, 
                                 SurfacePointParams &params,
                                 bool computeHessian) const;
      
private:
  float _radius;
};


/************************************************************************/
//                            TORUS                                     //
/************************************************************************/
class SimpleTorus : public Surface
{
public:
  SimpleTorus(float r, float R); 
  ~SimpleTorus(){};

  bool computeSurfacePointParams(const vector_type &pos, 
                                 SurfacePointParams &params,
                                 bool computeHessian) const;
      
private:
  float _r, _R;
};


/************************************************************************/
//                        HOLLOW CUBE                                   //
/************************************************************************/
class SimpleQuarticCube : public Surface
{
public:
  SimpleQuarticCube(float radius);
  ~SimpleQuarticCube(){};

  bool computeSurfacePointParams(const vector_type &pos, 
                                 SurfacePointParams &params,
                                 bool computeHessian) const;

  void minAndMaxCurvature(float &min, float &max) const 
  { min = -1.0; max = 5.0; }

private:
  float _radius;
};


/************************************************************************/
//                            ELLIPSE                                   //
/************************************************************************/
class SimpleEllipse : public Surface
{
public:
  SimpleEllipse(float a, float b, float c, float d, float e, float f); 
  ~SimpleEllipse(){};

  bool computeSurfacePointParams(const vector_type &pos, 
                                 SurfacePointParams &params,
                                 bool computeHessian) const;
      
private:
  float _a, _b, _c, _d, _e, _f;
};


/************************************************************************/
//                          CONSTANT                                    //
/************************************************************************/
class Constant : public Surface
{
public:
  Constant() {}; 
  ~Constant(){};

  bool computeSurfacePointParams(const vector_type &pos, 
                                 SurfacePointParams &params, 
                                 bool computeHessian) const;
      
private:
};


} // namespace particle_sys

#endif // __SIMPLE_SURFACES_H__
