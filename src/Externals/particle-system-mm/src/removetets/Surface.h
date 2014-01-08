#ifndef __SURFACE_H__
#define __SURFACE_H__

#include "mtxlib.h"
#include "defines.h"
#include "SurfaceParameters.h"


/************************************************************************/
//                     The Surface Base Class                           //
/************************************************************************/
class Surface
{
public:
  Surface();
  virtual ~Surface(){};

  //---------------------------------------------//
  //           virtual functions                 //
  //---------------------------------------------//
     
  //--------------------------------------------------------------------
  // functions for updating surface values at particle locations
  virtual bool computeSurfacePointParams(const vector_type &pos,
                                         SurfacePointParams &params,
                                         bool computeGradient=true,
                                         bool computeHessian=false) 
    const = 0;

  virtual void bounds(vector_type &s, vector_type &e) const
  { s = _d_start; e = _d_end; };

  //--------------------------------------------------------------------
  // a render function that can be filled in if there are specific
  //   characteristics of the surface to be rendered
  virtual void render() const {};

  //--------------------------------------------------------------------
  // a select function that can be filled in to specifically select
  //   something about the surface
  virtual void select(int x, int y) {};

  //--------------------------------------------------------------------
  // for debugging
  virtual void minAndMaxCurvature(float &min, float &max) const {};


  //---------------------------------------------//
  //           surface functions                 //
  //---------------------------------------------//

  // the isovalue for the implicit function
  virtual float isovalue() const { return _isovalue; };
  virtual void  isovalue(float i) { _isovalue = i; };

  // compute the curvature at a particle based on gordon's
  //   vis03 paper
  float computeCurvature(SurfacePointParams &params) const;

  // the same as above, only signed
  float computeSignedCurvature(SurfacePointParams &params) const;

  inline float min_isovalue() const { return _min_isovalue; };
  inline float max_isovalue() const { return _max_isovalue; };

  virtual bool projectionVector() const { return false; };
  virtual void getProjectionVector(const vector_type &pos, 
                                   vector_type &proj) const 
  { proj = 0.0; };

  virtual float maxSF() const { return MAX_VALUE; };

  inline bool inBounds(const vector_type &pos) const
  { if ( (pos(0) > _d_end(0)) ||
         (pos(1) > _d_end(1)) ||
         (pos(2) > _d_end(2)) ||
         (pos(0) < _d_start(0)) ||
         (pos(1) < _d_start(1)) ||
         (pos(2) < _d_start(2)) ) return false;
  else return true;
  };

protected:
  float _isovalue;
  vector_type _d_start, _d_end;

  float _min_isovalue, _max_isovalue;

  // overloaded helper functions to deal with 2D or 3D
  float computeCurvatureMagnitude(const matrix<2,2> &G,
                                  SurfacePointParams &params) const;
  float computeCurvatureMagnitude(const matrix<3,3> &G,
                                  SurfacePointParams &params) const;

  float computeSignedCurvatureMagnitude(const matrix<2,2> &G,
                                        SurfacePointParams &params) const;
  float computeSignedCurvatureMagnitude(const matrix<3,3> &G,
                                        SurfacePointParams &params) const;

  void computeEigen(const matrix<2,2> &G,
                    SurfacePointParams &params) const;
  void computeEigen(const matrix<3,3> &G,
                    SurfacePointParams &params) const;

};

#endif // __SURFACE_H__
