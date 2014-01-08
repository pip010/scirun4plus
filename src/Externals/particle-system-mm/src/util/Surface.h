#ifndef __SURFACE_H__
#define __SURFACE_H__

#include <utilExports.h>

#include <mtxlib.h>
#include <SurfaceParameters.h>
#include <constants.h>

class MM_util_SHARE Surface
{
public:
  Surface() {}
  virtual ~Surface() {}

  //---------------------------------------------//
  //           virtual functions                 //
  //---------------------------------------------//
     
  //--------------------------------------------------------------------
  // functions for updating surface values at particle locations
  virtual bool computeSurfacePointParams(const vec<3> &pos,
                                         SurfacePointParams &params,
                                         bool computeGradient=false,
                                         bool computeHessian=false) = 0;

  float computeCurvature(SurfacePointParams &params) const;

  // the isovalue for the implicit function
  virtual float isovalue() const { return _isovalue; };
  virtual void  isovalue(float i) { _isovalue = i; };

  inline void domain(vec<3> &start, vec<3> &end)
  { start = _start; end = _end; };

  bool inBounds(const vec<3> &pos) const;

  virtual int xdim() const { return 0; };
  virtual int ydim() const { return 0; };
  virtual int zdim() const { return 0; };
  
  virtual DataCenter xcenter() const { return UNKNOWN; }
  virtual DataCenter ycenter() const { return UNKNOWN; }
  virtual DataCenter zcenter() const { return UNKNOWN; }

protected:
  float _isovalue;
  vec<3> _start, _end;

  float computeCurvatureMagnitude(const matrix<3,3> &G,
                                  SurfacePointParams &params) const;

};

#endif // __SURFACE_H__
