//----------------------------------------------------------------------//
// FILE        : Intersection.h                                                
// DESCRIPTION : This surface reads in two distance field volumes and the 
//               first derivatives, along with a sizing field.
//----------------------------------------------------------------------//

#ifndef __INTERSECTION_H__
#define __INTERSECTION_H__

#include <system/systemExports.h>

#include <system/domain/Surface.h>
#include <system/domain/IOScalarField.h>
#include <features/mtxlib.h>
#include <system/defines.h>

namespace particle_sys 
{
  class System_SHARE Intersection : public Surface
  {
  public:
    Intersection(IOScalarField *io, int num_indicators,
                 int *indicators);
    ~Intersection();

    //--------------------------------------------------------------------
    // evaluate the surface parameters at this position
    bool computeSurfacePointParams(const vector_type &pos, 
                                   SurfacePointParams &params,
                                   bool computeHessian) const;

    bool projectOntoSurface(DynamicSurfacePoint *point,
                            vector_type &proj) const; 

    void projectMotion(const vector_type &pos,
                       const vector_type &n,
                       vector_type &m ) const;

    inline float maxSF() const
    { return _io_field->maxSF(); };

  private:
    IOScalarField *_io_field;
    int _num_indicators;
    int *_indicators;

  };

} // namespace particle_sys




#endif // __INTERSECTION_H__
