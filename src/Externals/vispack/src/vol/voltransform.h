 
//
// volxforminterp.h 
//
// This class takes two volumes and calculates the transformation
// that approximately maps one into the other.
//

#include "util/geometry.h"
#include "matrix.h"
#include "volume.h"
#include "nr/nr.h"
extern void tred2(float** a, int n, float* d, float* e);
extern void tqli(float* d, float* e, int n, float** z);


class GfxVolTransformInterp
{

  protected:
    GfxVolume<float> _vol0;
    GfxVolume<float> _vol1;

    boolean	_vol0_empty;
    boolean	_vol1_empty;
    boolean	_dirty;

// transformation info for the two volumes
    Quaternion  _quat0;
    Quaternion  _quat1;
    Vector3         _trans0;
    Vector3         _trans1;
    Vector3         _scale0;
    Vector3         _scale1;

    void moments();
		    
    
  public:
    GfxVolTransformInterp(GfxVolume<float>, GfxVolume<float>);
    GfxVolTransformInterp();

    void volumes(GfxVolume<float>, GfxVolume<float>);
    void volume0(GfxVolume<float>);
    void volume1(GfxVolume<float>);

    GfxTransform interp(float);

};
