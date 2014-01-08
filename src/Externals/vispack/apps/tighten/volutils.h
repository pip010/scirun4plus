#ifndef VOLUTILS_H
#define VOLUTILS_H 1

#include <volume.h>
#include <matrix.h>
//#include <vector.h>
#include <teem/nrrd.h>
#include <string>

typedef VISVolume<VISMatrix>  VolumeTensor;
typedef VISVolume<VISVector>  VolumeVector;
typedef VISVolume<float>      VolumeScalar;

class VolScale {
  public:
    VolScale() :
      spacing_x(1.0),
      spacing_y(1.0),
      spacing_z(1.0),
      center_x(nrrdCenterUnknown),
      center_y(nrrdCenterUnknown),
      center_z(nrrdCenterUnknown) {}
    
    double spacing_x;
    double spacing_y;
    double spacing_z;
    int center_x;
    int center_y;
    int center_z;
};


VolumeTensor readTensorVolumeFile(char* filename);
boolean writeTensorVolumeFile(char* filename, const VolumeTensor& vol_tensor);
boolean writeTensorVolumeFile(char* filename, const VolumeTensor& vol_tensor, const VolumeScalar& vol_scalar);
VolumeVector readVectorVolumeFile(char* filename);
boolean writeVectorVolumeFile(char* filename, const VolumeVector& vol_vector);

VolumeScalar readScalarVolumeFile(char* filename);
boolean writeScalarVolumeFile(char* filename, const VolumeScalar& vol_scalar);

VolumeScalar readScalarVolumeFile(char* filename,VolScale& scale);
boolean writeScalarVolumeFile(char* filename, const VolumeScalar& vol_scalar, VolScale& scale);

VolumeVector volumeFirstDerivatives(const VolumeScalar& vol);
VolumeVector volumeFirstDerivativesUpWind(const VolumeScalar& vol);
VolumeTensor volumeSecondDerivatives(const VolumeScalar& vol);
//VISVector vectorVolumeInterpolate(const VolumeVector& vol, const VISVector& vector);

#endif
