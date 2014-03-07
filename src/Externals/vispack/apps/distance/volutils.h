#include <volume.h>
#include <matrix.h>
//#include <vector.h>
#include <teem/nrrd.h>
#include <string>

typedef VISVolume<VISMatrix>  VolumeTensor;
typedef VISVolume<VISVector> VolumeVector;
typedef VISVolume<float> VolumeScalar;

VolumeTensor readTensorVolumeFile(char* filename);
boolean writeTensorVolumeFile(char* filename, const VolumeTensor& vol_tensor);
boolean writeTensorVolumeFile(char* filename, const VolumeTensor& vol_tensor, const VolumeScalar& vol_scalar);
VolumeVector readVectorVolumeFile(char* filename);
boolean writeVectorVolumeFile(char* filename, const VolumeVector& vol_vector);
VolumeScalar readScalarVolumeFile(char* filename);
boolean writeScalarVolumeFile(char* filename, const VolumeScalar& vol_scalar);
VolumeVector volumeFirstDerivatives(const VolumeScalar& vol);
VolumeVector volumeFirstDerivativesUpWind(const VolumeScalar& vol);
VolumeTensor volumeSecondDerivatives(const VolumeScalar& vol);
//VISVector vectorVolumeInterpolate(const VolumeVector& vol, const VISVector& vector);
