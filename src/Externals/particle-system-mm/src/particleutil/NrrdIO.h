#ifndef __PARTICLEUTIL_NRRDIO_H__
#define __PARTICLEUTIL_NRRDIO_H__

#include <particleutilExports.h>
#include <constants.h>

#include <multiDarrays.h>

#include <vector>

class NrrdWrapperRead;
class NrrdWrapperWrite;

class Particleutil_SHARE NrrdIO {
public:
  static void read_nrrd(const char *filename, array3D<float>& field,
                        float& sx, float& sy, float& sz,
                        int& xdim, int& ydim, int& zdim);

  static void read_nrrd(const char *filename, array3D<float>& field,
                        float& sx, float& sy, float& sz,
                        int& xdim, int& ydim, int& zdim,
                        DataCenter& center_x, DataCenter& center_y, DataCenter& center_z);

  // get size from field
  static void write_nrrd(const char* filename, const array3D<float>& field,
                         const float sx, const float sy, const float sz,
                         const DataCenter center_x = UNKNOWN, const DataCenter center_y = UNKNOWN, const DataCenter center_z = UNKNOWN);

  static void write_nrrd(const char* filename, const array3D<bool>& field,
                         const float sx, const float sy, const float sz,
                         const DataCenter center_x = UNKNOWN, const DataCenter center_y = UNKNOWN, const DataCenter center_z = UNKNOWN);        
                         
private:
  static void get_header(NrrdWrapperRead& nw, const unsigned int dim,
                         float& sx, float& sy, float& sz,
                         int& xdim, int& ydim, int& zdim,
                         DataCenter& center_x, DataCenter& center_y, DataCenter& center_z);

  static void get_data(NrrdWrapperRead& nw, array3D<float>& field,
                       int& xdim, int& ydim, int& zdim);

  static void set_header(NrrdWrapperWrite& nw, int nrrdType, const unsigned int dim,
                         const float sx, const float sy, const float sz,
                         const int xdim, const int ydim, const int zdim,
                         const DataCenter center_x, const DataCenter center_y, const DataCenter center_z);
};

#endif // __PARTICLEUTIL_NRRDIO_H__
