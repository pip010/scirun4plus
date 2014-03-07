#include <NrrdIO.h>

#include <teem/nrrd.h>

#include <iostream>
#include <vector>

class NrrdWrapperRead {
public:
  NrrdWrapperRead() : nrrd( nrrdNew() ) {}
  ~NrrdWrapperRead() { nrrdNuke(nrrd); }

  Nrrd* nrrd;
};

class NrrdWrapperWrite : public NrrdWrapperRead {
public:
  NrrdWrapperWrite() : NrrdWrapperRead(), nio(nrrdIoStateNew()) {}
  ~NrrdWrapperWrite() { nrrdIoStateNix(nio); }

  NrrdIoState* nio;
};

void NrrdIO::read_nrrd(const char *filename, array3D<float>& field,
                       float& sx, float& sy, float& sz,
                       int& xdim, int& ydim, int& zdim)
{
  DataCenter center_x, center_y, center_z;
  read_nrrd(filename, field, sx, sy, sz, xdim, ydim, zdim, center_x, center_y, center_z);
}

void NrrdIO::read_nrrd(const char *filename, array3D<float>& field,
                       float& sx, float& sy, float& sz,
                       int& xdim, int& ydim, int& zdim,
                       DataCenter& center_x, DataCenter& center_y, DataCenter& center_z)
{
  NrrdWrapperRead wrapper;
  Nrrd* nrrd = wrapper.nrrd;
  if ( nrrdLoad(nrrd, airStrdup(filename), 0) )
  {
    char *err = biffGetDone(NRRD);
    std::cerr << "NRRD error: nrrdLoad failed with error " << err << std::endl;
    free(err);
    return;
  }
  
  unsigned int dim = nrrd->dim;
  if (dim != 3)
  {
    std::cerr << "NRRD error: need 3D nrrd" << std::endl;
    // abort
    return;
  }
  
  // this is specific to BioMesh3D
  if (nrrd->type != nrrdTypeFloat) {
    std::cerr << "NRRD error: nrrd data must be of type float" << std::endl;
  }

  get_header(wrapper, dim, sx, sy, sz, xdim, ydim, zdim, center_x, center_y, center_z);
  get_data(wrapper, field, xdim, ydim, zdim);
}


void NrrdIO::write_nrrd(const char *filename, const array3D<float>& field,
                        const float sx, const float sy, const float sz,
                        const DataCenter center_x, const DataCenter center_y, const DataCenter center_z)
{
  NrrdWrapperWrite wrapper;
  Nrrd* nrrd = wrapper.nrrd;
  NrrdIoState* nio = wrapper.nio;
  set_header(wrapper, nrrdTypeFloat, 3, sx, sy, sz, field._idim, field._jdim, field._kdim, center_x, center_y, center_z);

  // write out data values
  float* dataptr = static_cast<float*>(nrrd->data);
  
  for (int k = 0; k < field._kdim; ++k)
  {
    for (int j = 0; j < field._jdim; ++j)
    {
      for (int i = 0; i < field._idim; ++i)
      {
        *dataptr = field(i,j,k);
        ++dataptr;
      }
    }
  }

  if (nrrdSave( airStrdup(filename), nrrd, nio)) 
  {
    char *err = biffGetDone(NRRD);
    std::cerr << "NRRD Error: nrrdSave for '" << filename << "' failed with error " << err << std::endl;
    free(err);
    
    return;
  }
}

void NrrdIO::write_nrrd(const char *filename, const array3D<bool>& field,
                        const float sx, const float sy, const float sz,
                        const DataCenter center_x, const DataCenter center_y, const DataCenter center_z)
{
  NrrdWrapperWrite wrapper;
  Nrrd* nrrd = wrapper.nrrd;
  NrrdIoState* nio = wrapper.nio;
  set_header(wrapper, nrrdTypeFloat, 3, sx, sy, sz, field._idim, field._jdim, field._kdim, center_x, center_y, center_z);
 
  // write out data values
  float* dataptr = static_cast<float*>(nrrd->data);
  
  for (int k = 0; k < field._kdim; ++k)
  {
    for (int j = 0; j < field._jdim; ++j)
    {
      for (int i = 0; i < field._idim; ++i)
      {
		// nrrd only supports c types
        *dataptr = static_cast<float>(field(i,j,k));
        ++dataptr;
      }
    }
  }

  if (nrrdSave( airStrdup(filename), nrrd, nio)) 
  {
    char *err = biffGetDone(NRRD);
    std::cerr << "NRRD Error: nrrdSave for '" << filename << "' failed with error " << err << std::endl;
    free(err);
    
    return;
  }
}


void NrrdIO::get_data(NrrdWrapperRead& nw, array3D<float>& field,
                      int& xdim, int& ydim, int& zdim)
{
  Nrrd* nrrd = nw.nrrd;
  // resize the field
  field.resize(xdim, ydim, zdim);
  // read in data values
  float* dataptr = static_cast<float*>(nrrd->data);

  for (int k = 0; k < zdim; ++k)
  {
    for (int j = 0; j < ydim; ++j)
    {
      for (int i = 0; i < xdim; ++i)
      {
        field(i,j,k) = *dataptr;
        ++dataptr;
      }
    }
  }
}

void NrrdIO::get_header(NrrdWrapperRead& nw, const unsigned int dim,
                        float& sx, float& sy, float& sz,
                        int& xdim, int& ydim, int& zdim,
                        DataCenter& center_x, DataCenter& center_y, DataCenter& center_z)
{
  Nrrd* nrrd = nw.nrrd;
  std::vector<int> size(dim);
  std::vector<double> spacing(dim);
   
  for (unsigned int i = 0; i < dim; i++)
  {
    size[i] = static_cast<int>(nrrd->axis[i].size);
    if ( airExists(nrrd->axis[i].spaceDirection[i]) )
    {
      // works, but in principle, should get vector magnitude
      spacing[i] = nrrd->axis[i].spaceDirection[i];
    }
    else
    {
      std::cerr << "NRRD error: invalid nrrd, spaceDirection must be set" << std::endl;
      return;
    }
  }

  if (nrrdCenterCell == nrrd->axis[0].center)
    center_x = ELEMENT;
  else if (nrrdCenterNode == nrrd->axis[0].center)
    center_x = NODE;
  else
    center_x = UNKNOWN;

  if (nrrdCenterCell == nrrd->axis[1].center)
    center_y = ELEMENT;
  else if (nrrdCenterNode == nrrd->axis[1].center)
    center_y = NODE;
  else
    center_y = UNKNOWN;

  if (nrrdCenterCell == nrrd->axis[2].center)
    center_z = ELEMENT;
  else if (nrrdCenterNode == nrrd->axis[2].center)
    center_z = NODE;
  else
    center_z = UNKNOWN;

  sx = static_cast<float>(spacing[0]);
  sy = static_cast<float>(spacing[1]);
  sz = static_cast<float>(spacing[2]);

  xdim = size[0];
  ydim = size[1];
  zdim = size[2];
}


void NrrdIO::set_header(NrrdWrapperWrite& nw, int nrrdType, const unsigned int dim,
                        const float sx, const float sy, const float sz,
                        const int xdim, const int ydim, const int zdim,
                        const DataCenter center_x, const DataCenter center_y, const DataCenter center_z)
{
  Nrrd* nrrd = nw.nrrd;
  NrrdIoState* nio = nw.nio;
  size_t nrrddims[3];
  nrrddims[0] = static_cast<size_t>(xdim);
  nrrddims[1] = static_cast<size_t>(ydim);
  nrrddims[2] = static_cast<size_t>(zdim);
  if ( nrrdAlloc_nva(nrrd, nrrdType, dim, nrrddims) != 0 )
  {
    char *err = biffGetDone(NRRD);
    std::cerr << "NRRD error: nrrdAlloc_nva failed with error " << err << std::endl;
    free(err);
    exit(1);
  }

  nrrd->spaceDim = dim;

  nrrd->spaceOrigin[0] = 0.0;
  nrrd->spaceOrigin[1] = 0.0;
  nrrd->spaceOrigin[2] = 0.0;

  nrrd->axis[0].size = xdim;
  nrrd->axis[1].size = ydim;
  nrrd->axis[2].size = zdim;
  
  nrrd->axis[0].spaceDirection[0] = sx;
  nrrd->axis[0].spaceDirection[1] = 0.0;
  nrrd->axis[0].spaceDirection[2] = 0.0;

  nrrd->axis[1].spaceDirection[0] = 0.0;
  nrrd->axis[1].spaceDirection[1] = sy;
  nrrd->axis[1].spaceDirection[2] = 0.0;

  nrrd->axis[2].spaceDirection[0] = 0.0;
  nrrd->axis[2].spaceDirection[1] = 0.0;
  nrrd->axis[2].spaceDirection[2] = sz;


  if (ELEMENT == center_x)
    nrrd->axis[0].center = nrrdCenterCell;
  else if (NODE == center_x)
    nrrd->axis[0].center = nrrdCenterNode;
  else
    nrrd->axis[0].center = nrrdCenterUnknown;

  if (ELEMENT == center_y)
    nrrd->axis[1].center = nrrdCenterCell;
  else if (NODE == center_y)
    nrrd->axis[1].center = nrrdCenterNode;
  else
    nrrd->axis[1].center = nrrdCenterUnknown;

  if (ELEMENT == center_z)
    nrrd->axis[2].center = nrrdCenterCell;
  else if (NODE == center_z)
    nrrd->axis[2].center = nrrdCenterNode;
  else
    nrrd->axis[2].center = nrrdCenterUnknown;

  // set encoding to be gzip
  // TODO: check if gzip is available?  If not, use raw
  nio->encoding = nrrdEncodingArray[nrrdEncodingTypeGzip];
  // set format to be nrrd
  nio->format = nrrdFormatArray[nrrdFormatTypeNRRD];
  // set endian to be endian of machine
  nio->endian = AIR_ENDIAN;
  nio->zlibLevel = 6;
}
