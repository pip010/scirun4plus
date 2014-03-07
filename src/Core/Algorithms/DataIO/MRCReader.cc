/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2009 Scientific Computing and Imaging Institute,
 University of Utah.
 
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the "Software"),
 to deal in the Software without restriction, including without limitation
 the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
 */
 
/*
 * MRC file format reader: http://www2.mrc-lmb.cam.ac.uk/image2000.html
 * Implementation follows EMAN2 and Chimera.
 */

#include <boost/shared_array.hpp>
#include <Core/Algorithms/DataIO/MRCReader.h>

#include <Core/ImportExport/Nrrd/NrrdIEPlugin.h>
#include <Core/Thread/Guard.h>
#include <Core/Util/Endian.h>
#include <Core/Datatypes/NrrdData.h>
#include <sci_debug.h>

#include <iostream>
#include <fstream>
#include <sstream>

namespace SCIRunAlgo {



MRCReader::MRCReader(ProgressReporter* pr)
  : AlgoLibrary(pr),
    host_is_big_endian_(isBigEndian()), swap_endian_(false),
    MASK_UPPER_(0xFF000000), MASK_LOWER_(0x000000FF)
{
  header_.machinestamp = 0;
  header_.map[3] = '\0';
}

MRCReader::~MRCReader() {}


bool MRCReader::read(const std::string& filename, NrrdDataHandle& nrrd_handle)
{
  try
  {
    std::ifstream::pos_type size = 0;

    std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
    if (! in)
    {
      error("Failed to open file " + filename);
      return false;
    }

    in.seekg(0, std::ios::end);
    size = in.tellg();
    in.seekg(0, std::ios::beg);
    std::ifstream::pos_type data_block_size = size - static_cast<std::ifstream::pos_type>(MRC_HEADER_LENGTH);

    in.read(reinterpret_cast<char*>(&header_), MRC_HEADER_LENGTH);
    process_header(&header_, MRC_HEADER_LENGTH);

    if (header_.mode == MRC_SHORT_COMPLEX || header_.mode == MRC_FLOAT_COMPLEX)
    {
      error("Complex mode is not supported.");
      return false;
    }
    // read in rest of file

    // allocate NRRD
    nrrd_handle = new NrrdData;
    Nrrd *n = nrrd_handle->nrrd_;
    if (nrrd_handle.get_rep() == 0)
    {
      error("Could not allocate nrrd.");
      return false;
    }

    Guard g(&(nrrd_handle->lock));
    // set NRRD type based on mode
    // Seg3D supports floats...

    // set NRRD header information
    const int NDIMS = 3;
    n->spaceDim = NDIMS;
    size_t sz[NDIMS];

    for (int i = 0; i < 3; ++i)
    {
      n->axis[i].spaceDirection[0] = 0.0;
      n->axis[i].spaceDirection[1] = 0.0;
      n->axis[i].spaceDirection[2] = 0.0;
    }

    bool useNewOrigin = true;
    if ((0 == header_.xorigin || 0 == header_.yorigin || 0 == header_.zorigin) &&
        (0 != header_.nxstart || 0 != header_.nystart || 0 != header_.nzstart))
    {
      // use n[x|y|z]start
      useNewOrigin = false;
    }

    // X=1, Y=2 and Z=3
    // axis corresponding to columns
    // fastest moving axis
    switch (header_.mapc) {
    case 1:
      n->axis[0].size = header_.nx;
      n->axis[0].spaceDirection[0] = header_.mx;
      if (useNewOrigin)
      {
        n->spaceOrigin[0] = header_.xorigin;
      }
      else
      {
        n->spaceOrigin[0] = header_.nxstart;
      }
      sz[0] = header_.nx;
      break;
    case 2:
      n->axis[0].size = header_.ny;
      n->axis[0].spaceDirection[0] = header_.my;
      if (useNewOrigin)
      {
        n->spaceOrigin[0] = header_.yorigin;
      }
      else
      {
        n->spaceOrigin[0] = header_.nystart;
      }
      sz[0] = header_.ny;
      break;
    case 3:
      n->axis[0].size = header_.nz;
      n->axis[0].spaceDirection[0] = header_.mz;
      if (useNewOrigin)
      {
        n->spaceOrigin[0] = header_.zorigin;
      }
      else
      {
        n->spaceOrigin[0] = header_.nzstart;
      }
      sz[0] = header_.nz;
      break;
    default:
      // TODO: error?
      break;
    }

    // axis corresponding to rows
    // slowest moving axis
    switch (header_.mapr) {
    case 1:
      n->axis[1].size = header_.nx;
      n->axis[1].spaceDirection[1] = header_.mx;
      if (useNewOrigin)
      {
        n->spaceOrigin[1] = header_.xorigin;
      }
      else
      {
        n->spaceOrigin[1] = header_.nxstart;
      }
      sz[1] = header_.nx;
      break;
    case 2:
      n->axis[1].size = header_.ny;
      n->axis[1].spaceDirection[1] = header_.my;
      if (useNewOrigin)
      {
        n->spaceOrigin[1] = header_.yorigin;
      }
      else
      {
        n->spaceOrigin[1] = header_.nystart;
      }
      sz[1] = header_.ny;
      break;
    case 3:
      n->axis[1].size = header_.nz;
      n->axis[1].spaceDirection[1] = header_.mz;
      if (useNewOrigin)
      {
        n->spaceOrigin[1] = header_.zorigin;
      }
      else
      {
        n->spaceOrigin[1] = header_.nzstart;
      }
      sz[1] = header_.nz;
      break;
    default:
      // TODO: error?
      break;
    }

    // axis corresponding to sections
    switch (header_.maps) {
    case 1:
      n->axis[2].size = header_.nx;
      n->axis[2].spaceDirection[2] = header_.mx;
      if (useNewOrigin)
      {
        n->spaceOrigin[2] = header_.xorigin;
      }
      else
      {
        n->spaceOrigin[2] = header_.nxstart;
      }
      sz[2] = header_.nx;
      break;
    case 2:
      n->axis[2].size = header_.ny;
      n->axis[2].spaceDirection[2] = header_.my;
      if (useNewOrigin)
      {
        n->spaceOrigin[2] = header_.yorigin;
      }
      else
      {
        n->spaceOrigin[2] = header_.nystart;
      }
      sz[2] = header_.ny;
      break;
    case 3:
      n->axis[2].size = header_.nz;
      n->axis[2].spaceDirection[2] = header_.mz;
      if (useNewOrigin)
      {
        n->spaceOrigin[2] = header_.zorigin;
      }
      else
      {
        n->spaceOrigin[2] = header_.nzstart;
      }
      sz[2] = header_.nz;
      break;
    default:
      // TODO: error?
      break;
    }

    n->axis[0].center = nrrdCenterUnknown;
    n->axis[1].center = nrrdCenterUnknown;
    n->axis[2].center = nrrdCenterUnknown;

    // cache alternate origin
    if (useNewOrigin)
    {
      nrrd_handle->set_property("nxstart", header_.nxstart, false);
      nrrd_handle->set_property("nystart", header_.nystart, false);
      nrrd_handle->set_property("nzstart", header_.nzstart, false);
    }
    else
    {
      nrrd_handle->set_property("xorigin", header_.xorigin, false);
      nrrd_handle->set_property("yorigin", header_.yorigin, false);
      nrrd_handle->set_property("zorigin", header_.zorigin, false);
    }

    nrrd_handle->set_property("cell_xlen", header_.xlen, false);
    nrrd_handle->set_property("cell_ylen", header_.ylen, false);
    nrrd_handle->set_property("cell_zlen", header_.zlen, false);

    nrrd_handle->set_property("cell_alpha", header_.alpha, false);
    nrrd_handle->set_property("cell_beta", header_.beta, false);
    nrrd_handle->set_property("cell_gamma", header_.gamma, false);
    nrrd_handle->set_property("dmin", header_.dmin, false);
    nrrd_handle->set_property("dmax", header_.dmax, false);
    nrrd_handle->set_property("dmean", header_.dmean, false);
    nrrd_handle->set_property("ispg", header_.ispg, false);
    nrrd_handle->set_property("nsymbt", header_.nsymbt, false);
    nrrd_handle->set_property("rms", header_.rms, false);

    boost::shared_array<char> data(new char[data_block_size]);
    in.read(data.get(), data_block_size);

    unsigned int nrrd_type;
    switch (header_.mode)
    {
    case MRC_CHAR:
      nrrd_type = get_nrrd_type<char>();
      break;
    case MRC_SHORT:
      nrrd_type = get_nrrd_type<short>();
      break;
    case MRC_FLOAT:
      nrrd_type = get_nrrd_type<float>();
      break;
    default:
      error("Unsupported MRC format.");
      return false;
    }

    // set NRRD data
    if ( nrrdWrap_nva(n, data.get(), nrrd_type, NDIMS, sz) != 0 )
    {
      char *err = biffGetDone(NRRD);
      error(std::string("NRRD error: nrrdAlloc_nva failed with error ") + err);
      free(err);
      data.reset();
      return false;
    }
  }
  // TODO: ios specific exceptions
  catch (...)
  {
    error("failed...");
    return false;
  }

  return true;
}

bool MRCReader::process_header(void* buffer, int buffer_len)
{
  int* long_word_buffer = reinterpret_cast<int*>(buffer);
  const int machine_stamp = long_word_buffer[53];

  // file endianness vs. host endianness
  // N.B. machine stamp field is not always implemented reliably
  if (host_is_big_endian_)
  {
    // big endian host, data read on little endian machine
    if (( (machine_stamp & MASK_LOWER_) == 0x44 ) || ( (machine_stamp & MASK_UPPER_) == 0x44000000 ))
    {
      swap_endian_ = true;
    }
    // big endian host, data read on big endian machine
    else if ( (machine_stamp &  MASK_LOWER_) == 0x11)
    {
      swap_endian_ = false;
    }
    else
    {
      // Guess endianness of data by checking upper bytes of nx:
      // assumes nx < 2^16
      char* char_buffer = reinterpret_cast<char*>(buffer);
      swap_endian_ = true;
      for (size_t j = MRC_LONG_WORD / 2; j < MRC_LONG_WORD; ++j) {
        std::cout << static_cast<int>(char_buffer[j]) << std::endl;
        if (char_buffer[j] != 0) {
          swap_endian_ = false;
          break;
        }
      }
    }
  }
  else // little endian (no other architectures supported by Core/Endian)
  {
    // little endian host, data read on little endian machine
    if (( (machine_stamp & MASK_LOWER_) == 0x44 ) || ( (machine_stamp & MASK_UPPER_) == 0x44000000 ))
    {
      swap_endian_ = false;
    }
    // little endian host, data read on big endian machine
    else if ( (machine_stamp &  MASK_LOWER_) == 0x11)
    {
      swap_endian_ = true;
    }
    else
    {
      // Guess endianness of data by checking lower bytes of nx:
      // assumes nx < 2^16
      char* char_buffer = reinterpret_cast<char*>(buffer);
      swap_endian_ = true;
      for (size_t j = 0; j < MRC_LONG_WORD / 2; ++j) {
        std::cout << static_cast<int>(char_buffer[j]) << std::endl;
        if (char_buffer[j] != 0) {
          swap_endian_ = false;
          break;
        }
      }
    }
  }

  if (swap_endian_)
  {
    // TODO: check this!!!
    for(int i = 0; i <  MRC_HEADER_LENGTH_LWORDS; ++i)
    {
      swapbytes(long_word_buffer[i]);
    }
  }

  MRCHeader *h = reinterpret_cast<MRCHeader*>(long_word_buffer);
  int temp = 0;

  temp = long_word_buffer[10];
  h->xlen = static_cast<float>(temp);
  temp = long_word_buffer[11];
  h->ylen = static_cast<float>(temp);
  temp = long_word_buffer[12];
  h->zlen = static_cast<float>(temp);

  temp = long_word_buffer[13];
  h->alpha = static_cast<float>(temp);
  temp = long_word_buffer[14];
  h->beta = static_cast<float>(temp);
  temp = long_word_buffer[15];
  h->gamma = static_cast<float>(temp);

  temp = long_word_buffer[19];
  h->dmin = static_cast<float>(temp);
  temp = long_word_buffer[20];
  h->dmax = static_cast<float>(temp);
  temp = long_word_buffer[21];
  h->dmean = static_cast<float>(temp);
  temp = long_word_buffer[22];
  h->ispg = static_cast<int>(temp);
  temp = long_word_buffer[23];
  h->nsymbt = static_cast<int>(temp);

  temp = long_word_buffer[49];
  h->xorigin = static_cast<float>(temp);
  temp = long_word_buffer[50];
  h->yorigin = static_cast<float>(temp);
  temp = long_word_buffer[51];
  h->zorigin = static_cast<float>(temp);
  // Force a sentinel null char
  h->map[MRC_LONG_WORD-1] = '\0';

  temp = long_word_buffer[54];
  h->rms = static_cast<float>(temp);

  return true;
}

  
}
