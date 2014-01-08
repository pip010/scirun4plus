#include <volutils.h>

VolumeTensor readTensorVolumeFile(char* filename)
{
  Nrrd *nin;
  char me[]="readTensorFile", *err;
  /* create a nrrd; at this point this is just an empty container */
  nin = nrrdNew();

  /* read in the nrrd from file */
  if (nrrdLoad(nin, filename, NULL)) 
    {
      err = biffGetDone(NRRD);
      fprintf(stderr, "%s: trouble reading \"%s\":\n%s", me, filename, err);
      free(err);
      return VolumeTensor();
    }

  //  nrrdKind3DSymMatrix,       /* 28: Mxx Mxy Mxz Myy Myz Mzz */
  //  nrrdKind3DMaskedSymMatrix, /* 29: mask Mxx Mxy Mxz Myy Myz Mzz */
  //  nrrdTypeFloat,         /*  9:          4-byte floating point */
  //  nrrdTypeDouble,        /* 10:          8-byte floating point */

  /* check dimension */
  int num_elements_per_voxel = (int)(nin->axis[0].size);
  if (!((nin->dim == 4)&&((num_elements_per_voxel == 7)||(num_elements_per_voxel == 6))
	&&((nin->type == nrrdTypeDouble)||(nin->type == nrrdTypeFloat))
	))
    {
      printf("%s: Not suitable tensor volume - \"%s\" is a %d-dimensional nrrd of type %d (%s)\n", 
	     me, filename, nin->dim, nin->type,
	     airEnumStr(nrrdType, nin->type));
      return VolumeTensor();
    }

  /* say something about the array */
  printf("%s: \"%s\" is a %d-dimensional nrrd of type %d (%s)\n", 
         me, filename, nin->dim, nin->type,
         airEnumStr(nrrdType, nin->type));
  printf("%s: the array contains %d elements, each %d bytes in size\n",
         me, (int)nrrdElementNumber(nin), (int)nrrdElementSize(nin));
  
  int w, h, d;
  VolumeTensor vol_return(w = nin->axis[1].size, 
			  h = nin->axis[2].size, 
			  d = nin->axis[3].size); 
  VISMatrix tensor(3,3);
  float *nrrd_data_float = (float*)nin->data;
  double *nrrd_data_double = (double*)nin->data;

  if (nin->type == nrrdTypeDouble)
  {
    for (int k = 0; k < d; k++)
      for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
          if (num_elements_per_voxel == 7) *(nrrd_data_double++);
          tensor.poke(0, 0) = *(nrrd_data_double++);
          tensor.poke(0, 1) = tensor.poke(1, 0) = *(nrrd_data_double++);
          tensor.poke(0, 2) = tensor.poke(2, 0) = *(nrrd_data_double++);
          tensor.poke(1, 1) = *(nrrd_data_double++);
          tensor.poke(1, 2) = tensor.poke(2, 1) = *(nrrd_data_double++);
          tensor.poke(2, 2)  = *(nrrd_data_double++);
          vol_return.poke(i, j, k) = tensor;
        }
  }
  else 
  {
    for (int k = 0; k < d; k++)
      for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
          if (num_elements_per_voxel == 7) *(nrrd_data_float++);
          tensor.poke(0, 0) = *(nrrd_data_float++);
          tensor.poke(0, 1) = tensor.poke(1, 0) = *(nrrd_data_float++);
          tensor.poke(0, 2) = tensor.poke(2, 0) = *(nrrd_data_float++);
          tensor.poke(1, 1) = *(nrrd_data_float++);
          tensor.poke(1, 2) = tensor.poke(2, 1) = *(nrrd_data_float++);
          tensor.poke(2, 2)  = *(nrrd_data_float++);
          vol_return.poke(i, j, k) = tensor;
        }
  }
  nrrdNuke(nin);
  return vol_return;
}


boolean writeTensorVolumeFile(char* filename, const VolumeTensor& vol_tensor)
{
  Nrrd *nout;
  char me[]="readTensorFile", *err;
  int w = vol_tensor.width(), h = vol_tensor.height(), d = vol_tensor.depth();
  /* create a nrrd; at this point this is just an empty container */
  nout = nrrdNew();
  nrrdAlloc_va(nout, nrrdTypeFloat, 4, 6, w, h, d);

  VISMatrix tensor(3,3);
  float *nrrd_data_float = (float*)nout->data;
  for (int k = 0; k < d; k++)
    for (int j = 0; j < h; j++)
      for (int i = 0; i < w; i++)
      {
        tensor = vol_tensor.peek(i, j, k);
        *(nrrd_data_float++) = tensor.peek(0, 0);
        *(nrrd_data_float++) = tensor.peek(0, 1); // assumes = tensor.poke(1, 0)
        *(nrrd_data_float++) = tensor.peek(0, 2); // assumes = tensor.poke(2, 0)
        *(nrrd_data_float++) = tensor.peek(1, 1);
        *(nrrd_data_float++) = tensor.peek(1, 2); // assumes = tensor.poke(2, 1)
        *(nrrd_data_float++) = tensor.peek(2, 2);
      }
  if (nrrdSave(filename, nout, NULL)) {
    err = biffGetDone(NRRD);
    fprintf(stderr, "%s: trouble writing \"%s\":\n%s", me, filename, err);
    free(err);
    return false;
  }
  nrrdNuke(nout);  
  return true;
}

boolean writeTensorVolumeFile(char* filename, const VolumeTensor& vol_tensor, const VolumeScalar& vol_scalar)
{
  Nrrd *nout;
  char me[]="readTensorFile", *err;
  int w = vol_tensor.width(), h = vol_tensor.height(), d = vol_tensor.depth();
  int i, j, k;
  /* create a nrrd; at this point this is just an empty container */
  nout = nrrdNew();
  nrrdAlloc_va(nout, nrrdTypeFloat, 4, 7, w, h, d);

  VISMatrix tensor(3,3);
  float *nrrd_data_float = (float*)nout->data;
  for (k = 0; k < d; k++)
    for (j = 0; j < h; j++)
      for (i = 0; i < w; i++)
	{
	  tensor = vol_tensor.peek(i, j, k);
	  *(nrrd_data_float++) = vol_scalar.peek(i, j, k);
	  *(nrrd_data_float++) = tensor.peek(0, 0);
	  *(nrrd_data_float++) = tensor.peek(0, 1); // assumes = tensor.poke(1, 0)
	  *(nrrd_data_float++) = tensor.peek(0, 2); // assumes = tensor.poke(2, 0)
	  *(nrrd_data_float++) = tensor.peek(1, 1);
	  *(nrrd_data_float++) = tensor.peek(1, 2); // assumes = tensor.poke(2, 1)
	  *(nrrd_data_float++) = tensor.peek(2, 2);
	}
  if (nrrdSave(filename, nout, NULL)) {
    err = biffGetDone(NRRD);
    fprintf(stderr, "%s: trouble writing \"%s\":\n%s", me, filename, err);
    free(err);
    return false;
  }
  nrrdNuke(nout);  
  return true;
}


VolumeVector volumeFirstDerivatives(const VolumeScalar& vol)
{
  int w = vol.width(), h = vol.height(), d = vol.depth();
  VISVolume<float> dx = vol.dx(), dy = vol.dy(), dz = vol.dz();
  VolumeVector ret(w, h, d);
  int i, j, k;
  for (k = 0; k < d; k++)
    for (j = 0; j < h; j++)
      for (i = 0; i < w; i++)
	ret.poke(i, j, k) = VISVector(dx.peek(i, j, k), dy.peek(i, j, k), dz.peek(i, j, k));
  return ret;
}

float resolveGrad(float grad_for, float grad_back, int up_or_down=-1)
{
  if (up_or_down < 0)
    {
      grad_for = VISmin(grad_for, 0.0f);
      grad_back = VISmax(grad_back, 0.0f);
      if (-1.0f*grad_for > grad_back)
	return(grad_for);
      else
	return(grad_back);
    }
  else
    {
      grad_for = VISmax(grad_for, 0.0f);
      grad_back = VISmin(grad_back, 0.0f);
      if (grad_for > -1.0f*grad_back)
	return(grad_for);
      else
	return(grad_back);
    }
}

VolumeVector volumeFirstDerivativesUpWind(const VolumeScalar& vol)
{
  int w = vol.width(), h = vol.height(), d = vol.depth();
  VISVolume<float> dx = vol.dx(), dy = vol.dy(), dz = vol.dz();
  VolumeVector ret(w, h, d);
  float grad_for, grad_back, grad_x, grad_y, grad_z;
  int i, j, k;
  for (k = 0; k < d; k++)
    for (j = 0; j < h; j++)
      for (i = 0; i < w; i++)
	{
	  if (i == (w-1))
	    grad_x = vol.peek(i, j, k) - vol.peek(i - 1, j, k);
	  else if (i == 0)
	    grad_x = vol.peek(i + 1, j, k) - vol.peek(i, j, k);
	  else
	    {
	      grad_back = vol.peek(i, j, k) - vol.peek(i - 1, j, k);
	      grad_for = vol.peek(i + 1, j, k) - vol.peek(i, j, k);
	      grad_x = resolveGrad(grad_for, grad_back, -1);
	    }
	  if (j == (h-1))
	    grad_y = vol.peek(i, j, k) - vol.peek(i, j - 1, k);
	  else if (j == 0)
	    grad_y = vol.peek(i, j + 1, k) - vol.peek(i, j, k);
	  else
	    {
	      grad_back = vol.peek(i, j, k) - vol.peek(i, j - 1, k);
	      grad_for = vol.peek(i, j + 1, k) - vol.peek(i, j, k);
	      grad_y = resolveGrad(grad_for, grad_back, -1);
	    }
	  if (k == (d - 1))
	    grad_z = vol.peek(i, j, k) - vol.peek(i, j, k - 1);
	  else if (k == 0)
	    grad_z = vol.peek(i, j, k + 1) - vol.peek(i, j, k);
	  else
	    {
	      grad_back = vol.peek(i, j, k) - vol.peek(i, j, k - 1);
	      grad_for = vol.peek(i, j, k + 1) - vol.peek(i, j, k);
	      grad_z = resolveGrad(grad_for, grad_back, -1);
	    }
	  ret.poke(i, j, k) = VISVector(grad_x, grad_y, grad_z);
	}
  return ret;
}



// VISVector vectorVolumeInterpolate(const VolumeVector& vol, const VISVector& vector)
// {
//   if (!vol_term.checkBounds((float)position[0], (float)position[1], (float)position[2]))
//     {
//       //         cout << "out of bounds" << endl;
//       return(VISVector());
//     }
//   int x_hi, x_lo, y_hi, y_lo, z_hi, z_lo;
//   int v_x_hi, v_x_lo, v_y_hi, v_y_lo, v_z_hi, v_z_lo;
//   x_lo = (int) floor(position[0]);
//   x_hi = (int) ceil(position[0]);
//   y_lo = (int) floor(position[1]);
//   y_hi = (int) ceil(position[1]);
//   z_lo = (int) floor(position[2]);
//   z_hi = (int) ceil(position[2]);
//   float vx, vy, vz;
// }

VolumeTensor volumeSecondDerivatives(const VolumeScalar& vol)
{
  int w = vol.width(), h = vol.height(), d = vol.depth();
  VISVolume<float> dxx = vol.dx(2), dyy = vol.dy(2), dzz = vol.dz(2), 
    dxy = vol.derivative(1, 1, 0), dxz = vol.derivative(1, 0, 1), 
    dyz = vol.derivative(0, 1, 1);
  VolumeTensor ret(w, h, d);
  VISMatrix tensor(3, 3);
  int i, j, k;
  for (k = 0; k < d; k++)
    for (j = 0; j < h; j++)
      for (i = 0; i < w; i++)
	{
	  tensor.poke(0, 0) = dxx.peek(i, j, k);
	  tensor.poke(1, 1) = dyy.peek(i, j, k);
	  tensor.poke(2, 2) = dzz.peek(i, j, k);
	  tensor.poke(1, 0) = tensor.poke(0, 1) = dxy.peek(i, j, k);
	  tensor.poke(2, 0) = tensor.poke(0, 2) = dxz.peek(i, j, k);
	  tensor.poke(2, 1) = tensor.poke(1, 2) = dyz.peek(i, j, k);
	  ret.poke(i, j, k) = tensor;
	}
  return ret;
}


VolumeVector readVectorVolumeFile(char* filename)
{
  Nrrd *nin;
  char me[]="readVectorFile", *err;
  /* create a nrrd; at this point this is just an empty container */
  nin = nrrdNew();

  /* read in the nrrd from file */
  if (nrrdLoad(nin, filename, NULL)) 
    {
      err = biffGetDone(NRRD);
      fprintf(stderr, "%s: trouble reading \"%s\":\n%s", me, filename, err);
      free(err);
      return VolumeVector();
    }

  //  nrrdKind3DSymMatrix,       /* 28: Mxx Mxy Mxz Myy Myz Mzz */
  //  nrrdKind3DMaskedSymMatrix, /* 29: mask Mxx Mxy Mxz Myy Myz Mzz */
  //  nrrdTypeFloat,         /*  9:          4-byte floating point */
  //  nrrdTypeDouble,        /* 10:          8-byte floating point */

  /* check dimension */
  int num_elements_per_voxel = (int)(nin->axis[0].size);
  if (!((nin->dim == 4)&&(num_elements_per_voxel == 3)
	&&((nin->type == nrrdTypeDouble)||(nin->type == nrrdTypeFloat))
	))
    {
      printf("%s: Not suitable vector volume - \"%s\" is a %d-dimensional nrrd of type %d (%s)\n", 
	     me, filename, nin->dim, nin->type,
	     airEnumStr(nrrdType, nin->type));
      return VolumeVector();
    }

  /* say something about the array */
  printf("%s: \"%s\" is a %d-dimensional nrrd of type %d (%s)\n", 
         me, filename, nin->dim, nin->type,
         airEnumStr(nrrdType, nin->type));
  printf("%s: the array contains %d elements, each %d bytes in size\n",
         me, (int)nrrdElementNumber(nin), (int)nrrdElementSize(nin));
  
  int w, h, d, i, j, k;
  VolumeVector vol_return(w = nin->axis[1].size, 
			  h = nin->axis[2].size, 
			  d = nin->axis[3].size); 
  VISVector vector(3);
  float *nrrd_data_float = (float*)nin->data;
  double *nrrd_data_double = (double*)nin->data;

  if (nin->type == nrrdTypeDouble)
    {
      for (k = 0; k < d; k++)
	for (j = 0; j < h; j++)
	  for (i = 0; i < w; i++)
	    {
	      vector.poke(0) = *(nrrd_data_double++);
	      vector.poke(1) = *(nrrd_data_double++);
	      vector.poke(2) = *(nrrd_data_double++);
	      vol_return.poke(i, j, k) = vector;
	    }
    }
  else 
    {
      for (k = 0; k < d; k++)
	for (j = 0; j < h; j++)
	  for (i = 0; i < w; i++)
	    {
	      vector.poke(0) = *(nrrd_data_float++);
	      vector.poke(1) = *(nrrd_data_float++);
	      vector.poke(2) = *(nrrd_data_float++);
	      vol_return.poke(i, j, k) = vector;
	    }
    }
  return vol_return;
}

VolumeScalar readScalarVolumeFile(char* filename)
{
  VolScale scale;
  return(readScalarVolumeFile(filename,scale));
}


VolumeScalar readScalarVolumeFile(char* filename, VolScale& scale)
{
  Nrrd *nin;
  char me[]="readScalarFile", *err;
  /* create a nrrd; at this point this is just an empty container */
  nin = nrrdNew();

  /* read in the nrrd from file */
  if (nrrdLoad(nin, filename, NULL)) 
    {
      err = biffGetDone(NRRD);
      fprintf(stderr, "%s: trouble reading \"%s\":\n%s", me, filename, err);
      free(err);
      return VolumeScalar();
    }

  //  nrrdKind3DSymMatrix,       /* 28: Mxx Mxy Mxz Myy Myz Mzz */
  //  nrrdKind3DMaskedSymMatrix, /* 29: mask Mxx Mxy Mxz Myy Myz Mzz */
  //  nrrdTypeFloat,         /*  9:          4-byte floating point */
  //  nrrdTypeDouble,        /* 10:          8-byte floating point */

  /* check dimension */
  int num_elements_per_voxel = (int)(nin->axis[0].size);
  if (!((nin->dim == 3)
    &&((nin->type == nrrdTypeDouble)
	   ||(nin->type == nrrdTypeFloat)
	   ||(nin->type == nrrdTypeUInt)
	   ||(nin->type == nrrdTypeUChar)
	   )
	))
    {
      printf("%s: Not suitable scalar volume - \"%s\" is a %d-dimensional nrrd of type %d (%s)\n", 
	     me, filename, nin->dim, nin->type,
	     airEnumStr(nrrdType, nin->type));
     return VolumeScalar();
    }
   
   scale.center_x = nin->axis[0].center;
   scale.center_y = nin->axis[1].center;
   scale.center_z = nin->axis[2].center;

  /* say something about the array */
  //  printf("%s: \"%s\" is a %d-dimensional nrrd of type %d (%s)\n", 
  //         me, filename, nin->dim, nin->type,
  //         airEnumStr(nrrdType, nin->type));
  //  printf("%s: the array contains %d elements, each %d bytes in size\n",
  //         me, (int)nrrdElementNumber(nin), (int)nrrdElementSize(nin));

  // TODO: temp fix to make sure we read nrrds with scaling
  if (nin->spaceDim == 3)
  {
    scale.spacing_x = nin->axis[0].spaceDirection[0];
    scale.spacing_y = nin->axis[1].spaceDirection[1];
    scale.spacing_z = nin->axis[2].spaceDirection[2];
  }
  else
  {
    scale.spacing_x = nin->axis[0].spacing;
    scale.spacing_y = nin->axis[1].spacing;
    scale.spacing_z = nin->axis[2].spacing;
  }

  int w = static_cast<int>(nin->axis[0].size);
  int h = static_cast<int>(nin->axis[1].size);
  int d = static_cast<int>(nin->axis[2].size);

  VolumeScalar vol_return(w, h, d); 

  float *nrrd_data_float = (float*)nin->data;
  double *nrrd_data_double = (double*)nin->data;
  unsigned int *nrrd_data_uint = (unsigned int*)nin->data;
  unsigned char *nrrd_data_uchar = (unsigned char*)nin->data;

  if (nin->type == nrrdTypeDouble)
  {
      for (int k = 0; k < d; k++)
	      for (int j = 0; j < h; j++)
	        for (int i = 0; i < w; i++)
	        {
	           vol_return.poke(i, j, k) = *(nrrd_data_double++);
	        }
  }
  else if (nin->type == nrrdTypeFloat)
  {
      for (int k = 0; k < d; k++)
	for (int j = 0; j < h; j++)
	  for (int i = 0; i < w; i++)
	    {
	      vol_return.poke(i, j, k) = *(nrrd_data_float++);
	    }
    }
  else if (nin->type == nrrdTypeUInt)
    {
      for (int k = 0; k < d; k++)
	for (int j = 0; j < h; j++)
	  for (int i = 0; i < w; i++)
	    {
	      vol_return.poke(i, j, k) = *(nrrd_data_uint++);
	    }
    }
  else if (nin->type == nrrdTypeUChar)
    {
      for (int k = 0; k < d; k++)
	for (int j = 0; j < h; j++)
	  for (int i = 0; i < w; i++)
	    {
	      vol_return.poke(i, j, k) = *(nrrd_data_uchar++);
	    }
    }
  else 
    {
      cout << "bad nrrd type in VolumeScalarRead " << endl;
      return(VISVolume<float>());
    }
  return vol_return;
}



boolean writeVectorVolumeFile(char* filename, const VolumeVector& vol_vector)
{
  Nrrd *nout;
  char me[]="readTensorFile", *err;
  int w = vol_vector.width(), h = vol_vector.height(), d = vol_vector.depth();

  /* create a nrrd; at this point this is just an empty container */
  nout = nrrdNew();
  nrrdAlloc_va(nout, nrrdTypeFloat, 4, 3, w, h, d);

  VISVector vector(3);
  float *nrrd_data_float = (float*)nout->data;
  for (int k = 0; k < d; k++)
    for (int j = 0; j < h; j++)
      for (int i = 0; i < w; i++)
	{
	  vector = vol_vector.peek(i, j, k);
	  *(nrrd_data_float++) = vector.peek(0);
	  *(nrrd_data_float++) = vector.peek(1); 
	  *(nrrd_data_float++) = vector.peek(2); 
	}
  
  NrrdIoState *nio = nrrdIoStateNew();
  // set encoding to be gzip
  nio->encoding = nrrdEncodingArray[nrrdEncodingTypeGzip];
  // set format to be nrrd
  nio->format = nrrdFormatArray[nrrdFormatTypeNRRD];
  // set endian to be endian of machine
  nio->endian = AIR_ENDIAN;
  nio->zlibLevel = 6;
  
  if (nrrdSave(filename, nout, nio)) {
    err = biffGetDone(NRRD);
    fprintf(stderr, "%s: trouble writing \"%s\":\n%s", me, filename, err);
    free(err);
    return false;
  }
  nrrdNuke(nout);  
  return true;
}


boolean writeScalarVolumeFile(char* filename, const VolumeScalar& vol_scalar)
{
  VolScale scale;
  return writeScalarVolumeFile(filename,vol_scalar,scale);
}

boolean writeScalarVolumeFile(char* filename, const VolumeScalar& vol_scalar, VolScale& scale)
{
  Nrrd *nout;
  char me[]="readTensorFile", *err;
  size_t w = vol_scalar.width(), h = vol_scalar.height(), d = vol_scalar.depth();
  int i, j, k;
  /* create a nrrd; at this point this is just an empty container */
  nout = nrrdNew();
  nrrdAlloc_va(nout, nrrdTypeFloat, 3, w, h, d);

  nout->spaceDim = 3;
  nout->spaceOrigin[0] = 0.0;
  nout->spaceOrigin[1] = 0.0;
  nout->spaceOrigin[2] = 0.0;

  nout->axis[0].spaceDirection[0] = scale.spacing_x;
  nout->axis[0].spaceDirection[1] = 0.0;
  nout->axis[0].spaceDirection[2] = 0.0;
  nout->axis[1].spaceDirection[0] = 0.0;
  nout->axis[1].spaceDirection[1] = scale.spacing_y;
  nout->axis[1].spaceDirection[2] = 0.0;
  nout->axis[2].spaceDirection[0] = 0.0;
  nout->axis[2].spaceDirection[1] = 0.0;
  nout->axis[2].spaceDirection[2] = scale.spacing_z;
  
  nout->axis[0].center = scale.center_x;
  nout->axis[1].center = scale.center_y;
  nout->axis[2].center = scale.center_z;
 
  float *nrrd_data_float = (float*)nout->data;

  for (k = 0; k < d; k++)
    for (j = 0; j < h; j++)
      for (i = 0; i < w; i++)
	{
	  *(nrrd_data_float++) = vol_scalar.peek(i, j, k);
	}
  NrrdIoState *nio = nrrdIoStateNew();
  // set encoding to be gzip
  nio->encoding = nrrdEncodingArray[nrrdEncodingTypeGzip];
  // set format to be nrrd
  nio->format = nrrdFormatArray[nrrdFormatTypeNRRD];
  // set endian to be endian of machine
  nio->endian = AIR_ENDIAN;
  nio->zlibLevel = 6;
  
  if (nrrdSave(filename, nout, nio)) {
    err = biffGetDone(NRRD);
    fprintf(stderr, "%s: trouble writing \"%s\":\n%s", me, filename, err);
    free(err);
    return false;
  }

  nrrdNuke(nout);  
  return true;
}

