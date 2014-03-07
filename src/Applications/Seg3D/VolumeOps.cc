//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//    File   : VolumeOps.cc
//    Author : McKay Davis
//    Date   : Tue Oct 31 17:42:53 2006


#include <Applications/Seg3D/VolumeOps.h>
#include <Core/Util/StringUtil.h>
#include <Core/Math/MinMax.h>

#include <algorithm>

namespace SCIRun {

size_t
VolumeOps::nrrd_type_size(NrrdDataHandle &nrrdh)
{
  switch (nrrdh->nrrd_->type) {
  case nrrdTypeChar: return sizeof(char); break;
  case nrrdTypeUChar: return sizeof(unsigned char); break;
  case nrrdTypeShort: return sizeof(short); break;
  case nrrdTypeUShort: return sizeof(unsigned short); break;
  case nrrdTypeInt: return sizeof(int); break;
  case nrrdTypeUInt: return sizeof(unsigned int); break;
  case nrrdTypeLLong: return sizeof(signed long long); break;
  case nrrdTypeULLong: return sizeof(unsigned long long); break;
  case nrrdTypeFloat: return sizeof(float); break;
  case nrrdTypeDouble: return sizeof(double); break;
  default:
    throw "Unsupported data type: "+to_string(nrrdh->nrrd_->type); break;
  }
  return 0;
}


size_t
VolumeOps::nrrd_elem_count(NrrdDataHandle &nrrdh)
{
  if (!nrrdh->nrrd_->dim) return 0;
  size_t size = nrrdh->nrrd_->axis[0].size;
  for (unsigned int a = 1; a < nrrdh->nrrd_->dim; ++a)
    size *= nrrdh->nrrd_->axis[a].size;
  return size;
}


size_t
VolumeOps::nrrd_data_size(NrrdDataHandle &nrrdh)
{
  return nrrd_elem_count(nrrdh) * nrrd_type_size(nrrdh);
}



NrrdDataHandle
VolumeOps::float_to_bit(NrrdDataHandle &inh, float ref, unsigned int mask)
{
  NrrdDataHandle nrrdh = VolumeOps::create_clear_nrrd(inh, LabelNrrdType);
  const size_t n = nrrd_elem_count(inh);
  float *src = (float *)inh->nrrd_->data;
  label_type *dst = (label_type *)nrrdh->nrrd_->data;
  for (size_t i = 0; i < n; ++i, ++src, ++dst) {
    if (*src >= ref) {
      *dst |= mask;
    }
  }

  return nrrdh;
}


void
VolumeOps::bit_set(NrrdDataHandle &dnrrd, unsigned int dbit)
{
  const size_t dsize = nrrd_elem_count(dnrrd);
  label_type *dst = (label_type *)dnrrd->nrrd_->data;
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);

  for (size_t i = 0; i < dsize; ++i, ++dst)
  {
    *dst |= dbit;
  }
}


void
VolumeOps::bit_clear(NrrdDataHandle &dnrrd, unsigned int dbit)
{
  dbit = ~dbit;
  const size_t dsize = nrrd_elem_count(dnrrd);
  label_type *dst = (label_type *)dnrrd->nrrd_->data;
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);

  for (size_t i = 0; i < dsize; ++i, ++dst)
  {
    *dst &= dbit;
  }
}


void
VolumeOps::bit_invert(NrrdDataHandle &dnrrd, unsigned int dbit,
                      NrrdDataHandle &snrrd, unsigned int sbit)
{
  const size_t ssize = nrrd_elem_count(snrrd);
  const size_t dsize = nrrd_elem_count(dnrrd);
  ASSERT(ssize == dsize);
  ASSERT(snrrd->nrrd_->type == LabelNrrdType);
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);

  label_type *src = (label_type *)snrrd->nrrd_->data;
  label_type *dst = (label_type *)dnrrd->nrrd_->data;
  for (size_t i = 0; i < dsize; ++i, ++dst, ++src) {
    if ((*src) & sbit) { *dst &= ~dbit; }
    else { *dst |= dbit; }
  }
}


void
VolumeOps::bit_copy(NrrdDataHandle &dnrrd, unsigned int dbit,
                    NrrdDataHandle &snrrd, unsigned int sbit)
{
  const size_t ssize = nrrd_elem_count(snrrd);
  const size_t dsize = nrrd_elem_count(dnrrd);
  ASSERT(ssize == dsize);
  ASSERT(snrrd->nrrd_->type == LabelNrrdType);
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);

  label_type *src = (label_type *)snrrd->nrrd_->data;
  label_type *dst = (label_type *)dnrrd->nrrd_->data;
  for (size_t i = 0; i < dsize; ++i, ++dst, ++src) {
    if ((*src) & sbit) { *dst |= dbit; }
    else { *dst &= ~dbit; }
  }
}


void
VolumeOps::bit_and(NrrdDataHandle &dnrrd, unsigned int dbit,
                   NrrdDataHandle &snrrd1, unsigned int sbit1,
                   NrrdDataHandle &snrrd2, unsigned int sbit2)
{
  const size_t ssize1 = nrrd_elem_count(snrrd1);
  const size_t ssize2 = nrrd_elem_count(snrrd2);
  const size_t dsize = nrrd_elem_count(dnrrd);
  ASSERT(ssize1 == dsize && ssize1 == ssize2);
  ASSERT(snrrd1->nrrd_->type == LabelNrrdType);
  ASSERT(snrrd2->nrrd_->type == LabelNrrdType);
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);

  label_type *src1 = (label_type *)snrrd1->nrrd_->data;
  label_type *src2 = (label_type *)snrrd2->nrrd_->data;
  label_type *dst = (label_type *)dnrrd->nrrd_->data;
  for (size_t i = 0; i < dsize; ++i, ++dst, ++src1, ++src2) {
    if ((*src1) & sbit1 && (*src2) & sbit2) { *dst |= dbit; }
    else { *dst &= ~dbit; }
  }
}


void
VolumeOps::bit_or(NrrdDataHandle &dnrrd, unsigned int dbit,
                  NrrdDataHandle &snrrd1, unsigned int sbit1,
                  NrrdDataHandle &snrrd2, unsigned int sbit2)
{
  const size_t ssize1 = nrrd_elem_count(snrrd1);
  const size_t ssize2 = nrrd_elem_count(snrrd2);
  const size_t dsize = nrrd_elem_count(dnrrd);
  ASSERT(ssize1 == dsize && ssize1 == ssize2);
  ASSERT(snrrd1->nrrd_->type == LabelNrrdType);
  ASSERT(snrrd2->nrrd_->type == LabelNrrdType);
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);

  label_type *src1 = (label_type *)snrrd1->nrrd_->data;
  label_type *src2 = (label_type *)snrrd2->nrrd_->data;
  label_type *dst = (label_type *)dnrrd->nrrd_->data;
  for (size_t i = 0; i < dsize; ++i, ++dst, ++src1, ++src2) {
    if (((*src1) & sbit1) || ((*src2) & sbit2)) { *dst |= dbit; }
    else { *dst &= ~dbit; }
  }
}


void
VolumeOps::bit_xor(NrrdDataHandle &dnrrd, unsigned int dbit,
                   NrrdDataHandle &snrrd1, unsigned int sbit1,
                   NrrdDataHandle &snrrd2, unsigned int sbit2)
{
  const size_t ssize1 = nrrd_elem_count(snrrd1);
  const size_t ssize2 = nrrd_elem_count(snrrd2);
  const size_t dsize = nrrd_elem_count(dnrrd);
  ASSERT(ssize1 == dsize && ssize1 == ssize2);
  ASSERT(snrrd1->nrrd_->type == LabelNrrdType);
  ASSERT(snrrd2->nrrd_->type == LabelNrrdType);
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);

  label_type *src1 = (label_type *)snrrd1->nrrd_->data;
  label_type *src2 = (label_type *)snrrd2->nrrd_->data;
  label_type *dst = (label_type *)dnrrd->nrrd_->data;
  for (size_t i = 0; i < dsize; ++i, ++dst, ++src1, ++src2) {
    if ((((*src1) & sbit1) || ((*src2) & sbit2))
        && !(((*src1) & sbit1) && ((*src2) & sbit2))) { *dst |= dbit; }
    else { *dst &= ~dbit; }
  }
}


// Remove s2 from s1 (dst <= s2 && !s1)
void
VolumeOps::bit_andnot(NrrdDataHandle &dnrrd, unsigned int dbit,
                      NrrdDataHandle &snrrd1, unsigned int sbit1,
                      NrrdDataHandle &snrrd2, unsigned int sbit2)
{
  const size_t ssize1 = nrrd_elem_count(snrrd1);
  const size_t ssize2 = nrrd_elem_count(snrrd2);
  const size_t dsize = nrrd_elem_count(dnrrd);
  ASSERT(ssize1 == dsize && ssize1 == ssize2);
  ASSERT(snrrd1->nrrd_->type == LabelNrrdType);
  ASSERT(snrrd2->nrrd_->type == LabelNrrdType);
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);

  label_type *src1 = (label_type *)snrrd1->nrrd_->data;
  label_type *src2 = (label_type *)snrrd2->nrrd_->data;
  label_type *dst = (label_type *)dnrrd->nrrd_->data;
  for (size_t i = 0; i < dsize; ++i, ++dst, ++src1, ++src2) {
    if ((*src1) & sbit1 && !((*src2) & sbit2)) { *dst |= dbit; }
    else { *dst &= ~dbit; }
  }
}


size_t
VolumeOps::bit_count(NrrdDataHandle &snrrd, unsigned int sbit)
{
  const size_t ssize = nrrd_elem_count(snrrd);
  ASSERT(snrrd->nrrd_->type == LabelNrrdType);

  label_type *src = (label_type *)snrrd->nrrd_->data;

  size_t count = 0;
  for (size_t i = 0; i < ssize; ++i, ++src) {
    if ((*src) & sbit) { count++; }
  }
  return count;
}


void
VolumeOps::merge_label_into_segmentation(NrrdDataHandle &seg,
                                         unsigned int svalue,
                                         NrrdDataHandle &label,
                                         unsigned int lbit)
{
  const size_t ssize = nrrd_elem_count(seg);
  const size_t lsize = nrrd_elem_count(label);
  ASSERT(ssize == lsize);
  ASSERT(label->nrrd_->type == LabelNrrdType);
  ASSERT(seg->nrrd_->type == nrrdTypeUChar);

  label_type *ldat = (label_type *)label->nrrd_->data;
  unsigned char *sdat = (unsigned char *)seg->nrrd_->data;
  for (size_t i = 0; i < ssize; ++i, ++ldat, ++sdat)
  {
    if ((*ldat) & lbit) { *sdat = svalue; }
  }
}


void
VolumeOps::extract_label_from_segmentation(NrrdDataHandle &label,
                                           unsigned int lbit,
                                           NrrdDataHandle &seg,
                                           unsigned int svalue)
{
  const size_t ssize = nrrd_elem_count(seg);
  const size_t lsize = nrrd_elem_count(label);
  ASSERT(ssize == lsize);
  ASSERT(label->nrrd_->type == LabelNrrdType);
  ASSERT(seg->nrrd_->type == nrrdTypeUInt);

  label_type *ldat = (label_type *)label->nrrd_->data;
  unsigned int *sdat = (unsigned int *)seg->nrrd_->data;
  for (size_t i = 0; i < ssize; ++i, ++ldat, ++sdat)
  {
    if ((*sdat) == svalue) { *ldat |= lbit; }
    else { *ldat &= ~lbit; }
  }
}


template<class T>
static void
templated_threshold(T *src, label_type *dst, size_t size,
                    unsigned int dbit, float minval, float maxval)
{
  for (size_t i = 0; i < size; ++i, ++dst, ++src)
  {
    if (*src >= minval && *src <= maxval) { *dst |= dbit; }
    else { *dst &= ~dbit; }
  }
}


void
VolumeOps::threshold(NrrdDataHandle &dnrrd, unsigned int dbit,
                     NrrdDataHandle &snrrd, float minval, float maxval)
{
  const size_t ssize = nrrd_elem_count(snrrd);
  const size_t dsize = nrrd_elem_count(dnrrd);
  ASSERT(ssize == dsize);
  ASSERT(dnrrd->nrrd_->type == LabelNrrdType);
  label_type *dst = (label_type *)dnrrd->nrrd_->data;

  switch (snrrd->nrrd_->type) {
  case nrrdTypeChar:
    templated_threshold((char *)snrrd->nrrd_->data, dst, ssize,
                        dbit, minval, maxval);
    break;
  case nrrdTypeUChar:
    templated_threshold((unsigned char *)snrrd->nrrd_->data, dst, ssize,
                        dbit, minval, maxval);
    break;
  case nrrdTypeShort:
    templated_threshold((short *)snrrd->nrrd_->data, dst, ssize,
                        dbit, minval, maxval);
    break;
  case nrrdTypeUShort:
    templated_threshold((unsigned short *)snrrd->nrrd_->data, dst, ssize,
                        dbit, minval, maxval);
    break;
  case nrrdTypeInt:
    templated_threshold((int *)snrrd->nrrd_->data, dst, ssize,
                        dbit, minval, maxval);
    break;
  case nrrdTypeUInt:
    templated_threshold((unsigned int *)snrrd->nrrd_->data, dst, ssize,
                        dbit, minval, maxval);
    break;
  case nrrdTypeFloat:
    templated_threshold((float *)snrrd->nrrd_->data, dst, ssize,
                        dbit, minval, maxval);
    break;
  case nrrdTypeDouble:
    templated_threshold((double *)snrrd->nrrd_->data, dst, ssize,
                        dbit, minval, maxval);
    break;
  default:
    ASSERTFAIL("Unsupported data type.");
  }
}


template<class T>
static void
templated_mask_data(T *dst, T *src, label_type *mask, size_t size,
                    unsigned int mbit, T newval)
{
  for (size_t i = 0; i < size; ++i)
  {
    dst[i] = (mask[i] & mbit) ? src[i] : newval;
  }
}


void
VolumeOps::mask_data(NrrdDataHandle &dnrrd, NrrdDataHandle &snrrd,
                     NrrdDataHandle &mnrrd, unsigned int mbit, float newval)
{
  const size_t dsize = VolumeOps::nrrd_elem_count(dnrrd);
  const size_t ssize = VolumeOps::nrrd_elem_count(snrrd);
  const size_t msize = VolumeOps::nrrd_elem_count(mnrrd);
  ASSERT(dsize == ssize && dsize == msize);
  ASSERT(dnrrd->nrrd_->type == snrrd->nrrd_->type);
  ASSERT(mnrrd->nrrd_->type == LabelNrrdType);

  label_type *mask = (label_type *)mnrrd->nrrd_->data;

  switch (snrrd->nrrd_->type) {
  case nrrdTypeChar:
    templated_mask_data((char *)dnrrd->nrrd_->data,
                        (char *)snrrd->nrrd_->data,
                        mask, ssize, mbit,
                        (char)newval);
    break;
  case nrrdTypeUChar:
    templated_mask_data((unsigned char *)dnrrd->nrrd_->data,
                        (unsigned char *)snrrd->nrrd_->data,
                        mask, ssize, mbit,
                        (unsigned char)newval);
    break;
  case nrrdTypeShort:
    templated_mask_data((short *)dnrrd->nrrd_->data,
                        (short *)snrrd->nrrd_->data,
                        mask, ssize, mbit,
                        (short)newval);
    break;
  case nrrdTypeUShort:
    templated_mask_data((unsigned short *)dnrrd->nrrd_->data,
                        (unsigned short *)snrrd->nrrd_->data,
                        mask, ssize, mbit,
                        (unsigned short)newval);
    break;
  case nrrdTypeInt:
    templated_mask_data((int *)dnrrd->nrrd_->data,
                        (int *)snrrd->nrrd_->data,
                        mask, ssize, mbit,
                        (int)newval);
    break;
  case nrrdTypeUInt:
    templated_mask_data((unsigned int *)dnrrd->nrrd_->data,
                        (unsigned int *)snrrd->nrrd_->data,
                        mask, ssize, mbit,
                        (unsigned int)newval);
    break;
  case nrrdTypeFloat:
    templated_mask_data((float *)dnrrd->nrrd_->data,
                        (float *)snrrd->nrrd_->data,
                        mask, ssize, mbit,
                        (float)newval);
    break;
  case nrrdTypeDouble:
    templated_mask_data((double *)dnrrd->nrrd_->data,
                        (double *)snrrd->nrrd_->data,
                        mask, ssize, mbit,
                        (double)newval);
    break;
  default:
    ASSERTFAIL("Unsupported data type.");
  }
}


NrrdDataHandle
VolumeOps::clear_nrrd(NrrdDataHandle &inh)
{
  memset(inh->nrrd_->data, 0, nrrd_data_size(inh));
  return inh;
}


NrrdDataHandle
VolumeOps::create_clear_nrrd(NrrdDataHandle &inh, unsigned int type)
{
  NrrdDataHandle nout = create_nrrd(inh, type);
  return clear_nrrd(nout);
}


NrrdDataHandle
VolumeOps::create_nrrd(NrrdDataHandle &inh, unsigned int type)
{
  Nrrd *src = inh->nrrd_;

  NrrdDataHandle nrrdh = new NrrdData();
  Nrrd *dst = nrrdh->nrrd_;

  nrrdBasicInfoCopy(dst, src, 0);
  nrrdAxisInfoCopy(dst, src, 0, 0);

  const size_t num = nrrd_elem_count(inh);

  if (!type) type = src->type;
  dst->type = type;

  switch (type) {
  case nrrdTypeChar: dst->data = new signed char[num]; break;
  case nrrdTypeUChar: dst->data = new unsigned char[num]; break;
  case nrrdTypeShort: dst->data = new signed short[num]; break;
  case nrrdTypeUShort: dst->data = new unsigned short[num]; break;
  case nrrdTypeInt: dst->data = new signed int[num]; break;
  case nrrdTypeUInt: dst->data = new unsigned int[num]; break;
  case nrrdTypeLLong: dst->data = new signed long long[num]; break;
  case nrrdTypeULLong: dst->data = new unsigned long long[num]; break;
  case nrrdTypeFloat: dst->data = new float[num]; break;
  case nrrdTypeDouble: dst->data = new double[num]; break;
  default: throw "unhandled type in NrrdVolume::get_clear_nrrd"; break;
  }

  return nrrdh;
}


bool
VolumeOps::compare_nrrd_info(NrrdDataHandle &nda, NrrdDataHandle &ndb)
{
  Nrrd *a = nda->nrrd_;
  Nrrd *b = ndb->nrrd_;
  
  if (a->dim != b->dim) return false;
  for (size_t i = 0; i < a->dim; i++)
  {
    if (a->axis[i].size != b->axis[i].size) return false;
  }

  // TODO: Compare transform also.  Just size will work but give wrong
  // results if the two nrrds are the same size but occupy different
  // spaces.

  return true;
}


NrrdDataHandle
VolumeOps::convert_to(NrrdDataHandle &ninh,
		      unsigned int type)
{
  if( ninh->nrrd_->type == type )
    return ninh;
  
  Nrrd *src = ninh->nrrd_;
  NrrdDataHandle nrrdh = create_nrrd(ninh, nrrdTypeFloat);
  Nrrd *dst = nrrdh->nrrd_;


  if (nrrdConvert(dst, src, type)) {
    char *err = biffGetDone(NRRD);
    cerr << string("Trouble resampling: ") + err
	 << "  input Nrrd: src->dim="<<src->dim<<"\n";
    free(err);
  }

  return nrrdh;
}


NrrdDataHandle
VolumeOps::bit_to_float(NrrdDataHandle &ninh,
                        unsigned int mask, float value)
{
  Nrrd *src = ninh->nrrd_;
  NrrdDataHandle nrrdh = create_nrrd(ninh, nrrdTypeFloat);
  Nrrd *dst = nrrdh->nrrd_;

  label_type *srcdata = (label_type *)src->data;
  float *dstdata = (float *)dst->data;
  const size_t num = nrrd_elem_count(ninh);
  for (size_t i = 0; i < num; ++i, ++srcdata, ++dstdata) {
    *dstdata = (*srcdata & mask) ? value : 0.0f;
  }

  return nrrdh;
}

NrrdStats
VolumeOps::nrrd_statistics(NrrdDataHandle &nrrdh, NrrdDataHandle &maskh, unsigned char bitmask)
{
  Nrrd *nrrd = nrrdh->nrrd_;
  Nrrd *mask = maskh->nrrd_;

  ASSERT(nrrd->dim > 3 && mask->dim > 3);
  ASSERT(nrrd->axis[0].size == mask->axis[0].size &&
         nrrd->axis[1].size == mask->axis[1].size &&
         nrrd->axis[2].size == mask->axis[2].size &&
         nrrd->axis[3].size == mask->axis[3].size);
  ASSERT(nrrd->type == nrrdTypeFloat);
  ASSERT(mask->type == nrrdTypeUChar);

  float *src = (float *)nrrd->data;
  unsigned char *test = (unsigned char *)mask->data;

  double mean = 0;
  double squared = 0;
  unsigned int n = 0;
  float min = AIR_POS_INF;
  float max = AIR_NEG_INF;

  const size_t size = nrrd_elem_count(nrrdh);
  for (size_t i = 0; i < size; ++i)
  {
    if (test[i]&bitmask) {
      mean += src[i];
      squared += src[i] * src[i];
      min = Min(min, src[i]);
      max = Max(max, src[i]);

      ++n;
    }
  }

  std::vector<float> v(n);
  for (size_t i = 0; i < size; ++i) {
    if (bitmask&test[i]) {
      v.push_back(src[i]);
    }
  }

  // Effective STL, item 31
  vector<float>::iterator begin(v.begin());
  vector<float>::iterator end(v.end());
  vector<float>::iterator goalPosition;
  
  goalPosition = begin + v.size() / 2;
  std::nth_element(begin, goalPosition, end);

  double spacingCalc[3], vecCalc[NRRD_SPACE_DIM_MAX];

  for (size_t j = 1; j <= 3; ++j) { 
    nrrdSpacingCalculate(nrrd, j, spacingCalc + (j-1), vecCalc); 
  }

  double volume = n * spacingCalc[0] * spacingCalc[1] * spacingCalc[2];

  mean = mean / n;
  const double deviation = sqrt(squared / n - mean * mean);
  return NrrdStats(min, max, mean, deviation, *goalPosition, n, volume);
}


void
VolumeOps::dice_32_into_8(NrrdDataHandle &src,
                          vector<NrrdDataHandle> &dst)
{
  size_t i;

  unsigned char *dstptr[4];
  unsigned int  *srcptr = (unsigned int *)(src->nrrd_->data);
  const size_t num = nrrd_elem_count(src);

  dst.clear();
  for (i = 0; i < 4; i++)
  {
    dst.push_back(create_nrrd(src, nrrdTypeUChar));
    dstptr[i] = (unsigned char *)(dst[i]->nrrd_->data);
  }

  for (i = 0; i < num; i++)
  {
    dstptr[0][i] = (unsigned char)(srcptr[i] >> 0) & 0xff;
    dstptr[1][i] = (unsigned char)(srcptr[i] >> 8) & 0xff;
    dstptr[2][i] = (unsigned char)(srcptr[i] >> 16) & 0xff;
    dstptr[3][i] = (unsigned char)(srcptr[i] >> 24) & 0xff;
  }
}

}
