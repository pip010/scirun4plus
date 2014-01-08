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
 *  ConvertFieldDataFromTensorsToIndices: Change a Field of indices (ints) into a Field or Tensors,
 *                      where the Tensor values are looked up in the
 *                      conductivity_table for each index
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   December 2002
 *
 */


#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Geometry/Tensor.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>


#include <vector>

namespace BioPSE {

using namespace SCIRun;

class ConvertFieldDataFromTensorsToIndices : public Module {

public:
  ConvertFieldDataFromTensorsToIndices(GuiContext *context);
  virtual ~ConvertFieldDataFromTensorsToIndices() {}
  virtual void execute();

#ifdef HAVE_HASH_MAP
  struct tensorhash
  {
    unsigned int operator()(const Tensor &t) const
    {
      unsigned char *s = (unsigned char *)(t.mat_);
      unsigned int tmp = 0;
#if defined(__ECC) || defined(_MSC_VER)
      hash_compare<unsigned int> h;
#else
      hash<unsigned int> h;
#endif // defined(__ECC) || defined(_MSC_VER)
      for (unsigned int i = 0; i < sizeof(double)*9; i++)
      {
        tmp = ( tmp << 5 ) - tmp + s[i];
      }
      return h(tmp);
    }

#if defined(__ECC) || defined(_MSC_VER)
    // These are particularly needed by ICC's hash stuff
    static const size_t bucket_size = 4;
    static const size_t min_buckets = 8;

    bool operator()(const Tensor &a, const Tensor &b) const
    {
      return (a < b);
    }
#endif // defined(__ECC) || defined(_MSC_VER)
  };
  struct tensorequal
  {
    //bool operator()(const Tensor &a, const Tensor &b) const
    bool operator()(const Tensor &a, const Tensor &b) const
    {
      return (a == b);
    }
  };
#if defined(__ECC) || defined(_MSC_VER)
  typedef hash_map<Tensor, int, tensorhash> tensor_map_type;
#else
  typedef hash_map<Tensor, int, tensorhash, tensorequal> tensor_map_type;
#endif // defined(__ECC) || defined(_MSC_VER)

#else
  // TODO: map is not a replacement for hash_map.  It's not-found
  // semantics are different.  No equal test.  Just implement our own
  // hash_map class if it's not there.
  struct tensorless
  {
    bool operator()(const Tensor &a, const Tensor &b) const
    {
      for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
          if( a.mat_[i][j] >= b.mat_[i][j])
            return false;
      return true;
    }
  };
  typedef std::map<Tensor, int, tensorless> tensor_map_type;
#endif

};


DECLARE_MAKER(ConvertFieldDataFromTensorsToIndices)


ConvertFieldDataFromTensorsToIndices::ConvertFieldDataFromTensorsToIndices(GuiContext *context)
  : Module("ConvertFieldDataFromTensorsToIndices", context, Filter, "Forward", "BioPSE")
{
}

void
ConvertFieldDataFromTensorsToIndices::execute()
{
  FieldHandle ifieldH;
  get_input_handle("TensorField", ifieldH);

  FieldInformation fi(ifieldH);
  fi.make_int();
  FieldHandle ofieldH = CreateField(fi,ifieldH->mesh());
  
  VField* src = ifieldH->vfield();
  VField* dst = ofieldH->vfield();
  dst->resize_values();
  
  std::vector<std::pair<std::string, Tensor> > conds;
  tensor_map_type tmap;
    
  size_type num_values = src->num_values();
  for (index_type idx=0; idx< num_values; idx++)
  {
    Tensor tensor;
    int    index;
    
    src->get_value(tensor,idx);
    const tensor_map_type::iterator loc = tmap.find(tensor);
    if (loc != tmap.end())
    {
      index = (*loc).second;
    }
    else
    {
      conds.push_back(std::pair<std::string, Tensor>("unknown", tensor));
      tmap[tensor] = static_cast<int>(conds.size() - 1);
      index = static_cast<int>(conds.size() - 1);
    }
  
    dst->set_value(index,idx);
  }
  
  dst->set_property("conductivity_table", conds, false);

  send_output_handle("IndexField", ofieldH);
}

} // End namespace BioPSE
