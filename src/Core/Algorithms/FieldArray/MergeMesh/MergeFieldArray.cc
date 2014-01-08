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

//! Class definition of this one
#include <Core/Algorithms/Fields/MergeFields/JoinFields.h>
#include <Core/Algorithms/FieldArray/MergeMesh/MergeFieldArray.h>

//! Need to find out what type of field we are dealing with
#include <Core/Datatypes/FieldInformation.h>

//! SearchGridT
#include <Core/Containers/SearchGridT.h>

namespace SCIRunAlgo {


bool 
MergeFieldArrayAlgo::run(FieldArrayHandle input, FieldHandle& output)
{
  // We do not do progress reporting, hence the lower level algorithm will do this
  algo_start("MergeFieldArray");
  
  // Get array of fields
  std::vector<FieldHandle>& inputs = input->array();

  JoinFieldsAlgo algo;
  // Copy parameters
  algo.set_progress_reporter(get_progress_reporter());
  algo.set_bool("merge_nodes",get_bool("merge_nodes"));
  algo.set_scalar("tolerance",get_scalar("tolerance"));
  algo.set_bool("match_node_values",get_bool("match_node_values"));
  algo.set_bool("make_no_data",get_bool("make_no_data"));

  // run field joiner
  bool ret = algo.run(inputs,output);
  algo_end();
  
  return (ret);
}

} // end namespace SCIRunAlgo
