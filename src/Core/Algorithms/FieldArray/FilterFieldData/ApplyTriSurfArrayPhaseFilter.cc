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

#include <Core/Algorithms/Fields/FilterFieldData/TriSurfPhaseFilter.h>
#include <Core/Algorithms/Fields/MergeFields/JoinFields.h>

#include <Core/Algorithms/FieldArray/FilterFieldData/ApplyTriSurfArrayPhaseFilter.h>

namespace SCIRunAlgo {

using namespace SCIRun;

bool 
ApplyTriSurfArrayPhaseFilterAlgo::
run(FieldArrayHandle& input, 
    FieldArrayHandle& output, 
    FieldArrayHandle& phaseline, 
    FieldArrayHandle& phasepoint)
{
  algo_start("ApplyTriSurfArrayPhaseFilter");
  
  std::vector<FieldHandle>& inputs = input->array();
  output = new FieldArray;
  phaseline = new FieldArray;
  phasepoint = new FieldArray;
  
  std::vector<FieldHandle>& outputs = output->array(); 
  std::vector<FieldHandle>& phaselines = phaseline->array(); 
  std::vector<FieldHandle>& phasepoints = phasepoint->array(); 

  size_t size = inputs.size();
  outputs.resize(size);
  phaselines.resize(size);
  phasepoints.resize(size);

  TriSurfPhaseFilterAlgo algo;
  JoinFieldsAlgo jalgo;
  algo.set_progress_reporter(get_progress_reporter());
  jalgo.set_progress_reporter(get_progress_reporter());
  
  jalgo.set_bool("merge_nodes",true);
  
  for (size_t j=0; j<size;j++)
  {
    if(!(algo.run(inputs[j],outputs[j],phaselines[j],phasepoints[j])))
    {
      algo_end();
      return (false);
    }
    std::vector<FieldHandle> handles(1);
    handles[0] = phaselines[j];
    jalgo.run(handles,phaselines[j]);
    update_progress(j,size);
  }
  
  algo_end();
  return (true);
}

}
