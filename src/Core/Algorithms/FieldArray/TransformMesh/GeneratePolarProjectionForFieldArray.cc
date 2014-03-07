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

#include <Core/Algorithms/Fields/TransformMesh/GeneratePolarProjection.h>
#include <Core/Algorithms/FieldArray/TransformMesh/GeneratePolarProjectionForFieldArray.h>

namespace SCIRunAlgo {

using namespace SCIRun;

bool 
GeneratePolarProjectionForFieldArrayAlgo::
run(FieldArrayHandle input, 
    FieldArrayHandle& output)
{
  algo_start("GeneratePolarProjectionForFieldArray");
  
  std::vector<FieldHandle>& inputs = input->array();
  output = new FieldArray;
  
  std::vector<FieldHandle>& outputs = output->array(); 

  size_t size = inputs.size();
  outputs.resize(size);

  GeneratePolarProjectionAlgo algo;
  algo.set_progress_reporter(get_progress_reporter());
  double spacing = get_scalar("spacing");
  Vector direction = get_vector("direction");
  Point origin = get_point("origin");
  algo.set_vector("direction",direction);
  algo.set_point("origin",origin);
  
  direction.normalize();
  
  for (size_t j=0; j<size;j++)
  {
    if(!(algo.run(inputs[j],outputs[j])))
    {
      algo_end();
      return (false);
    }
    
    // Translate field
    Vector trans = j*spacing*direction;
    VMesh* mesh = outputs[j]->vmesh();
    VMesh::size_type num_nodes = mesh->num_nodes();
    for (VMesh::Node::index_type idx; idx < num_nodes; idx++)
    {
      Point p;
      mesh->get_point(p,idx);
      p += trans;
      mesh->set_point(p,idx);
    }
    
    update_progress(j,size);
  }
  
  algo_end();
  return (true);
}

}
