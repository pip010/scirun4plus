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

#include <Core/Algorithms/Math/MappingMatrix/ConvertMappingOrderIntoMappingMatrix.h>

#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>

namespace SCIRunAlgo {

using namespace SCIRun;
 
bool
ConvertMappingOrderIntoMappingMatrixAlgo::
run(MatrixHandle input, 
    MatrixHandle& output, 
    size_type ncols)
{
  algo_start("ConvertMappingOrderIntoMappingMatrix");

  if (input.get_rep() == 0)
  {
    error("No input matrix");
    algo_end(); return (false);
  }  

  size_type size = input->get_data_size();
  double*           data = input->get_data_pointer(); 

  std::vector<index_type> order;
  try
  {
    order.resize(size);
    for (index_type p=0; p<size; p++) order[p] = static_cast<index_type>(data[p]);
  }
  catch (...)
  {
    error("Could not allocate output matrix");
    algo_end(); return (false);    
  }
 
  algo_end(); return(run(order,output,ncols));
}


bool
ConvertMappingOrderIntoMappingMatrixAlgo::
run(std::vector<index_type>& input, 
    MatrixHandle& output, 
    size_type ncols)
{
  algo_start("ConvertMappingOrderIntoMappingMatrix");

  bool transpose = get_bool("transpose");

  size_type nrows = input.size();  
  if (ncols == 0) ncols = nrows;

  SparseRowMatrix::Data outputData(nrows+1, nrows);
  const SparseRowMatrix::Rows& rr = outputData.rows();
  const SparseRowMatrix::Columns& cc = outputData.columns();
  const SparseRowMatrix::Storage& vv = outputData.data();

  if (!outputData.allocated())
  {
    error("Could not allocate output matrix");
    algo_end(); return (false);      
  }

  for (index_type r=0; r<nrows+1; r++) 
    rr[r] = r;
  for (index_type c=0; c<nrows; c++)
  {
    cc[c] = input[c];
    vv[c] = 1.0;
  }

  output = new SparseRowMatrix(nrows,ncols,outputData,nrows);
  if (output.get_rep() == 0)
  {
    error("Could not allocate output matrix");
    algo_end(); return (false);  
  }

  if (transpose)
  {
    MatrixHandle tmp = output; 
    output = tmp->make_transpose();
  }
  
  algo_end(); return (true);
} 
 
 
} // end namespace SCIRunAlgo


