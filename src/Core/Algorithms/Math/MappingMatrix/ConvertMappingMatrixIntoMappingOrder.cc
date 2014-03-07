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

#include <Core/Algorithms/Math/MappingMatrix/ConvertMappingMatrixIntoMappingOrder.h>

#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

namespace SCIRunAlgo {

using namespace SCIRun;

bool
ConvertMappingMatrixIntoMappingOrderAlgo::
run(MatrixHandle input, std::vector<index_type>& output)
{
  algo_start("ConvertMappingMatrixIntoMappingOrder");
  
  if (input.get_rep() == 0)
  {
    error("No input matrix");
    algo_end(); return (false);
  }
  
  MatrixHandle sparse = input->sparse();
  if (sparse.get_rep() == 0)
  {
    error("Could not convert mapping matrix into sparse matrix");
    algo_end(); return (false);
  }
  
  SparseRowMatrix* spr = matrix_cast::as_sparse(sparse);
  index_type* rr = spr->get_rows();
  index_type* cc = spr->get_cols();
  size_type nrows = spr->nrows();
  size_type ncols = spr->ncols();
  
  output.clear();
  output.reserve(nrows);
  for (index_type r=0; r<nrows; r++)
  {
    if (rr[r] != rr[r+1])
    {
      index_type c = cc[rr[r]];
      if (c < 0 || c >= ncols)
      {
        error("Mapping matrix was invalid");
        algo_end(); return (false);     
      }
      output.push_back(static_cast<unsigned int>(c)); 
    }
  }
  
  algo_end(); return (true);
}


bool
ConvertMappingMatrixIntoMappingOrderAlgo::
run(MatrixHandle input, MatrixHandle& output)
{
  algo_start("ConvertMappingMatrixIntoMappingOrder");
  
  if (input.get_rep() == 0)
  {
    error("No input matrix");
    algo_end(); return (false);
  }
  
  MatrixHandle sparse = input->sparse();
  if (sparse.get_rep() == 0)
  {
    error("Could not convert mapping matrix into sparse matrix");
    algo_end(); return (false);
  }
  
  SparseRowMatrix* spr = matrix_cast::as_sparse(sparse);
  index_type* rr = spr->get_rows();
  index_type* cc = spr->get_cols();
  size_type nrows = spr->nrows();
  size_type ncols = spr->ncols();
  
  output = new DenseMatrix(nrows,1);
  if (output.get_rep() == 0)
  {
    error("Could not allocate output matrix");
    algo_end(); return (false);
  }  

  double* data = output->get_data_pointer();
    
  for (index_type r=0; r<nrows; r++)
  {
    data[r]=0.0; 
    if (rr[r] != rr[r+1])
    {
      index_type c = cc[rr[r]];
      if (c < 0 || c >= ncols)
      {
        error("Mapping matrix was invalid");
        algo_end(); return (false);     
      }
      data[r]=static_cast<double>(c); 
    }
  }
  
  algo_end(); return (true);
}


} // end namespace SCIRunAlgo


