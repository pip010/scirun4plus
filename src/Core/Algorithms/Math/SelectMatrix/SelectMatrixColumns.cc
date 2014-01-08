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

#include <Core/Algorithms/Math/SelectMatrix/SelectMatrixColumns.h>

#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <algorithm>
#include <vector>

namespace SCIRunAlgo {

using namespace SCIRun;

bool
SelectMatrixColumnsAlgo::
run(MatrixHandle input, MatrixHandle& output, MatrixHandle columns)
{
  algo_start("SelectMatrixColumns");

  size_type num_sel_columns = columns->get_data_size();
  double* sel_columns_ptr = columns->get_data_pointer();

  std::vector<index_type> sel_columns;
  try
  {
    sel_columns.resize(num_sel_columns);
    for (index_type p=0; p< num_sel_columns; p++) 
      sel_columns[p] = static_cast<index_type>(sel_columns_ptr[p]);
  }
  catch (...)
  {
    error("Could not allocate enough memory");
    algo_end(); return (false); 
  }  
  
  bool ret = run(input,output,sel_columns);
  algo_end(); return(ret);
}

bool
SelectMatrixColumnsAlgo::
run(MatrixHandle input, MatrixHandle& output, 
    std::vector<index_type>& columns)
{
  algo_start("SelectMatrixColumns");

  if (input.get_rep() == 0)
  {
    error("No input matrix");
    algo_end(); return (false);
  }

  if (columns.size() == 0)
  {
    error("No row indices given");
    algo_end(); return (false);  
  }
  
  size_type m = input->nrows();
  size_type n = input->ncols();
    
  for (size_t r=0; r<columns.size(); r++)
  {
    if (columns[r] >= static_cast<index_type>(n))
    {
      error("Selected column exceeds matrix dimensions");
      algo_end(); return (false);
    }
  }
  
  SparseRowMatrix* sparse = matrix_cast::as_sparse(input);
  if (sparse)
  {
    if (!sparse->is_valid())
    {
      error("Sparse matrix is invalid");
      algo_end(); return (false);      
    }
   
    index_type *rr = sparse->get_rows();
    index_type *cc = sparse->get_cols();
    double *vv = sparse->get_vals();
    std::vector<size_t> s(n, n);
    for (size_t r=0; r< columns.size(); r++) s[columns[r]] = r;
      
    index_type k =0;
    for (index_type r=0; r<m; r++)
    {
      for (index_type q=rr[r]; q<rr[r+1]; q++)
      {
        if (s[cc[q]] < (size_t)n) k++;
      }
    }
    
    SparseRowMatrix::Data outputData(m+1, k);
    if (!outputData.allocated())
    {
      error("Could not allocate output matrix");
      algo_end(); return (false);      
    }    
    const SparseRowMatrix::Rows& nrr = outputData.rows();
    const SparseRowMatrix::Columns& ncc = outputData.columns();
    const SparseRowMatrix::Storage& nvv = outputData.data();

    k = 0;
    for (int r=0; r<m; r++)
    {
      nrr[r] =k;
      for (index_type q=rr[r]; q<rr[r+1]; q++)
      {
        if (s[cc[q]] < (size_t)n)
        {
          ncc[k] = s[cc[q]];
          nvv[k] = vv[q];
          k++;
        }
      }
    }    
    nrr[m] = k;

    output = new SparseRowMatrix(m, columns.size(), outputData, k);
    if (output.get_rep() == 0)
    {
      error("MatrixSelectRows: Could not allocate output matrix");
      algo_end(); return (false);          
    }
		        
    algo_end(); return (true);
  }
  else
  {
    MatrixHandle mat = input->dense();
    
    try
    {
      output = new DenseMatrix(m,columns.size());
    }
    catch (...)
    {
      error("Could not allocate output matrix");
      algo_end(); return (false);
    }
      
    if (mat.get_rep() == 0)
    {
      error("Could not convert matrix into dense matrix");
      algo_end(); return (false);    
    }
    
    double* src = mat->get_data_pointer();
    double* dst = output->get_data_pointer(); 
    
    if (dst==0 || src == 0)
    {
      error("Could not allocate output matrix");
      algo_end(); return (false);
    }
    
    size_type nc = static_cast<size_type>(columns.size());
    
    for (index_type p=0;p<m;p++)
    {
      for (index_type q=0;q<nc;q++)
      {
        dst[p*nc + q] = src[p*n +columns[q]];
      }
    }
			
    algo_end(); return (true);
  }
}


} // end namespace SCIRunAlgo

