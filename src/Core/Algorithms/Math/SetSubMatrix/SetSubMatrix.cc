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

#include <Core/Algorithms/Math/SetSubMatrix/SetSubMatrix.h>

#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <algorithm>
#include <vector>

namespace SCIRunAlgo {

using namespace SCIRun;

bool
SetSubMatrixAlgo::
run(MatrixHandle input, 
    MatrixHandle data,
    MatrixHandle& output, 
    MatrixHandle rows, 
    MatrixHandle columns)
{
  algo_start("SetSubMatrix");

  size_type num_sel_rows = rows->get_data_size();
  double* sel_rows_ptr = rows->get_data_pointer();

  size_type num_sel_columns = columns->get_data_size();
  double* sel_columns_ptr = columns->get_data_pointer();
  
  std::vector<index_type> sel_rows;
  std::vector<index_type> sel_columns;
  try
  {
    sel_rows.resize(num_sel_rows);
    for (int p=0; p< num_sel_rows; p++) 
      sel_rows[p] = static_cast<index_type>(sel_rows_ptr[p]);

    sel_columns.resize(num_sel_columns);
    for (int p=0; p< num_sel_columns; p++) 
      sel_columns[p] = static_cast<index_type>(sel_columns_ptr[p]);
  }
  catch (...)
  {
    error("Could not allocate enough memory");
    algo_end(); return (false); 
  }
  
  bool ret = run(input,data,output,sel_rows,sel_columns);
  algo_end(); return(ret);
}

bool
SetSubMatrixAlgo::
run(MatrixHandle input, 
    MatrixHandle data,
    MatrixHandle& output, 
    std::vector<index_type>& rows, 
    std::vector<index_type>& columns)
{
  algo_start("SelectSubMatrixAlgo");

  if (input.get_rep() == 0)
  {
    error("No input matrix");
    algo_end(); return (false);
  }

  if (data.get_rep() == 0)
  {
    error("No data matrix");
    algo_end(); return (false);
  }


  if (rows.size() == 0)
  {
    error("No row indices given");
    algo_end(); return (false);  
  }

  if (columns.size() == 0)
  {
    error("No column indices given");
    algo_end(); return (false);  
  }

  size_type m = input->nrows();
  size_type n = input->ncols();
    
  for (size_t r=0; r<rows.size(); r++)
  {
    if (rows[r] >= static_cast<index_type>(m))
    {
      error("Selected row exceeds matrix dimensions");
      algo_end(); return (false);
    }
  }

  for (size_t r=0; r<columns.size(); r++)
  {
    if (columns[r] >= static_cast<index_type>(n))
    {
      error("Selected column exceeds matrix dimensions");
      algo_end(); return (false);
    }
  }

  if (rows.size() != data->nrows())
  {
    error("The number of rows in the data matrix does not match the number of indices.");
    algo_end(); return (false);
  }

  if (columns.size() != data->ncols())
  {
    error("The number of columns in the data matrix does not match the number of indices.");
    algo_end(); return (false);
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
    // Column transformation matrix
    std::vector<index_type> s(n,n);
    for (index_type r=0;
              r< static_cast<index_type>(columns.size()); r++) 
    {
      s[columns[r]] = r;
    }
    
    // Row transformation matrix
    std::vector<index_type> t(m,m);
    for (index_type r=0;
              r< static_cast<index_type>(rows.size()); r++) 
    {
      t[rows[r]] = r;
    }

    // Insert sparse into sparse
    SparseRowMatrix* sparseData = matrix_cast::as_sparse(data);
    if (sparseData)
    {
      index_type *drr = sparseData->get_rows();
      index_type *dcc = sparseData->get_cols();
      double *dvv = sparseData->get_vals();
    
      SparseRowMatrix::Builder outputData;
      const SparseRowMatrix::Rows& nrows = outputData.allocate_rows(input->nrows() + 1);

      size_type k = 0;
      nrows[0] = 0;
      
      // First determine the nnz
      for (index_type r=0; r<m; r++)
      {
        index_type cidx = rr[r];
        index_type cidx_end = rr[r+1];
        
        if (t[r] < m)
        {
          index_type didx = drr[rows[t[r]]];
          index_type didx_end = drr[rows[t[r]]+1];
          
          // determine size
          for (size_t j = cidx; j<cidx_end; j++)
          {
            if (s[cc[j]] < n) continue;
            k++;
          }
        
          k+= didx_end-didx;
        }
        else
        {
          k+=cidx_end-cidx;
        }

        nrows[r+1] = k;
      }

      const SparseRowMatrix::Columns& ncolumns = outputData.allocate_columns(k);
      const SparseRowMatrix::Storage& nvalues = outputData.allocate_data(k);
      
      k = 0;
      for (index_type r=0; r<m; r++)
      {
        index_type cidx = rr[r];
        index_type cidx_end = rr[r+1];

        if (t[r] < m)
        {
          index_type didx = drr[rows[t[r]]];
          index_type didx_end = drr[rows[t[r]]+1];
        
          // determine size
          for (size_t j = cidx; j<cidx_end; j++)
          {
            if (s[cc[j]] < n) continue;
            ncolumns[k] = cc[j];
            nvalues[k] = vv[j];
            k++;
          }
        
          for (size_t j = didx; j<didx_end;j++)
          { 
            ncolumns[k] = s[dcc[j]];
            nvalues[k] = dvv[j];
            k++;
          }
        }
        else
        {
          for (size_t j = cidx; j<cidx_end; j++)
          {
            ncolumns[k] = cc[j];
            nvalues[k] = vv[j];
            k++;
          }        
        }
      }      
    
      output = new SparseRowMatrix(m, n, outputData.build(), k, true);
      if (output.get_rep() == 0)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);          
      }
      
      algo_end(); return (true);
    }
    else
    {
      double *dvv = data->get_data_pointer();
    
      SparseRowMatrix::Builder outputData;
      const SparseRowMatrix::Rows& nrows = outputData.allocate_rows(input->nrows() + 1);
      size_type k = 0;
      nrows[0] = 0;
      
      size_type dm = data->nrows();
      size_type dn = data->ncols();
      
      // First determine the nnz
      for (index_type r=0; r<m; r++)
      {
        index_type cidx = rr[r];
        index_type cidx_end = rr[r+1];
        
        if (t[r] < m)
        {          
          // determine size
          for (size_t j = cidx; j<cidx_end; j++)
          {
            if (s[cc[j]] < n) continue;
            k++;
          }
        
          index_type offset = t[r]*n;
          for (size_t j =0; j<dn; j++)
          {
            if (dvv[j+offset]) k++;
          }
        }
        else
        {
          k+=cidx_end-cidx;
        }

        nrows[r+1] = k;
      }

      const SparseRowMatrix::Columns& ncolumns = outputData.allocate_columns(k);
      const SparseRowMatrix::Storage& nvalues = outputData.allocate_data(k);
      
      k = 0;
      for (index_type r=0; r<m; r++)
      {
        index_type cidx = rr[r];
        index_type cidx_end = rr[r+1];

        if (t[r] < m)
        {        
          // determine size
          for (size_t j = cidx; j<cidx_end; j++)
          {
            if (s[cc[j]] < n) continue;
            ncolumns[k] = cc[j];
            nvalues[k] = vv[j];
            k++;
          }
        
          index_type offset = t[r]*n;
          for (size_t j =0; j<dn; j++)
          {
            if (dvv[j+offset]) 
            {
              ncolumns[k] = columns[j];
              nvalues[k] = dvv[j+offset];
              k++;
            }
          }
        }
        else
        {
          for (size_t j = cidx; j<cidx_end; j++)
          {
            ncolumns[k] = cc[j];
            nvalues[k] = vv[j];
            k++;
          }        
        }
      }      
    
      output = new SparseRowMatrix(m, n, outputData.build(), k, true);
      if (output.get_rep() == 0)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);          
      }
      
      algo_end(); return (true);    
    
    }

  }
  else
  {
    MatrixHandle mat = input->dense();
    
    if (mat.get_rep() == 0)
    {
      error("Could not convert matrix into dense matrix");
      algo_end(); return (false);    
    }
    
    
    output = mat->clone();
    
    if (output.get_rep() == 0)
    {
      error("Could not allocate output matrix");
      algo_end(); return (false);
    }
    
    double* src = data->get_data_pointer();
    double* dst = output->get_data_pointer(); 
    
    if (dst==0 || src == 0)
    {
      error("Could not allocate output matrix");
      algo_end(); return (false);
    }
    
    size_type m = output->nrows();
    size_type n = output->ncols();
    
    size_type nr = static_cast<size_type>(rows.size());
    size_type nc = static_cast<size_type>(columns.size());
    
    if (matrix_is::dense(data) || matrix_is::column(data))
    {
      for (index_type p=0;p<nr;p++)
      {
        for (index_type q=0;q<nc;q++)
        {
          dst[rows[p]*n + columns[q]] = src[p*nc + q];
        }
      }
    }
    else
    {
      for (index_type p=0;p<nr;p++)
      {
        for (index_type q=0;q<nc;q++)
        {
          dst[rows[p]*n + columns[q]] = data->get(p,q);
        }
      }    
    }

    algo_end(); return (true);
  }
}

} // end namespace SCIRunAlgo
