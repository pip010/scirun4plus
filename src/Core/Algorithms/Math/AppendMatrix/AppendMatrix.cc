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

#include <Core/Algorithms/Math/AppendMatrix/AppendMatrix.h>

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

namespace SCIRunAlgo {

using namespace SCIRun;

bool
AppendMatrixAlgo::
run(MatrixHandle input,MatrixHandle& output,
    MatrixHandle append)
{
  std::vector<index_type> dummy;
  return(run(input,output,append,dummy));
}

bool
AppendMatrixAlgo::
run(MatrixHandle input,MatrixHandle& output,
    MatrixHandle append,std::vector<index_type>& indices)
{
  algo_start("AppendMatrix");
  std::string method = get_option("method");
 
  if (method == "append_rows")
  {
    if (input.get_rep() == 0)
    {
      if (append.get_rep() == 0)
      {
        output = 0;
        warning("Base matrix and matrix to append are empty");
        algo_end(); return (true);
      }
    
      output = append;
      size_type m = append->nrows();
      indices.resize(m);
      for (index_type r=0; r<m; r++) indices[r] = r;
      algo_end(); return (true);
    }

    if (append.get_rep() == 0)
    {
      output = input;
      indices.clear();
      warning("Matrix to append is empty");
      algo_end(); return (true);
    } 
    
    size_type m = input->nrows();
    size_type n = input->ncols();
    
    size_type am = append->nrows();
    size_type an = append->ncols();
          
    if (an != n)
    {
      error("The number of columns in input matrix is not equal to number of columns in row matrix");
      algo_end(); return (false);
    }
    
    index_type newm = m+am;
    index_type newn = n;
    
    SparseRowMatrix* sparse = matrix_cast::as_sparse(input);
    if (sparse)
    {
      SparseRowMatrixHandle a = append->sparse();
      if (a.get_rep() == 0)
      {
        error("Could not convert matrix to sparse matrix");
        algo_end(); return (false);      
      }
    
      index_type *rr = sparse->get_rows();
      index_type *cc = sparse->get_cols();
      double *vv = sparse->get_vals();
    
      index_type *arr = a->get_rows();
      index_type *acc = a->get_cols();
      double *avv = a->get_vals();
   
      if (!sparse->is_valid() || !a->is_valid())
      {
        error("Sparse matrix is invalid");
        algo_end(); return (false);      
      }   
    
      size_type newnnz = rr[m]+arr[am];
      SparseRowMatrix::Data outputData(newm+1, newnnz);
      const SparseRowMatrix::Rows& nrr = outputData.rows();
      const SparseRowMatrix::Columns& ncc = outputData.columns();
      const SparseRowMatrix::Storage& nvv = outputData.data();
      
      if (!outputData.allocated())
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);      
      }

      for (index_type r=0;r<m;r++) nrr[r] = rr[r];
      for (index_type r=0;r<am;r++) nrr[m+r] = arr[r]+rr[m];
      nrr[newm] = newnnz;
      
      size_type nnz = rr[m];
      for (index_type r=0;r<nnz;r++) { ncc[r] = cc[r]; nvv[r] = vv[r]; }
      size_type annz = arr[am];
      for (index_type r=0;r<annz;r++) { ncc[r+nnz] = acc[r]; nvv[r+nnz] = avv[r]; } 
      
      indices.resize(newm-m);
      for (index_type r=m;r<newm;r++) indices[r-m]=r;
      
      output = new SparseRowMatrix(newm, newn, outputData, newnnz);
      if (output.get_rep() == 0)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);          
      }
      
      algo_end(); return (true);
    }
    else
    {
      MatrixHandle i = input->dense();
      MatrixHandle a = append->dense();
      if (a.get_rep() == 0 || i.get_rep() == 0)
      {
        error("Could not convert matrix to dense matrix");
        algo_end(); return (false);      
      }    
    
      output = new DenseMatrix(newm,newn);
      
      double *iptr = i->get_data_pointer();
      double *aptr = a->get_data_pointer();
      double *optr = output->get_data_pointer();
      
      if (aptr == 0 || iptr == 0 ||optr == 0)
      {
        error("Could not convert matrix to dense matrix");
        algo_end(); return (false);      
      }    
      
      index_type mn = m*n;
      for (index_type r=0; r<mn;r++) optr[r] = iptr[r];
      optr += mn;
      mn = am*an;
      for (index_type r=0; r<mn;r++) optr[r] = aptr[r]; 
    
      indices.resize(newm-m);
      for (index_type r=m;r<newm;r++) indices[r-m]=r;
      
      algo_end(); return (true);  
    }
  
  
  }
  else
  {
    if (input.get_rep() == 0)
    {
      if (append.get_rep() == 0)
      {
        output = 0;
        warning("MatrixAppendColumns: Base matrix and matrix to append are empty");
        return (true);
      }
    
      output = append;
      size_type n = append->ncols();
      indices.resize(n);
      for (index_type r=0; r<n; r++) indices[r] = r;
      algo_end(); return (true);
    }


    if (append.get_rep() == 0)
    {
      output = input;
      indices.clear();
      algo_end(); return (true);
    }  
    
    size_type m = input->nrows();
    size_type n = input->ncols();
    
    size_type am = append->nrows();
    size_type an = append->ncols();
    
    if (am != m)
    {
      error("MatrixAppendColumns: The number of rows in input matrix is not equal to number of rows in column matrix");
      algo_end(); return (false);
    }
    
    size_type newm = m;
    size_type newn = n+an;
    
    SparseRowMatrix* sparse = matrix_cast::as_sparse(input);
    if (sparse)
    {
      SparseRowMatrixHandle a = append->sparse();
      if (a.get_rep() == 0)
      {
        error("MatrixAppendColumns: Could not convert matrix to sparse matrix");
        algo_end(); return (false);      
      }
    
      index_type *rr = sparse->get_rows();
      index_type *cc = sparse->get_cols();
      double *vv = sparse->get_vals();
    
      index_type *arr = a->get_rows();
      index_type *acc = a->get_cols();
      double *avv = a->get_vals();
   
      if (!sparse->is_valid() || !a->is_valid())
      {
        error("Sparse matrix is invalid");
        algo_end(); return (false);      
      }   
    
      size_type newnnz = rr[m]+arr[am];
      SparseRowMatrix::Data outputData(newm+1, newnnz);
      const SparseRowMatrix::Rows& nrr = outputData.rows();
      const SparseRowMatrix::Columns& ncc = outputData.columns();
      const SparseRowMatrix::Storage& nvv = outputData.data();

      if (!outputData.allocated())
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);      
      }

      index_type k = 0;
      for (index_type r=0;r<m;r++) 
      {
        nrr[r] = k;
        for (index_type q=rr[r];q<rr[r+1]; q++)
        {
          ncc[k] = cc[q];
          nvv[k] = vv[q];
          k++;
        }
        for (index_type q=arr[r];q<arr[r+1]; q++)
        {
          // *** bug fix *** - Peter Johnston 18/12/12
          // Offset was missing in the column labelling of the sparse matrix format
          ncc[k] = acc[q] + n;
          nvv[k] = avv[q];
          k++;
        }    
      }
      nrr[m] = k;
      
      indices.resize(newn-n);
      for (index_type r=n;r<newn;r++) indices[r-n] = r;

      output = new SparseRowMatrix(newm, newn, outputData, newnnz);
      if (output.get_rep() == 0)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);          
      }
          
      algo_end(); return (true);
    }
    else
    {
      MatrixHandle i = input->dense();
      MatrixHandle a = append->dense();
      if (a.get_rep() == 0 || i.get_rep() == 0)
      {
        error("Could not convert matrix to dense matrix");
        algo_end(); return (false);      
      }    
    
      output = new DenseMatrix(newm,newn);
      
      double *iptr = i->get_data_pointer();
      double *aptr = a->get_data_pointer();
      double *optr = output->get_data_pointer();
      
      if (aptr == 0 || iptr == 0 ||optr == 0)
      {
        error("Could not convert matrix to dense matrix");
        algo_end(); return (false);      
      }    
      
      for (index_type r=0; r<m; r++)
      {
        for (index_type q=0;q<n;q++) 
        { 
          optr[q] = iptr[q];
        }
        optr += n;
        iptr += n;
        for (index_type q=0;q<an;q++) 
        { 
          optr[q] = aptr[q];
        }
        optr += an;
        aptr += an;
      }
      
      indices.resize(newn-n);
      for (index_type r=n;r<newn;r++) indices[r-n] = r;
      
      algo_end(); return (true);
    }
  }
}


} // end namespace SCIRunAlgo
