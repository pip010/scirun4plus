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

#include <Core/Algorithms/Math/AppendMatrix/AppendMatrices.h>

#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

namespace SCIRunAlgo {

using namespace SCIRun;

bool
AppendMatricesAlgo::
run(std::vector<MatrixHandle>& inputs,MatrixHandle& output)
{
  algo_start("AppendMatrices");
  std::string method = get_option("method");
  
  // check whether we have any matrices
  if (inputs.empty())
  {
    error("No input matrices");
    algo_end(); return (false);
  }
  
  // If there is one inout, nothing needs to be appended
  if (inputs.size() == 1)
  {
    output = inputs[0];
    algo_end(); return (true);
  }
  
  // separate algorithms for row and column appending
  if (method == "append_rows")
  {
    // Get size of matrix 1
    size_type m =inputs[0]->nrows();
    size_type n =inputs[0]->ncols();
    
    // check whether all matrices have the same number of columns
    for (size_t j=1; j<inputs.size();j++)
    {
      if (inputs[j]->ncols() != n)
      {
        error("A matrix with a different number of columns cannot be appended row wise");
        algo_end(); return (false);
      }
      m += inputs[j]->nrows();
    }

    // n and m will have new dimensions
    if (matrix_is::sparse(inputs[0]))
    {
      // sparse version
      SparseRowMatrixHandle srm = inputs[0]->sparse();
      size_type nnz = srm->get_nnz();
      
      // Figure out the total number of non zeros
      for (size_t j=1; j<inputs.size(); j++)
      {
        //Only merge sparse matrices with sparse matrices
        if (!matrix_is::sparse(inputs[j]))
        {
          error("The algorithm does not support appending matrices of different types");
          algo_end(); return (false);
        }
        srm = inputs[j]->sparse();
        nnz += srm->get_nnz();
      }

      // Setup new matrix, reserve its memory
      SparseRowMatrix::Data outputData(m+1, nnz);
    
      if (!outputData.allocated())
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);
      }
      
      index_type ms = 0;
      index_type me = 0;
      index_type nns = 0;
      index_type nne = 0;
      const SparseRowMatrix::Rows& rr = outputData.rows();
      const SparseRowMatrix::Columns& cc = outputData.columns();
      const SparseRowMatrix::Storage& vv = outputData.data();
      
      for (size_t j=0; j<inputs.size(); j++)
      {
        srm = inputs[j]->sparse();
        
        index_type* rows = srm->get_rows();
        index_type* columns = srm->get_cols();
        double* a = srm->get_vals();

        me = ms + srm->nrows();
                
        for (index_type r= ms; r<me; r++)
        {
          rr[r] = rows[r-ms]+nns;
        }
        ms = me;
        
        nne = nns+srm->get_nnz();
        for (index_type c=nns; c<nne; c++)
        {
          cc[c] = columns[c-nns];
          vv[c] = a[c-nns];
        }
        nns = nne;
      }
      
      rr[m] = nnz;
      
      output = new SparseRowMatrix(m, n, outputData, nnz);
      if (output.get_rep() == 0)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);          
      }
    }
    else if (matrix_is::dense(inputs[0]))
    {
      for (size_t j=1; j<inputs.size(); j++)
      {
        if (!matrix_is::dense(inputs[j]))
        {
          error("The algorithm does not support appending matrices of different types");
          algo_end(); return (false);
        }
      }

      // Create new matrix
      output = new DenseMatrix(m,n);
      if (output.get_rep() == 0)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);          
      }
     
      // copy all data sequentially into new matrix             
      size_type offset = 0;
      double* ndata = output->get_data_pointer();
      
      for (size_t j=0; j<inputs.size(); j++)
      {      
        double* data = inputs[j]->get_data_pointer();
        size_type size = inputs[j]->get_data_size();
        for (index_type j=0; j<size; j++)
        {
          ndata[offset+j] = data[j];
        }
        offset += size;
      }
    }
    else if (matrix_is::column(inputs[0]))
    {
      for (size_t j=1; j<inputs.size(); j++)
      {
        if (!(matrix_is::column(inputs[j])))
        {
          error("The algorithm does not support appending matrices of different types");
          algo_end(); return (false);
        }
      }

      // Create new matrix
      output = new ColumnMatrix(m);
      if (output.get_rep() == 0)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);          
      }
     
      // copy all data sequentially into new matrix             
      size_type offset = 0;
      double* ndata = output->get_data_pointer();
      
      for (size_t j=0; j<inputs.size(); j++)
      {      
        double* data = inputs[j]->get_data_pointer();
        size_type size = inputs[j]->get_data_size();
        for (index_type j=0; j<size; j++)
        {
          ndata[offset+j] = data[j];
        }
        offset += size;
      }
    }
    else
    {
      error("This algorithm is not support for this type");
      algo_end(); return (false);                
    }
    
    algo_end(); return (true);
  }
  else
  {
    // Get size of matrix 1
    size_type m = inputs[0]->nrows();
    size_type n = inputs[0]->ncols();
    
    // check whether all matrices have the same number of columns
    for (size_t j=1; j<inputs.size();j++)
    {
      if (inputs[j]->nrows() != m)
      {
        error("A Matrix with a different number of rows cannot be appended column wise");
        algo_end(); return (false);
      }
      n += inputs[0]->ncols();
    }
  

    // n and m will have new dimensions
    if (matrix_is::sparse(inputs[0]))
    {
      // sparse version
      SparseRowMatrix* srm = inputs[0]->sparse();
      size_type nnz = srm->get_nnz();

      SparseRowMatrix::Builder outputData;
      const SparseRowMatrix::Rows& rr = outputData.allocate_rows(m+1);
      
      if (!rr)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);
      }
            
      for (size_type j=0; j<m+1; j++) rr[j] = 0;
      
      std::vector<SparseRowMatrixHandle> srms(inputs.size());
      // Figure out the total number of non zeros
      for (size_t j=1; j<inputs.size(); j++)
      {
        //Only merge sparse matrices with sparse matrices
        if (!matrix_is::sparse(inputs[j]))
        {
          error("The algorithm does not support appending matrices of different types");
          algo_end(); return (false);
        }
        srms[j] = inputs[j]->sparse();
        nnz += srms[j]->get_nnz();
        
        index_type* rows = srms[j]->get_rows();
        for (size_type r=1; r<m+1; r++)
        {
          rr[r] += rows[r]-rows[r-1];
        }
      }
    
      size_type s = 0;
      for (size_type r=1;r<m+1;r++)
      {
        s += rr[r]; rr[r] = s;
      }
       
      // Setup new matrix, reserve its memory
      const SparseRowMatrix::Columns& cc = outputData.allocate_columns(nnz);
      const SparseRowMatrix::Storage& vv = outputData.allocate_data(nnz);
    
      if (!cc || !vv)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);
      }

      index_type k = 0;
      for (size_type r=0;r<m;r++)
      {
        for (size_t j=0; j < srms.size(); j++)
        {
          index_type* rows = srms[j]->get_rows();
          index_type* columns = srms[j]->get_cols();
          double* values = srms[j]->get_vals();
          index_type start = rows[r];
          index_type end = rows[r+1];
          for(index_type i= start; i< end; i++, k++)
          {
            cc[k] = columns[i];
            vv[k] = values[i];
          }
        }
      }
      
      output = new SparseRowMatrix(m, n, outputData.build(), nnz);
      if (output.get_rep() == 0)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);          
      }
    }
    else if (matrix_is::dense(inputs[0]))
    {
      for (size_t j=1; j<inputs.size(); j++)
      {
        if (!(matrix_is::dense(inputs[j])))
        {
          error("The algorithm does not support appending matrices of different types");
          algo_end(); return (false);
        }
      }

      // Create new matrix
      output = new DenseMatrix(m,n);
      if (output.get_rep() == 0)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);          
      }
     
      // copy all data sequentially into new matrix             
      size_type offset = 0;
      double* ndata = output->get_data_pointer();
      
      for (size_t j=0; j<inputs.size(); j++)
      {      
        double* data = inputs[j]->get_data_pointer();
        size_type mm = inputs[j]->ncols();
        index_type k = 0;
        for (index_type p=0; p<n; p++)
          for(index_type q=0; q<mm; q++, k++)
          {
            ndata[p*m+q+offset] = data[k];
          }
        offset += mm;  
      }
    }
    else if (matrix_is::column(inputs[0]))
    {
      for (size_t j=1; j<inputs.size(); j++)
      {
        if (!(matrix_is::column(inputs[j])))
        {
          error("The algorithm does not support appending matrices of different types");
          algo_end(); return (false);
        }
      }
      
      // Create new matrix
      output = new DenseMatrix(m,n);
      if (output.get_rep() == 0)
      {
        error("Could not allocate output matrix");
        algo_end(); return (false);          
      }
     
      // copy all data sequentially into new matrix             
      size_type offset = 0;
      double* ndata = output->get_data_pointer();
      
      for (size_t j=0; j<inputs.size(); j++)
      {      
        double* data = inputs[j]->get_data_pointer();
        index_type k = 0;
        for (index_type p=0; p<n; p++, k++)
        {
          ndata[p*m+offset] = data[k];
        }
        offset++;  
      }
    }
    else
    {
      error("This algorithm is not support for this type");
      algo_end(); return (false);                
    }
  }
  return(false);
}


} // end namespace SCIRunAlgo
