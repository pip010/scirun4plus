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


#include <Core/Algorithms/Math/SortMatrix/SortMatrixColumns.h>
#include <Core/Algorithms/Math/SelectMatrix/SelectMatrixColumns.h>

#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <algorithm>
#include <vector>

namespace SCIRunAlgo {

using namespace SCIRun;


class SortAscDenseMatrix : public std::binary_function<index_type,index_type,bool>
{
  public:
    SortAscDenseMatrix(DenseMatrix* mat)
    {
      m_ = static_cast<index_type>(mat->nrows());
      n_ = static_cast<index_type>(mat->ncols());
      data_ = mat->get_data_pointer();
    }
    
    bool operator()(index_type i1, index_type i2)
    {
      for (index_type p=0; p<m_; p++)
      {
        double val1 = data_[i1+p*n_];
        double val2 = data_[i2+p*n_];
        if (val1 == val2) continue;
        return (val1 < val2);
      }
      return (false);        
    }

  private:
    index_type m_;
    index_type n_;
    double*      data_;
};

class SortDescDenseMatrix : public std::binary_function<index_type,index_type,bool>
{
  public:
    SortDescDenseMatrix(DenseMatrix* mat)
    {
      m_ = static_cast<index_type>(mat->nrows());
      n_ = static_cast<index_type>(mat->ncols());
      data_ = mat->get_data_pointer();
    }
    
    bool operator()(index_type i1, index_type i2)
    {
      for (index_type p=0; p<m_; p++)
      {
        double val1 = data_[i1+p*n_];
        double val2 = data_[i2+p*n_];
        if (val1 == val2) continue;
        return (val1 > val2);
      }
      return (false);        
    }

  private:
    index_type m_;
    index_type n_;
    double*      data_;
};
 


bool 
SortMatrixColumnsAlgo::
get_sortorder(MatrixHandle input, std::vector<index_type>& order)
{
  algo_start("SortMatrixColumns");
  bool asc = get_bool("ascending");

  if (input.get_rep() == 0)
  {
    error("No input matrix");
    algo_end(); return (false); 
  }

  if (!(matrix_is::dense(input)||matrix_is::column(input)))
  {
    error("This function has only been implemented for DenseMatrix and ColumnMatrix");
    algo_end(); return (false);
  }

  size_type ncols = input->ncols();
  order.resize(ncols);
  for (index_type p=0; p<ncols; p++) order[p] = p;
  
  if (matrix_is::dense(input))
  {
    DenseMatrix* dm = matrix_cast::as_dense(input);
    if (asc == true)
      std::sort(order.begin(),order.end(),SortAscDenseMatrix(dm));
    else
      std::sort(order.begin(),order.end(),SortDescDenseMatrix(dm));      
  }

  if (matrix_is::column(input))
  {
    // only one column, so this one is obvious
    order.resize(1); order[0] = 0;
  }

  algo_end(); return (true);
}


bool 
SortMatrixColumnsAlgo::
get_sortorder(MatrixHandle input, MatrixHandle& output)
{
  algo_start("SortMatrixColumns");

  if (input.get_rep() == 0)
  {
    error("No input matrix");
    algo_end(); return (false); 
  }
  
  size_type ncols = input->ncols();
  std::vector<index_type> order;
  if (!(get_sortorder(input,order))) 
  { return (false); }

  output = new DenseMatrix(ncols,1);
  if (output.get_rep() == 0)
  {
    error("Could not allocate output matrix");
    algo_end(); return (false);
  }

  double* data = output->get_data_pointer();
  for (index_type p=0; p<ncols; p++) data[p] = static_cast<double>(order[p]);
  
  algo_end(); return (true);
}


bool 
SortMatrixColumnsAlgo::
run(MatrixHandle input, MatrixHandle& output)
{
  algo_start("SortMatrixColumns");

  std::vector<index_type> order;
  if(!(get_sortorder(input,order)))
  {
    algo_end(); return (false);
  }  
  
  SelectMatrixColumnsAlgo algo;
  algo.set_progress_reporter(this);
  if(!(algo.run(input,output,order)))
  {
    algo_end(); return (false);
  }

  algo_end(); return(true);
}

} // end namespace SCIRunAlgo

