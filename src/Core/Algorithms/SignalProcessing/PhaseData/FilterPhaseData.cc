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

//! Datatypes that the algorithm uses
#include <Core/Algorithms/SignalProcessing/PhaseData/FilterPhaseData.h>
#include <Core/Datatypes/MatrixTypeConverter.h>
#include <Core/Datatypes/DenseMatrix.h>

namespace SCIRunAlgo {

using namespace SCIRun;

bool
FilterPhaseDataAlgo::run(MatrixHandle input, MatrixHandle& output)
{
  algo_start("FilterPhaseData");
  
  if (!(matrix_is::dense(input)))
  {
    error("The input matrix needs to be a dense matrix");
    algo_end(); return (false);
  }
  
  int filter_length = get_int("filter_length");
  
  size_type n,m;
  m = input->nrows();
  n = input->ncols();
  
  output = new DenseMatrix(m,n);
  if (output.get_rep() == 0)
  {
    error("Could not allocate output matrix");
    algo_end(); return (false);
  }
  
  double* data_in = input->get_data_pointer();
  double* data_out = output->get_data_pointer();

  size_type flength = ((filter_length-1)/2);
  double val;
  double mphase;
  
  int cnt =0;
  
  for (index_type row=0; row < m; row++)
  {
    double* din = data_in + n*row;
    double* dout = data_out + n*row;
    for (index_type col=0; col<n; col++)
    {
      index_type start = col-flength;
      index_type end = (col+flength)+1;
      if (start < 0) start = 0;
      if (end > n) end = n;
      
      val = 0.0;
      mphase = din[col]-M_PI;
      double dphase;
      for (index_type k=start;k<end;k++)
      {
        dphase = fmod(din[k]-mphase,2*M_PI);
        while (dphase < 0) dphase += 2*M_PI;
        val += dphase;
      }
      val = val/static_cast<double>(end-start);
      val += mphase;
      
      dout[col] = val;
    }
    cnt++; if(cnt == 10) { cnt = 0; update_progress(row,m); }
  }
  
  algo_end(); return (true);
}
    
} // end namespace SCIRunAlgo
