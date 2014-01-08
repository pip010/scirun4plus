/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2012 Scientific Computing and Imaging Institute,
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

#include <Core/Algorithms/Math/RootMeanSquareError/RootMeanSquareError.h>

#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixOperations.h>

namespace SCIRunAlgo {

bool
RootMeanSquareErrorAlgo::run(const MatrixHandle& inputHandle, const MatrixHandle& estimateHandle,
                             MatrixHandle& outputHandle)
{
  algo_start("RootMeanSquareError");
  std::string mode;
  get_option("mode", mode);

  double sum = 0;

  MatrixHandle diffMatrixHandle = inputHandle - estimateHandle;
  DenseMatrix* diffMatrix = diffMatrixHandle->dense();

  const index_type M = diffMatrix->nrows();
  const index_type N = diffMatrix->ncols();
  index_type size = M * N;

  // TODO: need better way to traverse rows & columns...
  // Use ParallelLinearAlgebra::norm?

  if (mode == "vector")
  {
    // column matrix
    for (index_type r = 0; r < M; ++r)
    {
      sum += diffMatrix->get(r, 0) * diffMatrix->get(r, 0);
    }      
  }
  else if (mode == "roworder")
  {
    for (index_type r = 0; r < M; ++r)
    {
      for (index_type c = 0; c < N; ++c)
      {
        sum += diffMatrix->get(r, c) * diffMatrix->get(r, c);
      }
    }
  }
  else if (mode == "columnorder")
  {
    for (index_type c = 0; c < N; ++c)
    {
      for (index_type r = 0; r < M; ++r)
      {
        sum += diffMatrix->get(r, c) * diffMatrix->get(r, c);
      }
    }
  }
  else
  {
    error("Invalid algorithm option.");
    return false;
  }

  sum = sqrt(sum/size);
  outputHandle = new DenseMatrix(1, 1);
  outputHandle->put(0, 0, sum);

  algo_end();
  return true;
}

}
