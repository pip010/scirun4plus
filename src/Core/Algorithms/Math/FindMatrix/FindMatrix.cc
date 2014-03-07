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

#include <Core/Algorithms/Math/FindMatrix/FindMatrix.h>
#include <Core/Algorithms/Math/FindMatrix/FindZeroMatrixRows.h>
#include <Core/Algorithms/Math/FindMatrix/FindZeroMatrixColumns.h>
#include <Core/Algorithms/Math/FindMatrix/FindNonZeroMatrixRows.h>
#include <Core/Algorithms/Math/FindMatrix/FindNonZeroMatrixColumns.h>

namespace SCIRunAlgo {

using namespace SCIRun;

bool 
FindMatrixAlgo::
run(MatrixHandle input, std::vector<index_type>& output)
{
  std::string method = get_option("method");
  if (method == "zero_rows")
  {
    FindZeroMatrixRowsAlgo algo;
    algo.set_progress_reporter(this);
    return(algo.run(input,output));
  }
  if (method == "zero_columns")
  {
    FindZeroMatrixColumnsAlgo algo;
    algo.set_progress_reporter(this);
    return(algo.run(input,output));
  }
  if (method == "nonzero_rows")
  {
    FindNonZeroMatrixRowsAlgo algo;
    algo.set_progress_reporter(this);
    return(algo.run(input,output));
  }
  if (method == "nonzero_columns")
  {
    FindNonZeroMatrixColumnsAlgo algo;
    algo.set_progress_reporter(this);
    return(algo.run(input,output));
  }
  
  return (false);
}


} // end namespace SCIRunAlgo
