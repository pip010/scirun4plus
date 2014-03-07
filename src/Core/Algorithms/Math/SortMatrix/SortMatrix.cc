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

#include <Core/Algorithms/Math/SortMatrix/SortMatrix.h>
#include <Core/Algorithms/Math/SortMatrix/SortMatrixRows.h>
#include <Core/Algorithms/Math/SortMatrix/SortMatrixColumns.h>

// STL Fucntions we need
#include <algorithm>
#include <vector>

namespace SCIRunAlgo {

using namespace SCIRun;

bool 
SortMatrixAlgo::
run(MatrixHandle input, std::vector<index_type>& output)
{
  std::string method = get_option("method");
  if (method == "sortorder_rows")
  {
    SortMatrixRowsAlgo algo;
    algo.set_progress_reporter(this);
    algo.set_bool("ascending",get_bool("ascending"));
    return(algo.get_sortorder(input,output));
  }
  else if (method == "sortorder_columns")
  {
    SortMatrixColumnsAlgo algo;
    algo.set_progress_reporter(this);
    algo.set_bool("ascending",get_bool("ascending"));
    return(algo.get_sortorder(input,output));
  }

  return (false);
}

bool 
SortMatrixAlgo::
run(MatrixHandle input, MatrixHandle& output)
{
  std::string method = get_option("method");
  if (method == "sort_rows")
  {
    SortMatrixRowsAlgo algo;
    algo.set_progress_reporter(this);
    algo.set_bool("ascending",get_bool("ascending"));
    return(algo.run(input,output));
  }
  else if (method == "sort_columns")
  {
    SortMatrixColumnsAlgo algo;
    algo.set_progress_reporter(this);
    algo.set_bool("ascending",get_bool("ascending"));
    return(algo.run(input,output));
  }
  else if (method == "sortorder_rows")
  {
    SortMatrixRowsAlgo algo;
    algo.set_progress_reporter(this);
    algo.set_bool("ascending",get_bool("ascending"));
    return(algo.get_sortorder(input,output));
  }
  else if (method == "sortorder_columns")
  {
    SortMatrixColumnsAlgo algo;
    algo.set_progress_reporter(this);
    algo.set_bool("ascending",get_bool("ascending"));
    return(algo.get_sortorder(input,output));
  }

  return (false);
}

} // end namespace SCIRunAlgo

