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

#include <gtest/gtest.h>

#include <boost/timer.hpp>

#include <Core/Datatypes/MatrixFwd.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/MatrixOperations.h>

#include <Testing/Util/MatrixTestUtilities.h>

using namespace boost::assign; 
using namespace SCIRun;
using namespace SCIRun::TestUtils;

#define TIME_MATRIX_ACCESS_USING(access, code) void TimeMatrixAccess##access(int size, int times)\
{\
  DenseMatrixHandle mH = make_random(size,size);\
  DenseMatrix& m = *mH;\
  boost::timer t;\
  {\
    for (int n = 0; n < times; ++n)\
      for (index_type i = 0; i < m.nrows(); ++i)\
        for (index_type j = 0; j < m.ncols(); ++j)\
          code *= 2;\
  }\
  double time = t.elapsed();\
  std::cout << "Took " << time << " seconds for " << size << " square matrix, " << times << " trials, using " << #access << std::endl;\
}\

TIME_MATRIX_ACCESS_USING(Pointer, m[i][j])
TIME_MATRIX_ACCESS_USING(Operator, m(i,j))

TEST(MatrixElementAccessorPerformanceTest, DISABLED_CompareAccessFunctions)
{
  TimeMatrixAccessPointer(5000, 100);
  TimeMatrixAccessOperator(5000, 100);
}