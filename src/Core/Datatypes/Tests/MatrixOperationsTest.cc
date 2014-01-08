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

#include <sci_defs/lapack_defs.h>

#include <Core/Datatypes/MatrixFwd.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/MatrixOperations.h>
#include <Core/Datatypes/MatrixTypeConverter.h>
#include <Core/Exceptions/DimensionMismatch.h>

#include <Testing/Util/MatrixTestUtilities.h>

using namespace boost::assign; 
using namespace SCIRun;
using namespace SCIRun::TestUtils;

TEST(MatrixOperationsTest, ZeroMatrixDense)
{
  DenseMatrix dm1(60, 500);
  dm1.zero();
  EXPECT_EQ(30000, count_single_value(dm1, 0));

  DenseMatrix dm2(400, 1);
  dm2.zero();
  EXPECT_EQ(400, count_single_value(dm2, 0));
  
  DenseMatrix dm3(40, 10);
  dm3.zero();
  EXPECT_EQ(400, count_single_value(dm3, 0));
}

TEST(MatrixOperationsTest, ZeroMatrixSparseRow)
{
  SparseRowMatrix::Data data1(501, 20);
  SparseRowMatrix srm1(500, 20, data1, 20);
  srm1.zero();
  EXPECT_EQ(20, count_single_value(srm1, 0));
}

TEST(MatrixOperationsTest, AddSparse)
{
  SparseRowMatrixHandle srm1 = MAKE_SPARSE_ROW_MATRIX_POINTER((1,0,0)(0,0,0)(0,0,0));
  EXPECT_TRUE(matrix_is::sparse(srm1));

  SparseRowMatrixHandle srm2 = MAKE_SPARSE_ROW_MATRIX_POINTER((2,0,0)(0,0,0)(0,0,0));
  EXPECT_TRUE(matrix_is::sparse(srm2));

  MatrixHandle m = srm1 + srm2;
  EXPECT_TRUE(matrix_is::sparse(m));
  EXPECT_MATRIX_EQ_TO(*m, (3,0,0)(0,0,0)(0,0,0));
}


TEST(MatrixOperationsTest, AddSparseConverted)
{
  DenseMatrixHandle dm1 = MAKE_DENSE_MATRIX_POINTER((1,0,0)(0,0,0)(0,0,0)(0,0,0));
  SparseRowMatrixHandle srm1 = dm1->sparse();
  EXPECT_TRUE(matrix_is::sparse(srm1));
  
  DenseMatrixHandle dm2 = MAKE_DENSE_MATRIX_POINTER((2,0,0)(0,0,0)(0,0,0)(0,0,0));
  SparseRowMatrixHandle srm2 = dm2->sparse();
  EXPECT_TRUE(matrix_is::sparse(srm2));
  
  MatrixHandle m1 = dm1 + dm2;
  EXPECT_TRUE(matrix_is::dense(m1));
  EXPECT_MATRIX_EQ_TO(*m1, (3,0,0)(0,0,0)(0,0,0)(0,0,0));
  EXPECT_TRUE(matrix_is::dense(m1));
  
  MatrixHandle m2 = srm1 + srm2;
  EXPECT_TRUE(matrix_is::sparse(m2));
  EXPECT_SPARSE_ROW_MATRIX_EQ_TO(*m2, (3,0,0)(0,0,0)(0,0,0)(0,0,0));
}

TEST(MatrixOperationsTest, MultiplyDenseValues)
{
  DenseMatrix m1value = MAKE_DENSE_MATRIX((1,2)(3,4));
  DenseMatrix m2value = MAKE_DENSE_MATRIX((-1,0)(0,1));

  DenseMatrix prod1 = m1value * m2value;

  EXPECT_MATRIX_EQ_TO(prod1, (-1,2)(-3,4));
}

TEST(MatrixOperationsTest, DenseToSparse_SparseMultiply)
{
  DenseMatrix d = MAKE_DENSE_MATRIX((1,0,0)(0,0,2)(1,0,-1));
  
  SparseRowMatrixHandle s = d.sparse();
  DenseMatrixHandle dCopy = s->dense();
  
  EXPECT_MATRIX_EQ(*dCopy, d);
  
  EXPECT_TRUE(matrix_is::sparse(s));
  
  MatrixHandle d2 = dCopy * dCopy;
  MatrixHandle s2 = s * s;
  EXPECT_TRUE(matrix_is::sparse(s2));
  
  DenseMatrixHandle s2Dense = s2->dense();
  EXPECT_MATRIX_EQ(*s2Dense, *d2);
}

TEST(MatrixOperationsTest, MultiplyDensePointers)
{
  DenseMatrixHandle m1ptr = MAKE_DENSE_MATRIX_POINTER((1,2)(3,4));
  DenseMatrixHandle m2ptr = MAKE_DENSE_MATRIX_POINTER((-1,0)(0,1));

  MatrixHandle prod1 = m1ptr * m2ptr;

  EXPECT_MATRIX_EQ_TO(*prod1, (-1,2)(-3,4));
}

TEST(MatrixOperationsTest, MultiplyDenseByColumnValues)
{
  DenseMatrix m1value = MAKE_DENSE_MATRIX((1,2)(3,4));
  ColumnMatrix m2value = MAKE_COLUMN_MATRIX((-1)(1));

  MatrixHandle prod1T = (m1value * m2value).make_transpose();

  EXPECT_MATRIX_EQ_TO(*prod1T, (1,1));
}

TEST(MatrixOperationsTest, MultiplyDenseByColumnPointers)
{
  DenseMatrixHandle m1ptr = MAKE_DENSE_MATRIX_POINTER((1,2)(3,4));
  ColumnMatrixHandle m2ptr = MAKE_COLUMN_MATRIX_POINTER((-1)(1));

  MatrixHandle prod1T = (m1ptr * m2ptr)->make_transpose();

  EXPECT_MATRIX_EQ_TO(*prod1T, (1,1));
}

namespace
{
  DenseMatrixHandle Matrix1()
  {
    return DenseMatrix::identity(2);
  }
}

TEST(MatrixOperationsTest, MultiplyDensePointersToMatrix)
{
  DenseMatrixHandle m1ptr = MAKE_DENSE_MATRIX_POINTER((1,2)(3,4));
  DenseMatrixHandle m2ptr = MAKE_DENSE_MATRIX_POINTER((-1,0)(0,1));

  MatrixHandle prod1 = m1ptr * m2ptr;

  EXPECT_MATRIX_EQ_TO(*prod1, (-1,2)(-3,4));
  EXPECT_TRUE(matrix_is::dense(prod1));

  MatrixHandle m1ptrBase = m1ptr;
  MatrixHandle m2ptrBase = m2ptr;

  prod1 = m1ptrBase * m2ptrBase;

  EXPECT_MATRIX_EQ_TO(*prod1, (-1,2)(-3,4));
  EXPECT_TRUE(matrix_is::dense(prod1));

  DenseMatrixHandle m1ptrCopy = m1ptr;
  DenseMatrixHandle m2ptrCopy = Matrix1();
  MatrixHandle m3ptrCopy = Matrix1();  //should use the new template copy ctor

  prod1 = m1ptrCopy * m2ptrCopy * m3ptrCopy;

  EXPECT_MATRIX_EQ(*prod1, *m1ptrCopy);
  EXPECT_TRUE(matrix_is::dense(prod1));
}

TEST(MatrixOperationsTest, AddScalarAndMatrixFunction)
{
  const DenseMatrix m = MAKE_DENSE_MATRIX((1,2)(3,4));
  
  DenseMatrix result(2, 2, 0);
  Add(1, result, m);
  EXPECT_MATRIX_EQ_TO(result, (1,2)(3,4));

  result.zero();
  EXPECT_MATRIX_EQ_TO(result, (0,0)(0,0));
  Add(-1, result, m);
  EXPECT_MATRIX_EQ_TO(result, (-1,-2)(-3,-4));

  result.zero();
  Add(2, result, m);
  EXPECT_MATRIX_EQ_TO(result, (2,4)(6,8));

  result.zero();
  Add(-2, result, m);
  EXPECT_MATRIX_EQ_TO(result, (-2,-4)(-6,-8));

  result.zero();
  Add(0, result, m);
  EXPECT_MATRIX_EQ_TO(result, (0,0)(0,0));

  result.zero();
  Add(1, result, m);
  EXPECT_MATRIX_EQ_TO(result, (1,2)(3,4));
  Add(-1, result, m);
  EXPECT_MATRIX_EQ_TO(result, (-2,-4)(-6,-8));
}

TEST(MatrixOperationsTest, DenseScalarMultiply)
{
  DenseMatrix m = MAKE_DENSE_MATRIX((1.0,2.0)(3.0,4.0));
  m.scalar_multiply(2.0);
  EXPECT_MATRIX_EQ_TO(m, (2.0,4.0)(6.0,8.0));
  m.scalar_multiply(-0.5);
  EXPECT_MATRIX_EQ_TO(m, (-1.0,-2.0)(-3.0,-4.0));
}

TEST(MatrixOperationsTest, SparseScalarMultiply)
{
  SparseRowMatrix m = MAKE_SPARSE_ROW_MATRIX((1.0,0.0,0.0)(0.0,2.0,0.0)(3.0,0.0,4.0));

  m.scalar_multiply(2);
  EXPECT_SPARSE_ROW_MATRIX_EQ_TO(m, (2.0,0.0,0.0)(0.0,4.0,0.0)(6.0,0.0,8.0));

  m.scalar_multiply(-0.5);
  EXPECT_SPARSE_ROW_MATRIX_EQ_TO(m, (-1.0,0.0,0.0)(0.0,-2.0,0.0)(-3.0,0.0,-4.0));
  
  m = m * -1;
  EXPECT_SPARSE_ROW_MATRIX_EQ_TO(m, (1.0,0.0,0.0)(0.0,2.0,0.0)(3.0,0.0,4.0));  
}

TEST(MatrixOperationsTest, FindInfiniteComponent)
{
  ColumnMatrix x(2);
  x[0] = std::numeric_limits<double>::infinity();
  x[1] = 1;
  EXPECT_TRUE(std::find_if(x.begin(), x.end(), IsInfinite) != x.end());
  EXPECT_TRUE(std::find_if(x.begin(), x.end(), IsFinite) != x.end());
  x[1] = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(std::find_if(x.begin(), x.end(), IsFinite) != x.end());
  x.zero();
  EXPECT_FALSE(std::find_if(x.begin(), x.end(), IsInfinite) != x.end());
}

TEST(MatrixOperationsTest, IsZeroTest)
{
  DenseMatrix zero(2,3,0);
  EXPECT_TRUE(zero.is_zero());
  SparseRowMatrixHandle sparse = zero.sparse();
  EXPECT_TRUE(sparse->is_zero());  
  DenseColMajMatrixHandle colMaj = zero.dense_col_maj();
  EXPECT_TRUE(colMaj->is_zero());
  ColumnMatrixHandle column = zero.column();
  EXPECT_TRUE(column->is_zero());
  zero[1][2] = 1;
  EXPECT_FALSE(zero.is_zero());
  sparse = zero.sparse();
  EXPECT_FALSE(sparse->is_zero());  
  colMaj = zero.dense_col_maj();
  EXPECT_FALSE(colMaj->is_zero());
  column = zero.column();
  (*column)[0] = 1;
  EXPECT_FALSE(column->is_zero());
}

TEST(MatrixOperationsTest, DISABLED_MEMORYIsZeroPerformanceTestDenseMatrix)
{
  int size = sqrt(2e9/8);
  std::cout << "DenseMatrix num of elements: " << size*size << std::endl;
  boost::timer t;
  DenseMatrix zero(size,size,0);
  double time = t.elapsed();
  std::cout << "time for allocation: " << time << " s" << std::endl;
  t.restart();
  EXPECT_TRUE(zero.is_zero());
  time = t.elapsed();
  std::cout << "time for is_zero true: " << time << " s" << std::endl;
  zero[size-1][size-1] = 1.0;
  t.restart();
  EXPECT_FALSE(zero.is_zero());
  time = t.elapsed();
  std::cout << "time for is_zero false worse case: " << time << " s" << std::endl;
  zero[0][0] = 1.0;
  t.restart();
  EXPECT_FALSE(zero.is_zero());
  time = t.elapsed();
  std::cout << "time for is_zero false best case: " << time << " s" << std::endl;
}

TEST(MatrixOperationsTest, DISABLED_MEMORYIsZeroPerformanceTestColumnMatrix)
{
  int size = 2e9/8;
  std::cout << "ColumnMatrix num of elements: " << size << std::endl;
  boost::timer t;
  ColumnMatrix zero(size);
  zero.zero();
  double time = t.elapsed();
  std::cout << "time for allocation and zeroing: " << time << " s" << std::endl;
  t.restart();
  EXPECT_TRUE(zero.is_zero());
  time = t.elapsed();
  std::cout << "time for is_zero true: " << time << " s" << std::endl;
  zero[size-1] = 1.0;
  t.restart();
  EXPECT_FALSE(zero.is_zero());
  time = t.elapsed();
  std::cout << "time for is_zero false worse case: " << time << " s" << std::endl;
  zero[0] = 1.0;
  t.restart();
  EXPECT_FALSE(zero.is_zero());
  time = t.elapsed();
  std::cout << "time for is_zero false best case: " << time << " s" << std::endl;
}

#if defined(HAVE_LAPACK)

TEST(MatrixOperationsTest, TestSolveLapack)
{
  DenseMatrix M = MAKE_DENSE_MATRIX(
                                    (2, 4, 7, 9, 8)
                                    (6, 9, 2, 5, 2)
                                    (6, 3, 5, 1, 8)
                                    (1, 5, 6, 1, 2)
                                    (1, 2, 8, 2, 9));
  
  ColumnMatrix rhs = MAKE_COLUMN_MATRIX((107) (60) (71) (43) (82));
  ColumnMatrix output(5);

  EXPECT_NO_THROW(LinearAlgebra::solve_lapack(M, rhs, output));

  EXPECT_COLUMN_MATRIX_EQ_TO(output,
                             (1.0)
                             (2.0)
                             (3.0)
                             (4.0)
                             (5.0));
}

TEST(MatrixOperationsTest, TestSolveLapackDimensionMismatchException)
{
  DenseMatrix M = MAKE_DENSE_MATRIX(
                                    (2, 4, 7, 9)
                                    (6, 9, 2, 5)
                                    (6, 3, 5, 1)
                                    (1, 5, 6, 1)
                                    (1, 2, 8, 2));
  
  ColumnMatrix rhs = MAKE_COLUMN_MATRIX((107) (60) (71) (43) (82));
  ColumnMatrix output(5);

  EXPECT_THROW(LinearAlgebra::solve_lapack(M, rhs, output), DimensionMismatch);
}

#endif

