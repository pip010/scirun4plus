/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2010 Scientific Computing and Imaging Institute,
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
#include <stdexcept>

#include <Core/Datatypes/BlockMatrix.h>

#include <Testing/Util/MatrixTestUtilities.h>

using namespace SCIRun;
using namespace SCIRun::TestUtils;
using namespace boost::assign;

TEST(BlockMatrixTests, CanBlockMultiplySquare)
{
  BlockMatrix m1(2,2), m2(2,2);

  m1(0,0) = DenseMatrix(2,2,1.0);
  m1(0,1) = DenseMatrix(2,2,2.0);
  m1(1,0) = DenseMatrix(2,2,3.0);
  m1(1,1) = DenseMatrix(2,2,4.0);

  DenseMatrixHandle d1(m1.to_dense());

  EXPECT_MATRIX_EQ_TO(*d1, 
    (1,1,2,2)
    (1,1,2,2)
    (3,3,4,4)
    (3,3,4,4));

  m2(0,0) = DenseMatrix(2,2,-1.0);
  m2(0,1) = DenseMatrix(2,2,2);
  m2(1,0) = DenseMatrix(2,2,0);
  m2(1,1) = DenseMatrix(2,2,-3.0);

  d1 = m2.to_dense();

  EXPECT_MATRIX_EQ_TO(*d1, 
    (-1,-1,2,2)
    (-1,-1,2,2)
    (0,0,-3,-3)
    (0,0,-3,-3));

  BlockMatrix product = m1 * m2;

  DenseMatrixHandle productAsDense = product.to_dense();

  EXPECT_MATRIX_EQ_TO(*productAsDense, 
    (-2,-2,-8,-8)
    (-2,-2,-8,-8)
    (-6,-6,-12,-12)
    (-6,-6,-12,-12));

  DenseMatrixHandle product01 = product(0,1);

  EXPECT_MATRIX_EQ_TO(*product01,
    (-8,-8)
    (-8,-8));
}

TEST(BlockMatrixTests, CanBlockAdd)
{
  BlockMatrix m1(2,3), m2(2,3);

  m1(0,0) = DenseMatrix(2,2,1.0);
  m1(0,1) = DenseMatrix(2,2,2.0);
  m1(0,2) = DenseMatrix(2,1,3.0);
  m1(1,0) = DenseMatrix(2,2,4.0);
  m1(1,1) = DenseMatrix(2,2,5.0);
  m1(1,2) = DenseMatrix(2,1,6.0);

  DenseMatrixHandle d1(m1.to_dense());

  EXPECT_MATRIX_EQ_TO(*d1, 
    (1,1,2,2,3)
    (1,1,2,2,3)
    (4,4,5,5,6)
    (4,4,5,5,6));

  m2(0,0) = DenseMatrix(2,2,-1);
  m2(0,1) = DenseMatrix(2,2,-2);
  m2(0,2) = DenseMatrix(2,1,-3);
  m2(1,0) = DenseMatrix(2,2,-4);
  m2(1,1) = DenseMatrix(2,2,-5);
  m2(1,2) = DenseMatrix(2,1,-6);

  d1 = m2.to_dense();

  EXPECT_MATRIX_EQ_TO(*d1, 
    (-1,-1,-2,-2,-3)
    (-1,-1,-2,-2,-3)
    (-4,-4,-5,-5,-6)
    (-4,-4,-5,-5,-6));

  BlockMatrix sum = m1 + m2;

  DenseMatrixHandle sumAsDense = sum.to_dense();

  EXPECT_MATRIX_EQ_TO(*sumAsDense, 
    (0,0,0,0,0)
    (0,0,0,0,0)
    (0,0,0,0,0)
    (0,0,0,0,0));
}

TEST(BlockMatrixTests, CanBlockMultiplyNonsquarePartitionsOfNonsquareBlocks)
{
  BlockMatrix m1(2, 3);

  m1(0,0) = MAKE_DENSE_MATRIX(
    (1,2,3)
    (6,5,4));
  m1(0,1) = MAKE_DENSE_MATRIX(
    (4,3,2,1)
    (8,7,6,5));
  m1(0,2) = MAKE_DENSE_MATRIX(
    (1,2)
    (3,4));
  m1(1,0) = MAKE_DENSE_MATRIX(
    (-1,-2,-3)
    (-4,-5,-6)
    (-7,-8,-9));
  m1(1,1) = MAKE_DENSE_MATRIX(
    (12,11,10,9)
    (8,7,6,5)
    (4,3,2,1));
  m1(1,2) = MAKE_DENSE_MATRIX(
    (-5,-6)
    (-3,-4)
    (-1,-2));

  DenseMatrixHandle dense(m1.to_dense());
  DenseMatrixHandle transpose(dense->make_transpose());
  EXPECT_MATRIX_EQ_TO(*transpose, 
    (1, 6, -1, -4, -7)
    (2, 5, -2, -5, -8)
    (3, 4, -3, -6, -9)
    (4, 8, 12,  8,  4)
    (3, 7, 11,  7,  3)
    (2, 6, 10,  6,  2)
    (1, 5,  9,  5,  1)
    (1, 3, -5, -3, -1)
    (2, 4, -6, -4, -2));

  BlockMatrix m2(3, 2);
  m2(0,0) = DenseMatrix(3,1,1);
  m2(0,1) = DenseMatrix(3,2,-1);
  m2(1,0) = DenseMatrix(4,1,-1);
  m2(1,1) = DenseMatrix(4,2,1);
  m2(2,0) = DenseMatrix(2,1,1);
  m2(2,1) = DenseMatrix(2,2,-1);

  dense = m2.to_dense();
  EXPECT_MATRIX_EQ_TO(*dense, 
    (1, -1, -1)
    (1, -1, -1)
    (1, -1, -1)
    (-1, 1,  1)
    (-1, 1,  1)
    (-1, 1,  1)
    (-1, 1,  1)
    (1, -1, -1)
    (1, -1, -1));

  BlockMatrix product = m1 * m2;

  DenseMatrixHandle productAsDense = product.to_dense();

  EXPECT_MATRIX_EQ_TO(*productAsDense, 
    (-1, 1, 1)
    (-4, 4, 4)
    (-59, 59, 59)
    (-48, 48, 48)
    (-37, 37, 37));
}

TEST(BlockMatrixTests, ThrowsWhenBlocksHaveMismatchedSize)
{
  EXPECT_THROW({
  BlockMatrix m(2, 3);

  m(0,0) = MAKE_DENSE_MATRIX(
                           (1,2,3)
                           (6,5,4));
  m(0,1) = MAKE_DENSE_MATRIX(
                           (4,3,2,1)
                           (8,7,6,5)
                           (1,1,1,1));  // # rows don't match in this row partition
  }, std::invalid_argument);

  EXPECT_THROW({
    BlockMatrix m(2, 3);

    m(0,0) = MAKE_DENSE_MATRIX(
      (1,2,3)
      (6,5,4));
    m(1,0) = MAKE_DENSE_MATRIX(
      (4,3,2,1)
      (8,7,6,5)
      );  // # cols don't match in this col partition
  }, std::invalid_argument);
}

TEST(BlockMatrixTests, CanDeepCopyBlocksBetweenMatrices)
{
  BlockMatrix m1(2,3), m2(2,3);

  m1(0,0) = DenseMatrix(2,2,1.0);
  m1(0,1) = DenseMatrix(2,2,2.0);
  m1(0,2) = DenseMatrix(2,1,3.0);
  m1(1,0) = DenseMatrix(2,2,4.0);
  m1(1,1) = DenseMatrix(2,2,5.0);
  m1(1,2) = DenseMatrix(2,1,6.0);

  DenseMatrixHandle d1(m1.to_dense());

  EXPECT_MATRIX_EQ_TO(*d1, 
    (1,1,2,2,3)
    (1,1,2,2,3)
    (4,4,5,5,6)
    (4,4,5,5,6));

  m2(0,0) = m1(0,0);
  m2(0,1) = m1(0,1);
  m2(0,2) = m1(0,2);
  m2(1,0) = m1(1,0);
  m2(1,1) = m1(1,1);
  m2(1,2) = m1(1,2);

  DenseMatrixHandle d2 = m2.to_dense();

  EXPECT_MATRIX_EQ(*d1, *d2);

  m1(0,0) = DenseMatrix(2,2,-1.0);
  d1 = m1.to_dense();

  EXPECT_MATRIX_NOT_EQ(*d1, *d2);
}

TEST(BlockMatrixTests, CanCopyAndAssign)
{
  BlockMatrix m1(2,3);

  m1(0,0) = DenseMatrix(2,2,1.0);
  m1(0,1) = DenseMatrix(2,2,2.0);
  m1(0,2) = DenseMatrix(2,1,3.0);
  m1(1,0) = DenseMatrix(2,2,4.0);
  m1(1,1) = DenseMatrix(2,2,5.0);
  m1(1,2) = DenseMatrix(2,1,6.0);

  DenseMatrixHandle d1(m1.to_dense());

  EXPECT_MATRIX_EQ_TO(*d1, 
    (1,1,2,2,3)
    (1,1,2,2,3)
    (4,4,5,5,6)
    (4,4,5,5,6));

  BlockMatrix m2 = m1;
  
  DenseMatrixHandle d2 = m2.to_dense();

  EXPECT_MATRIX_EQ(*d1, *d2);

  m1(0,0) = DenseMatrix(2,2,-1.0);
  d1 = m1.to_dense();

  EXPECT_MATRIX_NOT_EQ(*d1, *d2);

  m1 = m2;

  d1 = m1.to_dense();
  d2 = m2.to_dense();

  EXPECT_MATRIX_EQ(*d1, *d2);

  m2(0,0) = DenseMatrix(2,2,-17.0);
  d1 = m1.to_dense();
  d2 = m2.to_dense();

  EXPECT_MATRIX_NOT_EQ(*d1, *d2);
}