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

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/ColumnMatrix.h>

#include <Testing/Util/MatrixTestUtilities.h>

using namespace boost::assign; 
using namespace SCIRun;
using namespace SCIRun::TestUtils;

TEST(MatrixInitializerTest, Basic)
{
  std::vector<boost::tuple<double,double,double> >
    matrixData = 
    tuple_list_of
    (1,2,3)
    (4,5,6)
    (7,8,9)
    (10,11,12);

  EXPECT_EQ(4, matrixData.size());
  EXPECT_EQ(3, matrixData[0].get<2>());

  DenseMatrix matrix = convertDataToMatrix(matrixData);

  EXPECT_EQ(4, matrix.nrows());
  EXPECT_EQ(3, matrix.ncols());
  
  EXPECT_EQ( 1, matrix[0][0] );
  EXPECT_EQ( 5, matrix[1][1] );
  EXPECT_EQ( 12, matrix[3][2] );

  DenseMatrix almostAwesome = convertDataToMatrix(to_vector(tuple_list_of
  (1, 0, 0, 0, 0)
  (0, 0.1, 3, 0, 0)
  (0, 0, 0, 0, 0)
  (0, 4, 0, 0, 0)));

  EXPECT_EQ(4, almostAwesome.nrows());
  EXPECT_EQ(5, almostAwesome.ncols());
  EXPECT_EQ( 1, almostAwesome[0][0] );
  EXPECT_EQ( 0, almostAwesome[1][1] );
  EXPECT_EQ( 3, almostAwesome[1][2] );

  DenseMatrix exampleIntMatrix = MAKE_DENSE_MATRIX(
    (1, 0, 0)
    (0, 0, 3)
    (0, 0, 0)
    (0, 4, 0));

  EXPECT_EQ(4, exampleIntMatrix.nrows());
  EXPECT_EQ(3, exampleIntMatrix.ncols());
  EXPECT_EQ(1, exampleIntMatrix[0][0]);
  EXPECT_EQ(0, exampleIntMatrix[1][1] );
  EXPECT_EQ(3, exampleIntMatrix[1][2]);

  DenseMatrix exampleDoubleMatrix = MAKE_DENSE_MATRIX(
    (4.0, 0.0)
    (3, 0.0)
    (2.23607, 0.0)
    (0, 0.0));

  EXPECT_EQ(4, exampleDoubleMatrix.nrows());
  EXPECT_EQ(2, exampleDoubleMatrix.ncols());
  EXPECT_EQ(4, exampleDoubleMatrix[0][0]);
  EXPECT_EQ(3, exampleDoubleMatrix[1][0] );
  EXPECT_EQ(2.23607, exampleDoubleMatrix[2][0]);

  ColumnMatrix exampleColumn = MAKE_COLUMN_MATRIX(
    (4.0)
    (3)
    (2.23607)
    (0));

  EXPECT_EQ(4, exampleColumn.nrows());
  EXPECT_EQ(1, exampleColumn.ncols());
  EXPECT_EQ(4, exampleColumn[0]);
  EXPECT_EQ(3, exampleColumn[1]);
  EXPECT_EQ(2.23607, exampleColumn[2]);
}

TEST(MatrixInitializerTest, VectorApproach)
{
  std::vector<std::vector<double> > v;
  v += list_of(1.0).repeat(9,0.0),
    list_of(0.0)(1).repeat(8,0.0);

  EXPECT_EQ(2, v.size());
  EXPECT_EQ(10, v[0].size());
  EXPECT_EQ(10, v[1].size());

  DenseMatrix matrix = convertDataToMatrix(v);
  EXPECT_EQ(10, matrix.ncols());
  EXPECT_EQ(2, matrix.nrows());
}

TEST(MatrixInitializerTest, DenseMatrixConstructorWithValue)
{
  double test_val = 3.14;
  DenseMatrix dm1(2, 2, test_val);
  EXPECT_EQ(4, count_single_value(dm1, test_val));

  test_val = -0.0039;
  DenseMatrix dm2(100, 100, test_val);
  EXPECT_EQ(10000, count_single_value(dm2, test_val));
}
