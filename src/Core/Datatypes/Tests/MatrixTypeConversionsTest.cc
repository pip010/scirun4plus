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

#include <Core/Datatypes/MatrixFwd.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/MatrixOperations.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <Testing/Util/MatrixTestUtilities.h>

using namespace boost::assign; 
using namespace SCIRun;
using namespace SCIRun::TestUtils;

TEST(MatrixTypeConversionsTest, DenseToSparse)
{
  DenseMatrixHandle dm = MAKE_DENSE_MATRIX_POINTER((0.02,0.0,0.0)(0.0,100.2,0.0)(0.0,0.0,0.0)(0.0,0.0,0.75));
  EXPECT_TRUE(matrix_is::dense(dm));

  SparseRowMatrixHandle srm = dm->sparse();
  EXPECT_TRUE(matrix_is::sparse(srm));

  SparseRowMatrixHandle srm2 = MAKE_SPARSE_ROW_MATRIX_POINTER((0,0,0)(0,0,0.1));
  EXPECT_TRUE(matrix_is::sparse(srm2));
}

TEST(MatrixTypeConversionsTest, ColumnToDense)
{
  ColumnMatrix c = MAKE_COLUMN_MATRIX((1)(2)(3));

  DenseMatrixHandle d = DenseMatrixHandle(c.dense())->make_transpose();

  EXPECT_MATRIX_EQ_TO(*d, (1,2,3));

  ColumnMatrixHandle c2(c.clone());
  DenseMatrixHandle d2 = matrix_cast::as_dense(c2);
  EXPECT_EQ(0, d2.get_rep());
  EXPECT_EQ(d.get_rep(), matrix_cast::as_dense(d));
}

TEST(TypeMapTests, CanFindDerived)
{
  PersistentTypeIDPtr t = Persistent::find_derived("Matrix", "MatrixBase");
  EXPECT_TRUE(t.get() != 0);

  EXPECT_TRUE(Persistent::is_base_of(t->parent, t->type));
  EXPECT_TRUE(Persistent::is_base_of("MatrixBase", "ColumnMatrix"));
}
