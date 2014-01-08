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

#include <stdexcept>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <gtest/gtest.h>

#include <Packages/BioPSE/Dataflow/Modules/Inverse/SolveInverseProblemWithTikhonovSVD.h>
#include <Packages/BioPSE/Dataflow/Modules/Inverse/Tests/TikhonovTestDataSource.h>

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Exceptions/InvalidState.h>
#include <Core/Datatypes/MatrixOperations.h>
#include <Core/Datatypes/ColumnMatrixFunctions.h>

#include <Testing/Util/MatrixTestUtilities.h>
#include <Testing/Util//MatlabMatrixReader.h>

using namespace SCIRun;
using namespace BioPSE;
using namespace boost::assign;

using namespace SCIRun::TestUtils;

/*
* NOTE
* Most of these tests are broken until matrix dimension requirements are solidified.
*/

bool ScopedMatrixPrinting::print_matrix = false;

TEST(TikhonovSVDTests, DISABLED_CtorThrowsForSVDCaseWhenSandVHaveUnequalRows)
{
  // M = 2
  // N = 3
	// number of singular values must == nrows of V.  This seems to be restricting matrix input to square?
	DenseMatrix U = MAKE_DENSE_MATRIX(
		(1,0,0)
		(0,1,0));
	ColumnMatrix S = MAKE_COLUMN_MATRIX(
		(1)
		(0));
	DenseMatrix V = MAKE_DENSE_MATRIX(
		(1,2,3)
		(0,0,1)
		(0,1,0));
	ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
		(1)
		(-1));

	DenseMatrixHandle Sdense = S.dense();
	ASSERT_THROW(TikhonovSVDAlgorithm(rhs, U, *Sdense, V, 0, TikhonovSVDAlgorithm::SVD), InvalidState);
}

TEST(TikhonovSVDTests, DISABLED_CtorThrowsForSVDCaseWhenVRowsDoesNotEqualUColumns)
{
  //  (m x n)(n x 1) = (m x 1)
  //  (m x m)(min(n,m) singular values, m x n expanded)(n x n) (n x 1) = (m x 1)
  // M = 3
  // N = 2
	DenseMatrix U = MAKE_DENSE_MATRIX(
		(1,0,1)
		(0,1,1)
		(0,0,1));
	ColumnMatrix S = MAKE_COLUMN_MATRIX(
		(1)
		(0));
	DenseMatrix V = MAKE_DENSE_MATRIX(
		(1,2)
		(2,1));
	ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
		(1)
		(-1)
    (0));

	DenseMatrixHandle Sdense = S.dense();
	ASSERT_THROW(TikhonovSVDAlgorithm(rhs, U, *Sdense, V, 0, TikhonovSVDAlgorithm::SVD), InvalidState);
}

TEST(TikhonovSVDTests, CtorThrowsForSVDCaseWhenRhsRowsDoesNotEqualURows)
{
  //this one makes sense: (_m_ x n)(n x 1) = (_m_ x 1)
  DenseMatrix U = MAKE_DENSE_MATRIX(
    (1,0)
    (0,1));
  ColumnMatrix S = MAKE_COLUMN_MATRIX(
    (1)
    (0));
  DenseMatrix V = MAKE_DENSE_MATRIX(
    (1,2)
    (2,1));
  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (1)
    (-1)
    (0));

  DenseMatrixHandle Sdense = S.dense();
  ASSERT_THROW(TikhonovSVDAlgorithm(rhs, U, *Sdense, V, 0, TikhonovSVDAlgorithm::SVD), InvalidState);
}

TEST(TikhonovSVDTests, CtorThrowsForGSVDCaseWhenXMatrixIsNull)
{
  DenseMatrix U = MAKE_DENSE_MATRIX(
    (1,0)
    (0,1));
  ColumnMatrix S = MAKE_COLUMN_MATRIX(
    (1)
    (0));
  DenseMatrix V = MAKE_DENSE_MATRIX(
    (1,2)
    (2,1));
  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (1)
    (-1)
    (0));

  DenseMatrixHandle Sdense = S.dense();
  ASSERT_THROW(TikhonovSVDAlgorithm(rhs, U, *Sdense, V, 0, TikhonovSVDAlgorithm::GSVD), InvalidState);
}

TEST(TikhonovSVDTests, DISABLED_CtorThrowsForGSVDCaseWhenURowsIsLessThanXRows)
{
  DenseMatrix U = MAKE_DENSE_MATRIX(
    (1,0,0)
    (0,1,0));
  ColumnMatrix S = MAKE_COLUMN_MATRIX(
    (1)
    (0));
  DenseMatrix V = MAKE_DENSE_MATRIX(
    (1,2)
    (2,1));
  DenseMatrix X = MAKE_DENSE_MATRIX(
    (1,2,0)
    (2,3,0)
    (3,4,0));
  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (1)
    (-1));

  DenseMatrixHandle Sdense = S.dense();
  ASSERT_THROW(TikhonovSVDAlgorithm(rhs, U, *Sdense, V, &X, TikhonovSVDAlgorithm::GSVD), InvalidState);
}

TEST(TikhonovSVDTests, DISABLED_CtorThrowsForGSVDCaseWhenXIsNotSquare)
{
  // M = 2
  // N = 3
  // P = 2
  DenseMatrix U = MAKE_DENSE_MATRIX(
    (1,0,1)
    (0,1,1));
  ColumnMatrix S = MAKE_COLUMN_MATRIX(
    (1)
    (0));
  DenseMatrix V = MAKE_DENSE_MATRIX(
    (1,2)
    (2,1));
  DenseMatrix X = MAKE_DENSE_MATRIX(
    (1,2)
    (1,3)
    (1,4));
  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (1)
    (-1));

  DenseMatrixHandle Sdense = S.dense();
  ASSERT_THROW(TikhonovSVDAlgorithm(rhs, U, *Sdense, V, &X, TikhonovSVDAlgorithm::GSVD), InvalidState);
}

TEST(TikhonovSVDTests, CtorThrowsForGSVDCaseWhenSRowsDoesNotEqualVRows)
{
  DenseMatrix U = MAKE_DENSE_MATRIX(
    (1,0)
    (0,1));
  ColumnMatrix S = MAKE_COLUMN_MATRIX(
    (1)
    (0));
  DenseMatrix V = MAKE_DENSE_MATRIX(
    (1,2)
    (2,1)
    (0,1));
  DenseMatrix X = MAKE_DENSE_MATRIX(
    (1,2)
    (0,1));
  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (1)
    (-1));

  DenseMatrixHandle Sdense = S.dense();
  ASSERT_THROW(TikhonovSVDAlgorithm(rhs, U, *Sdense, V, &X, TikhonovSVDAlgorithm::GSVD), InvalidState);
}

TEST(TikhonovSVDTests, DISABLED_CtorThrowsForGSVDCaseWhenUColumnsDoesNotEqualXRows)
{
  DenseMatrix U = MAKE_DENSE_MATRIX(
    (1,0,0)
    (0,1,1));
  ColumnMatrix S = MAKE_COLUMN_MATRIX(
    (1)
    (0));
  DenseMatrix V = MAKE_DENSE_MATRIX(
    (1,2)
    (2,1));
  DenseMatrix X = MAKE_DENSE_MATRIX(
    (1,2)
    (1,3));
  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (1)
    (-1));

  DenseMatrixHandle Sdense = S.dense();
  ASSERT_THROW(TikhonovSVDAlgorithm(rhs, U, *Sdense, V, &X, TikhonovSVDAlgorithm::GSVD), InvalidState);
}

TEST(TikhonovSVDTests, CtorThrowsForGSVDCaseWhenVRowsIsGreaterThanXRows)
{
  DenseMatrix U = MAKE_DENSE_MATRIX(
    (1,0)
    (0,1));
  ColumnMatrix S = MAKE_COLUMN_MATRIX(
    (1)
    (0));
  DenseMatrix V = MAKE_DENSE_MATRIX(
    (1,2)
    (2,1)
    (0,1));
  DenseMatrix X = MAKE_DENSE_MATRIX(
    (1,2)
    (0,1));
  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (1)
    (-1));

  DenseMatrixHandle Sdense = S.dense();
  ASSERT_THROW(TikhonovSVDAlgorithm(rhs, U, *Sdense, V, &X, TikhonovSVDAlgorithm::GSVD), InvalidState);
}

TEST(TikhonovSVDTests, CtorThrowsForGSVDCaseWhenRhsRowsDoesNotEqualURows)
{
  DenseMatrix U = MAKE_DENSE_MATRIX(
    (1,0)
    (0,1));
  ColumnMatrix S = MAKE_COLUMN_MATRIX(
    (1)
    (0));
  DenseMatrix V = MAKE_DENSE_MATRIX(
    (1,2)
    (2,1));
  DenseMatrix X = MAKE_DENSE_MATRIX(
    (1,2)
    (0,1));
  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (1)
    (-1)
    (0));

  DenseMatrixHandle Sdense = S.dense();
  ASSERT_THROW(TikhonovSVDAlgorithm(rhs, U, *Sdense, V, &X, TikhonovSVDAlgorithm::GSVD), InvalidState);
}

//TODO: make these value-parameterized Google Tests
TEST(TikhonovSVDTests, CanRegularizeSquareInvertibleCaseUsingGSVD)
{
  TestTikhonovSVDData data(underDeterminedData);
  TikhonovSVDAlgorithm algo(*data.matrixMeasDatD->column(), *data.matrixU->dense(), *data.matrixS->dense(), *data.matrixV->dense(), data.matrixX->dense(), TikhonovSVDAlgorithm::GSVD);

  algo.tikhonov_fun(0.0001);

  ColumnMatrixHandle solution = algo.get_solution();

  EXPECT_MATRIX_EQ_TOLERANCE(*data.expectedSolutionGSVD, *solution, .1);
}

TEST(TikhonovSVDTests, CanRegularizeSquareInvertibleCaseUsingGSVDFirstTest)
{
  TestTikhonovSVDData data(squareData);
  TikhonovSVDAlgorithm algo(*data.matrixMeasDatD->column(), *data.matrixU->dense(), *data.matrixS->dense(), *data.matrixV->dense(), data.matrixX->dense(), TikhonovSVDAlgorithm::GSVD);

  algo.tikhonov_fun(0.0001);

  ColumnMatrixHandle solution = algo.get_solution();

  EXPECT_MATRIX_EQ_TOLERANCE(*data.expectedSolutionGSVD, *solution, .1);
}

TEST(TikhonovSVDTests, CanRegularizeUnderdeterminedCaseUsingSVD)
{
  TestTikhonovSVDData data(underDeterminedData);
  TikhonovSVDAlgorithm algo(*data.matrixMeasDatD->column(), *data.aU->dense(), *data.aS->dense(), *data.aV->dense(), 0, TikhonovSVDAlgorithm::SVD);

  algo.tikhonov_fun(0.0001);

  ColumnMatrixHandle solution = algo.get_solution();

  EXPECT_MATRIX_EQ_TOLERANCE(*data.expectedSolutionSVD, *solution, .1);
}

TEST(NeedMoreExamples, DISABLED_CopyPastedCodeForSVDs)
{
  SvdForTesting(MAKE_DENSE_MATRIX(
    (1,2,3)
    (1,4,5)
    (5,6,7)
    (4,0,0)
    ));
 
  SvdForTesting(MAKE_DENSE_MATRIX(
    (1,2,3,4)
    (1,4,5,0)
    (5,6,7,0)
    (4,0,0,0)
    ));

  SvdForTesting(MAKE_DENSE_MATRIX(
    (1,2,3,4,1)
    (1,4,5,0,0)
    (5,6,7,0,0)
    (4,0,0,0,0)
    ));
}

TEST(TikhonovSVDTests, CanRegularizeUsingSVDOverdetermined)
{
  SvdForTesting over(MAKE_DENSE_MATRIX(
    (1,2,3)
    (1,4,5)
    (5,6,7)
    (4,0,0)
    ));

  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (7)
    (13)
    (15)
    (-4)
    );

  TikhonovSVDAlgorithm algo(rhs, *over.U_->dense(), *over.S_->dense(), *over.V_->dense(), 0, TikhonovSVDAlgorithm::SVD);

  algo.tikhonov_fun(0.01);

  ColumnMatrixHandle solution = algo.get_solution();

  //std::cout << "solution:\n" << matrix_to_string(*solution) << std::endl;
  //std::cout << "inverse:\n" << matrix_to_string(*algo.get_inverse_matrix()) << std::endl;

  //EXPECT_MATRIX_EQ_TOLERANCE(*data.expectedSolutionSVD, *solution, .1);
}


TEST(TikhonovSVDTests, CanRegularizeUsingSVDUnderdetermined)
{
  SvdForTesting under(MAKE_DENSE_MATRIX(
    (1,2,3,4,1)
    (1,4,5,0,0)
    (5,6,7,0,0)
    (4,0,0,0,0)
    ));

  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (1)
    (-2)
    (-6)
    (-4)
    );

  TikhonovSVDAlgorithm algo(rhs, *under.U_->dense(), *under.S_->dense(), *under.V_->dense(), 0, TikhonovSVDAlgorithm::SVD);

  algo.tikhonov_fun(0.0001);

  ColumnMatrixHandle solution = algo.get_solution();

  //std::cout << "solution:\n" << matrix_to_string(*solution) << std::endl;
  //std::cout << "inverse:\n" << matrix_to_string(*algo.get_inverse_matrix()) << std::endl;

  //EXPECT_MATRIX_EQ_TOLERANCE(*data.expectedSolutionSVD, *solution, .1);
}


TEST(TikhonovSVDTests, CanRegularizeUsingSVDSquare)
{
  SvdForTesting under(MAKE_DENSE_MATRIX(
    (1,2,3,4)
    (1,4,5,0)
    (5,6,7,0)
    (4,0,0,0)
    ));

  ColumnMatrix rhs = MAKE_COLUMN_MATRIX(
    (1)
    (-2)
    (-6)
    (-4)
    );

  TikhonovSVDAlgorithm algo(rhs, *under.U_->dense(), *under.S_->dense(), *under.V_->dense(), 0, TikhonovSVDAlgorithm::SVD);

  algo.tikhonov_fun(0.0001);

  ColumnMatrixHandle solution = algo.get_solution();

  //std::cout << "solution:\n" << matrix_to_string(*solution) << std::endl;
  //std::cout << "inverse:\n" << matrix_to_string(*algo.get_inverse_matrix()) << std::endl;

  //EXPECT_MATRIX_EQ_TOLERANCE(*data.expectedSolutionSVD, *solution, .1);
}

void TikhonovTestHelper::DafangTestCode_SVDCase(int rows, int cols)
{
  ForwardProblem p = MakeForwardProblem(rows, cols);
  DenseMatrix& A = *p.first;
  ColumnMatrix& rhs = *p.second;

  SvdForTesting svd(A);

  TikhonovWithLambdaFunc t(svd, rhs);

  std::vector<int> lambdas;
  lambdas += 1,2,3,4,5,6,7,8,9,10;
  BOOST_FOREACH(int lambda, lambdas)
  {
    t(lambda);
  }
}

std::pair<DenseMatrixHandle, ColumnMatrixHandle> TikhonovTestHelper::MakeForwardProblem(int rows, int cols) 
{
  DenseMatrixHandle A = make_random(rows, cols, -5, 5);
  ColumnMatrix x = column_with_integer_sequence(1, cols);
  ColumnMatrix b = *A * x;

  PrintMatrix("A", *A);
  {
    //ScopedMatrixPrinting q;
    PrintMatrix("x", x);
  }

  PrintMatrix("b", b);

  ColumnMatrix e = 0.001 * b.vector_norm() * make_random_column(cols);
  PrintMatrix("e", e);

  ColumnMatrix rhs = b + e;
  PrintMatrix("rhs", rhs);

  ColumnMatrixHandle rhsH(rhs.clone());
  return std::make_pair(A, rhsH);
}



TEST(TikhonovSVDTests, DafangSquareCase)
{
  const int rows = 10;
  const int cols = 10;
  
  TikhonovTestHelper::DafangTestCode_SVDCase(rows, cols);
}

TEST(TikhonovSVDTests, DafangUnderdeterminedCase)
{
  const int rows = 5;
  const int cols = 10;

  //ASSERT_THROW(
    TikhonovTestHelper::DafangTestCode_SVDCase(rows, cols)
      ;
    //, InvalidState);
}

TEST(TikhonovSVDTests, DafangOverdeterminedCase)
{
  const int rows = 10;
  const int cols = 5;

  //ASSERT_THROW(
    TikhonovTestHelper::DafangTestCode_SVDCase(rows, cols)
      ;
    //, InvalidState);
}

TEST(TikhonovSVDTests, DISABLED_CanRegularizeUnderdeterminedCaseUsingSVDAnnotated)
{
  TestTikhonovSVDData data(underDeterminedData);
  ScopedMatrixPrinting s;
  PrintMatrix("A", *data.A);
  PrintMatrix("aU", *data.aU);
  PrintMatrix("aS", *data.aS);
  PrintMatrix("aV", *data.aV);

  TikhonovSVDAlgorithm algo(*data.matrixMeasDatD->column(), *data.aU->dense(), *data.aS->dense(), *data.aV->dense(), 0, TikhonovSVDAlgorithm::SVD);

  const double lambda = 0.0001;
  std::cout << "Lambda: " << lambda << std::endl;
  algo.tikhonov_fun(lambda);

  ColumnMatrixHandle solution = algo.get_solution();
  PrintMatrix("solution", *solution);
  PrintMatrix("expected solution", *data.expectedSolutionSVD);

  EXPECT_MATRIX_EQ_TOLERANCE(*data.expectedSolutionSVD, *solution, .1);
}