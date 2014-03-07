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
#include <gtest/gtest.h>

#include <Packages/BioPSE/Dataflow/Modules/Inverse/SolveInverseProblemWithTikhonov.h>
#include <Packages/BioPSE/Dataflow/Modules/Inverse/Tests/TikhonovTestDataSource.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/SparseRowMatrixFromMap.h>

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <Testing/Util/MatrixTestUtilities.h>
#include <Testing/Util/MatlabMatrixReader.h>

using namespace SCIRun;
using namespace BioPSE;
using namespace boost::assign;
using namespace SCIRun::TestUtils;

// TODO: need tests for new combinations of automatic, underdetermined and
// overdetermined subcases

// TODO: SolveLinearSystemWithLapackTest_simpletest, SolveLinearSystemWithLapackTest_regularized
// and SolveLinearSystemWithLapackTest_unregularized_will_fail test LinearAlgebra::solve_lapack,
// not the module algorithm. Move to Core/Datasets tests???

TEST(SolveLinearSystemWithLapackTest_simpletest, TestSolveLinSysTestLapack_simpletest)
{

  //simple test:  Id * x = b , x = b
  MatrixHandle measuredData = underInput->read_matrix("data")->column();
  EXPECT_EQ(118, forwardMatrixUnder->nrows());
  EXPECT_EQ(11172, forwardMatrixUnder->ncols());
  EXPECT_EQ(118, measuredData->nrows());
  EXPECT_EQ(1, measuredData->ncols());
  
  MatrixHandle Identity = DenseMatrix::identity(118);

  ColumnMatrixHandle solution = new ColumnMatrix(118);
  
  EXPECT_NO_THROW(LinearAlgebra::solve_lapack(*(Identity->dense()), *(matrix_cast::as_column(measuredData)),  *(matrix_cast::as_column(solution))));
  
  EXPECT_MATRIX_EQ(*solution, *measuredData);
  
}


TEST(SolveLinearSystemWithLapackTest_regularized, TestSolveLinSysTestLapack_regularized)
{
  // regularized example
  MatrixHandle measuredData = underInput->read_matrix("data")->column();
  EXPECT_EQ(118, forwardMatrixUnder->nrows());
  EXPECT_EQ(11172, forwardMatrixUnder->ncols());
  EXPECT_EQ(118, measuredData->nrows());
  EXPECT_EQ(1, measuredData->ncols());
  
  
  MatrixHandle input = open_matlab_file(linsolvetest_input.string())->read_matrix("input");
  
  ColumnMatrixHandle solution = new ColumnMatrix(118);
  
  EXPECT_NO_THROW(LinearAlgebra::solve_lapack(*(input->dense()), *(matrix_cast::as_column(measuredData)),  *(matrix_cast::as_column(solution))));
  
  MatrixHandle expectedOutput = open_matlab_file(linsolvetest_expected_output.string())->read_matrix("res");
  
  EXPECT_MATRIX_EQ(*solution, *expectedOutput);
}


TEST(SolveLinearSystemWithLapackTest_unregularized, TestSolveLinSysTestLapack_unregularized)
{
  // example not regularized
  //
  // in SCIRun 4.x and LAPACK, the result from solve_lapack does not match the result from Matlab
  MatrixHandle measuredData = underInput->read_matrix("data")->column();
  EXPECT_EQ(118, forwardMatrixUnder->nrows());
  EXPECT_EQ(11172, forwardMatrixUnder->ncols());
  EXPECT_EQ(118, measuredData->nrows());
  EXPECT_EQ(1, measuredData->ncols());
  
  
  MatrixHandle input = open_matlab_file(linsolvetest_input1.string())->read_matrix("input");
  
  ColumnMatrixHandle solution = new ColumnMatrix(118);
  
  EXPECT_NO_THROW(LinearAlgebra::solve_lapack(*(input->dense()), *(matrix_cast::as_column(measuredData)),  *(matrix_cast::as_column(solution))));
  
  MatrixHandle expectedOutput = open_matlab_file(linsolvetest_expected_output1.string())->read_matrix("res");
  
  EXPECT_MATRIX_NOT_EQ(*solution, *expectedOutput);
}


TEST(TikhonovAlgorithmTest, TestNet1Underdetermined)
{
  MatrixHandle measuredData = underInput->read_matrix("data")->column();
  EXPECT_EQ(118, forwardMatrixUnder->nrows());
  EXPECT_EQ(11172, forwardMatrixUnder->ncols());
  EXPECT_EQ(118, measuredData->nrows());
  EXPECT_EQ(1, measuredData->ncols());

  TikhonovAlgorithmImpl algo(forwardMatrixUnder, measuredData);

  EXPECT_NO_THROW(algo.run(TikhonovAlgorithmImpl::Input("lcurve", 0, 0, 10, 50, 1000000, 0)));

  MatrixHandle expectedOutput = open_matlab_file(net1output.string())->read_matrix("Matrix1");
  EXPECT_EQ(11172, expectedOutput->nrows());
  EXPECT_EQ(1, expectedOutput->ncols());

  MatrixHandle inverseSolution = algo.get_inverse_solution();

  ASSERT_TRUE(inverseSolution.get_rep());

  EXPECT_MATRIX_EQ(*inverseSolution, *expectedOutput);
}


TEST(TikhonovAlgorithmTest, TestNet2UnderdeterminedWithWeighting)
{
  MatrixHandle measuredData = forwardMatrixUnder->submatrix(0, 0, forwardMatrixUnder->nrows() - 1, 0)->column();
  EXPECT_EQ(118, forwardMatrixUnder->nrows());
  EXPECT_EQ(11172, forwardMatrixUnder->ncols());
  EXPECT_EQ(118, measuredData->nrows());
  EXPECT_EQ(1, measuredData->ncols());

  MatrixHandle W = open_matlab_file(Wpath.string())->read_matrix("W");
  EXPECT_TRUE(matrix_is::sparse(W));
  EXPECT_EQ(11172, W->nrows());
  EXPECT_EQ(11172, W->ncols());

  MatrixHandle C = open_matlab_file(Cpath.string())->read_matrix("C");
  EXPECT_EQ(118, C->nrows());
  EXPECT_EQ(118, C->ncols());

  TikhonovAlgorithmImpl algo(forwardMatrixUnder, measuredData,
                             TikhonovAlgorithmImpl::automatic, TikhonovAlgorithmImpl::solution_constrained, TikhonovAlgorithmImpl::residual_constrained,
                             W, C);

  algo.run(TikhonovAlgorithmImpl::Input("lcurve", 0, 0, 5, 50, 1000000, 0));

  MatrixHandle expectedOutput = open_matlab_file(net2output.string())->read_matrix("Matrix1");
  EXPECT_EQ(11172, expectedOutput->nrows());
  EXPECT_EQ(1, expectedOutput->ncols());

  MatrixHandle inverseSolution = algo.get_inverse_solution();
  ASSERT_TRUE(inverseSolution.get_rep());

  EXPECT_MATRIX_EQ(*inverseSolution, *expectedOutput);
}

TEST(TikhonovAlgorithmTest, TestNet3OverdeterminedCase)
{
  EasyMatlabMatrixReader input(overDeterminedForwardMatrixFile.string());
  MatrixHandle forwardMatrixOver = input.read_matrix("lead");
  MatrixHandle measuredData = open_matlab_file(measuredDataFile.string())->read_matrix("data")->column();
  EXPECT_EQ(118, forwardMatrixOver->nrows());
  EXPECT_EQ(63, forwardMatrixOver->ncols());
  EXPECT_EQ(118, measuredData->nrows());
  EXPECT_EQ(1, measuredData->ncols());

  TikhonovAlgorithmImpl algo(forwardMatrixOver, measuredData);

  algo.run(TikhonovAlgorithmImpl::Input("lcurve", 0, 0, 40, 50, 1000000, 0));

  MatrixHandle expectedOutput = open_matlab_file(net3output.string())->read_matrix("Matrix1");
  EXPECT_EQ(63, expectedOutput->nrows());
  EXPECT_EQ(1, expectedOutput->ncols());

  MatrixHandle inverseSolution = algo.get_inverse_solution();
  ASSERT_TRUE(inverseSolution.get_rep());

  EXPECT_MATRIX_EQ(*inverseSolution, *expectedOutput);
}
