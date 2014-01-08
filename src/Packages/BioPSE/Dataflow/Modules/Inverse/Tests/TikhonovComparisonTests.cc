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

#include <stdexcept>
#include <boost/filesystem.hpp>
#include <boost/timer.hpp>
#include <gtest/gtest.h>

#include <Packages/BioPSE/Dataflow/Modules/Inverse/SolveInverseProblemWithTikhonov.h>
#include <Packages/BioPSE/Dataflow/Modules/Inverse/SolveInverseProblemWithTikhonovSVD.h>
#include <Packages/BioPSE/Dataflow/Modules/Inverse/Tests/TikhonovTestDataSource.h>

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixOperations.h>

#include <Testing/Util/MatrixTestUtilities.h>
#include <Testing/Util/MatlabMatrixReader.h>

using namespace SCIRun;
using namespace BioPSE;
using namespace boost::assign;
using namespace SCIRun::TestUtils;

struct ScopedTimer
{
  explicit ScopedTimer(const std::string& name) : name_(name) 
  {
    std::cout << "Starting timer " << name_ << std::endl;
    t_.restart();
  }

  ~ScopedTimer()
  {
    double elapsed = t_.elapsed();
    std::cout << "Timer " << name_ << " stopped at " << elapsed << " seconds." << std::endl;
  }

  std::string name_;
  boost::timer t_;
};

void CompareTikhonovMethods(MatrixHandle measuredData, MatrixHandle forwardMatrix, double lambda, double error = 1e-15, double lambdaDelta = 0, int lambdaLoops = 1)
{
  SvdForTesting svd(*forwardMatrix->dense());

  TikhonovAlgorithmImpl algoTikhonov(forwardMatrix, measuredData);
  TikhonovSVDAlgorithm algoTikhonovSVD(*measuredData->column(), *svd.U_, *svd.S_->dense(), *svd.V_->dense(), 0, TikhonovSVDAlgorithm::SVD);

  double testLambda = lambda;
  for (int i = 0; i < lambdaLoops; ++i)
  {
    std::cout << "Comparing algorithms with lambda = " << testLambda << std::endl;
    
    algoTikhonov.run(TikhonovAlgorithmImpl::Input("single", testLambda, 0, 0, 0, 0, 0));

    MatrixHandle expectedOutput = open_matlab_file(net1output.string())->read_matrix("Matrix1");
    ColumnMatrixHandle inverseSolution = algoTikhonov.get_inverse_solution()->column();

    ASSERT_TRUE(inverseSolution.get_rep());

    algoTikhonovSVD.tikhonov_fun(testLambda);

    ColumnMatrixHandle inverseSolutionSVD = algoTikhonovSVD.get_solution();

    EXPECT_COLUMN_MATRIX_EQ_BY_TWO_NORM(*inverseSolution, *inverseSolutionSVD, error);
    testLambda += lambdaDelta;
  }
}

TEST(TikhonovComparisonTest, CompareTestNet1Underdetermined)
{
  MatrixHandle measuredData = underInput->read_matrix("data")->column();
  //needs a larger lambda due to high noise
  CompareTikhonovMethods(measuredData, forwardMatrixUnder, 0.16, 2e-8);
}

TEST(TikhonovComparisonTest, CompareTestSmallSquareData)
{
  TestTikhonovSVDData data(squareData);
  CompareTikhonovMethods(data.matrixMeasDatD, data.A, 0.0001, 1e-12);
}

TEST(TikhonovComparisonTest, CompareTestSmallUnderdeterminedData)
{
  TestTikhonovSVDData data(underDeterminedData);
  CompareTikhonovMethods(data.matrixMeasDatD, data.A, 0.0001, 1e-13);
}

TEST(TikhonovComparisonTest, CompareTestNet3OverdeterminedCase)
{
  EasyMatlabMatrixReader input(overDeterminedForwardMatrixFile.string());
  MatrixHandle forwardMatrixOver = input.read_matrix("lead");
  MatrixHandle measuredData = open_matlab_file(measuredDataFile.string())->read_matrix("data")->column();

  CompareTikhonovMethods(measuredData, forwardMatrixOver, 0.0001, 2e-7);
}

TEST(TikhonovComparisonTest, DISABLED_LotsOfUnderdeterminedRandomProblems)
{
  std::ofstream file("C:\\Dev\\Dropbox\\TikhonovSVD\\compare\\undertest100.txt");
  for (int rows = 100; rows < 150; rows += 10)
  {
    for (int cols = rows + 1; cols < rows*100; cols += 200)
    {
      std::ostringstream ostr;
      ostr << "SOLVING UNDERDETERMINED PROBLEM, rows = " << rows << ", cols = " << cols << std::endl;
      std::cout << ostr.str();
      file << ostr.str();

      TikhonovTestHelper::ForwardProblem p = TikhonovTestHelper::MakeForwardProblem(rows, cols);
      DenseMatrixHandle& A = p.first;
      ColumnMatrixHandle& rhs = p.second;
      CompareTikhonovMethods(rhs, A, 0.0001);
    }
  }
}