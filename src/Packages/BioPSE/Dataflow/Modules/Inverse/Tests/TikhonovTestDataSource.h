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

#ifndef TIKHONOV_TEST_DATA_SOURCE_H__
#define TIKHONOV_TEST_DATA_SOURCE_H__

#include <boost/filesystem.hpp>

#include <Core/Datatypes/Matrix.h>
#include <Testing/Util/MatlabMatrixReader.h>
#include <Packages/BioPSE/Dataflow/Modules/Inverse/SolveInverseProblemWithTikhonovSVD.h>

namespace SCIRun
{
  namespace TestUtils
  {
    namespace 
    {
      // Visual Studio's debugger use a different build output than bin\Debug, so one needs to hard code for debugging purposes.
      // TODO: switch testing path root to variable (should be enabled in SCIRun 5)
#if (DEBUG && _WIN32)
      boost::filesystem::path testDataRoot = boost::filesystem::path("C:\\Dev\\trunk_ref\\SCIRunUnitTestData");
#elif (_WIN32) // Windows Release builds?
      boost::filesystem::path testDataRoot = boost::filesystem::path("..") / ".." / ".." / "SCIRunUnitTestData";
#else // *nix platforms
      boost::filesystem::path testDataRoot = boost::filesystem::path("..") / ".." / "SCIRunUnitTestData";
#endif
      boost::filesystem::path tikhonovRoot = testDataRoot / "Tikhonov";
      boost::filesystem::path tikhonovSVDRoot = testDataRoot / "TikhonovSVD";

      boost::filesystem::path underDeterminedForwardMatrixFile = tikhonovRoot / "lead_eeg_fem.mat";
      boost::filesystem::path overDeterminedForwardMatrixFile = tikhonovRoot / "lead_eeg_fem_small.mat";
      boost::filesystem::path measuredDataFile = tikhonovRoot / "data_set1.mat";
      boost::filesystem::path Wpath = tikhonovRoot / "W.mat";
      boost::filesystem::path Cpath = tikhonovRoot / "C.mat";
      boost::filesystem::path net1output = tikhonovRoot / "net1_output.mat";
      
      boost::filesystem::path linsolvetest_input = tikhonovRoot / "linsolvetest_input.mat";
      boost::filesystem::path linsolvetest_expected_output = tikhonovRoot / "linsolvetest_expected_output.mat";
      
      boost::filesystem::path linsolvetest_input1 = tikhonovRoot / "linsolvetest_input1.mat";
      boost::filesystem::path linsolvetest_expected_output1 = tikhonovRoot / "linsolvetest_expected_output1.mat";
      
      boost::filesystem::path net2output = tikhonovRoot / "net2_output.mat";
      boost::filesystem::path net3output = tikhonovRoot / "net3_output.mat";

      MatlabMatrixReaderHandle underInput = open_matlab_file(measuredDataFile.string());
      MatrixHandle forwardMatrixUnder = underInput->read_matrix("lead");

      boost::filesystem::path underDeterminedData = tikhonovSVDRoot / "TknvTest_underDetermine1.mat";
      boost::filesystem::path squareData = tikhonovSVDRoot / "TknvTest1.mat";
    }

    //TODO MAJOR: extract ComputeSVD algorithm header rather than copy-pasting here.
    struct SvdForTesting
    {
      SvdForTesting(const DenseMatrix& input)
      {
        MatrixHandle imH(input.clone());

        // The input matrix is an m-by-n matrix, find and store those sizes
        int m = imH->nrows();
        int n = imH->ncols();
        int numsingularvals = std::min(m,n);

        // We don't care what type of matrix the input is, just obtain a dense matrix
        // (because that's the only type that has a built-in svd member function at the moment)
        DenseMatrix *denseinput = imH->dense();

        // Create the matrices for the singular value decomposition
        S_ = new ColumnMatrix(numsingularvals);
        ColumnMatrix *S = S_->column();

        U_ = new DenseMatrix(m,m);
        DenseMatrix *U = U_->dense();

        V_ = new DenseMatrix(n,n);
        DenseMatrix *V = V_->dense();

        // Compute the SVD of the input matrix

        // Since this code used to throw an exception, let the upstream caller handle
        // any exceptions
        LinearAlgebra::svd(*denseinput, *U, *S, *V);
        V_ = V_->make_transpose();
      }

      ColumnMatrixHandle S_;
      DenseMatrixHandle U_,V_;
    };

    inline void printBeginning(const Matrix<double>& m, size_t elements = 10)
    {
      std::copy(m.get_data_pointer(), m.get_data_pointer() + elements, std::ostream_iterator<double>(std::cout, "\n"));
    }

    struct TestTikhonovSVDData
    {
      TestTikhonovSVDData(const boost::filesystem::path& file)
      {
        EasyMatlabMatrixReader reader(file.string());
        matrixU = reader.read_matrix("U");
        matrixS = reader.read_matrix("sm");
        matrixV = reader.read_matrix("V");
        matrixX = reader.read_matrix("X");
        matrixMeasDatD = reader.read_matrix("rhs");
        expectedSolutionGSVD = reader.read_matrix("x2");
        expectedSolutionSVD = reader.read_matrix("x3");
        aU = reader.read_matrix("aU");
        aS = reader.read_matrix("aS");
        aV = reader.read_matrix("aV");
        A = reader.read_matrix("A");
        columns = matrixU->nrows();
        solutionRows = matrixV->nrows();
      }

      MatrixHandle matrixMeasDatD, matrixU, matrixS, matrixV, matrixX, expectedSolutionGSVD, expectedSolutionSVD, aU, aS, aV, A;
      int columns, solutionRows;
    };

     //TODO move to Matrix test helper header in unit test util library
    struct ScopedMatrixPrinting
    {
      static bool print_matrix;
      ScopedMatrixPrinting()
      {
        print_matrix = true;
      }
      ~ScopedMatrixPrinting()
      {
        print_matrix = false;
      }
    };

    inline void PrintMatrix(const char* name, const Matrix<double>& m, bool print = ScopedMatrixPrinting::print_matrix)
    {
      if (print)
      {
        std::cout << name << "\n" << 
          m.nrows() << "x" << m.ncols() << "\n" <<
          matrix_to_string(m) 
          << std::endl;
      }
    }

    inline void read_and_print(EasyMatlabMatrixReader& r, const char* name)
    {
      MatrixHandle m = r.read_matrix(name);
      PrintMatrix(name, *m);
    }

    

    class TikhonovTestHelper
    {
    public:
      static void DafangTestCode_SVDCase(int rows, int cols);
      typedef std::pair<DenseMatrixHandle, ColumnMatrixHandle> ForwardProblem;
      static ForwardProblem MakeForwardProblem(int rows, int cols);
    private:
      class TikhonovWithLambdaFunc
      {
      public:
        TikhonovWithLambdaFunc(const SvdForTesting& svd, const ColumnMatrix& rhs) :
            Sdense_(svd.S_->dense()),
              algo_(rhs, *svd.U_, *Sdense_, *svd.V_, 0, BioPSE::TikhonovSVDAlgorithm::SVD)
            {
            }

            void operator()(double lambda) 
            {
              algo_.tikhonov_fun(lambda);
              PrintMatrix("solution", *algo_.get_solution());
            }
      private:
        DenseMatrixHandle Sdense_;
        BioPSE::TikhonovSVDAlgorithm algo_;
      };
    };
  }
}



#endif
