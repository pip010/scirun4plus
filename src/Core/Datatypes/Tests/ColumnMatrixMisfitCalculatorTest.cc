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
#include <Testing/Util/MatrixTestUtilities.h>

#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Algorithms/Math/ColumnMisfitCalculator/ColumnMatrixMisfitCalculator.h>

using namespace boost::assign; 
using namespace SCIRun;
using namespace SCIRun::TestUtils;

void printStats(const ColumnMatrixMisfitCalculator& calc)
{
  std::cout << std::setprecision(15);
  std::cout << "p-value: " << calc.getPValue() << std::endl;
  std::cout << "CorrelationCoefficient: " << calc.getCorrelationCoefficient() << std::endl;
  std::cout << "InverseCorrelationCoefficient: " << calc.getInverseCorrelationCoefficient() << std::endl;
  std::cout << "RelativeRMS: " << calc.getRelativeRMS() << std::endl;
  std::cout << "RMS: " << calc.getRMS() << std::endl;
}

void checkStats(const ColumnMatrixMisfitCalculator& calc, double expectedCC, double expectedInvCC, double expectedRelRMS, double expectedRMS)
{
  const double error = 1e-12;
  EXPECT_NEAR(expectedCC, calc.getCorrelationCoefficient(), error);
  EXPECT_NEAR(expectedInvCC, calc.getInverseCorrelationCoefficient(), error);
  EXPECT_NEAR(expectedRelRMS, calc.getRelativeRMS(), 1e-8);
  EXPECT_NEAR(expectedRMS, calc.getRMS(), error);
}

TEST(ColumnMatrixMisfitCalculatorTest, TestOnArbitraryData)
{
  ColumnMatrix x = MAKE_COLUMN_MATRIX((1)(2)(3)(-3)(-2)(-1));
  ColumnMatrix y = MAKE_COLUMN_MATRIX((1.1)(2.1)(3.1)(-3.1)(-2.1)(-1.1));

  ColumnMatrixMisfitCalculator calc1(x, y, 1);
  checkStats(calc1, 0.999859290353657, 0.00014070964634260719, 0.3, 0.1);

  ColumnMatrixMisfitCalculator calc2(x, y, 2);
  checkStats(calc2, 0.999859290353657, 0.00014070964634260719, 0.0271979048246, 0.1);

  ColumnMatrixMisfitCalculator calc3(x, y, 3);
  checkStats(calc3, 0.999859290353657, 0.00014070964634260719, 0.0026207413942089, 0.1);
}

TEST(ColumnMatrixMisfitCalculatorTest, InfiniteVectors)
{
  ColumnMatrix x = MAKE_COLUMN_MATRIX((1)(2));
  ColumnMatrix y = MAKE_COLUMN_MATRIX((1.1)(std::numeric_limits<double>::infinity()));

  ColumnMatrixMisfitCalculator calc(x, y, 2);
  EXPECT_DOUBLE_EQ(1e6, calc.getCorrelationCoefficient());
  EXPECT_DOUBLE_EQ(1e6, calc.getInverseCorrelationCoefficient());
  EXPECT_DOUBLE_EQ(1e6, calc.getRelativeRMS());
  EXPECT_DOUBLE_EQ(std::numeric_limits<double>::infinity(), calc.getRMS());
}
