/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2013 Scientific Computing and Imaging Institute,
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

#include <Core/Math/Mat.h>


TEST(Matsolve3by3Test, IdentityTrivialCase)
{
  double A[3][3];
  double rhs[3];

  // initialize
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      A[i][j] = 0;
    }
  }
  A[0][0] = 1;
  A[1][1] = 1;
  A[2][2] = 1;
  
  rhs[0] = 1;
  rhs[1] = 1;
  rhs[2] = 1;
  
  matsolve3by3(A, rhs);
  std::cout << "rhs result: [" << rhs[0] << " " << rhs[1] << " " << rhs[2] << "]" << std::endl;
  EXPECT_EQ(rhs[0], 1);
  EXPECT_EQ(rhs[1], 1);
  EXPECT_EQ(rhs[2], 1);
}

TEST(Matsolve3by3Test, NonTrivialCase)
{
//  double A[3][3];
//  double rhs[3];
//  
//  // initialize
//  for (int i = 0; i < 3; ++i)
//  {
//    for (int j = 0; j < 3; ++j)
//    {
//      A[i][j] = 0;
//    }
//  }
//  A[0][0] = 1;
//  A[1][1] = 1;
//  A[2][2] = 1;
//  
//  rhs[0] = 1;
//  rhs[1] = 1;
//  rhs[2] = 1;
//  
//  matsolve3by3(A, rhs);
  FAIL() << "TODO";
}

TEST(MinNormLeastSquare3Test, IdentityTrivialCase)
{
  double *A[3];
  A[0] = new double[3];
  A[1] = new double[3];
  A[2] = new double[3];
  double *b = new double[3];
  double *x = new double[3];
  double *bprime = new double[3];
  
  // initialize
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      A[i][j] = 0;
    }
  }
  A[0][0] = 1;
  A[1][1] = 1;
  A[2][2] = 1;
  
  b[0] = 1;
  b[1] = 1;
  b[2] = 1;
  
  min_norm_least_sq_3(A, b, x, bprime, 3);
  std::cout << "bprime result: [" << bprime[0] << " " << bprime[1] << " " << bprime[2] << "]" << std::endl;
  EXPECT_EQ(bprime[0], 1);
  EXPECT_EQ(bprime[1], 1);
  EXPECT_EQ(bprime[2], 1);
}

TEST(MinNormLeastSquare3Test, NonTrivialCase)
{
  FAIL() << "TODO";
}

TEST(MinNormLeastSquare3Test, NonInvertibleInputCase)
{
  double A[5][3];
  // initialize
  for (int i = 0; i < 5; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      A[i][j] = 0;
    }
  }
  A[0][0] = 1;
  A[1][1] = 1;
  A[2][2] = 1;
  
  FAIL() << "TODO";
}