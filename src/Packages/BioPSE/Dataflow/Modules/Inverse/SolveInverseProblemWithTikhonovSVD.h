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

//    File       : SolveInverseProblemWithTikhonovSVD.h
//    Author     : Yesim Serinagaoglu & Alireza Ghodrati
//    Date       : 07 Aug. 2001
//    Last update: Dec 2011

#ifndef Packages_BioPSE_Dataflow_Modules_Inverse_SolveInverseProblemWithTikhonovSVD_h
#define Packages_BioPSE_Dataflow_Modules_Inverse_SolveInverseProblemWithTikhonovSVD_h

#include <vector>
#include <boost/utility.hpp>
#include <Core/Datatypes/MatrixFwd.h>
#include <Packages/BioPSE/Dataflow/Modules/Inverse/share.h>

namespace BioPSE 
{
  //TODO: move to Algorithm library
  //TODO: unit test
  class SCISHARE TikhonovSVDAlgorithm : boost::noncopyable
  {
  public:
    enum Method {SVD, GSVD};

    TikhonovSVDAlgorithm(const SCIRun::ColumnMatrix& matrixMeasDatRHS, const SCIRun::DenseMatrix& matrixU, const SCIRun::DenseMatrix& matrixS,
      const SCIRun::DenseMatrix& matrixV, const SCIRun::DenseMatrix* matrixX, Method method);
    
    ~TikhonovSVDAlgorithm();
    SCIRun::ColumnMatrixHandle get_solution() const;
    SCIRun::DenseMatrixHandle get_inverse_matrix() const;
    
    void tikhonov_fun(double lam) const;
    void prep_lcurve_data(double& rho, double& eta, double lam);
    double find_corner(const std::vector<double>& rho, const std::vector<double>& eta, const std::vector<double>& lambdaArray, 
      SCIRun::ColumnMatrix& kapa, int& lambda_index, int nLambda);

  private:
    const SCIRun::ColumnMatrix& matrixMeasDatRHS_;
    const SCIRun::DenseMatrix& matrixU_;
    const SCIRun::DenseMatrix& matrixS_;	
    const SCIRun::DenseMatrix& matrixV_;
    const SCIRun::DenseMatrix* matrixX_;

    SCIRun::ColumnMatrixHandle Uy_;
    SCIRun::DenseMatrixHandle inverseMat_;
    SCIRun::ColumnMatrixHandle solution_;

    static double inner_product(const SCIRun::DenseMatrix& A, SCIRun::size_type col_num, const SCIRun::ColumnMatrix& w);
  };
}

#endif