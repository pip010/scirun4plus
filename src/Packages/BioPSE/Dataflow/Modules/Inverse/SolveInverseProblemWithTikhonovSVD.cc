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

//    File       : SolveInverseProblemWithTikhonovSVD.cc
//    Author     : Yesim Serinagaoglu & Alireza Ghodrati
//    Date       : 07 Aug. 2001
//    Last update: Dec 2011

#include <Packages/BioPSE/Dataflow/Modules/Inverse/SolveInverseProblemWithTikhonovSVD.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Core/Exceptions/InvalidState.h>
#include <cmath>
#include <iostream>
#include <sstream>

namespace BioPSE 
{

using namespace SCIRun;

class SolveInverseProblemWithTikhonovSVD : public Module 
{
private:
  GuiDouble 	lambda_fix_;
  GuiDouble 	lambda_sld_;
  GuiInt    	have_ui_;
  GuiString 	reg_method_;
  GuiDouble	  lambda_min_;
  GuiDouble	  lambda_max_;
  GuiInt		  lambda_num_;
  GuiDouble 	tex_var_;	

  

public:
  explicit SolveInverseProblemWithTikhonovSVD(GuiContext *context);
  virtual void execute();
private:
  double find_lambda(TikhonovSVDAlgorithm& algo);
  void plot_lcurve_graph(int nLambda, const std::vector<double>& rho, const std::vector<double>& eta, int lambda_index, double lower_y, double lambda);  
};
	
// MODULE MAKER
DECLARE_MAKER(SolveInverseProblemWithTikhonovSVD)

  SolveInverseProblemWithTikhonovSVD::SolveInverseProblemWithTikhonovSVD(GuiContext *context)
    : Module("SolveInverseProblemWithTikhonovSVD", context, Source, "Inverse", "BioPSE"),
      lambda_fix_(context->subVar("lambda_fix")),
      lambda_sld_(context->subVar("lambda_sld")),
      have_ui_(context->subVar("have_ui")),
      reg_method_(context->subVar("reg_method")),
      lambda_min_(context->subVar("lambda_min")),
      lambda_max_(context->subVar("lambda_max")),
      lambda_num_(context->subVar("lambda_num")),
      tex_var_(context->subVar("tex_var"))
{
}

  TikhonovSVDAlgorithm::TikhonovSVDAlgorithm(const ColumnMatrix& matrixMeasDatRHS, const DenseMatrix& matrixU, const DenseMatrix& matrixS, const DenseMatrix& matrixV, const DenseMatrix* matrixX, Method method)
    : matrixMeasDatRHS_(matrixMeasDatRHS), 
    matrixU_(matrixU), 
    matrixS_(matrixS), 
    matrixV_(matrixV), 
    matrixX_(matrixX)
  {
    const int M = matrixU.nrows();
    int N;

    if (method == SVD)
    {
      N = matrixV.nrows();
      //TODO: this code is fine, but ComputeSVD module is sending in extra 0 columns that fail the middle test below.
      // need to standardize the output of ComputeSVD module to match this.
      if (
        //matrixS.nrows() != N //TODO--VERIFY REQUIREMENTS
        //|| 
        //matrixU.ncols() != N //TODO--VERIFY REQUIREMENTS
        //|| 
        matrixMeasDatRHS_.nrows() != M )
      {
        throw InvalidState("Input matrix dimensions incorrect.");
      }
    }
    else
    {
      if (!matrixX_)
        throw InvalidState("X matrix not specified.");

      N = matrixX->nrows();
      const int P = matrixV.nrows();
      if (M < N)
      {
        //throw InvalidState("The forward matrix should be overdetermined.");  //TODO--VERIFY REQUIREMENTS
      }
      if (
        //matrixX->ncols() != N     //TODO--VERIFY REQUIREMENTS
        //|| 
        //matrixS.nrows() != P   //TODO--VERIFY REQUIREMENTS
        //|| 
        //matrixU.ncols() != N   //TODO--VERIFY REQUIREMENTS
        //|| 
        P > N 
        || matrixMeasDatRHS_.nrows() != M)
      {
        throw InvalidState("Input matrix dimensions incorrect.");
      }
    }

    int columns = M, solutionRows = N;
    Uy_ = new ColumnMatrix(matrixU.ncols());
    inverseMat_ = new DenseMatrix(solutionRows, columns),
    solution_ = new ColumnMatrix(solutionRows);

    for (size_type i = 0; i < matrixU_.ncols(); i++)
      (*Uy_)[i] = inner_product(matrixU_, i, matrixMeasDatRHS_);
  }

  TikhonovSVDAlgorithm::~TikhonovSVDAlgorithm() {}
  ColumnMatrixHandle TikhonovSVDAlgorithm::get_solution() const { return solution_; }
  DenseMatrixHandle TikhonovSVDAlgorithm::get_inverse_matrix() const { return inverseMat_; } 
  
////////////////////////////////////////////////////////////////////////////
// THIS FUNCTION returns the inner product of one column of matrix A
// and w , B=A(:,i)'*w
///////////////////////////////////////////////////////////////////////////
double
TikhonovSVDAlgorithm::inner_product(const DenseMatrix& A, size_type col_num, const ColumnMatrix& w)
{
  int nRows = A.nrows();
  double B = 0;
  for (int i = 0; i < nRows; i++)
    B += A[i][col_num] * w[i];
  return B;
}

//////////////////////////////////////////////////////////////////////
// THIS FUNCTION returns regularized solution by tikhonov method
//////////////////////////////////////////////////////////////////////

namespace TikhonovSVDAlgorithmDetail
{
  struct IsAlmostZero : public std::unary_function<double, bool>
  {
    bool operator()(double d) const
    {
      return fabs(d) < 1e-14;
    }
  };

  size_type count_non_zero_entries_in_column(const DenseMatrix& S, size_t column)
  {
    return std::count_if(S[column], S[column] + S.nrows(), std::not1(IsAlmostZero()));
  }
}

void
TikhonovSVDAlgorithm::tikhonov_fun(double lambda) const
{
  ColumnMatrix& X_reg = *solution_;
  DenseMatrix& inverseMatrix = *inverseMat_;
  const ColumnMatrix& Uy = *Uy_;
  const DenseMatrix& U = matrixU_;
  const DenseMatrix& S = matrixS_;
  const DenseMatrix& V = matrixV_;

  using namespace TikhonovSVDAlgorithmDetail;
  if (S.ncols() == 1)
  { // SVD case
    int rank = count_non_zero_entries_in_column(S, 0);
    const size_type v_rows = V.nrows();
    X_reg.zero();
    for (int i = 0; i < rank; i++)
    {
      double filterFactor_i = S[i][0] / (lambda*lambda + S[i][0] * S[i][0]) * Uy[i];
      for (int j = 0; j < v_rows; j++)
      {
        X_reg[j] += filterFactor_i * V[j][i];
      }
    }

    //Finding Regularized Inverse Matrix
    if (V.ncols() == U.ncols()) //TODO--VERIFY REQUIREMENTS
    {
      DenseMatrix Mat_temp(V.nrows(), V.ncols());
      for (int i = 0; i < rank; i++)
      {
        double temp = S[i][0] / (lambda * lambda + S[i][0] * S[i][0]);
        for (int j = 0; j < V.nrows(); j++)
        {
          Mat_temp[j][i] = temp * V[j][i];
        }
      }
      DenseMatrixHandle Utranspose(U.make_transpose());
      Mult(inverseMatrix, Mat_temp, *Utranspose);
    }
    else
      inverseMatrix.zero();
  }
  else
  { //GSVD case
    const DenseMatrix& X = *matrixX_;

    int rank0 = count_non_zero_entries_in_column(S, 0); 
    int rank1 = count_non_zero_entries_in_column(S, 1); 
    if (rank0 != rank1)
    {
      throw InvalidState("singular value vectors do not have same rank");
    }
    int rank = rank0;

    X_reg.zero();

    ////TODO: make member functions for scalar*row/column, return new column vs in-place

    for (int i = 0; i < rank; i++)
    {
      double filterFactor_i = S[i][0] / (lambda * lambda * S[i][1] * S[i][1] + S[i][0] * S[i][0]) * Uy[i];
      for (int j = 0; j < X.nrows(); j++)
      {
        X_reg[j] += filterFactor_i * X[j][i];
      }
    }

    int minDimension = std::min(U.nrows(), U.ncols());
    for (int i = rank; i < minDimension; i++)
    {
      for (int j = 0; j < X.nrows(); j++)
      {
        X_reg[j] += Uy[i] * X[j][i];
      }
    }
		
    //Finding Regularized Inverse Matrix
    if (V.ncols() == U.ncols()) //TODO--VERIFY REQUIREMENTS
    {
      DenseMatrix Mat_temp = X;
      for (int i = 0; i < rank; i++)
      {
        double temp = S[i][0] / (lambda*lambda*S[i][1] * S[i][1] + S[i][0] * S[i][0]);
        for (int j = 0; j < X.nrows(); j++)
        {
          Mat_temp[j][i] = temp * X[j][i];
        }
      }
      DenseMatrixHandle Utranspose(U.make_transpose());
      Mult(inverseMatrix, Mat_temp, *Utranspose);
    }
    else
      inverseMatrix.zero();
  }
}



/////////////////////////////////////////////////////////////	
// THIS FUNCTION Calculate ro and eta for lcurve
/////////////////////////////////////////////////////////////
void
TikhonovSVDAlgorithm::prep_lcurve_data(double& rho, double& eta, double lam)
{
  const ColumnMatrix& Uy = *Uy_;
  const DenseMatrix& U = matrixU_;
  const DenseMatrix& S = matrixS_;
  const DenseMatrix& V = matrixV_;
  const ColumnMatrix& y = matrixMeasDatRHS_;

  ColumnMatrix AX_reg(U.nrows());
  
  int rank = S.nrows();
  if (S.ncols() == 1)
  {
    ColumnMatrix X_reg(V.nrows());		
    double temp = S[0][0]/(lam*lam + S[0][0] * S[0][0]) * Uy[0];
    for (int j = 0; j < V.nrows(); j++)
    {
      X_reg[j] = temp * V[j][0];
    }
    for (int j = 0; j < U.nrows(); j++)
    {
      AX_reg[j] = temp * S[0][0] * U[j][0];
    }
    for (int i = 1; i < rank; i++)
    {
      temp = S[i][0] / (lam*lam + S[i][0] * S[i][0]) * Uy[i];
      for (int j = 0; j < V.nrows(); j++)
      {
        X_reg[j] = X_reg[j] + temp * V[j][i];
      }
      temp *= S[i][0];
      for (int j = 0; j < U.nrows(); j++)
      {
        AX_reg[j] += temp * U[j][i];
      }
    }
    // Calculate the norm of Ax-b and Rx  
    rho = 0;
    eta = 0;
    for (int i = 0; i < U.nrows(); i++)
    {
      AX_reg[i] -= y[i];
      rho += AX_reg[i] * AX_reg[i];
    }
    for (int i = 0; i < V.nrows(); i++)
    {
      eta += X_reg[i] * X_reg[i];
    }
    rho = sqrt(rho);
    eta = sqrt(eta);
  }
  else
  {
    ColumnMatrix RX_reg(V.nrows());
    double temp = S[0][0]/(lam*lam*S[0][1] * S[0][1] + S[0][0] * S[0][0]) * Uy[0];

    for (int j = 0; j < U.nrows(); j++)
    {
      AX_reg[j] = temp * S[0][0] * U[j][0];
    }
    for (int j = 0; j < V.nrows(); j++)
    {
      RX_reg[j] = temp * S[0][1] * V[j][0];
    }
    for (int i = 1; i < rank; i++)
    {
      double temp = S[i][0] / (lam*lam*S[i][1] * S[i][1] + S[i][0] * S[i][0]) * Uy[i];
      double temp1 = temp * S[i][0];
      for (int j = 0; j < U.nrows(); j++)
      {
        AX_reg[j] += temp1 * U[j][i];
      }
      temp1 = temp * S[i][1];
      for (int j = 0; j < V.nrows(); j++)
      {
        RX_reg[j] += temp1 * V[j][i];
      }
    }
    const DenseMatrix& X = *matrixX_;
    for (int i = rank; i < X.ncols(); i++)
    {
      for (int j = 0; j < U.nrows(); j++)
      {
        AX_reg[j] += Uy[i] * U[j][i];
      }
    }		

    // Calculate the norm of Ax-b and Rx  
    rho = 0;
    eta = 0;
    for (int i = 0; i < U.nrows(); i++)
    {
      AX_reg[i] -= y[i];
      rho += AX_reg[i] * AX_reg[i];
    }
    for (int i = 0; i < V.nrows(); i++)
    {
      eta += RX_reg[i] * RX_reg[i];
    }
    rho = sqrt(rho);
    eta = sqrt(eta);
  }
}

////////////////////////////////////////////
// FIND CORNER
////////////////////////////////////////////
double
TikhonovSVDAlgorithm::find_corner(const std::vector<double> &rho, const std::vector<double> &eta,
                        const std::vector<double> &lambdaArray, ColumnMatrix& kapa,
                        int& lambda_index, int nLambda)
{
  std::vector<double> deta(nLambda), ddeta(nLambda), drho(nLambda), ddrho(nLambda), lrho(nLambda), leta(nLambda);

  for (int i = 0; i < nLambda; i++)
  {
    lrho[i] = log10(rho[i]);
    leta[i] = log10(eta[i]);	
    if (i > 0)
    {
      deta[i] = (leta[i]-leta[i-1])/(lambdaArray[i]-lambdaArray[i-1]);
      drho[i] = (lrho[i]-lrho[i-1])/(lambdaArray[i]-lambdaArray[i-1]);
    }
    if (i > 1)
    {
      ddeta[i] = (deta[i]-deta[i-1])/(lambdaArray[i]-lambdaArray[i-1]);
      ddrho[i] = (drho[i]-drho[i-1])/(lambdaArray[i]-lambdaArray[i-1]);
    }
  }
  drho[0] = drho[1];
  deta[0] = deta[1];
  ddrho[0] = ddrho[2];
  ddrho[1] = ddrho[2];
  ddeta[0] = ddeta[2];
  ddeta[1] = ddeta[2];	

  double maxKapa = -1e10;
  //TODO: use std::max_element here
  lambda_index = 0;
  for (int i = 0; i < nLambda; i++)
  {
    kapa[i] = std::abs(drho[i]*ddeta[i] - ddrho[i]*deta[i]) / 
      sqrt(pow((deta[i]*deta[i]+drho[i]*drho[i]),3));  
    if (kapa[i] > maxKapa)
    {	
      maxKapa = kapa[i];
      lambda_index = i;
    }
  }
  double lambda_cor = lambdaArray[lambda_index];
  return lambda_cor;
}



/////////////////////////////////////////
// MODULE EXECUTION
/////////////////////////////////////////
void
SolveInverseProblemWithTikhonovSVD::execute()
{
  // DEFINE MATRIX HANDLES FOR INPUT/OUTPUT PORTS
  MatrixHandle hMatrixMeasDat, hMatrixU, hMatrixS, hMatrixV ,hMatrixX;
    	    
  if (! get_input_handle("U", hMatrixU)) return;
  if (! get_input_handle("S", hMatrixS)) return;
  if (! get_input_handle("V", hMatrixV)) return;

  if (!get_input_handle("MeasuredPots", hMatrixMeasDat)) return;

  // TYPE CHECK

  // TODO: need to discuss lifetime/memory management of these
  // potential unmanaged pointers here--if input is not of the desired type, new matrices are allocated.
  // to be handled by future Matrix memory refactoring.
  
  // TODO: if the conversion to ColumnMatrix is ever removed,
  // check hMatrixMeasDat dimensions (should be Mx1)
  const ColumnMatrix* matrixMeasDatD = hMatrixMeasDat->column();

  const DenseMatrix* matrixU = hMatrixU->dense();
  const DenseMatrix* matrixS = hMatrixS->dense();	
  const DenseMatrix* matrixV = hMatrixV->dense();
  const DenseMatrix* matrixX = 0;
	
  //TODO: move this verification code to algo class ctor.

  //IDEA: defer X parameter with boost::function, then ctor can hold this logic.
  TikhonovSVDAlgorithm::Method method;

  if (matrixS->ncols() == 1)
    method = TikhonovSVDAlgorithm::SVD;
  else if (matrixS->ncols() == 2)
  {
    method = TikhonovSVDAlgorithm::GSVD;

    if (!get_input_handle("X", hMatrixX)) 
      return;
    matrixX = hMatrixX->dense();
  }
  else
  {
    error("S matrix dimensions incorrect.");
    return;	
  }

  try
  {
	  TikhonovSVDAlgorithm algo(*matrixMeasDatD, *matrixU, *matrixS, *matrixV, matrixX, method);
		
	  //find_lambda is tightly coupled with Module calls, needs more work.
	  double lambda = find_lambda(algo);
	
	  algo.tikhonov_fun(lambda);
	
	  //...........................................................
	  // SEND RESULTS TO THE OUTPUT PORTS
	  ColumnMatrixHandle solution = algo.get_solution();
	  send_output_handle("InverseSoln", solution);
	
	  ColumnMatrixHandle regParameter = new ColumnMatrix(1);
	  (*regParameter)[0] = lambda;
	  send_output_handle("RegParam", regParameter);
	
	  DenseMatrixHandle inverse = algo.get_inverse_matrix();
	  send_output_handle("RegInverseMat", inverse);
  }
  catch (InvalidState& e)
  {
  	error(e.message());
    return;
  }
}

double SolveInverseProblemWithTikhonovSVD::find_lambda(TikhonovSVDAlgorithm& algo)
{
  double lambda = 0;
  
  if (reg_method_.get() == "single")
  {
    // Use single fixed lambda value, entered in UI
    lambda = lambda_fix_.get();
    msg_stream_ << "  method = " << reg_method_.get() << "\n";//DISCARD
  }
  else if (reg_method_.get() == "slider")
  {
    // Use single fixed lambda value, select via slider
    lambda = tex_var_.get(); 
    msg_stream_ << "  method = " << reg_method_.get() << "\n";//DISCARD
  }
  else if (reg_method_.get() == "lcurve")
  {
    // Use L-curve, lambda from corner of the L-curve
    msg_stream_ << "method = " << reg_method_.get() << "\n";//DISCARD

    //TODO: move all this state into algo methods--but then what about plot_lcurve_graph?  too many return values...

    int nLambda = lambda_num_.get();
    ColumnMatrix kapa(nLambda);
    std::vector<double> lambdaArray(nLambda), rho(nLambda), eta(nLambda);

    lambdaArray[0] = lambda_min_.get();
    double lam_step = pow(10.0, log10(lambda_max_.get()/lambda_min_.get())/(nLambda-1));

    for (int j = 0; j < nLambda; j++)
    {
      if (j) 
        lambdaArray[j] = lambdaArray[j-1] * lam_step;

      double rhotemp, etatemp;
      algo.prep_lcurve_data(rhotemp, etatemp, lambdaArray[j]);
      rho[j] = rhotemp;
      eta[j] = etatemp;
    }

    int lambda_index = 0;
    lambda = algo.find_corner(rho, eta, lambdaArray, kapa, lambda_index, nLambda);

    double lower_y = eta[0]/10.;
    if (eta[nLambda-1] < lower_y)  
      lower_y = eta[nLambda-1];

    plot_lcurve_graph(nLambda, rho, eta, lambda_index, lower_y, lambda);
  } 

  return lambda;
}

void SolveInverseProblemWithTikhonovSVD::plot_lcurve_graph(int nLambda, const std::vector<double>& rho, const std::vector<double>& eta, int lambda_index, double lower_y, double lambda)
{
  std::ostringstream str;
  str << get_id() << " plot_graph \" ";
  for (int i = 0; i < nLambda; i++)
    str << log10( rho[i] ) << " " << log10( eta[i] ) << " ";
  str << "\" \" " << log10( rho[0]/10 ) << " " << log10( eta[lambda_index] ) << " ";
  str << log10( rho[lambda_index] ) << " " << log10( eta[lambda_index] ) << " ";
  str << log10( rho[lambda_index] ) << " " << log10( lower_y ) << " \" ";
  str << lambda << " ; update idletasks";
  TCLInterface::execute(str.str());
}


} // End namespace BioPSE
