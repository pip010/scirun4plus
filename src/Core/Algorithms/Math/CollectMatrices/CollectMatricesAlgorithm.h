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

#ifndef CORE_ALGORITHMS_MATH_COLLECTMATRICES_H
#define CORE_ALGORITHMS_MATH_COLLECTMATRICES_H

#include <boost/noncopyable.hpp>
#include <Core/Algorithms/Util/AlgoBase.h>
#include <Core/Datatypes/MatrixFwd.h>
#include <Core/Datatypes/SparseRowMatrixFromMap.h>  //TODO: remove
#include <Core/Algorithms/Math/share.h>

namespace SCIRun {
  namespace Algo {

    class SCISHARE CollectMatricesAlgorithmBase : public SCIRunAlgo::AlgoBase, boost::noncopyable
    {
    public:
      virtual MatrixHandle concat_cols(MatrixHandle m1H, MatrixHandle m2H) const = 0;
      virtual MatrixHandle concat_rows(MatrixHandle m1H, MatrixHandle m2H) const = 0;
    };

    class SCISHARE CollectDenseMatricesAlgorithm : public CollectMatricesAlgorithmBase
    {
    public:
      virtual MatrixHandle concat_cols(MatrixHandle m1H, MatrixHandle m2H) const;
      virtual MatrixHandle concat_rows(MatrixHandle m1H, MatrixHandle m2H) const;
    private:
      void copy_matrix(MatrixHandle mh, DenseMatrix& out) const;
    };

    class SCISHARE CollectSparseRowMatricesAlgorithm : public CollectMatricesAlgorithmBase
    {
    public:
      MatrixHandle concat_cols(MatrixHandle m1H, MatrixHandle m2H) const;
      MatrixHandle concat_rows(MatrixHandle m1H, MatrixHandle m2H) const;
    private:
      void check_args(MatrixHandle m1H, MatrixHandle m2H) const;
      void copy_shifted_contents(SparseRowMatrix* sparse, SparseRowMatrixFromMap::Values& shiftedValues, size_type rowShift, size_type columnShift) const;
    };
  }
}

#endif