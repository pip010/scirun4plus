/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2010 Scientific Computing and Imaging Institute,
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

#ifndef BLOCK_MATRIX_H__
#define BLOCK_MATRIX_H__

#include <boost/shared_ptr.hpp>
#include <iosfwd>

#include <Core/Datatypes/MatrixFwd.h>
#include <Core/Datatypes/share.h>

namespace SCIRun
{
  class BlockMatrixImpl;

  //************************************************************************
  // BlockMatrix: simple representation of block matrices with value type 
  // semantics.  Current interface includes addition, multiplication, and 
  // conversion to DenseMatrix when done processing, as well as block get/set
  // via function call syntax.  Current implementation needs to copy in
  // DenseMatrix blocks; future versions will attempt to avoid copying.
  //************************************************************************

  class SCISHARE BlockMatrix
  {
  public:
    BlockMatrix(size_type partitionRows, size_type partitionCols);
    BlockMatrix(const BlockMatrix& matrix);
    BlockMatrix& operator=(const BlockMatrix& matrix);
    ~BlockMatrix();

    // these enable (i,j) syntax while also adding block matrix size consistency checks.
    class SCISHARE ConstBlockMatrixElementProxy
    {
    public:
      ConstBlockMatrixElementProxy(BlockMatrixImpl& impl, size_type i, size_type j);
      operator const DenseMatrixHandle() const;
    protected:
      BlockMatrixImpl* impl_;
      size_type i_, j_;
    private:
      ConstBlockMatrixElementProxy& operator=(const ConstBlockMatrixElementProxy&);
    };

    class SCISHARE BlockMatrixElementProxy : public ConstBlockMatrixElementProxy
    {
    public:
      BlockMatrixElementProxy(BlockMatrixImpl& impl, size_type i, size_type j);
      BlockMatrixElementProxy& operator=(const DenseMatrix&);
      BlockMatrixElementProxy& operator=(const BlockMatrixElementProxy&);
    };

    BlockMatrixElementProxy operator()(size_type i, size_type j);
    ConstBlockMatrixElementProxy operator()(size_type i, size_type j) const;

    DenseMatrixHandle to_dense() const;

    friend SCISHARE BlockMatrix operator*(const BlockMatrix& lhs, const BlockMatrix& rhs);
    friend SCISHARE BlockMatrix operator+(const BlockMatrix& lhs, const BlockMatrix& rhs);
    friend SCISHARE std::ostream& operator<<(std::ostream& o, const BlockMatrix& m);

  private:
    explicit BlockMatrix(boost::shared_ptr<BlockMatrixImpl> impl);
    boost::shared_ptr<BlockMatrixImpl> impl_;
  };
}

#endif