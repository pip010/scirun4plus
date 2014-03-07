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

#include <Core/Datatypes/BlockMatrix.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <numeric>
#include <boost/lambda/core.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <stdexcept>

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/ColumnMatrix.h>

namespace SCIRun
{
namespace BlockMatrixDetail
{
  //Idea: DenseMatrix will inherit most of this functionality once we templatize the Matrix classes.
template <typename ElementType>
class BasicMatrix
{
public:
  typedef BasicMatrix<ElementType> my_type;
  typedef ElementType element_type;

  //NOTE: do not use boost::matrix types outside of this class!  We don't want to introduce yet another matrix type into the system.
  typedef boost::numeric::ublas::matrix<ElementType> ImplType;
  typedef typename ImplType::size_type size_type;

public:
  //dummy ctor to allow block multiply.
  BasicMatrix(size_type size) : empty_(true) {}

  BasicMatrix(const ImplType& data) : impl_(data), empty_(false) {}
  
  BasicMatrix() : empty_(true) {}
  BasicMatrix(size_type rows, size_type cols) : impl_(rows, cols), empty_(false) {}

  BasicMatrix(const DenseMatrix& dense) : impl_(dense.nrows(), dense.ncols())
  {
    // TODO: replace double when DenseMatrix is a template
    const bool isDouble = boost::is_same<ElementType, double>::value;
    BOOST_STATIC_ASSERT(isDouble);

    // faster way to do this?
    for (index_type i = 0; i < dense.nrows(); ++i)
      for (index_type j = 0; j < dense.ncols(); ++j)
        impl_(i,j) = dense[i][j];
  }

  size_type rows() const { return impl_.size1(); }
  size_type cols() const { return impl_.size2(); }

  my_type& operator+=(const my_type& rhs)
  {
    if (empty_)
    {
      impl_ = rhs.impl_;
      empty_ = false;
    }
    else
      impl_ += rhs.impl_;
    return *this;
  }

public:
  
  size_type totalRows() const
  {
    return std::accumulate(impl_.begin1(), impl_.end1(), static_cast<size_t>(0), RowSizeAccumulator());
  }

  size_type totalCols() const
  {
    return std::accumulate(impl_.begin2(), impl_.end2(), static_cast<size_t>(0), ColumnSizeAccumulator());
  }
  
  //TODO: there will need to be checks in the interface to BlockMatrix, that make sure dimensions of blocks match up with previously added blocks.
  void insert_block_element(size_type i, size_type j, const ElementType& element) 
  { 
    if (!numRowsEqualAcrossRowPartition(i, element.rows()))
      throw std::invalid_argument("Invalid block added, row count doesn't match this row partition.");
    if (!numColsEqualAcrossColPartition(j, element.cols()))
      throw std::invalid_argument("Invalid block added, column count doesn't match this column partition.");
    impl_.insert_element(i, j, element); 
  }

  const ElementType& operator()(size_type i, size_type j) const { return impl_(i,j); }
  ElementType& operator()(size_type i, size_type j) { return impl_(i,j); }


public:
  //TODO: this ought to belong in a template partial specialization.  
  friend class ConstBlockMatrixIterator;
  class ConstBlockMatrixIterator : public std::iterator<std::forward_iterator_tag, typename ElementType::element_type>
  {
  public:
    typedef ConstBlockMatrixIterator my_type;
    typedef BasicMatrix<ElementType> matrix_type;
    typedef typename std::iterator<std::forward_iterator_tag, typename ElementType::element_type> my_base;
    typedef typename my_base::value_type value_type;
    typedef typename my_base::pointer pointer;
    
  private:
    const matrix_type* blockMatrix_;
    size_type block_i_, block_j_, rowInBlock_;

    typedef typename ElementType::ImplType block_impl_type;
    typedef typename boost::numeric::ublas::matrix_row<const block_impl_type> block_row_type;
    typedef typename block_row_type::iterator block_row_iterator;

    boost::shared_ptr<block_row_type> currentRow_;
    block_row_iterator rowIter_;
    size_type numRowPartitions_, numColPartitions_;
    size_type currentBlockRows_;
    bool isEnd_;

  public:
    explicit ConstBlockMatrixIterator(const BasicMatrix<ElementType>* matrix = 0) : blockMatrix_(matrix),
      block_i_(0), block_j_(0), rowInBlock_(0)
    {
      if (matrix)
      {
        numRowPartitions_ = matrix->rows();
        numColPartitions_ = matrix->cols();
        if (0 == numColPartitions_ || 0 == numRowPartitions_)
          throw std::invalid_argument("Cannot iterate over degenerate matrix");

        setCurrentBlockRows();
        resetRow();
        isEnd_ = false;
      }
      else
      {
        isEnd_ = true;
      }
    }
  private:
    void resetRow()
    {
      currentRow_.reset(new block_row_type((*blockMatrix_)(block_i_, block_j_).impl_, rowInBlock_));
      rowIter_ = currentRow_->begin();
    }

    void setCurrentBlockRows()
    {
      currentBlockRows_ = (*blockMatrix_)(block_i_, block_j_).rows();
    }

  public:
    const value_type& operator*() const
    {
      if (isEnd_)
        throw std::out_of_range("Cannot dereference end iterator");

      return *rowIter_;
    }

    const pointer operator->() const
    {
      return &(**this);
    }

    my_type& operator++()
    {
      if (!isEnd_)
      {
        ++rowIter_;
        if (rowIter_ == currentRow_->end())
        {
          // move to next block.
          ++block_j_;
          if (block_j_ == numColPartitions_)
          {
            block_j_ = 0;
            ++rowInBlock_;
            
            if (rowInBlock_ == currentBlockRows_)
            {
              rowInBlock_ = 0;
              ++block_i_;
              if (block_i_ == numRowPartitions_)
                isEnd_ = true;
              if (!isEnd_)
                setCurrentBlockRows();
            }
          }
          if (!isEnd_)
          {
            resetRow();
          }
        }
      }
      
      return *this;
    }

    my_type operator++(int)
    {
      my_type orig(*this);
      ++(*this);
      return orig;
    }

    bool operator==(const my_type& rhs)
    {
      if (isEnd_ && rhs.isEnd_)
        return true;

      if (isEnd_ != rhs.isEnd_)
        return false;

      return blockMatrix_ == rhs.blockMatrix_
        && rowIter_ == rhs.rowIter_;
    }

    bool operator!=(const my_type& rhs)
    {
      return !(*this == rhs);
    }
  };


  ConstBlockMatrixIterator block_begin() const
  {
    return ConstBlockMatrixIterator(this);
  }

  ConstBlockMatrixIterator block_end() const
  {
    return ConstBlockMatrixIterator();
  }

public:
  //TODO: for some reason couldn't get friend iterator class to work with this as private.
  ImplType impl_;
private:
  bool empty_;
  
  struct RowSizeAccumulator 
  {
    size_type operator()(size_type sum, const ElementType& m) const
    {
      return sum + m.rows();   
    }
  };

  struct ColumnSizeAccumulator 
  {
    size_type operator()(size_type sum, const ElementType& m) const
    {
      return sum + m.cols();   
    }
  };

  //TODO: refactor duplication between these two helpers, without getting too unreadable. Macro?
  bool numRowsEqualAcrossRowPartition(size_type rowToCheck, size_type newNumRows) 
  {
    boost::numeric::ublas::matrix_row<ImplType> row(impl_, rowToCheck);
    BOOST_FOREACH(typename ImplType::value_type& block, row)
    {
      size_type existingBlockRows = block.rows();
      if (existingBlockRows != 0 && newNumRows != existingBlockRows)
        return false;
    }
    return true;

    //a bit too yucky, broke down and wrote the loop.  Need to study boost::filtered_range/iterator.
    //return row.end() == std::find_if(row.begin(), row.end(), boost::bind(std::not_equal_to<size_type>(), boost::bind(&element_type::rows, _1), newNumRows));
  }

  bool numColsEqualAcrossColPartition(size_type colToCheck, size_type newNumCols) 
  {
    boost::numeric::ublas::matrix_column<ImplType> col(impl_, colToCheck);
    BOOST_FOREACH(typename ImplType::value_type& block, col)
    {
      size_type existingBlockCols = block.cols();
      if (existingBlockCols != 0 && newNumCols != existingBlockCols)
        return false;
    }
    return true;
  }
};

template <typename T>
BasicMatrix<T> operator*(const BasicMatrix<T>& a, const BasicMatrix<T>& b)
{
  return BasicMatrix<T>(prod(a.impl_, b.impl_));
}

template <typename T>
BasicMatrix<T> operator+(const BasicMatrix<T>& a, const BasicMatrix<T>& b)
{
  return BasicMatrix<T>(a.impl_ + b.impl_);
}

template <typename T>
std::ostream& operator<<(std::ostream& o, const BasicMatrix<T>& m)
{
  return o << m.impl_;
}

}

class BlockMatrixImpl
{
public:
  typedef BlockMatrixDetail::BasicMatrix<BlockMatrixDetail::BasicMatrix<double> > matrix_impl_type;  //can remove this layer of indirection after DenseMatrix is templatized.
  BlockMatrixImpl(size_type rows, size_type cols) : matrix_(rows, cols) {}
  BlockMatrixImpl(const matrix_impl_type& matrix) : matrix_(matrix) {} //really want this to be move ctor

  void insert_block_element(size_type i, size_type j, const DenseMatrix& block)
  {
    matrix_.insert_block_element(i, j, block);
  }

  DenseMatrixHandle to_dense() const 
  {
    DenseMatrixHandle dense(new DenseMatrix(matrix_.totalRows(), matrix_.totalCols()));

    std::copy(matrix_.block_begin(), matrix_.block_end(), dense->get_data_pointer());

    return dense;
  }

  matrix_impl_type matrix_;
};

BlockMatrix::BlockMatrix(size_type rows, size_type cols)
  : impl_(new BlockMatrixImpl(rows, cols))
{
}

BlockMatrix::BlockMatrix(boost::shared_ptr<BlockMatrixImpl> impl)
: impl_(impl)
{
}

BlockMatrix::~BlockMatrix() {}

BlockMatrix::BlockMatrixElementProxy::BlockMatrixElementProxy(BlockMatrixImpl& impl, size_type i, size_type j) : ConstBlockMatrixElementProxy(impl, i, j) {}

BlockMatrix::ConstBlockMatrixElementProxy::ConstBlockMatrixElementProxy(BlockMatrixImpl& impl, size_type i, size_type j) : impl_(&impl), i_(i), j_(j) {}

BlockMatrix::BlockMatrixElementProxy& BlockMatrix::BlockMatrixElementProxy::operator=(const DenseMatrix& block)
{
  impl_->insert_block_element(i_, j_, block);
  return *this;
}

BlockMatrix::BlockMatrixElementProxy& BlockMatrix::BlockMatrixElementProxy::operator=(const BlockMatrixElementProxy& block)
{
  DenseMatrixHandle blockData(block);
  impl_->insert_block_element(i_, j_, *blockData);
  return *this;
}

BlockMatrix::ConstBlockMatrixElementProxy::operator const DenseMatrixHandle() const
{
  //this ugly copying will go away once DenseMatrix is templatized.
  if (i_ >= static_cast<size_type>( impl_->matrix_.rows() ) || j_ >= static_cast<size_type>( impl_->matrix_.cols() ))
    throw std::out_of_range("Out of matrix bounds");
  const BlockMatrixDetail::BasicMatrix<double>& blockData = impl_->matrix_(i_, j_);
  DenseMatrixHandle block(new DenseMatrix(blockData.rows(), blockData.cols()));
  for (index_type i = 0; i < block->nrows(); ++i)
    for (index_type j = 0; j < block->ncols(); ++j)
      (*block)[i][j] = blockData(i,j);
  return block;
}

BlockMatrix::BlockMatrixElementProxy BlockMatrix::operator()(size_type i, size_type j)
{
  return BlockMatrixElementProxy(*impl_, i, j);
}

BlockMatrix::ConstBlockMatrixElementProxy BlockMatrix::operator()(size_type i, size_type j) const
{
  return ConstBlockMatrixElementProxy(*impl_, i, j);
}

DenseMatrixHandle BlockMatrix::to_dense() const
{
  return impl_->to_dense();
}

BlockMatrix operator*(const BlockMatrix& lhs, const BlockMatrix& rhs)
{
  boost::shared_ptr<BlockMatrixImpl> impl(new BlockMatrixImpl(lhs.impl_->matrix_ * rhs.impl_->matrix_));
  return BlockMatrix(impl);
}

BlockMatrix operator+(const BlockMatrix& lhs, const BlockMatrix& rhs)
{
  boost::shared_ptr<BlockMatrixImpl> impl(new BlockMatrixImpl(lhs.impl_->matrix_ + rhs.impl_->matrix_));
  return BlockMatrix(impl);
}

std::ostream& operator<<(std::ostream& o, const BlockMatrix& m)
{
  return o << m.impl_->matrix_;
}

BlockMatrix::BlockMatrix(const BlockMatrix& matrix) : impl_(new BlockMatrixImpl(*matrix.impl_))
{
}

BlockMatrix& BlockMatrix::operator=(const BlockMatrix& matrix)
{
  impl_.reset(new BlockMatrixImpl(*matrix.impl_));
  return *this;
}

}
