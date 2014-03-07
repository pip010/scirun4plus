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

#ifndef TESTING_UTIL_MATRIXTESTUTILITIES
#define TESTING_UTIL_MATRIXTESTUTILITIES 1

#include <gtest/gtest.h>
#include <boost/test/floating_point_comparison.hpp>

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/MatrixOperations.h>
#include <Core/Math/MiscMath.h>

#include <boost/assign/list_of.hpp> 
#include <boost/assign/std/vector.hpp>
#include <boost/assert.hpp> 
#include <boost/type_traits.hpp>

#include <stdexcept>

namespace SCIRun 
{ 

namespace TestUtils
{

inline long long count_single_value(const Matrix<double>& m, const double val)
{
  return static_cast<long long>(std::count(m.begin(), m.end(), val));
}

inline bool equal_size(const Matrix<double>& m1, const Matrix<double>& m2)
{
  return m1.nrows() == m2.nrows() 
    && m1.ncols() == m2.ncols();
}

inline bool compare_exactly(const Matrix<double>& m1, const Matrix<double>& m2)
{
  return equal_size(m1, m2) &&
    std::equal(m1.begin(), m1.end(), m2.begin());
}

//grr--linkage? I'll use a macro
#define DEFAULT_MATRIX_PERCENT_TOLERANCE 1e-5

//TODO add tolerance parameter, improve error reporting
inline bool compare_with_tolerance(const Matrix<double>& m1, const Matrix<double>& m2, double percentTolerance = DEFAULT_MATRIX_PERCENT_TOLERANCE)
{
  using namespace boost::test_tools;
  return equal_size(m1, m2) &&
    std::equal(m1.begin(), m1.end(), m2.begin(),
    close_at_tolerance<double>(percent_tolerance(percentTolerance)));
}

inline ::testing::AssertionResult compare_with_tolerance_readable(const Matrix<double>& m1, const Matrix<double>& m2, double percentTolerance, int printSize = 50)
{
  return compare_with_tolerance(m1, m2, percentTolerance) ? 
    ::testing::AssertionSuccess() :
  ::testing::AssertionFailure() << "Matrix 1: \n"<< matrix_to_string(m1).substr(0, printSize) << "\nMatrix 2: \n" << matrix_to_string(m2).substr(0, printSize);
}

inline double relative_infinity_norm(const ColumnMatrix& x, const ColumnMatrix& xhat)
{
  ColumnMatrix diff = xhat - x;
  double num = diff.infinity_norm();
  double denom = xhat.infinity_norm();

  return num / denom; 
}

inline ::testing::AssertionResult compare_with_relative_infinity_norm(const ColumnMatrix& x, const ColumnMatrix& xhat, double relativeError = DEFAULT_MATRIX_PERCENT_TOLERANCE, int printSize = 50)
{
  return relative_infinity_norm(x, xhat) < relativeError ? 
    ::testing::AssertionSuccess() :
  ::testing::AssertionFailure() << "ColumnMatrix 1: \n"<< matrix_to_string(x).substr(0, printSize) << "ColumnMatrix 2: \n" << matrix_to_string(xhat).substr(0, printSize);
}

inline ::testing::AssertionResult compare_with_two_norm(const ColumnMatrix& x, const ColumnMatrix& xhat, double error = 1e-15, int printSize = 50)
{
  double delta = (x - xhat).vector_norm();
  return delta < error ? 
    ::testing::AssertionSuccess() :
  ::testing::AssertionFailure() <<
    "Vectors are " << delta << " apart, expect less than " << error << " distance apart.\n" <<
    "ColumnMatrix 1: \n"<< matrix_to_string(x).substr(0, printSize) << "ColumnMatrix 2: \n" << matrix_to_string(xhat).substr(0, printSize);
}

// Note: this approach is limited to matrices with <= 5 columns, due to the implementation of boost::assign::tuple_list_of.
template <class Tuple>
inline DenseMatrix convertDataToMatrix(const std::vector<Tuple>& matrixData) 
{
  typedef typename boost::tuples::element<0,Tuple>::type element0type;
  typedef boost::is_convertible<element0type, double> CanConvertToDouble;
  BOOST_STATIC_ASSERT(CanConvertToDouble::value);

  const int columns = boost::tuples::length<Tuple>::value;
  const bool tupleShouldBeAllOneType = columns * sizeof(element0type) == sizeof(Tuple);
  BOOST_STATIC_ASSERT(tupleShouldBeAllOneType);

  const size_t rows = matrixData.size();
  DenseMatrix matrix(rows, columns);

  for (size_t i = 0; i < rows; ++i)
  {
    const element0type* tupleData = reinterpret_cast<const element0type*>(&matrixData[i]);
    std::copy(tupleData, tupleData + columns, matrix[i]);
  }

  return matrix;
}
  
template <class Tuple>
inline SparseRowMatrix convertDataToSparseRowMatrix(const std::vector<Tuple>& matrixData)
{
  typedef typename boost::tuples::element<0,Tuple>::type element0type;
  typedef boost::is_convertible<element0type, double> CanConvertToDouble;
  BOOST_STATIC_ASSERT(CanConvertToDouble::value);
  
  const size_t columns = boost::tuples::length<Tuple>::value;
  const bool tupleShouldBeAllOneType = columns * sizeof(element0type) == sizeof(Tuple);
  BOOST_STATIC_ASSERT(tupleShouldBeAllOneType);
  if (0 == columns)
    throw std::invalid_argument("input is empty");
  
  const size_t rows = matrixData.size();
  if (0 == rows)
    throw std::invalid_argument("input is empty");
  
  SparseRowMatrix::Builder builder;
  const SparseRowMatrix::Rows& rr = builder.allocate_rows(rows + 1);

  index_type nnz = 0;
  for (size_t r = 0; r < rows; ++r)
  {
    const element0type* tupleData = reinterpret_cast<const element0type*>(&matrixData[r]);
    nnz += std::count_if(tupleData, tupleData + columns, nonzero<double>);
  }
  const SparseRowMatrix::Columns& cc = builder.allocate_columns(nnz);
  const SparseRowMatrix::Storage& d = builder.allocate_data(nnz);

  index_type counter = 0;
  for (size_t r = 0; r < rows; ++r)
  {
    rr[r] = counter;
    // TODO: is there a better way to iterate over contents of tuple?
    const element0type* tupleData = reinterpret_cast<const element0type*>(&matrixData[r]);
    for (size_t c = 0; c < columns; ++c)
    {
      if ( nonzero<double>( *(tupleData + c) ) )
      {
        cc[counter] = c;
        d[counter++] = *(tupleData + c);
      }
    }
  }
  rr[rows] = counter;

  SparseRowMatrix matrix(rows, columns, builder.build(), nnz);
  return matrix;
}
  
inline ColumnMatrix convertColumnDataToMatrix(const std::vector<double>& columnData)
{
  ColumnMatrix matrix(columnData.size());
  std::copy(columnData.begin(), columnData.end(), matrix.get_data_pointer());
  return matrix;
}

template <class Cont>
inline std::vector<typename Cont::value_type> to_vector(const Cont& cont)
{
  if (cont.empty())
    throw std::invalid_argument("Matrix<double> construction will not work without any data");

  std::vector<typename Cont::value_type> v;
  return cont.to_container(v);
}

#define MAKE_DENSE_MATRIX(x) (convertDataToMatrix(to_vector(tuple_list_of x)))

#define MAKE_DENSE_MATRIX_POINTER(x) (convertDataToMatrix(to_vector(tuple_list_of x))).clone()

#define MAKE_COLUMN_MATRIX(x) (convertColumnDataToMatrix(list_of x))

#define MAKE_COLUMN_MATRIX_POINTER(x) (convertColumnDataToMatrix(list_of x)).clone()

#define MAKE_SPARSE_ROW_MATRIX(x) (convertDataToSparseRowMatrix(to_vector(tuple_list_of x)))

#define MAKE_SPARSE_ROW_MATRIX_POINTER(x) (convertDataToSparseRowMatrix(to_vector(tuple_list_of x))).clone()

//TODO improve failure reporting
#define EXPECT_MATRIX_EQ_TOLERANCE(actual, expected, tolerance) EXPECT_TRUE(compare_with_tolerance_readable((actual), (expected), (tolerance)))
#define EXPECT_MATRIX_EQ(actual, expected) EXPECT_MATRIX_EQ_TOLERANCE((actual), (expected), DEFAULT_MATRIX_PERCENT_TOLERANCE)
#define EXPECT_MATRIX_EQ_TO(actual, expected) EXPECT_MATRIX_EQ((actual), MAKE_DENSE_MATRIX(expected))

#define EXPECT_SPARSE_ROW_MATRIX_EQ_TO(actual, expected) EXPECT_MATRIX_EQ((actual), MAKE_SPARSE_ROW_MATRIX(expected))

#define EXPECT_COLUMN_MATRIX_EQ_TO(actual, expected) EXPECT_MATRIX_EQ((actual), MAKE_COLUMN_MATRIX(expected))

#define EXPECT_COLUMN_MATRIX_EQ_BY_RELATIVE_INFINITY_NORM(actual, expected, error) EXPECT_TRUE(compare_with_relative_infinity_norm((actual), (expected), (error)))
#define EXPECT_COLUMN_MATRIX_EQ_BY_TWO_NORM(actual, expected, error) EXPECT_TRUE(compare_with_two_norm((actual), (expected), (error)))
#define EXPECT_COLUMN_MATRIX_EQ_BY_TWO_NORM_WITHIN_MACHINE_EPSILON(actual, expected) EXPECT_TRUE(compare_with_two_norm((actual), (expected)))

#define EXPECT_MATRIX_NOT_EQ_TOLERANCE(actual, expected, tolerance) EXPECT_FALSE(compare_with_tolerance_readable((actual), (expected), (tolerance)))
#define EXPECT_MATRIX_NOT_EQ(actual, expected) EXPECT_MATRIX_NOT_EQ_TOLERANCE((actual), (expected), DEFAULT_MATRIX_PERCENT_TOLERANCE)

// vector of vector approach.  This ought to be an official DenseMatrix constructor.
inline DenseMatrix convertDataToMatrix(const std::vector<std::vector<double> >& matrixData)
{
  const size_t rows = matrixData.size();
  if (0 == rows)
    throw std::invalid_argument("input is empty");
  
  const size_t columns = matrixData[0].size();
  if (0 == columns)
    throw std::invalid_argument("input is empty");

  DenseMatrix matrix(rows, columns);

  for (size_t i = 0; i < rows; ++i)
  {
    if (matrixData[i].size() != columns)
      throw std::invalid_argument("incomplete row");
    std::copy(matrixData[i].begin(), matrixData[i].end(), matrix[i]);
  }
  return matrix;
}

inline SparseRowMatrix convertDataToSparseRowMatrix(const std::vector<std::vector<double> >& matrixData)
{
  const size_t rows = matrixData.size();
  if (0 == rows)
    throw std::invalid_argument("input is empty");
  
  const size_t columns = matrixData[0].size();
  if (0 == columns)
    throw std::invalid_argument("input is empty");
  
  SparseRowMatrix::Builder builder;
  const SparseRowMatrix::Rows& rr = builder.allocate_rows(rows + 1);

  index_type nnz = 0;
  for (size_t r = 0; r < rows; ++r)
  {
    nnz += std::count_if((matrixData[r]).begin(), (matrixData[r]).end(), nonzero<double>);
  }

  const SparseRowMatrix::Columns& cc = builder.allocate_columns(nnz);
  const SparseRowMatrix::Storage& d = builder.allocate_data(nnz);

  index_type counter = 0;
  for (size_t r = 0; r < rows; ++r)
  {
    rr[r] = counter;
    for (size_t c = 0; c < columns; ++c)
    {
      if ( nonzero<double>((matrixData[r])[c]) )
      {
        cc[counter] = c;
        d[counter++] = (matrixData[r])[c];
      }
    }
  }
  rr[rows] = counter;

  SparseRowMatrix matrix(rows, columns, builder.build(), nnz);
  return matrix;
}

// this duplication indicates we should be able to make move ctors to convert between various matrix types on the stack.
// this one fn can just call the other and return a quick shallow copy, only the ptr is copied.
inline DenseMatrixHandle make_random(size_t rows, size_t cols, double randMin = 0, double randMax = 1)
{
  DenseMatrixHandle m(new DenseMatrix(rows, cols));
  const double range = randMax - randMin;
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j)
      (*m)[i][j] = (static_cast<double>(rand()) / RAND_MAX) * range + randMin;
  return m;
}

inline ColumnMatrix make_random_column(size_t rows, double randMin = 0, double randMax = 1)
{
  ColumnMatrix m(rows);
  const double range = randMax - randMin;
  for (size_t i = 0; i < rows; ++i)
    m[i] = (static_cast<double>(rand()) / RAND_MAX) * range + randMin;
  return m;
}

inline ColumnMatrix column_with_integer_sequence(int start, size_t count)
{
  ColumnMatrix c(count);
  for (size_t i = 0; i < count; ++i)
    c[i] = start++;
  return c;
}

}

}

#endif
