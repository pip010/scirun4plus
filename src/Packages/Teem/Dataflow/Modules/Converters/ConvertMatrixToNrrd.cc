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


/*
 *  ConvertMatrixToNrrd.cc: Converts a SCIRun Matrix to Nrrd(s).  
 *
 *  Written by:
 *   Darby Van Uitert
 *   April 2004
 *
 */

#include <Dataflow/Network/Module.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <Dataflow/Network/Ports/NrrdPort.h>

#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>

#include <Dataflow/Network/Ports/MatrixPort.h>

namespace SCITeem {

using namespace SCIRun;

class ConvertMatrixToNrrd : public Module {
public:
  int          matrix_generation_;

  ConvertMatrixToNrrd(GuiContext*);

  virtual ~ConvertMatrixToNrrd();

  virtual void execute();

  void create_and_send_column_matrix_nrrd(MatrixHandle mat);
  void create_and_send_dense_matrix_nrrd(MatrixHandle mat);
  void create_and_send_sparse_matrix_nrrd(MatrixHandle mat);

};


DECLARE_MAKER(ConvertMatrixToNrrd)
ConvertMatrixToNrrd::ConvertMatrixToNrrd(GuiContext* ctx)
  : Module("ConvertMatrixToNrrd", ctx, Source, "Converters", "Teem"),
    matrix_generation_(-1)
{
}


ConvertMatrixToNrrd::~ConvertMatrixToNrrd()
{
}


void
ConvertMatrixToNrrd::execute()
{
  MatrixHandle matH;
  if (!get_input_handle("Matrix", matH)) return;

  if (matrix_generation_ != matH->generation ||
      !oport_cached("Data") &&
      !oport_cached("Rows") &&
      !oport_cached("Columns"))
  {
    matrix_generation_ = matH->generation;
    Matrix<double>* matrix = matH.get_rep();
    
    if (matrix_is::column(matrix)) {
      create_and_send_column_matrix_nrrd(matH);
    } else if (matrix_is::dense(matrix)) {
      create_and_send_dense_matrix_nrrd(matH);
    } else {
      create_and_send_sparse_matrix_nrrd(matH);
    }
  }
}


void
ConvertMatrixToNrrd::create_and_send_column_matrix_nrrd(MatrixHandle matH)
{
  ColumnMatrix* matrix = matrix_cast::as_column(matH);

  size_t size[NRRD_DIM_MAX];
  size[0] = matrix->nrows();
  
  NrrdData *nd = new NrrdData();
  nrrdAlloc_nva(nd->nrrd_, nrrdTypeDouble, 1, size);

  nrrdAxisInfoSet_va(nd->nrrd_, nrrdAxisInfoKind, nrrdKindDomain);
  nrrdAxisInfoSet_va(nd->nrrd_, nrrdAxisInfoLabel, "column-data");

  double *val = (double*)nd->nrrd_->data;
  double *data = matrix->get_data_pointer();

  for(unsigned int i=0; i<size[0]; i++) {
    *val = *data;
    ++data;
    ++val;
  }

  // Send the data nrrd.
  NrrdDataHandle dataH(nd);
  send_output_handle("Data", dataH);
}


void
ConvertMatrixToNrrd::create_and_send_dense_matrix_nrrd(MatrixHandle matH)
{
  DenseMatrix* matrix = matrix_cast::as_dense(matH);

  size_type rows = matrix->nrows();
  size_type cols = matrix->ncols();
  
  NrrdData *nd = new NrrdData();
  size_t size[NRRD_DIM_MAX];
  size[0] = cols; size[1] = rows;
  nrrdAlloc_nva(nd->nrrd_, nrrdTypeDouble, 2, size);

  char *labels[NRRD_DIM_MAX];
  labels[0] = airStrdup("dense-columns");
  labels[1] = airStrdup("dense-rows");
  nrrdAxisInfoSet_nva(nd->nrrd_, nrrdAxisInfoLabel, labels);
  nd->nrrd_->axis[0].kind = nrrdKindDomain;
  nd->nrrd_->axis[1].kind = nrrdKindDomain;

  double *val = (double*)nd->nrrd_->data;
  double *data = matrix->get_data_pointer();

  for(index_type r=0; r<rows; r++) 
  {
    for(index_type c=0; c<cols; c++) 
    {
      *val = *data;
      ++data;
      ++val;
    }
  }
  // send the data nrrd
  NrrdDataHandle dataH(nd);
  send_output_handle("Data", dataH);
}


void
ConvertMatrixToNrrd::create_and_send_sparse_matrix_nrrd(MatrixHandle matH)
{
  SparseRowMatrix* matrix = matrix_cast::as_sparse(matH);

  size_t nnz[NRRD_DIM_MAX];
  nnz[0] = static_cast<size_t>(matrix->get_nnz());
  size_type rows = matrix->nrows();

  // create 3 nrrds (data, rows, cols)
  NrrdData *data_n = new NrrdData();
  nrrdAlloc_nva(data_n->nrrd_, nrrdTypeDouble, 1, nnz);
  nrrdAxisInfoSet_nva(data_n->nrrd_, nrrdAxisInfoLabel, "sparse-data");
  data_n->nrrd_->axis[0].kind = nrrdKindDomain;

  NrrdData *rows_n = new NrrdData();
  size_t sparse_size[NRRD_DIM_MAX];
  sparse_size[0] = rows + 1;
  
  if (sizeof(index_type) == 4)
  {
    nrrdAlloc_nva(rows_n->nrrd_, nrrdTypeInt, 1, sparse_size);
  }
  else
  {
    nrrdAlloc_nva(rows_n->nrrd_, nrrdTypeLLong, 1, sparse_size);
  }
  
  nrrdAxisInfoSet_nva(rows_n->nrrd_, nrrdAxisInfoLabel, "sparse-rows");
  rows_n->nrrd_->axis[0].kind = nrrdKindDomain;

  NrrdData *cols_n = new NrrdData();

  if (sizeof(index_type) == 4)
  {
    nrrdAlloc_nva(cols_n->nrrd_, nrrdTypeInt, 1, nnz);
  }
  else
  {
    nrrdAlloc_nva(cols_n->nrrd_, nrrdTypeLLong, 1, nnz);  
  }
  
  nrrdAxisInfoSet_nva(cols_n->nrrd_, nrrdAxisInfoLabel, "sparse-columns");
  cols_n->nrrd_->axis[0].kind = nrrdKindDomain;

  // pointers to nrrds
  double *data_p = (double*)data_n->nrrd_->data;
  index_type *rows_p = static_cast<index_type* >(rows_n->nrrd_->data);
  index_type *cols_p = static_cast<index_type* >(cols_n->nrrd_->data);

  // points to matrix arrays
  index_type *rr = matrix->get_rows();
  index_type *cc = matrix->get_cols();
  double *d = matrix->get_vals();

  // copy data and cols (size nnz)
  index_type i = 0;
  for (i=0; i<static_cast<size_type>(nnz[0]); i++) 
  {
    data_p[i] = d[i];
    cols_p[i] = cc[i];
  }

  // copy rows (size rows+1)
  for (i=0; i<rows+1; i++) 
  {
    rows_p[i] = rr[i];
  }
  
  // send nrrds
  NrrdDataHandle dataH(data_n);
  NrrdDataHandle rowsH(rows_n);
  NrrdDataHandle colsH(cols_n);
  
  send_output_handle("Data", dataH);
  send_output_handle("Rows", rowsH);
  send_output_handle("Columns", colsH);
}


} // End namespace Teem


