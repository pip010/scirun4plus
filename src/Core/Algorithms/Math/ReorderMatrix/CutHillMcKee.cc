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

#include <Core/Algorithms/Math/ReorderMatrix/CutHillMcKee.h>

#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

namespace SCIRunAlgo {

using namespace SCIRun;

bool 
CutHillMcKeeAlgo::run(MatrixHandle  input, 
                      MatrixHandle& output,
                      MatrixHandle& mapping)
{
  algo_start("CutHillMcKee");
  
  bool calcmapping = get_bool("build_mapping");

  index_type *rr, *cc;
  double *d;
  size_type n,m,nnz;

  if (input.get_rep() == 0)
  {
    error("No input matrix was found");
    algo_end(); return (false);
  }
  
  if (input->ncols() != input->nrows())
  {
    error("Matrix is not square");
    algo_end(); return (false);  
  }
  
  if (!matrix_is::sparse(input)) 
  {
    error("Matrix is not sparse");
    algo_end(); return (false);
  }
  
  SparseRowMatrix* sinput = matrix_cast::as_sparse(input);
  
  output = sinput->clone();
  if (output.get_rep() == 0)
  {
    error("Could not copy sparse matrix");
    algo_end(); return (false);  
  }
  
  m  = sinput->nrows();
  n  = sinput->ncols();
  nnz = sinput->get_nnz();
  rr = sinput->get_rows();
  cc = sinput->get_cols();
  d  = sinput->get_vals();
  
  index_type *drr, *dcc;
  double *dd;

  SparseRowMatrix* soutput = matrix_cast::as_sparse(output);  
  drr = soutput->get_rows();
  dcc = soutput->get_cols();
  dd  = soutput->get_vals();

  index_type *degree;
  index_type ns,nr,nq,nx,nss;
  
  // Count number of connections
  degree = drr;
  for (index_type p = 0;p<m;p++) degree[p] = (degree[p+1] - degree[p]);
  
  // clear the mask of already processed nodes

  std::vector<index_type> X(m);
  std::vector<index_type> Q(m);
  std::vector<index_type> S(m);
  std::vector<index_type> R(m);
    
  for (index_type p=0;p<m;p++) Q[p] = 0;
  
  index_type root = 0;
  for (index_type p=0;p<m;p++) if((degree[p] < degree[root])&&(Q[p] == 0)) root = p;

  nss = 0;
  nr = 0;
  ns = 0;
  nq = 0;
  nx = 0;
  
  R[nr++] = root;
  Q[root] = 1; nq++;
  X[nx++] = root;
  
  nr = 0;
  for (index_type p = rr[root];p<rr[root+1];p++) if (Q[cc[p]] == 0) { R[nr++] = cc[p]; Q[cc[p]] = 1; nq++; }
 
  while(true)
  {
    if(nr > 0) 
    { 
      for (index_type i =0;i<nr;i++)
      {
        index_type t;
        index_type k = i;
        // sort entries by order
        for (index_type j=i+1;j<nr;j++) { if (degree[R[k]] > degree[R[j]]) { k = j;}
        if (k > i)  {t = R[i]; R[i] = R[k]; R[k] = t;}  }
      }
      for (index_type j=0;j<nr;j++) { X[nx++] = R[j]; S[ns++] = R[j]; }
      nr  = 0;
    }
    if (nx == m) break;
    
    if (nss < ns)
    {
      root = S[nss++];
      nr = 0;
      for (index_type p = rr[root];p<rr[root+1];p++) if (Q[cc[p]] == 0) { R[nr++] = cc[p]; Q[cc[p]] = 1; nq++; }
    }
    else
    {
      for (index_type p=0;p<m;p++) if(Q[p] == 0) { root = p; break;}
      for (index_type p=0;p<m;p++) if((degree[p] < degree[root])&&(Q[p] == 0)) root = p;
      Q[root] = 1; nq++;
      X[nx++] = root;
      nr = 0;
      for (index_type p = rr[root];p<rr[root+1];p++) if (Q[cc[p]] == 0) { R[nr++] = cc[p]; Q[cc[p]] = 1; nq++; }
    }    
  }

  // Finish mapping matrix and inverse mapping matrix.
  if (calcmapping)
  {
    SparseRowMatrix::Data outputData(m+1, m);
    if (!outputData.allocated())
    {
      error("Could not reserve space for mapping matrix");    
      algo_end(); return (false);
    }    
    const SparseRowMatrix::Rows& mrr = outputData.rows();
    const SparseRowMatrix::Columns& mcc = outputData.columns();
    const SparseRowMatrix::Storage& md = outputData.data();
    
    for (index_type p = 0; p < m+1; p++) 
      mrr[p] = p; 
    for (index_type p = 0; p < m; p++) 
      md[p] = 1.0;
    for (index_type p = 0; p < m; p++) 
      mcc[X[p]] = p;
    mapping = new SparseRowMatrix(m,m,outputData,m);
  }
    
   for (index_type p=0;p<m;p++) { Q[X[p]] = p; }
  
  // finish the reordering new matrix
  index_type i,j;
  drr[0] = 0;
  i = 0;
  for (index_type p=0;p<m;p++)
  {
    j = Q[p];
    drr[p+1] = drr[p]+(rr[j+1]-rr[j]);
    for (index_type r=rr[j];r<rr[j+1];r++)
    {
      dcc[i] = X[cc[r]];
      dd[i] = d[r];
      i++;
    }
  }

  algo_end(); return (true);
}

} // end namespace SCIRunAlgo
