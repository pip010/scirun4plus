#include <lapack.h>
#include <matrix.h>
#include <array.h>
#include <mathutil.h>


void svdTallMatrix(const VISMatrix &A, VISMatrix &U, VISMatrix &W, VISMatrix &V)
{
  //A=U*W*transpose(V)
  //The following stuff are needed by SVD of LAPACK//
  int cols = A.c();
  int rows = A.r();
  integer  M = rows;
  integer  N = cols;
  integer S = VISmax(M, N);
  char  ju  = 'S';   //'All' columns of U are returned
  char  jvt = 'A';   //'All' columns of VT are returned
  integer  la = M;   //leading dimension of A: >=max(1,M)
  integer  lu = M;   //leading dimension of U: >= M
  integer  lvt = N;  //leading dimension of VT: >= N
  integer  lwk;      //dimension of the array work
  integer  ret;      //return value: 0=successful exit
  datatype *a, *work;     //input array and work space
  datatype *u, *w, *vt;   //output arrays
  // from the lapack documentation...
  lwk =VISmax(3*VISmin(M, N) + VISmax(M, N), 5*VISmin(M, N));
  work = new datatype [lwk];   
  a = A.createLAPK();
  w = new datatype [N]; 
  //    w = new datatype [S]; // this should be N but lapack seems to be broken
  vt = new datatype [N*N];
  //vt = new datatype [S*S]; // this should be NxN but lapack seems to be broken
  //  u = new datatype [S*S];
  //u = new datatype [N*N];
  u = new datatype [M*N];

  //printLAPK(a, M, N);
#ifdef VISMATRIXDOUBLE


  //    cout << "about to do svd" << endl;
  dgesvd_(&ju, &jvt, &M, &N, a, &la, w, u, &lu, vt, &lvt, work, &lwk, &ret);
#else
  sgesvd_(&ju, &jvt, &M, &N, a, &la, w, u, &lu, vt, &lvt, work, &lwk, &ret);
#endif
  //    cout << "done svd" << endl;

  VISMatrix um(M, N);
  um = 0.0f;

  um.becomeLAPK(u, M, N);  //u from lapack
  U = um;  //U is to be returned
  //printLAPK(u, M, N);
  VISMatrix vm(N, N);
  vm = 0.0f;

  vm.becomeLAPK(vt, N, N);  //NOTE: vt here!

  V = vm.t();  //V is to be returned

  //printLAPK(vt, N, N);
  VISMatrix wm(N, N);
  int j;
  wm = 0.0f;
  for (j=0; j<VISmin(N, M); j++)
    wm.poke(j, j) = w[j];  //w from lapack
  for (; j<cols; j++)
    wm.poke(j, j) = 0.0f;

  W = wm;  //W is to be returned
  //printLAPK(w, M, N);

  destroyLAPK(u);  destroyLAPK(vt);
  destroyLAPK(w);  destroyLAPK(work);
  destroyLAPK(a);

  if (ret == 0){
    //	cout << "SVD called successfully." << endl << endl;
  }
  else{
    if(ret < 0)
      cout << "Illegal value for " << -(ret) << "th argument!" << endl;
    else
      {
        cout << "Error: " << ret << " bidiagonals not zero!" << endl;
	cout << "input is " << A << endl;
      }
    //      exit(0);
    U = V = W = VISMatrix();
  }
}//svd



VISMatrix multiply(const VISMatrix& mat1, const VISMatrix& mat2)
{
  int J = mat1.r(), K = mat1.c(), L = mat2.c();
  int size = J*L;
  if (!(K == mat2.r()))
    return(VISMatrix());
  VISMatrix ret(J, L);
  datatype *ret_buf; 
  int j, k, l;
  const datatype *mat1_buf, *mat2_buf;
  
  ret_buf = ret.dataBufRef();
  for (j = 0; j < size; j++)
    *(ret_buf++) = 0.0f;

  ret_buf = ret.dataBufRef();
  mat1_buf = mat1.dataBuf();
  mat2_buf = mat2.dataBuf();
  for (l = 0; l < L; l++)
    {
      mat1_buf = mat1.dataBuf();
      for (k = 0; k < K; k++)      
	{
	  ret_buf = ret.dataBufRef() + l*J;
	  for (j = 0; j < J; j++)      
	    *(ret_buf++) += *(mat1_buf++)*(*mat2_buf);
	}
      mat2_buf++;
    }
  return(ret);
}

VISVector removeItem(const VISVector& vec, int item)
{
  int i, j, N = vec.n() - 1;
  VISVector ret(N);
  for (i = j = 0; i < N; i++, j++)
    {
      if (i == item) j++;
      ret.poke(i) = vec.peek(j);
    }
  return(ret);
}

VISVector removeItems(const VISVector& vec, const VISArray< boolean > &to_remove, int num_gone)
{
  int i, j, N = vec.n();
  VISVector ret(N - num_gone);
  for (i = j = 0; j < N; j++)
    {
      if (!to_remove[j]) 
	ret.poke(i++) = vec.peek(j);
    }
  return(ret);
}

VISMatrix removeColumns(const VISMatrix& mat, const VISArray< boolean > &to_remove, int num_gone)
{
  int i, j, k, M = mat.r(), N = mat.c();
  VISMatrix ret(M, N - num_gone);
  for (i = 0, k = 0; k < N; k++)
    {
      if (!to_remove[k]) 
	{
	  for (j = 0; j < M; j++)
	    {
	      ret.poke(j, i) = mat.peek(j, k);
	    }
	  i++;  
	}
    }
  return(ret);
}

VISMatrix removeColumn(const VISMatrix& mat, int col)
{
  int i, j, k, M = mat.r(), N = mat.c() -1;
  VISMatrix ret(M, N);
  for (i = 0, k = 0; i < N; i++, k++)
    {
      if (i == col) k++;
      for (j = 0; j < M; j++)
	{
	  ret.poke(j, i) = mat.peek(j, k);
	}
    }
  return(ret);
}

// assumes that cols(A) > rows(A)
VISVector solvLinearSystemUnderconstrainted(const VISMatrix& A, const VISVector& b)
{
  VISMatrix V, W, U;
  //  (A.t()).svd(U, W, V);
  svdTallMatrix(A.t(), U, W, V);
  int i;
  float value;
  int N = W.r();
  for (i = 0; i < N; i++)
    {
      value = W.peek(i, i);
      if (fabs(value) > 1.0e-10)
	W.poke(i, i) = 1.0f/value;
      else
	W.poke(i, i) = 0.0f;
    }
  return(U*W*(V.t())*b);
}


#define EPSILON (1.0e-6)

VISVector solvePosConstrainedUpdate3(const VISMatrix &A, const VISVector  &b, const VISVector  &solution, const VISVector &update, const VISMatrix &A_inv = VISMatrix())
{
  int iteration = 0;
  int N = solution.n();
  VISMatrix A_local;
  VISVector b_local;
  int i, j, k;
  VISVector adjustment, new_solution, first_solution, this_update, first_update, current_solution;
  VISArray <boolean> constrained(N), to_remove(N);
  float error;
  boolean done = true;
  int this_item, other_item;
  int num_constraints = 0;

  A_local = A;

  //  adjustment = solvLinearEqn(A_local, -1.0f*(A_local*update));
  if (A_inv.isValid())
    adjustment = -1.0*A_inv*(A_local*update);
  else
    adjustment = solvLinearSystemUnderconstrainted(A_local, -1.0f*(A_local*update));
  // project the update onto the tangent space of constraints
  first_update = update + adjustment;
  //  cout << "raw update " << update << endl;
  //  cout << "projected update " << first_update << endl;
  // find a solution that satifies the hard constraints but is, perhaps negative
  current_solution = solution + first_update;


  // check for negative  and zero entries
  for (i = 0; i < N; i++)
    {
      if (solution[i] < 0.0f)
	{
	  cout << "Warning: got negative input on solvePosConstrainedUpdate" << endl;
	  return(VISVector());
	}
      constrained.poke(i) = false;
    }

  //  cout << " error " << (A*first_solution - b).norm() << endl;
  //  VISMatrix U, V, W;

  new_solution = current_solution;

  //  cout << "new solution " << new_solution << endl;
  int num_removed = 0;
  for (i = 0, this_item = 0; i < N; i++)
    {
      if (new_solution[i] < 0.0)
	{
	  to_remove.poke(this_item) = true;
	  constrained.poke(i) = true;
	  num_constraints++;
	  num_removed++;
	    //	  cout << "adding constraint " << i << endl;
	  done = false;
	}
      else
	to_remove.poke(i) = false;
      this_item++;
    }
  current_solution = removeItems(current_solution, to_remove, num_removed);
  //  cout << "sol after" << current_solution << endl;
  //  cout << "A_local before" << A_local << endl;
  A_local = removeColumns(A_local, to_remove, num_removed);
  //  cout << "A_local after" << A_local << endl;

  //  cout << "got " << num_constraints << " constraints after first update" << endl;

  while (!done)
    {
      //      this_update = solvLinearEqn(A_local, b_local - A_local*first_solution);
      this_update = solvLinearSystemUnderconstrainted(A_local, b - A_local*current_solution);

      //      cout << "new_solution " << new_solution << endl;
      //      cout << "this_update " << this_update << endl;

      new_solution = current_solution + this_update;

      //      error = (A*new_solution - b).norm();
      //      cout << "error is " << error << endl;
      //
      //      error = (A_local*new_solution - b_local).norm();
      //      cout << "real error is " << error << endl;
      //      cout << "new_solution " << new_solution << endl;

      done = true;
      int num_removed = 0;
      for (i = 0, this_item = 0; i < N; i++)
	{
	  if (!constrained[i])
	    {
	      if (new_solution[this_item] < 0.0)
		{
		  to_remove.poke(this_item) = true;
		  constrained.poke(i) = true;
		  num_constraints++;
		  num_removed++;
		    //	  cout << "adding constraint " << i << endl;
		    done = false;
		}
	      else
		to_remove.poke(this_item) = false;
	      this_item++;
	    }
	}
      A_local = removeColumns(A_local, to_remove, num_removed);
      current_solution = removeItems(current_solution, to_remove, num_removed);

      //      cout << "got " << num_constraints << " constraints " << endl;
    }

  current_solution = VISVector(N);
  for (i = 0, this_item = 0; i < N; i++)
    if (!constrained[i])
      {
	current_solution.poke(i) = new_solution.peek(this_item);
	this_item++;
      }
    else
      current_solution.poke(i) = 0.0f;

  return(current_solution);
}




VISVector projectPosConstrainedSolution3(const VISMatrix &A, const VISVector  &x, const VISVector  &b, const VISMatrix &A_inv = VISMatrix())
{
  int iteration = 0;
  int N = x.n();
  VISMatrix A_local;
  int i, j, k;
  VISVector adjustment, new_solution, old_solution, this_update, first_update, current_solution;
  VISArray <boolean> constrained(N), to_remove(N);
  float error;
  boolean done = true;
  int this_item, other_item;
  int num_constraints = 0;

  old_solution = x;
  
  // check for negative  and zero entries
  for (i = 0; i < N; i++)
    {
      if (old_solution[i] < 0.0f)
	{
	  cout << "Warning: got negative input on solvePosConstrainedUpdate" << endl;
	  return(VISVector());
	}
      constrained.poke(i) = false;
    }

  //  cout << " error " << (A*old_solution - b).norm() << endl;
  //  VISMatrix U, V, W;

  A_local = A;

  if (A_inv.isValid())
    {
      //      cout << "about to do old multiply " << endl;
      first_update = A_inv*(b - A*old_solution);
      //      cout << "about to do new multiply " << endl;
      //      first_update = multiply(A_inv, (b - multiply(A, old_solution)));
      //      cout << "done" << endl;
    }
  else
    first_update = solvLinearSystemUnderconstrainted(A, b - A*old_solution);
  //  cout << "first update with under c is " << first_update << endl;
  //  first_update = solvLinearEqn(A, b - A*old_solution);
  //  cout << "first update with svd c is " << first_update << endl;
  //  exit(-1);

  current_solution = old_solution + first_update;

  new_solution = current_solution;

  //  cout << "new solution " << new_solution << endl;
  int num_removed = 0;
  for (i = 0, this_item = 0; i < N; i++)
    {
      if (new_solution[i] < 0.0)
	{
	  to_remove.poke(this_item) = true;
	  constrained.poke(i) = true;
	  num_constraints++;
	  num_removed++;
	  //	  cout << "adding constraint " << i << endl;
	  done = false;
	}
      else
	to_remove.poke(this_item) = false;
      this_item++;
    }
  //  cout << "num_removed " << num_removed << endl;
  //  cout << "num constraints " << num_removed << endl;
  //  cout << "to_remove " << to_remove << endl;
  current_solution = removeItems(current_solution, to_remove, num_removed);
  //  cout << "new solution " << current_solution << endl;
  //  cout << "A_local before" << A_local << endl;
  A_local = removeColumns(A_local, to_remove, num_removed);
  //  cout << "A_local after" << A_local << endl;

  //  cout << "got " << num_constraints << " constraints after first update" << endl;

  while (!done)
    {
      //      cout << "current_solution " << current_solution << endl;
      //      cout << "A_local " << A_local << endl;
      //      this_update = solvLinearEqn(A_local, b_local - A_local*current_solution);
      this_update = solvLinearSystemUnderconstrainted(A_local, b - A_local*current_solution);

      //      cout << "new_solution " << new_solution << endl;
      //      cout << "this_update " << this_update << endl;

      new_solution = current_solution + this_update;

      //      error = (A*new_solution - b).norm();
      //      cout << "error is " << error << endl;

      //      error = (A_local*new_solution - b_local).norm();
      //      cout << "real error is " << error << endl;

      //      cout << "new_solution " << new_solution << endl;

      done = true;
      num_removed = 0;
      for (i = 0, this_item = 0; i < N; i++)
	{
	  if (!constrained[i])
	    {
	      if (new_solution[this_item] < 0.0)
		{
		  to_remove.poke(this_item) = true;
		  constrained.poke(i) = true;
		  num_constraints++;
		  num_removed++;
		  //		  cout << "adding constraint " << i << " : " << this_item << endl;
		  done = false;
		}
	      else
		to_remove.poke(this_item) = false;
	      this_item++;
	    }
	}
      //      cout << "num_removed later " << num_removed << endl;
      A_local = removeColumns(A_local, to_remove, num_removed);
      current_solution = removeItems(current_solution, to_remove, num_removed);
      //      cout << new_solution << endl;
      //      cout << "got num constraints " << num_constraints << endl;
    }

  current_solution = VISVector(N);
  for (i = 0, this_item = 0; i < N; i++)
    if (!constrained[i])
      {
	current_solution.poke(i) = new_solution.peek(this_item);
	this_item++;
      }
    else
      current_solution.poke(i) = 0.0f;

  //  for (i = 0; i < N; i++)
  //    if (fabs(new_solution[i]) < EPSILON)
  //      new_solution.poke(i) = 0.0f;
  ///  cout << new_solution << endl;
  return(current_solution);
}




VISVector matrixToVector(const VISMatrix &mat)
{
  VISVector ret(mat.c()*mat.r());
  int i, j, k = 0;
  for (i = 0; i < mat.c(); i++)
    for (j = 0; j < mat.r(); j++)
    {
      ret.poke(k++) = mat.peek(j, i);
    }
  return(ret);
}

VISMatrix vectorToMatrix(const VISVector &vec)
{
  int M, N;
  M = (int)(sqrt(N = vec.n()));
  if (M*M != N)
    {
      cout << "got bad vector input on vector to matrix" << endl;
      return(VISMatrix());
    }
      
  VISMatrix ret(M, M);
  int i, j, k = 0;
  for (i = 0; i < M; i++)
    for (j = 0; j < M; j++)
    {
      ret.poke(j, i) = vec.peek(k++);
    }
  return(ret);
}


int main(int argc, char** argv)
{
  int dim = 1;
  int num_pts = 10;
  int num_objects = 2;
  int i, j, k, l;
  VISMatrix A(2*num_pts, num_pts*num_pts), update, A_inv;
  VISVector x, b(2*num_pts);
  b = 1.0f;
  A = 0.0f;
  VISMatrix *points;
  VISVector point_tmp(dim);
  VISMatrix **distances;
  VISMatrix *weights, *weights_tmp;
  VISVector vec_tmp;
  float domain_size = 1.0;
  float dt = -1.0/domain_size;
  int seed = 1;
  srand(seed);


  // test multiply

  // declare some memory for arrays...
  points = new VISMatrix[num_objects];
  distances = new VISMatrix*[num_objects];
  for (k = 0; k < num_objects; k++)
    distances[k] = new VISMatrix[num_objects];
  weights = new VISMatrix[num_objects];
  weights_tmp = new VISMatrix[num_objects];

  for (k = 0; k < num_objects; k++)
    {
      points[k] = VISMatrix(dim, num_pts);
      weights[k] = VISMatrix(num_pts, num_pts);
      weights_tmp[k] = VISMatrix(num_pts, num_pts);
      for (l = 0; l < k; l++)
	distances[k][l] = VISMatrix(num_pts, num_pts);
    }

  for (k = 0; k < num_objects; k++)
    {
      for (i = 0; i < num_pts; i++)
	{
	  for (l = 0; l < dim; l++)
	    point_tmp.poke(l) = domain_size*rand1();
	  (points[k]).pokeROI(0, i, point_tmp);
	}
      //      cout << "points " << k << points[k] << endl;
    }

  for (k = 0; k < num_objects; k++)
    for (l = 0; l < k; l++)
      {
	for (i = 0; i < num_pts; i++)
	  for (j = 0; j < num_pts; j++)
	    {
	      // the first index is always has it's points on the row and it's allways the greater of the two.
	      (distances[k][l]).poke(i, j) = ((points[k]).vec(i) - (points[l]).vec(j)).norm();
	      //	      (distances[k][l]).poke(i, j) = power(((points[k]).vec(i) - (points[l]).vec(j)).norm(), 2);
	    }
      }


  // constraints that rows and columns sum to one...  
  for (i = 0; i < num_pts; i++)
    for (j = 0; j < num_pts; j++)
      {
	A.poke(i, j + i*num_pts) = 1.0f;
      }
  for (i = 0; i < num_pts; i++)
    for (j = 0; j < num_pts; j++)
      {
	A.poke(i + num_pts, i + j*num_pts) = 1.0f;
      }

  VISMatrix V, W, U;
  //  (A.t()).svd(U, W, V);
  svdTallMatrix(A.t(), U, W, V);
  //  cout << U*W*V.t() << endl;
  float value;
  int N = W.r();
  //  cout << W << endl;
  for (i = 0; i < N; i++)
    {
      value = W.peek(i, i);
      if (fabs(value) > 1.0e-10)
	W.poke(i, i) = 1.0f/value;
      else
	W.poke(i, i) = 0.0f;
    }
  A_inv = U*W*(V.t());
  //  A_inv = V*W_inv*(U.t());
  //  A_inv = A.inverseSVD();
  //  cout << A*A_inv << endl;

  //    cout << A_inv*A << endl;
  //    cout << V*V.t() << endl;
  //    cout << V.t()*V << endl;
  //    cout << U.t()*U << endl;
  //    cout << W_inv*W << endl;
  //    cout << "difference is " << ((W_inv*(U.t()*U)*W) - VISIdentity(25)).norm() << endl;
  //    cout << (W_inv*(U.t()*U)*W) << endl;
    
  
  //  cout << A << endl;

  //  cout << "distances " << distances[1][0] << endl;


  // INITIALIZE weights to be random
  for (k = 0; k < num_objects; k++)
    {
      weights[k] = VISMatrix(num_pts, num_pts);
      for (i = 0; i < num_pts; i++)
	for (j = 0; j < num_pts; j++)
	  (weights[k]).poke(i, j) = 1.0 + 0.01*(rand1() - 0.5);
      // project them onto the constraints
      weights[k] *= 1.0/(float)num_pts;
      //      cout << "about to do constrained solution " << endl;
      weights[k] = vectorToMatrix(projectPosConstrainedSolution3(A, matrixToVector(weights[k]), b, A_inv));
      //      cout << "weights " << k << weights[k] << endl;
    }

  update = VISMatrix(num_pts, num_pts);
  int num_initial_iterations = 30;
  float energy;

  float start_energy = 0.0f;
  for (k = 0; k < num_objects; k++)
    for (l = 0; l < k; l++)
      start_energy += ((weights[k])*distances[k][l]*(weights[l]).t()).trace();
  start_energy *= (2.0/(num_objects*(num_objects-1)));
  cout << "Start energy " << start_energy << endl;

  for (int kk = 0; kk < 60; kk++)
    {
      for (k = 0; k < num_objects; k++)
	{
	  update = 0.0f;
	  for (l = 0; l < num_objects; l++)
	    if (k != l)
	      if (k > l)
		update += (dt/num_objects)*weights[l]*(distances[k][l]).t();
	      else 
		update += (dt/num_objects)*weights[l]*distances[l][k];
	  weights_tmp[k] = vectorToMatrix(solvePosConstrainedUpdate3(A, b, matrixToVector(weights[k]), matrixToVector(update), A_inv));
	}

      for (k = 0; k < num_objects; k++)
	{
	  weights[k] = weights_tmp[k];
	  //	  cout << "weights " << k << weights.peek(k) << endl;
	}

      if (kk > num_initial_iterations)
	for (k = 0; k < num_objects; k++)
	  {
	    weights[k] *= 1.5f;
	    weights[k] = vectorToMatrix(projectPosConstrainedSolution3(A, matrixToVector(weights[k]), b, A_inv));
	  }
      energy = 0.0f;
      for (k = 0; k < num_objects; k++)
	for (l = 0; l < k; l++)
	  energy += ((weights[k])*distances[k][l]*(weights[l]).t()).trace();
      energy *= (2.0/(num_objects*(num_objects-1)));
      cout << "Iteration " << kk << " energy " << energy << endl;

    }

  for (k = 0; k < num_objects; k++)
    cout << "weights " << k << ": " << weights[k] << endl;

  cout << "correspondences " << (weights[1].t())*weights[0];

  cout << "REGRESSION: " << endl;
  cout << "rand seed: " << seed << endl;
  cout << "num pts: " << num_pts << endl;
  cout << "num objects: " << num_objects << endl;
  cout << "dimension: " << dim << endl;
  cout << "start energy: " << start_energy << endl;
  cout << "final energy: " << energy << endl;
}




