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


// assumes that cols(A) > rows(A)
VISVector solvLinearSystemUnderconstrainted(const VISMatrix& A,const VISVector& b)
{
  VISMatrix V, W, U;
  //(A.t()).svd(U, W, V);
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

VISMatrix vectorToMatrix(const VISVector &vec, int num_rows = -1, int num_cols = -1)
{
  int M, N;
  if (num_rows  > 0)
    {
      M = num_rows;
      N = num_cols;
      if (M*N != vec.n())
	{
	  cout << "ERROR vectorToMatrix rows/cols/size mismatch " << endl;
	  return(VISMatrix());
	}
    }
  else
    {
      M = (int)(sqrt(N = vec.n()));
        if (M*M != N)
	  {
	    cout << "got bad vector input on vector to matrix" << endl;
	    return(VISMatrix());
	  }
    }
      
  VISMatrix ret(M, N);
  int i, j, k = 0;
  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++)
      {
	ret.poke(i, j) = vec.peek(k++);
      }
  return(ret);
}


int main(int argc, char** argv)
{
  int dim = 3;
  int num_pts;
  int num_labels;
  int num_objects;
  int i, j, k, l;
  VISMatrix A, update, A_inv;
  VISVector x, b;
  VISMatrix *points;
  VISVector point_tmp;
  VISMatrix **distances;
  VISMatrix *weights, *weights_tmp;
  VISVector vec_tmp;
  float domain_size = 1.0;
  //  float dt = -1.0/domain_size;
  float dt;
  int seed = 10;
  srand(seed);

  FILE* in = fopen(argv[1], "r" );
  FILE* in_pts;

  fscanf(in, "%d", &num_objects);
  fscanf(in, "%d", &dim);
  fscanf(in, "%d", &num_pts);
  fscanf(in, "%d", &num_labels);
  char filename[80];
  float flt_tmp;

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
      weights[k] = VISMatrix(num_labels, num_pts);
      weights_tmp[k] = VISMatrix(num_labels, num_pts);
      for (l = 0; l <= k; l++)
	distances[k][l] = VISMatrix(num_pts, num_pts);
    }

  point_tmp = VISVector(dim);
  for (k = 0; k < num_objects; k++)
    {
      fscanf(in, "%s", filename);
      //      cout << filename << endl;
      in_pts = fopen(filename, "r");
      for (i = 0; i < num_pts; i++)
	{
	  for (l = 0; l < dim; l++)
	    {
	      fscanf(in_pts, "%f", &flt_tmp);
	      point_tmp.poke(l) = flt_tmp;
	      }
	  (points[k]).pokeROI(0, i, point_tmp);
	}
      fclose(in_pts);
      //      cout << "points " << k << points[k] << endl;
    }
  fclose(in);


  for (k = 0; k < num_objects; k++)
    {
      for (l = 0; l < k; l++)
	{
	  for (i = 0; i < num_pts; i++)
	    for (j = 0; j < num_pts; j++)
	    {
	      // the first index is always has it's points on the row and it's allways the greater of the two.
	      (distances[k][l]).poke(i, j) = ((points[k]).vec(i) - (points[l]).vec(j)).norm();
	      //	      (distances[k][l]).poke(i, j) = power(((points[k]).vec(i) - (points[l]).vec(j)).norm(), 2);
	    }
	  flt_tmp = VISmax((distances[k][l]).max(), (float)flt_tmp);
      }
    }

  for (k = 0; k < num_objects; k++)
    for (i = 0; i < num_pts; i++)
      for (j = 0; j < num_pts; j++)
      {
	// the first index is always has it's points on the row and it's allways the greater of the two.
	(distances[k][k]).poke(i, j) = 1.0f/(((points[k]).vec(i) - (points[k]).vec(j)).norm() + 0.2/flt_tmp);
	//	(distances[k][k]).poke(i, j) = -1.0f*((points[k]).vec(i) - (points[k]).vec(j)).norm();
	//	      (distances[k][l]).poke(i, j) = power(((points[k]).vec(i) - (points[l]).vec(j)).norm(), 2);
      }

  //  dt = -1.0/flt_tmp;
  //  dt = -1.0/flt_tmp;
  dt = -0.2/flt_tmp;

  //  cout << "distances " << distances[1][0] << endl;

  A = VISMatrix(num_labels, num_labels*num_pts);
  A = 0.0f;
  b = VISVector(num_labels);
  b = 1.0f;

  // constraints that rows and columns sum to one...  
  for (i = 0; i < num_labels; i++)
    for (j = 0; j < num_pts; j++)
      {
	A.poke(i, j*num_labels + i) = 1.0f;
      }

  //  cout << "done matrix a" << endl;

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

  //cout << A_inv*A << endl;
  //    cout << V*V.t() << endl;
  //    cout << V.t()*V << endl;
  //    cout << U.t()*U << endl;
  //    cout << W_inv*W << endl;
  //    cout << "difference is " << ((W_inv*(U.t()*U)*W) - VISIdentity(25)).norm() << endl;
  //    cout << (W_inv*(U.t()*U)*W) << endl;
    
  
  //  cout << A << endl;



  //  cout << "about to do rand init" << endl;
  // INITIALIZE weights to be random
  for (k = 0; k < num_objects; k++)
    {
      weights[k] = VISMatrix(num_labels, num_pts);
      for (i = 0; i < num_labels; i++)
	for (j = 0; j < num_pts; j++)
	  (weights[k]).poke(i, j) = 1.0 + 0.1*(rand1() - 0.5);
      // project them onto the constraints
      weights[k] *= 1.0/(float)num_pts;
      //      cout << "about to do constrained solution " << endl;
      weights[k] = vectorToMatrix(projectPosConstrainedSolution3(A, matrixToVector(weights[k]), b, A_inv), 
				  num_labels, num_pts);
      //      cout << "weights " << k << weights[k] << endl;
    }

  //  cout << "done rand init" << endl;

  update = VISMatrix(num_labels, num_pts);
  float energy;
  VISMatrix matrix_tmp;
  float corres_weight = 0.5*(1.0/num_objects),  separation_weight = (1.0/num_labels);

  float start_energy = 0.0f;
  for (k = 0; k < num_objects; k++)
    {
      for (l = 0; l < k; l++)
	start_energy += corres_weight*((weights[k])*distances[k][l]*(weights[l]).t()).trace();
      matrix_tmp = weights[k]*distances[k][l]*weights[k].t();
      for (i = 0; i < num_labels; i++)
	matrix_tmp.poke(i, i) = 0.0f;
      start_energy += separation_weight*matrix_tmp.sum();
    }
  //  start_energy *= (2.0/(num_objects*(num_objects-1)));
  
  cout << "Start energy " << start_energy << endl;
  // pushes weights on the same object away from each other
  VISMatrix repulsion_mask(num_labels, num_labels);
  repulsion_mask = 1.0;
  repulsion_mask -= VISIdentity(num_labels);

      
  //  for (k = 0; k < num_objects; k++)
  //    cout << "weights " << k << weights[k] << endl;

  //  VISMatrix update_correspondence, update_repulsion, updates[num_objects];

  int num_initial_iterations = 200;
  for (int kk = 0; kk < 250; kk++)
    {
      for (k = 0; k < num_objects; k++)
	{
	  update = 0.0f;
	  //	  cout << "weights " << k << weights[k] << endl;
	  for (l = 0; l < num_objects; l++)
	    if (k > l)
	      update += corres_weight*weights[l]*(distances[k][l]).t();
	    else if (k < l)
	      update += corres_weight*weights[l]*distances[l][k];
	    else  if (k == l)
	      update += separation_weight*repulsion_mask*(weights[k]*distances[k][l]);
	  //	  cout << "update " << update << endl;
	  //	  cout << "weights " << k << weights[k]+update << endl;
	  weights_tmp[k] = vectorToMatrix(solvePosConstrainedUpdate3(A, b, matrixToVector(weights[k]), matrixToVector(dt*update), A_inv), 
					  num_labels, num_pts);
	}
      
      for (k = 0; k < num_objects; k++)
	{
	  weights[k] = weights_tmp[k];
	  //	  cout << "weights " << k << weights[k] << endl;
	}

       if (kk > num_initial_iterations)
 	for (k = 0; k < num_objects; k++)
 	  {
 	    weights[k] *= 1.5f;
 	    weights[k] = vectorToMatrix(projectPosConstrainedSolution3(A, matrixToVector(weights[k]), b, A_inv),
 					num_labels, num_pts);
 	  }

      energy = 0.0f;
      float separation_energy = 0.0;
      for (k = 0; k < num_objects; k++)
	{
	  for (l = 0; l < k; l++)
	    energy += corres_weight*((weights[k])*distances[k][l]*(weights[l]).t()).trace();
	  matrix_tmp = weights[k]*distances[k][l]*weights[k].t();
	  for (i = 0; i < num_labels; i++)
	    matrix_tmp.poke(i, i) = 0.0f;
	  separation_energy += separation_weight*matrix_tmp.sum();
	}



      //	for (l = 0; l < k; l++)
      //	  ;
      //	  energy += ((weights[k])*distances[k][l]*(weights[l]).t()).trace();
      //      energy *= (2.0/(num_objects*(num_objects-1)));
      cout << "Iteration " << kk << " energy " << energy << endl;
      cout << "Iteration " << kk << " separation_energy " << separation_energy << endl;
      cout << "Iteration " << kk << " total energy " << separation_energy + energy << endl;
    }

//       for (k = 0; k < num_objects; k++)
// 	for (l = 0; l < k; l++)
// 	  {
// 	    cout << ((weights[k])*distances[k][l]*(weights[l]).t()) << endl;
// 	    cout << "trace " << ((weights[k])*distances[k][l]*(weights[l]).t()).trace() << endl;
// 	  }

  float max_weight, max_index;
  for (k = 0; k < num_objects; k++)
    {
      cout << "object " << k << endl;
      for (j = 0; j < num_labels; j++)
	{
	  max_index = 0; max_weight = weights[k].peek(j, 0);
	  for (i = 1; i < num_pts; i++)
	      if (weights[k].peek(j, i) > max_weight)
		{
		  max_weight = weights[k].peek(j, i);
		  max_index = i;
		}
	  cout << points[k].peek(0, max_index) << " " << points[k].peek(1, max_index) << " "<< points[k].peek(2, max_index) << endl;
	}
      cout << "weights " << k << ": " << weights[k] << endl;
      cout << "weights means " << num_labels*(weights[k].meanOfCols()) << endl;
    }

  //  cout << "correspondences " << (weights[1].t())*weights[0];

  cout << "REGRESSION: " << endl;
  cout << "rand seed: " << seed << endl;
  cout << "num pts: " << num_pts << endl;
  cout << "num objects: " << num_objects << endl;
  cout << "dimension: " << dim << endl;
  cout << "start energy: " << start_energy << endl;
  cout << "final energy: " << energy << endl;

  for (k = 0; k < num_objects; k++)
    {
      for (l = 0; l < k; l++)
	energy += ((weights[k])*distances[k][l]*(weights[l]).t()).trace();
    }

  energy *= (2.0f)/(num_objects*(num_objects -1));
  energy *= (1.0f)/(num_labels);
  cout << "average distance to points" << energy << endl;

}




