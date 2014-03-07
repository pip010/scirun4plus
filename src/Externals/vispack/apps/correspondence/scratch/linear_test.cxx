
#include <matrix.h>
#include <array.h>
#include <mathutil.h>


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
VISVector solvLinearSystemUnderconstrainted(const VISMatrix& A,const VISVector& b)
{
  VISMatrix V, W, U;
  (A.t()).svd(U, W, V);
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

// // // solve a underconstrained/least-squares linear system with nonnegativity constraint...
// VISVector solvePosConstrainedUpdate(const VISMatrix &A, const VISVector  &solution, const VISVector &update)
// {
//   int iteration = 0;
//   int N = solution.n();
//   VISMatrix A_local;
//   int i, j, k;
//   VISVector adjustment, new_solution, old_solution;
//   VISVector new_constraint(N);
//   VISArray< boolean > constrained(N);
//   float remaining_update = 1.0f, update_amount;
//   boolean done = false;
  
//   // check for negative  and zero entries
//   for (i = 0; i < N; i++)
//     {
//       if (solution[i] < 0.0f)
// 	{
// 	  cout << "Warning: got negative input on solvePosConstrainedUpdate" << endl;
// 	  return(VISVector());
// 	}
//       constrained.poke(i) = false;
//     }

//   old_solution = solution;

//   A_local = A;
//   while (!done)
//     {
//       // project update onto the constraints
//       adjustment = solveLinearEqn(A_local, (-1.0f*remaining_update)*(A_local*update));
//       this_update = remaining_update + adjustment;
      
//       new_solution = old_solution + this_update;
      
//       min_index = 0;
//       min_solution = new_solution[min_index];
//       for (i = 0; i < N; i++)
// 	{
// 	  if (new_solution[i] < min_solution)
// 	    {
// 	      min_index = i;
// 	      min_solution = new_solution[i];
// 	    }
// 	}

//       update_amount = 1.0f;
//       if (min_solution < 0.0f)
// 	{
// 	  update_amount = (this_update[min_index] - min_solution)/(this_update[min_index]);
// 	  this_update *= update_amount;
// 	  new_solution = old_solution += this_update;
// 	}
//       // this is the amount you have left to update
//       remaining_update *= (1.0f - update_amount);
//       old_solution = new_solution;
//       if (fabs(remaining_update) < EPSILON) done = true;
//     }
// }

#define EPSILON (1.0e-6)

VISVector solvePosConstrainedUpdate(const VISMatrix &A, const VISVector  &b, const VISVector  &solution, const VISVector &update)
{
  int iteration = 0;
  int N = solution.n();
  VISMatrix A_local;
  VISVector b_local;
  int i, j, k;
  int min_index;
  float min_solution;
  VISVector adjustment, new_solution, first_solution, this_update, first_update, current_solution;
  VISVector new_constraint(N);
  VISArray <boolean> constrained(N);
  float error;
  VISVector vec_short(1);
  vec_short = 0.0;
  boolean done = true;


  A_local = A;
  b_local = b;


  //  adjustment = solvLinearEqn(A_local, -1.0f*(A_local*update));
  adjustment = solvLinearSystemUnderconstrainted(A_local, -1.0f*(A_local*update));
  // project the update onto the tangent space of constraints
  first_update = update + adjustment;
  //  cout << "raw update " << update << endl;
  //  cout << "projected update " << first_update << endl;
  // find a solution that satifies the hard constraints but is, perhaps negative
  first_solution = solution + first_update;

  

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

  new_solution = first_solution;

  for (i = 0; i < N; i++)
    {
      if (new_solution[i] < -EPSILON)
	{
	  new_constraint = 0.0f;
	  new_constraint.poke(i) = 1.0f;
	  A_local = A_local.concatRow(new_constraint.t());
	  //	  vec_short = -new_solution[i];
	  vec_short = 0.0f;
	  b_local = b_local.concatRow(vec_short);
	  //	  cout << "added constraint at " << i << endl;
	  constrained.poke(i) = true;
	  done = false;
	}
    }

  while (!done)
    {
      //      this_update = solvLinearEqn(A_local, b_local - A_local*first_solution);
      this_update = solvLinearSystemUnderconstrainted(A_local, b_local - A_local*first_solution);


      //      cout << "new_solution " << new_solution << endl;

      //      cout << "this_update " << this_update << endl;

      new_solution = first_solution + this_update;

      //      error = (A*new_solution - b).norm();
      //      cout << "error is " << error << endl;
      //
      //      error = (A_local*new_solution - b_local).norm();
      //      cout << "real error is " << error << endl;
      //      cout << "new_solution " << new_solution << endl;

      done = true;

      for (i = 0; i < N; i++)
	{
	  if (new_solution[i] < 0.0f)
	    if (!constrained[i])
	      {
		new_constraint = 0.0f;
		new_constraint.poke(i) = 1.0f;
		A_local = A_local.concatRow(new_constraint.t());
		//	  vec_short = -new_solution[i];
		vec_short = 0.0f;
		b_local = b_local.concatRow(vec_short);
		//		cout << "added constraint at " << i << endl;
		constrained.poke(i) = true;
		done = false;
	      }
	}
    }

  for (i = 0; i < N; i++)
    if (constrained[i])
      new_solution.poke(i) = 0.0f;

  //  for (i = 0; i < N; i++)
  //    if (fabs(new_solution[i]) < EPSILON)
  //      new_solution.poke(i) = 0.0f;

  return(new_solution);
}



VISVector projectPosConstrainedSolution2(const VISMatrix &A, 
					 const VISVector  &x, 
					 const VISVector  &b)
{
  int iteration = 0;
  int N = x.n();
  VISMatrix A_local;
  VISVector b_local;
  int i, j, k;
  int min_index;
  float min_solution;
  VISVector adjustment, new_solution, old_solution, this_update, first_update, current_solution;
  VISVector new_constraint(N);
  VISArray <boolean> constrained(N);
  float error;
  VISVector vec_short(1);
  vec_short = 0.0;
  boolean done = true;

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
  b_local = b;

  first_update = solvLinearSystemUnderconstrainted(A, b - A*old_solution);
  //  cout << "first update with under c is " << first_update << endl;
  //  first_update = solvLinearEqn(A, b - A*old_solution);
  //  cout << "first update with svd c is " << first_update << endl;
  //  exit(-1);

  current_solution = old_solution + first_update;

  new_solution = current_solution;

  for (i = 0; i < N; i++)
    {
      if (new_solution[i] < 0.0)
	{
	  new_constraint = 0.0f;
	  new_constraint.poke(i) = 1.0f;
	  A_local = A_local.concatRow(new_constraint.t());
	  //	  vec_short = -new_solution[i];
	  vec_short = 0.0f;
	  b_local = b_local.concatRow(vec_short);
	  //	  cout << "added constraint at " << i << endl;
	  constrained.poke(i) = true;
	  current_solution.poke(i) = 0.0f;
	  done = false;
	}
    }

  while (!done)
    {
      //      this_update = solvLinearEqn(A_local, b_local - A_local*current_solution);
      this_update = solvLinearSystemUnderconstrainted(A_local, b_local - A_local*current_solution);

      //      cout << "new_solution " << new_solution << endl;
      //      cout << "this_update " << this_update << endl;

      new_solution = current_solution + this_update;

      //      error = (A*new_solution - b).norm();
      //      cout << "error is " << error << endl;

      //      error = (A_local*new_solution - b_local).norm();
      //      cout << "real error is " << error << endl;

      //      cout << "new_solution " << new_solution << endl;

      done = true;

      for (i = 0; i < N; i++)
	{
	  if (new_solution[i] < 0.0)
	    if (!constrained[i])
	      {
		new_constraint = 0.0f;
		new_constraint.poke(i) = 1.0f;
		A_local = A_local.concatRow(new_constraint.t());
		//	  vec_short = -new_solution[i];
		vec_short = 0.0f;
		b_local = b_local.concatRow(vec_short);
		//		cout << "added constraint at " << i << endl;
		constrained.poke(i) = true;
		current_solution.poke(i) = 0.0f;
		done = false;
	      }
	}
      //      cout << new_solution << endl;
    }


  for (i = 0; i < N; i++)
    if (constrained[i])
      new_solution.poke(i) = 0.0f;

  //  for (i = 0; i < N; i++)
  //    if (fabs(new_solution[i]) < EPSILON)
  //      new_solution.poke(i) = 0.0f;
  ///  cout << new_solution << endl;
  return(new_solution);
}


VISVector solvePosConstrainedUpdate3(const VISMatrix &A, const VISVector  &b, const VISVector  &solution, const VISVector &update, const VISMatrix &A_inv = VISMatrix())
{
  int iteration = 0;
  int N = solution.n();
  VISMatrix A_local;
  VISVector b_local;
  int i, j, k;
  VISVector adjustment, new_solution, first_solution, this_update, first_update, current_solution;
  VISArray <boolean> constrained(N);
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


  for (i = 0, this_item = 0; i < N; i++)
    {
      if (new_solution[i] < 0.0)
	{
	  A_local = removeColumn(A_local, this_item);
	  current_solution = removeItem(current_solution, this_item);
	  constrained.poke(i) = true;
	  num_constraints++;
	  //	  cout << "adding constraint " << i << endl;
	  done = false;
	}
      else
	this_item++;
    }
  cout << "got " << num_constraints << " constraints after first update" << endl;

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

      for (i = 0, this_item = 0, other_item = 0; i < N; i++)
	{
	  if (!constrained[i])
	    {
	      if (new_solution[other_item] < 0.0)
		{
		  A_local = removeColumn(A_local, this_item);
		  current_solution = removeItem(current_solution, this_item);
		  constrained.poke(i) = true;
		  num_constraints++;
		  done = false;
		  //		  cout << "adding constraint " << i << endl;
		}
	      else
		this_item++;
	      other_item++;
	    }
	}
      cout << "num contraints " << num_constraints << endl;
      //      cout << new_solution << endl;
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
  VISArray <boolean> constrained(N);
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
    first_update = A_inv*(b - A*old_solution);
  else
    first_update = solvLinearSystemUnderconstrainted(A, b - A*old_solution);
  //  cout << "first update with under c is " << first_update << endl;
  //  first_update = solvLinearEqn(A, b - A*old_solution);
  //  cout << "first update with svd c is " << first_update << endl;
  //  exit(-1);

  current_solution = old_solution + first_update;

  new_solution = current_solution;

  //  cout << "current_solution " << current_solution << endl;
  for (i = 0, this_item = 0; i < N; i++)
    {
      if (new_solution[i] < 0.0)
	{
	  A_local = removeColumn(A_local, this_item);
	  current_solution = removeItem(current_solution, this_item);
	  constrained.poke(i) = true;
	  num_constraints++;
	  //	  cout << "adding constraint " << i << endl;
	  done = false;
	}
      else
	this_item++;
    }

  cout << "added " << num_constraints << " constraints after first update" << endl;

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

      // this is complicated, since you have 3 difference size arrays floating around....
      //      num_constraints = 0;
      for (i = 0, this_item = 0, other_item = 0; i < N; i++)
	{
	    {
	      if (new_solution[other_item] < 0.0)
		{
		  A_local = removeColumn(A_local, this_item);
		  current_solution = removeItem(current_solution, this_item);
		  constrained.poke(i) = true;
		  done = false;
		  //		  cout << "adding constraint " << i << endl;
		  num_constraints++;
		}
	      else
		this_item++;
	      other_item++;
	    }
	}
      cout << "num constraints " << num_constraints << endl;
      //      cout << new_solution << endl;
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
  int num_pts = 5;
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


  // INITIALIZE weights to be random
  for (k = 0; k < num_objects; k++)
    {
      weights[k] = VISMatrix(num_pts, num_pts);
      for (i = 0; i < num_pts; i++)
	for (j = 0; j < num_pts; j++)
	  (weights[k]).poke(i, j) = 1.0 + 0.01*(rand1() - 0.5);
      // project them onto the constraints
      weights[k] *= 1.0/(float)num_pts;
      weights[k] = vectorToMatrix(projectPosConstrainedSolution3(A, matrixToVector(weights[k]), b, A_inv));
      cout << "weights " << k << weights[k] << endl;
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
  (A.t()).svd(U, W, V);
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
  A_inv = U*W*(V.t());
  A_inv = A.inverseSVD();
  
  cout << A << endl;
  cout << A*A_inv << endl;


  cout << "distances " << distances[1][0] << endl;

  update = VISMatrix(num_pts, num_pts);

  int num_initial_iterations = 40;
  float energy;
  for (int kk = 0; kk < 80; kk++)
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
	    weights[k] *= 1.3f;
	    weights[k] = vectorToMatrix(projectPosConstrainedSolution3(A, matrixToVector(weights[k]), b, A_inv));
	  }
      energy = 0.0f;
      for (k = 0; k < num_objects; k++)
	for (l = 0; l < k; l++)
	  energy += ((weights[k])*distances[k][l]*(weights[l]).t()).trace();
      cout << "Energy " << energy*(2.0/(num_objects*(num_objects-1))) << endl;
    }

  for (k = 0; k < num_objects; k++)
    cout << "weights " << k << weights[k] << endl;

}




