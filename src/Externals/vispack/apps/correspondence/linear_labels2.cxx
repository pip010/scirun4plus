#include <lapack.h>
#include <matrix.h>
#include <array.h>
#include <mathutil.h>
#include <vector.h>


VISVector sumOfCols(const VISMatrix  &A)
{
  int i, j, k;
  int rows = A.r(), cols = A.c();
  float total;
  VISVector ret(rows);
  for (i = 0; i < rows; i++)
    {
      total = 0.0f;
      for (j = 0; j < cols; j++)
	total += A.peek(i, j);
      ret.poke(i) = total;
    }
  return(ret);
}

#define EPSILON (1.0e-6)

// this one only ensures that the rows sum to one...
VISMatrix projectUpdateOntoLinearConstraints(const VISMatrix  &update)
{
  int i, j, k;
  float adjust;
  int rows = update.r(), cols = update.c();
  VISMatrix ret(rows, cols);
  for (i = 0; i < rows; i++)
    {
      adjust = 0.0f;
      for (j = 0; j < cols; j++)
	adjust += update.peek(i, j);
      adjust = adjust/cols;
      for (j = 0; j < cols; j++)
	ret.poke(i, j) = update.peek(i, j) - adjust;
    }
  return(ret);
}

// this one only ensures that the rows sum to one...
VISMatrix projectSolutionOntoLinearConstraints(const VISMatrix  &solution)
{
  int i, j, k;
  float adjust;
  int rows = solution.r(), cols = solution.c();
  VISMatrix ret(rows, cols);
  for (i = 0; i < rows; i++)
    {
      adjust = -1.0f;
      for (j = 0; j < cols; j++)
	adjust += solution.peek(i, j);
      adjust = adjust/cols;
      for (j = 0; j < cols; j++)
	ret.poke(i, j) = solution.peek(i, j) - adjust;
    }
  return(ret);
}

VISMatrix projectOntoPositiveSolution(const VISMatrix  &solution)
{
  int i, j, k;
  float adjust, total, this_value;
  int num_items;
  int rows = solution.r(), cols = solution.c();
  boolean done, done2;
  std::vector<float> positive_values;
  std::vector<float>::reverse_iterator vector_iter;
  VISMatrix ret(rows, cols);

  for (i = 0; i < rows; i++)
    {
      done = true;
      for (j = 0; j < cols; j++)
	if (solution.peek(i, j) < 0.0f) done = false;
      if (!done)
	{
	  positive_values.clear();
	  for (j = 0; j < cols; j++)
	    if (solution.peek(i, j) > 0.0f)
	      positive_values.push_back(solution.peek(i, j));
	  positive_values.push_back(0.0);
	  std::sort(positive_values.begin(), positive_values.end());
	  num_items = 0;
	  total = 0.0f;
	  adjust = 0.0f;
	  //	  cout << "sorted values ";
	  done2 = false;
	  for (vector_iter = positive_values.rbegin(); ((vector_iter != positive_values.rend())&&(!done2)); vector_iter++)
	    {
	      this_value = *vector_iter;
	      if (adjust < this_value)
		{
		  num_items++;
		  total += this_value;
		  adjust = (total - 1.0)/(float)num_items;
		}
	      else (done2 = true);
	      //	      cout << this_value << " ";
	    }
	  //	  cout << " adjust is " << adjust << " " << endl;
	  for (j = 0; j < cols; j++)
	    {
	      this_value = solution.peek(i, j) - adjust;
	      ret.poke(i, j) = VISmax(this_value, 0.0f);
	    }
	}
      else 
	for (j = 0; j < cols; j++)
	  ret.poke(i, j) = solution.peek(i, j);
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

  cout << "done distances" << endl;

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
  dt = -0.4/flt_tmp;

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
      //      cout << "weights sum before" << k << sumOfCols(weights[k]) << endl;
      weights[k] = projectSolutionOntoLinearConstraints(weights[k]);
      //      cout << "weights sum middle" << k << sumOfCols(weights[k]) << endl;
      weights[k] = projectOntoPositiveSolution(weights[k]);
      //      cout << "weights " << k << weights[k] << endl;
      //      cout << "weights sum " << k << sumOfCols(weights[k]) << endl;
    }

  cout << "done rand init" << endl;

  update = VISMatrix(num_labels, num_pts);
  float energy;
  VISMatrix matrix_tmp;
  float corres_weight = 0.1*(1.0/num_objects), separation_weight = (1.0/num_labels);

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

  boolean done = false;

  int num_initial_iterations = 150;
  for (int kk = 0; ((kk < 200)&&!done); kk++)
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
	  //	  weights_tmp[k] = vectorToMatrix(solvePosConstrainedUpdate3(A, b, matrixToVector(weights[k]), matrixToVector(dt*update), A_inv), 
	  //					  num_labels, num_pts);

	  weights_tmp[k] = projectOntoPositiveSolution(weights[k] + projectUpdateOntoLinearConstraints(dt*update));
	  //	  cout << "updates sum " << k << " " << sumOfCols(projectUpdateOntoLinearConstraints(dt*update)) << endl;
	  //	  cout << "before pos weights sum " << k << " " << sumOfCols(weights[k] + projectUpdateOntoLinearConstraints(dt*update)) << endl;
	  //	  cout << "**********" << endl;
	  //	  cout << "weights sum " << k << " " << sumOfCols(weights_tmp[k]) << endl;
	}
      
       if (kk > num_initial_iterations)
 	for (k = 0; k < num_objects; k++)
 	  {
 	    weights_tmp[k] *= 1.5f;
	    weights_tmp[k] = projectSolutionOntoLinearConstraints(weights_tmp[k]);
	    weights_tmp[k] = projectOntoPositiveSolution(weights_tmp[k]);
 	  }

       float change = 0.0f;

       for (k = 0; k < num_objects; k++)
	 {
	   change += (weights[k] - weights_tmp[k]).norm();
	   weights[k] = weights_tmp[k];
	   //	  cout << "weights " << k << weights[k] << endl;
	 }

       if (change/num_objects < 1.0e-10)
	 done = true;

      energy = 0.0f;
      for (k = 0; k < num_objects; k++)
	{
	  for (l = 0; l < k; l++)
	    energy += corres_weight*((weights[k])*distances[k][l]*(weights[l]).t()).trace();
	  matrix_tmp = weights[k]*distances[k][l]*weights[k].t();
	  for (i = 0; i < num_labels; i++)
	    matrix_tmp.poke(i, i) = 0.0f;
	  energy += separation_weight*matrix_tmp.sum();
	}

      //	for (l = 0; l < k; l++)
      //	  ;
      //	  energy += ((weights[k])*distances[k][l]*(weights[l]).t()).trace();
      //      energy *= (2.0/(num_objects*(num_objects-1)));
      cout << "Iteration " << kk << " energy " << energy << endl;
    
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
      //      cout << "weights " << k << ": " << weights[k] << endl;
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




