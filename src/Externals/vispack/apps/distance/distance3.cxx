#include <volume.h> 
#include <volumefile.h> 
#include <volutils.h>
#include <ScalarField.h>
#include <AnalyticScalarFields.h>
#include <ReconstructionField.h>
#include <IOScalarField.h>
#include <octree.h>
#include <iostream>
#include <fstream>
//#include <visfield.h>

#define COMPUTE_HESSIAN (false)

#undef EPSILON
#define EPSILON (1.0e-20)

// KEY PARAMETERS....
//#define PROJECTION_TOLERANCE (1.0e-4)
#define PROJECTION_TOLERANCE (1.0e-3)
#define SHORTENING_TOLERANCE (1.0 - 5.0e-3)
#define CURVATURE_TOLERANCE (0.9)

typedef VISVolume< vec<3> > VolumeVec;

#define VOL_W 128
#define VOL_H 128
#define VOL_D 128
#define BOX_W (VOL_W/3.0)
#define BOX_H (VOL_H/3.0)
#define BOX_D (VOL_D/3.0)

double box(unsigned i, unsigned j, unsigned k) 
{
  double pos_x = fabs(i - VOL_W/2.0),
    pos_y = fabs(j - VOL_H/2.0),
    pos_z = fabs(k - VOL_D/2.0);
  pos_x -= BOX_W;
  pos_y -= BOX_H;
  pos_z -= BOX_D;

  if ((pos_x <= 0.0)&&(pos_y <= 0.0)&&(pos_z <= 0.0))
    {
      return(VISmax(pos_x, VISmax(pos_y, pos_z)));
    }

  double total = 0.0;
  if (pos_x > 0.0) total+= pos_x*pos_x;
  if (pos_y > 0.0) total+= pos_y*pos_y;
  if (pos_z > 0.0) total+= pos_z*pos_z;

  return(sqrt(total));
}

void gaussDiffuse(VolumeScalar &vol, double sigma) 
{
  double max_dt = 1.0/12.0, time = sigma*sigma/2.0, extra_time;
  int iterations = (int)floor(time/max_dt);
  extra_time = time-iterations*max_dt;
  int i, j, k, l;
  double laplacian;
  int w = vol.width(), h = vol.height(), d = vol.depth();
  VolumeScalar update(w, d, h);
  for (l = 0; l < iterations+1; l++)
    {
      cout << "gaussDiffuse about to do iteration " << l << endl;
      for (k = 0; k < d; k++)
	for (j = 0; j < h; j++) 
	  for (i = 0; i < w; i++) 
	    {
	      update.poke(i, j, k) = max_dt*-6.0*(vol.peek(i, j, k) 
						   + vol.peek(VISmax(i - 1, 0), j, k) + vol.peek(VISmin(i + 1, w - 1), j, k)
						   + vol.peek(i, VISmax(j - 1, 0), k) + vol.peek(i, VISmin(j + 1, h - 1), k)
						   + vol.peek(i, j, VISmax(k - 1, 0)) + vol.peek(i, j, VISmin(k + 1, d - 1)));
	    }
      if (l < iterations)
	vol += update;
      else
	vol += (float)(extra_time/max_dt)*update;
    }
}


// moves a point onto the zero set of a function defined by a analytic
// scalar field in the form of a ScalarField object
// MOVES only along the ray connecting the pixel position and the surface point
int constraintProject(vec<3> &ret, const vec<3> &pixel_pos, 
const vec<3> &surface_pos, ScalarField<3> *field, 
double tolerance, boolean debug = false)
{
  vec<3> gradient, new_position = surface_pos, update, direction; 
  double value;
  ScalarFieldParams<3> params;
  
  if (field->checkBounds(new_position))
    //    value = field->value(new_position);
    {
      field->computeScalarFieldParams(new_position, params, COMPUTE_HESSIAN);
      value = params._F;
      gradient = params._Fx;
    }
  else
    return(-1);
  double error = FLT_MAX;
  int iter = 0, max_iter = 1000; 

  direction = surface_pos - pixel_pos;
  direction = direction.normalize();

  //  double lamdba = 1.0;
  //  boolean done = false;

  if (debug)
    {
      cout << "starting value " << value << endl;
      cout << "starting position " << new_position << endl;
    }

  vec<3> old_position;
  double old_value;

  double lambda, min_lambda = 1.0e-6;
  lambda = min_lambda;
  while ((error > tolerance) &&(iter++ < max_iter)&&(lambda < 1.0e12))
    //while (error > 1.0e-4)
    {
      old_value = value;
      old_position = new_position;
      //      gradient = field->gradient(new_position);
      update = -1.0*value*direction/(DotProduct(gradient, direction) + lambda);
      new_position += update;

      if (field->checkBounds(new_position))
	{
	  field->computeScalarFieldParams(new_position, params,COMPUTE_HESSIAN);
	  value = params._F;
      //	  value = field->value(new_position);	
	  if (fabs(value) > fabs(old_value))
	    {
	      lambda *= 10.0;
	      value = old_value; 
	      new_position = old_position;
	    }
	  else
	    {
	      lambda = VISmax(lambda/10.0, (double)min_lambda);
	      gradient = params._Fx;
	      error = fabs(value/DotProduct(gradient, direction));
	    }
	}
      else
	{
	  if (debug)
	    {
	      cout << "out of bounds at " << new_position << endl;
	    }
	  lambda *= 10.0;
	  new_position = old_position;
	}
      
      if (debug)
	{
	  cout << "value is " << value << endl;
	  cout << "lambda is " << lambda << endl;
	  cout << "gradient " << gradient << endl;
	  cout << "new position " << new_position << endl;
	}
      //      error = update.length();
    }
  ret = new_position;
  if (iter >= max_iter) 
    {
      //      if (!debug)
      //	constraintProject(pixel_pos, surface_pos, field,
      //	tolerance, true);
      if (debug)
	{
	  cout <<  "warning on project - exceeded iterations with value of " 
	       << value << " and tolerance of " << tolerance << " position " << new_position << endl;
		
	}
      return(-1);
    }
  return(0);
}


// moves a point onto the zero set of a function defined by a analytic
// scalar field in the form of a ScalarField object
int constraintProject(vec<3> &ret, 
const vec<3> &surface_pos, ScalarField<3> *field, 
double tolerance, boolean debug = false)
{
  vec<3> gradient, new_position = surface_pos, update; 
  double value;
  ScalarFieldParams<3> params;
  
  if (field->checkBounds(new_position))
    //    value = field->value(new_position);
    {
      field->computeScalarFieldParams(new_position, params,COMPUTE_HESSIAN);
      value = params._F;
      gradient = params._Fx;
    }
  else
    return(-1);

  double error = FLT_MAX;
  int iter = 0, max_iter = 50; 

  //  double lamdba = 1.0;
  //  boolean done = false;

  if (debug)
    {
      cout << "starting value " << value << endl;
      cout << "starting position " << new_position << endl;
    }

  vec<3> old_position;
  double old_value;

  double lambda, min_lambda = 1.0e-6;
  lambda = min_lambda;
  while ((error > tolerance) &&((iter++) < max_iter)&&(lambda < 1.0e12))
    //while (error > 1.0e-4)
    {
      old_value = value;
      old_position = new_position;
      //      gradient = field->gradient(new_position);
      update = -1.0*value*gradient/(gradient.lengthSqr() + lambda);
      new_position += update;

      if (iter > 100)
	cout << "got a problem" << endl;

      if (field->checkBounds(new_position))
	{
	  field->computeScalarFieldParams(new_position, params);
	  value = params._F;
      //	  value = field->value(new_position);	
	  if (fabs(value) > fabs(old_value))
	    {
	      lambda *= 10.0;
	      value = old_value; 
	      new_position = old_position;
	      if (debug) cout << "lambda_up" << endl;
	      error = fabs(value/gradient.length());
	    }
	  else
	    {
	      lambda = VISmax(lambda/10.0, (double)min_lambda);
	      gradient = params._Fx;
	      error = fabs(value/gradient.length());
	      if (debug) 
		{
		  cout << "got update of difference " << (old_position - new_position) << endl;
		  cout << "update is " << update << endl;
		}
	      
	    }
	}
      else
	{
	  if (debug)
	    {
	      cout << "out of bounds at " << new_position << endl;
	    }
	  lambda *= 10.0;
	  new_position = old_position;
	}
      
      if (debug)
	{
	  //	  cout << "value is " << value << endl;
	  //	  cout << "updated value is " << params._F << endl;
	  //	  cout << "lambda is " << lambda << endl;
	  cout << "old position " << old_position << endl;
	  field->computeScalarFieldParams(old_position, params);
	  cout << "old value is " << params._F << endl;
	  cout << "old gradient " << params._Fx << endl;
	  cout << "new position " << new_position << endl;
	  field->computeScalarFieldParams(new_position, params);
	  cout << "new value is " << params._F << endl;
	  cout << "new gradient " << params._Fx << endl;
	  cout << "difference of old and new positions " << (old_position - new_position) << endl;
	  //	  cout << "error is " << error << endl;
	}
      //      error = update.length();
    }
  ret = new_position;
  if (iter >= max_iter) 
    {
      //      if (!debug)
      //	constraintProject(pixel_pos, surface_pos, field,
      //	tolerance, true);
      if (debug)
	{
	  cout <<  "warning on project - exceeded iterations with value of " 
	       << value << " and tolerance of " << tolerance << " position " << new_position << endl;
		
	}
      return(-1);
    }
  return(0);
}


// computes an upwind derivative for the field - "upwind" is defined as
// "toward the zero set"  
// Used to initialize surface positions from a distance field.
VolumeVec upwindGradient(const VolumeScalar &vol)
{
  int w, h, d; 
  int i, j, k;
  VolumeVec ret(w = vol.width(), h = vol.height(), d = vol.depth());
  double grad_x, grad_y, grad_z, grad_before, grad_after;
  double this_value;

  for (k = 0; k < d; k++) 
    for (j = 0; j < h; j++) 
      for (i = 0; i < w; i++) 
	{
	  this_value = vol.peek(i, j, k);
	  if (i < (w - 1)) 
	    grad_before = (vol.peek(i + 1, j, k) - vol.peek(i, j, k));
	  else 
	    grad_before = 0;
			
	  if (i > 0)
	    grad_after = (vol.peek(i, j, k) - vol.peek(i - 1, j, k));
	  else 
	    grad_after = 0.0;

	  // look back to the zero set...
	  if (this_value > 0.0)
	    {
	      grad_before = VISmin(grad_before, 0.0);
	      grad_after = VISmax(grad_after, 0.0);		       
	      grad_x = VISMaxAbs(grad_before, grad_after);
	    }
	  else
	    {
	      grad_before = VISmax(grad_before, 0.0);
	      grad_after = VISmin(grad_after, 0.0);		       
	      grad_x = VISMaxAbs(grad_before, grad_after);
	    }


	  if (j < (h - 1)) 
	    grad_before = (vol.peek(i, j + 1, k) 
			     - vol.peek(i, j, k));
	  else 
	    grad_before = 0.0;

	  if (j > 0)
	    grad_after = (vol.peek(i, j, k)
			    - vol.peek(i, j - 1, k) );
	  else
	    grad_after = 0.0;

	  // look back to the zero set...
	  if (this_value > 0.0)
	    {
	      grad_before = VISmin(grad_before, 0.0);
	      grad_after = VISmax(grad_after, 0.0);		       
	      grad_y = VISMaxAbs(grad_before, grad_after);
	    }
	  else
	    {
	      grad_before = VISmax(grad_before, 0.0);
	      grad_after = VISmin(grad_after, 0.0);		       
	      grad_y = VISMaxAbs(grad_before, grad_after);
	    }


	  if (k < (d - 1)) 
	    grad_before = (vol.peek(i, j, k + 1) 
			     - vol.peek(i, j, k));
	  else 
	    grad_before = 0.0;

	  if (k > 0)
	    grad_after = (vol.peek(i, j, k)
			    - vol.peek(i, j, k - 1) );
	  else
	    grad_after = 0.0;

	  // look back to the zero set...
	  if (this_value > 0.0)
	    {
	      grad_before = VISmin(grad_before, 0.0);
	      grad_after = VISmax(grad_after, 0.0);		       
	      grad_z = VISMaxAbs(grad_before, grad_after);
	    }
	  else
	    {
	      grad_before = VISmax(grad_before, 0.0);
	      grad_after = VISmin(grad_after, 0.0);		       
	      grad_z = VISMaxAbs(grad_before, grad_after);
	    }
	  
	  ret.poke(i, j,k) = vec<3> (grad_x, grad_y, grad_z);
	  
	}
  return(ret);
}



// // moves a point onto the zero set of a function defined by a analytic
// // scalar field in the form of a ScalarField object
// int constraintProject(vec<3> &ret, 
// 		      const vec<3> &pixel_pos, 
// 		      const vec<3> &surface_pos, 
// 		      const VolumeScalar &vol, 
// 		      const VolumeVec &grad_vol, 
// 		      double tolerance, 
// 		      boolean debug = false)
// {
//   vec<3> gradient, new_position = surface_pos, update; 
//   double x, y, z;
//   double value;
//   if (vol.checkBounds((float)(x = new_position[0]), 
// 		      (float)(y = new_position[1]), 
// 		      (float)(z = new_position[2])))
//     {
//       //      cout << "place 1" << endl;
//       value = vol.interp(x, y, z);
//     }
//   else
//     {
//       if (debug) 
// 	cout << "got interp out of bounds for constrain project" << endl;
//       cout << "new position" << new_position << endl;
//       return(-1);
//     }
//   double error = FLT_MAX;
//   int iter = 0, max_iter = 500; 

//   //  double lamdba = 1.0;
//   //  boolean done = false;

//   if (debug)
//     {
//       cout << "starting value " << value << endl;
//       cout << "starting position " << new_position << endl;
//     }

//   vec<3> old_position;
//   double old_value;

//   double lambda, min_lambda = 1.0e-6;
//   lambda = min_lambda;
//   while ((fabs(value) > tolerance) &&(iter++ < max_iter)&&(lambda < 1.0e12))
//     //while (error > 1.0e-4)
//     {
//       old_value = value;
//       old_position = new_position;
//       //      cout << "place 2" << endl;
//       gradient = grad_vol.interp(x, y, z);
//       update = -1.0*value*gradient/(gradient.lengthSqr() + lambda);
//       new_position += update;

//       if (vol.checkBounds((float)(x = new_position[0]), 
// 			  (float)(y = new_position[1]), 
// 			  (float)(z = new_position[2])))
// 	{
// 	  //	  cout << "place 3" << endl;
// 	  value = vol.interp(x, y, z);
// 	  if (fabs(value) > fabs(old_value))
// 	    {
// 	      lambda *= 10.0;
// 	      value = old_value; 
// 	      new_position = old_position;
// 	      x = new_position[0]; y = new_position[1]; z = new_position[2];
// 	    }
// 	  else
// 	    {
// 	      lambda = VISmax(lambda/10.0, (double)min_lambda);
// 	    }
// 	}
//       else
// 	{
// 	  if (debug)
// 	    {
// 	      cout << "out of bounds at " << new_position << endl;
// 	    }
// 	  lambda *= 10.0;
// 	  new_position = old_position;
// 	  x = new_position[0]; y = new_position[1]; z = new_position[2];
// 	}
//       if (debug)
// 	{
// 	  cout << "value is " << value << endl;
// 	  cout << "lambda is " << lambda << endl;
// 	  cout << "gradient " << gradient << endl;
// 	  cout << "new position " << new_position << endl;
// 	}
      
//     }
//   ret = new_position;
//   if (iter >= max_iter) 
//     {
//       //      if (!debug)
//       //	constraintProject(pixel_pos, surface_pos, field,
//       //	tolerance, true);
//       if (debug)
// 	{
// 	  cout <<  "warning on project - exceeded iterations with value of " 
// 	       << value << " and tolerance of " << tolerance << " position " << new_position << endl;
		
// 	}
//       return(-1);
//     }
//   return(0);
// }



// This moves the point along the surface to short the distance
// between surface_pos and pixel_pos.  This is done by simply
// projecting the difference vector onto the tangent plane of the
// surface. 
// This should return zero when the difference vector is perpendicular
// to the surface, which is, by definition, and extremal point of
// position. 
int shortenDistance(vec<3> &ret, const vec<3> &pixel_pos, const vec<3> &surface_pos, 
		    ScalarField<3> *field, double &dt, boolean debug=false)
{
  vec<3> gradient, difference_vector, new_position;
  vec<3>  update;
  matrix<3,3> tangent_plane;
  ScalarFieldParams<3> params;
  static matrix<3, 3> ident = Identity<3>();
  if (debug) cout << "pixel is " << pixel_pos << endl;
  difference_vector = pixel_pos - surface_pos;
  if (debug) cout << "surface_position is " << surface_pos << endl;
  if (debug) cout << "difference vector " << difference_vector << endl;
  if (field->checkBounds(surface_pos))
    {
      field->computeScalarFieldParams(surface_pos, params,COMPUTE_HESSIAN);
      gradient = params._Fx;
    }
  else return(1);
  gradient /= gradient.length(EPSILON);
  if (debug)  cout << "gradient " << gradient << endl;
  if (debug)  cout << "diff dot gradient " << DotProduct(gradient, difference_vector/difference_vector.length(1.0e-20)) << endl;
  tangent_plane =  ident - DirectProduct(gradient, gradient); 
  update = tangent_plane*difference_vector;
  if (debug)  cout << "update " << update << endl;
  new_position = surface_pos + dt*update;
  //  cout << "new position " << new_position << endl;
  ret = new_position;
  //  difference_vector = pixel_pos - new_position;
  //  field->computeScalarFieldParams(new_position, params);
  //  gradient = params._Fx;
  return(0);
  //  return(sqrt(1.0 - power(dotProduct(gradient, difference_vector.normalize(1.0e-20)))));
}

#define VOL_SIZE (128)
#define DT_GROWTH (10.0)

int minimizeDistance(vec<3> &ret, const vec<3> &pixel_pos, 
		     const vec<3> &surface_pos, ScalarField<3> *field, 
		     double tolerance, double tolerance_project, boolean debug = false)
{
  double error= FLT_MAX, old_distance, new_distance, tmp_distance;
  double dot_prod = 0.0;
  vec<3> surface_position = surface_pos, new_position, tmp_position, gradient, difference_vector;
  boolean done2, return_error_code, done = false;
  int iterations = 0;
  double dt = 0.1;
  ScalarFieldParams<3> params;
  double dt_min = 1.0e-20;
  int max_iter = 1000;
  double this_change;

  double original_distance = (surface_position - pixel_pos).length();

  difference_vector = pixel_pos - surface_position;
  difference_vector = difference_vector.normalize(1.0e-20);
  field->computeScalarFieldParams(surface_position, params);
  gradient = params._Fx;
  gradient = gradient.normalize(1.0e-20);
  if ((dot_prod  = fabs(DotProduct(difference_vector, gradient))) > tolerance)
    done = true;
  if (debug) cout << "diff dot gradient before min " << fabs(DotProduct(difference_vector, gradient)) << endl;
  
  //  while(error > tolerance)
  while (!done) // (dot_prod < tolerance)&&(iterations++ < max_iter)
    {
      //      if (iterations == 5000) debug = true; 
      //		if (iterations++ > 50) debug = true;
      old_distance = (surface_position - pixel_pos).length();
      // this moves the surface point in the tangent plane
      // to be closer to the pixel position
      return_error_code = shortenDistance(tmp_position, pixel_pos, surface_position, field, dt, debug);
      // put point back on surface
      //      tmp_distance = (pixel_pos - tmp_position).length();
      if (!return_error_code)
	{
	  //	  return_error_code = constraintProject(new_position, pixel_pos, tmp_position, field, tolerance_project /*, debug*/);
	  return_error_code = constraintProject(new_position, tmp_position, field, tolerance_project /*, debug*/);
	  //		    new_position = constraintProject(pixel_pos, new_position, field, tolerance,true);
	}
      if (!return_error_code)
	// compute new dot sign of angle
	{
	  new_distance = (pixel_pos - new_position).length();
	}
      else 
	{
	  //	  cout << "got error code on constraint project " << endl;
	  //	  cout << "surf pos " << surface_position;
	  //	  cout << "pixel pos " << pixel_pos;
	  //	  constraintProject(new_position, pixel_pos, tmp_position, field, tolerance, true);
	  new_distance = FLT_MAX;
	}

      // this is an adaptive time step scheme
      done2 = false;

      // code snippet to get data for a specific pixel -- for debugging
      //      if (pixel_pos == vec<3>(14.0, 9.0, 9.0)) 
      //	{
      //	  cout << "got target vector" << endl;
      //	  debug=true;
      //	}

      // if time step is too big, then you overshot, and you
      // need to shorten this update
      if ((new_distance > old_distance))
	while (!done2) 
	  {
	    //	    if (tmp_distance < old_distance)
	    //	      cout << "tmp distance is better but not projected" << endl;
	    //	   
	    //reduce time step
	    //	    dt = VISmax(dt/2.0, 1.0e-20);
	    dt = VISmax(dt/DT_GROWTH, (double)FLT_MIN);
	    if (debug) cout << "down dt is " << dt << endl;
	    return_error_code = shortenDistance(tmp_position, pixel_pos, surface_position, field, dt, debug);
	    // put point back on surface
	    if (!return_error_code)
	      {
		//		return_error_code = constraintProject(new_position, pixel_pos, tmp_position, field, tolerance_project/*, debug*/);
		return_error_code = constraintProject(new_position, tmp_position, field, tolerance_project/*, debug*/);
		//		    new_position = constraintProject(pixel_pos, new_position, field, tolerance,true);
	      }
	    //	    tmp_distance = (pixel_pos -tmp_position).length();
	    if (!return_error_code)
	      {
		// compute new distance
		new_distance = (pixel_pos - new_position).length();
		this_change = (surface_position - new_position).length()/old_distance;
	      }
	    else
	      {
		new_distance = FLT_MAX;
		//		cout << "WARNING - could not get valid output from constrain project" << endl;
		//		constraintProject(new_position, pixel_pos, tmp_position, field, tolerance, true);
	      }

	    if (new_distance <= old_distance) done2 = true;
	    //	    if (iterations++ > 2000) done2 = true;
	    if (this_change < 1.0e-4) {done2 = true; done = true;}

	    if (debug) cout << "pixel position" << pixel_pos << endl;
	    if (debug) cout << "surface position " << surface_position << endl;
	    if (debug) cout << "new distance " << new_distance << endl;
	    if (debug) cout << "new-old  distance " << (new_distance - old_distance) << endl;
	    if (debug) cout << "new-old  distance before project" << ((pixel_pos - tmp_position).length() - old_distance) << endl;
	    if (debug) cout << "this_change " << this_change << endl;
	    if (debug) cout << "iterations " << iterations << endl;
	    {
	      field->computeScalarFieldParams(pixel_pos, params,COMPUTE_HESSIAN);
	      if (debug) cout << "value " << params._F << endl;
	    }
	    {
	      field->computeScalarFieldParams(surface_pos, params,COMPUTE_HESSIAN);
	      if (debug) cout << "surface value " << params._F << endl;
	    }
	    if (iterations++ >= max_iter)
	      return(-1);
	  }
      // be optimistic.  if it worked, go for a bigger time
      // step next time around
      else 
	{
	  if (iterations++ >= max_iter)
	    return(-1);
	  dt*=DT_GROWTH;
	  if (debug) cout << "up dt is " << dt << endl;
	}
      //	  project points to surface
      //		error = (new_position -
      //		surface_positionition).length()/new_distance;
      if (!return_error_code)
	{
	  //	  error = (new_position - surface_position).length()/new_distance;
	  field->computeScalarFieldParams(new_position, params,COMPUTE_HESSIAN);
	  gradient = params._Fx;
	  gradient = gradient.normalize(1.0e-20);
	  difference_vector = pixel_pos - new_position;
	  difference_vector = difference_vector.normalize(1.0e-20);
	  
	  if (debug) 
	    {
	      //	      cout << "error is " << error << endl;
	      cout  << "old " << surface_position << " new " << new_position << endl;
	    }
	  surface_position = new_position;
	}
      else
	{
	  //	  error = 0.0;
	  done = true;
	  //	  dot_prod = 1.0;
	  //	      cout << "WARNING - could not get valid output from optmize" << endl;
	  //	      cout << "new position " << new_position << endl;
	  //	      cout << "old position " << surface_position << endl;
	  //	      cout << "dt " << dt << endl;
	}
      if ((dot_prod  = fabs(DotProduct(difference_vector, gradient))) > tolerance)
	done = true;
      if (debug) cout << "diff dot gradient in min " << fabs(DotProduct(difference_vector, gradient)) << endl;
	

    }
      
//       double bbox[3][2];
//       vec<3> node;
//       bbox[0][0] = floor(surface_pos[0]);
//       bbox[1][0] = floor(surface_pos[1]);
//       bbox[2][0] = floor(surface_pos[2]);
//       bbox[0][1] = ceil(surface_pos[0]);
//       bbox[1][1] = ceil(surface_pos[1]);
//       bbox[2][1] = ceil(surface_pos[2]);
//      cout << "got max iter for pixel" << pixel_pos << " and dot prod is " << dot_prod << endl;
//       cout << "origina dist " << original_distance << " new distance " << new_distance << endl;
//       cout << "surface position is " << surface_pos << endl;
//       cout << "dt is " << dt << endl;
//       field->computeScalarFieldParams(surface_pos, params);
//       cout << "surface value " << params._F << endl;
//       field->computeScalarFieldParams(surface_pos, params);
//       cout << "surface gradient " << params._Fx << endl;
//       field->computeScalarFieldParams(node = vec<3>(bbox[0][0], bbox[1][0], bbox[2][0]), params);
//       cout << "value at node " << node << " is " << params._F << endl;
//       field->computeScalarFieldParams(node = vec<3>(bbox[0][1], bbox[1][0], bbox[2][0]), params);
//       cout << "value at node " << node << " is " << params._F << endl;
//       field->computeScalarFieldParams(node = vec<3>(bbox[0][0], bbox[1][1], bbox[2][0]), params);
//       cout << "value at node " << node << " is " << params._F << endl;
//       field->computeScalarFieldParams(node = vec<3>(bbox[0][1], bbox[1][1], bbox[2][0]), params);
//       cout << "value at node " << node << " is " << params._F << endl;
//       field->computeScalarFieldParams(node = vec<3>(bbox[0][0], bbox[1][0], bbox[2][1]), params);
//       cout << "value at node " << node << " is " << params._F << endl;
//       field->computeScalarFieldParams(node = vec<3>(bbox[0][1], bbox[1][0], bbox[2][1]), params);
//       cout << "value at node " << node << " is " << params._F << endl;
//       field->computeScalarFieldParams(node = vec<3>(bbox[0][0], bbox[1][1], bbox[2][1]), params);
//       cout << "value at node " << node << " is " << params._F << endl;
//       field->computeScalarFieldParams(node = vec<3>(bbox[0][1], bbox[1][1], bbox[2][1]), params);
//       cout << "value at node " << node << " is " << params._F << endl;

  ret = surface_position;
  return(0);
}


// computes an upwind derivative for the field - "upwind" is defined as
// "toward the zero set"  
// Used to initialize surface positions from a distance field.
int upwindGradient(vec<3> &ret, const VolumeScalar &vol, int i, int j, int k)
{
  if (!vol.checkBounds(i, j, k))
    return(-1);
  
  double grad_x, grad_y, grad_z, grad_before, grad_after;
  double this_value;
  int w = vol.width(), h = vol.height(), d = vol.depth();

  {
    this_value = vol.peek(i, j, k);
    if (i < (w - 1)) 
      grad_before = (vol.peek(i + 1, j, k) - vol.peek(i, j, k));
    else 
      grad_before = 0;
			
    if (i > 0)
      grad_after = (vol.peek(i, j, k) - vol.peek(i - 1, j, k));
    else 
      grad_after = 0.0;

    // look back to the zero set...
    if (this_value > 0.0)
      {
	grad_before = VISmin(grad_before, 0.0);
	grad_after = VISmax(grad_after, 0.0);		       
	grad_x = VISMaxAbs(grad_before, grad_after);
      }
    else
      {
	grad_before = VISmax(grad_before, 0.0);
	grad_after = VISmin(grad_after, 0.0);		       
	grad_x = VISMaxAbs(grad_before, grad_after);
      }

    if (j < (h - 1)) 
      grad_before = (vol.peek(i, j + 1, k) 
		     - vol.peek(i, j, k));
    else 
      grad_before = 0.0;

    if (j > 0)
      grad_after = (vol.peek(i, j, k)
		    - vol.peek(i, j - 1, k) );
    else
      grad_after = 0.0;

    // look back to the zero set...
    if (this_value > 0.0)
      {
	grad_before = VISmin(grad_before, 0.0);
	grad_after = VISmax(grad_after, 0.0);		       
	grad_y = VISMaxAbs(grad_before, grad_after);
      }
    else
      {
	grad_before = VISmax(grad_before, 0.0);
	grad_after = VISmin(grad_after, 0.0);		       
	grad_y = VISMaxAbs(grad_before, grad_after);
      }

    if (k < (d - 1)) 
      grad_before = (vol.peek(i, j, k + 1) 
		     - vol.peek(i, j, k));
    else 
      grad_before = 0.0;

    if (k > 0)
      grad_after = (vol.peek(i, j, k)
		    - vol.peek(i, j, k - 1) );
    else
      grad_after = 0.0;

    // look back to the zero set...
    if (this_value > 0.0)
      {
	grad_before = VISmin(grad_before, 0.0);
	grad_after = VISmax(grad_after, 0.0);		       
	grad_z = VISMaxAbs(grad_before, grad_after);
      }
    else
      {
	grad_before = VISmax(grad_before, 0.0);
	grad_after = VISmin(grad_after, 0.0);		       
	grad_z = VISMaxAbs(grad_before, grad_after);
      }
  }
  ret = vec<3>(grad_x, grad_y, grad_z);
  return(0);
}


vec<3> findBoxIntersection(vec<3> pixel_pos, vec<3> gradient)
{
  if (fabs(gradient[0]) > fabs(gradient[1]))
    if (fabs(gradient[0]) > fabs(gradient[2]))
      {
	return(pixel_pos + gradient/fabs(gradient[0]));
      }
    else
      {
	return(pixel_pos + gradient/fabs(gradient[2]));
      }
  else 
    if (fabs(gradient[1]) > fabs(gradient[2]))
      {
	return(pixel_pos + gradient/fabs(gradient[1]));
      }
    else
      {
	return(pixel_pos + gradient/fabs(gradient[2]));
      }
}


int main(int argc, char** argv)
{
  VISVolumeFile volume_file;
  ScalarField<3> *field;

  // to restrict compution to a bbox, e.g. for debugging...
  int x_lo, x_hi, y_lo, y_hi, z_lo, z_hi;

  //  if (argc > 1)
  
    
  int i, j, k, l;

  // This is a GUESS at how much room miriah's recon stuff needs in order not to crash.  Who knows...
  int recon_border = 0;

  // This just sets up the field and the distance transform. 
  // This needs to be hacked to however the application needs to do it.
  //  if (!vol.isValid())
  // SETTING UP DATA OBJECTS...
  //     {
  //       //      vol = VolumeScalar(VOL_W, VOL_H, VOL_D);
  //       vol = VolumeScalar(VOL_SIZE, VOL_SIZE, VOL_SIZE);
  //       read_raw("dt_cube.raw", vol);
  //       //      cout << "vol max " << vol.max() << " and min " << vol.min() << endl;
  //       vol_grad = upwindGradient(vol);
  //       //      field = new Ellipse(BOX_W, BOX_H, BOX_D);
  //       field = new SplineCube(VOL_W, VOL_H, VOL_D, (0.2f/2.0)*VOL_D);
  //     }

  //reconstruction field on the volume
  //  vol = VolumeScalar(VOL_SIZE, VOL_SIZE, VOL_SIZE);
  //  vol = VolumeScalar(271, 390, 282);

  //  read_raw("dt_cube.raw", vol);
  //  vol = VolumeScalar(volume_file.read_double("dt_pelvis.vol"));
  //  vol = VolumeScalar(volume_file.read_double("dt_cube.vol"));

  VISImageFile im_file;
  ScalarFieldParams<3> params;
  int material;
  float epsilon;
  int args = 1;

  //  vol = VolumeScalar(volume_file.read_double("smth_dt_brain_flux.vol"));
  //  im_file.write(((vol.image()).becomeFlat()).abs(), "smth_dt_brain_flux.fts");  
  //  exit(-1);

  char infile[80], outname_pts[80], outname_vol[80];

  if (argc > args)
    sscanf(argv[args++], "%s", infile);
  else
    sscanf("test.vol", "%s", infile);

  if (argc > args)
    sscanf(argv[args++], "%d", &material);
  else
    sscanf("0", "%d", &material);

  if (argc > args)
    sscanf(argv[args++], "%f", &epsilon);
  else
    sscanf("-1.0f", "%f", &epsilon);


  
  //  field = new VolumetricField("cube.vol");
  //   field = new VolumetricField("dt_pelvis.vol");
  //  field = new VolumetricField(fname1);
  field = new IOScalarField<3>(infile, IOScalarField<3>::INTERPOLATING);
  ((IOScalarField<3>*)field)->setMainIndicatorFunction(material);
  if (epsilon > 0.0f)
    ((IOScalarField<3>*)field)->setEpsilon(epsilon);

  cout << "done field" << endl;
  double spacing[3];
  vec<3> field_size;

  cout << "about to do the bbh" << endl;
  field_size = field->boundingBoxHigh();
  int w = field->xdim(), h = field->ydim(), d = field->zdim();
  spacing[0] = field_size(0)/(double)(w-1);
  spacing[1] = field_size(1)/(double)(h-1);
  spacing[2] = field_size(2)/(double)(d-1);
  cout << "spacing = " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << endl;
  cout << "dimensions = " << w << ", " << h << ", " << d << endl;
  vec<3> gradient, difference_vector, new_position, tmp_position,
    pixel_position, surface_position, bad_position;
  double value, old_distance, new_distance;

  std::filebuf *fb;
  //  std::ostream os_pts;
  if (argc > args)
    {
      sscanf(argv[args++], "%s", outname_pts);
      //      os_pts.open(outname_pts);
      fb = new std::filebuf();
      fb->open(outname_pts, ios::out);
      //      os_pts = ostream();
      cout << "opened outname " << outname_pts << endl;
    }
  else
    {
      //      os_pts = std::cout;
      std::ostream os_pts(std::cout.rdbuf());
      //      os_pts = ostream(std::cout.rdbuf());
    }
  
  std::ostream os_pts(fb);


  if (argc > args)
    {
      sscanf(argv[args++], "%s", outname_vol);
    }
  else
    {
      sscanf("vol_ma.fts", "%s", outname_vol);
    }

  if (argc > args)
    {
      sscanf(argv[args++], "%d", &x_lo);
      sscanf(argv[args++], "%d", &x_hi);
      sscanf(argv[args++], "%d", &y_lo);
      sscanf(argv[args++], "%d", &y_hi);
      sscanf(argv[args++], "%d", &z_lo);
      sscanf(argv[args++], "%d", &z_hi);
    }
  else
    {
      x_lo = y_lo = z_lo = 0;
      x_hi = w; y_hi = h; z_hi = d;
    }



  
  VolumeVec vol_position(w, h, d);
  VolumeScalar vol(w, h, d);
  vec<3> this_grad;
  double this_value;
  double dt, error;
  boolean done, done2;

  bad_position = -FLT_MAX;

  // INITIALIZE THINGS USING DISTANCE FIELD IF AVAILABLE...


  double max_value = -FLT_MAX, min_value = FLT_MAX;

  VolumeScalar vol_out(w, h, d);
  //   double max_value = -FLT_MAX, min_value = FLT_MAX;
  //   // find min and max value of field in order to set tolerances
  cout << "about to sample field"  << endl;
  for (k = recon_border; k < (d - recon_border); k++)
    {
      cout << "sameple slice " << k << " of " << d << endl;
      for (j = recon_border; j < (h - recon_border); j++)
	for (i = recon_border; i < (w - recon_border); i++)
	  {
	    //pixel_position = vec<3>((double)i, (double)j, (double)k);
            pixel_position = vec<3>((double)i*spacing[0],
                                    (double)j*spacing[1],
                                    (double)k*spacing[2]);
	    field->computeScalarFieldParams(pixel_position, params);
	    this_value = params._F;
	    max_value = VISmax(max_value, this_value);
	    min_value = VISmin(min_value, this_value);
	    vol.poke(i, j, k) = this_value;
	  }
    }


  int ll;
    int N[6][3];
  N[0][0] = 1;
  N[0][1] = 0;
  N[0][2] = 0;
  N[1][0] = -1;
  N[1][1] = 0;
  N[1][2] = 0;
  N[2][0] = 0;
  N[2][1] = 1;
  N[2][2] = 0;
  N[3][0] = 0;
  N[3][1] = -1;
  N[3][2] = 0;
  N[4][0] = 0;
  N[4][1] = 0;
  N[4][2] = 1;
  N[5][0] = 0;
  N[5][1] = 0;
  N[5][2] = -1;


  // set the border to some clearly outside value
  //  vol.setBorder(100.0, recon_border);
  double next_value;
  vec<3> this_position, next_position, position;
  Point *this_point, point_tmp;
  double scaling = 1.0;

  OctreeNode<Point> octree;
  vec<3> start(0.0);
  octree.setDimensions( start, field_size);
  cout << "field size " << field_size << endl;

  int num_pts = 0;
  for (k = 0; k < (d - 1); k++)
    {
      //      cout << "zeros slice " << k << endl;
      for (j = 0; j < (h - 1); j++)
	for (i = 0; i < (w - 1); i++)
	  {
	    for (ll = 0; ll < 6; ll+=2)
	      {
		this_value = vol.peek(i, j , k);
		next_value = vol.peek(i + N[ll][0], j + N[ll][1], k+ N[ll][2]);
		if (this_value*next_value <= 0.0)
		  {
		    this_position  = vec<3>(spacing[0]*(double)i,
                                            spacing[1]*(double)j,
                                            spacing[2]*(double)k);
		    next_position  = vec<3>(spacing[0]*(double)(i + N[ll][0]),
                                            spacing[1]*(double)(j + N[ll][1]),
                                            spacing[2]*(double)(k + N[ll][2]));
		    this_point =
                      new Point((this_position + (next_position - this_position)*((this_value)/(this_value  - next_value))));
		    //		    this_point = new Point(this_position);
		  //		  cout << "this point is " << *this_point << endl;
		    //		    cout << "this " << this_position << " next " << next_position << " value " << ((this_value)/(this_value  - next_value))  << endl;
		      //		    cout << this_point->position()[0] << " " << this_point->position()[1] << " " << this_point->position()[2] << " 0.0 0.0 0.0 0.1" <<endl;
		    //		    cout << this_point->position()[0] << ", " << this_point->position()[1]<< ", " << this_point->position()[2] << endl;
		    octree.addPoint(this_point);
		    num_pts++;
		  }
	      }
	  }
    }

  cout << "there are " << num_pts << " pts " << endl;


//   for (k = 0; k < d ; k++)
//     {
//       cout << "distance slice " << k << endl;
//     for (j = 0; j < h; j++)
//       for (i = 0; i < w; i++)
// 	{
// 	  point_tmp = Point(vec<3>(scaling*i, scaling*j, scaling*k));
// 	  //	  cout << "this point is " << point_tmp << endl;
// 	  surface_position = (octree.closestPoint(&point_tmp)).position(); 

// 	  //	  cout << "surface point is " << surface_position << endl;
// 	  vol.poke(i, j, k) = (point_tmp- surface_position).length();
// 	  //	  vol_out.poke(i, j, k) = (octree.distanceToClosestPoint(&point_tmp)); 
// 	}
//     }

  int iterations = 0;
  boolean debug, return_error_code;
  double mc_distance, grad_des_distance;
  //  double tolerance_project, tolerance = (1.0 - 5.0e-4);//1.0e-4;
  double tolerance_project, tolerance = SHORTENING_TOLERANCE;
  //  tolerance_project = 1.0e-5;   // this works ok
  tolerance_project = PROJECTION_TOLERANCE;  // is this enough

  for (k = z_lo; k < z_hi; k++)
    {
      cout << "about to do slice " << k << endl;	    
      for (j = y_lo; j < y_hi; j++)
	{
	  //	  cout << "about to do row " << j << endl;	    
	  for (i = x_lo; i < x_hi; i++)
	    {
	      //	      if (field->value(pixel_position) > -1.0) 
		{
		  iterations = 0;
		  //	    if ((j == 172)&&(i == 9)) debug = true; else 
		  debug = false;
		  pixel_position = vec<3>(spacing[0]*(double)i,
                                          spacing[1]*(double)j,
                                          spacing[2]*(double)k);

		  // 	  point_tmp 
		  // 	  //	  cout << "this point is " << point_tmp << endl;
		  // 	  surface_position = 

		  done = false;
		  error = 1.0;
		  dt = 0.1;
		  //pixel_position = vec<3>((double)i, (double)j, (double)k);

		  //	  cout << "this point is " << point_tmp << endl;
		  surface_position = (octree.closestPoint(&Point(pixel_position))).position(); 
		  //		  return_error_code = constraintProject(new_position, pixel_position,
                  //                                                      surface_position, field, tolerance_project);
		  return_error_code = constraintProject(new_position, surface_position, field, tolerance_project);
		  //	      cout << "befor min dist " << surface_position << pixel_position << endl;
		  if (!return_error_code)
		    return_error_code = minimizeDistance(new_position, pixel_position, new_position,
                                                         field, tolerance, tolerance_project);
		  //		  mc_distance = (new_position - pixel_position).length();
		  // 	      if (mc_distance < 2.0)
		  // 		{
		  // 		  return_error_code = constraintProject(tmp_position, pixel_position, pixel_position, field, tolerance_project);
		  // 		  return_error_code = minimizeDistance(tmp_position, pixel_position, tmp_position, field, tolerance, tolerance_project);
		  // 		  grad_des_distance = (tmp_position - pixel_position).length();
		  // 		  if (mc_distance > grad_des_distance)
		  // 		    {
		  // 		      //		      cout << "got a modification for distance" << endl;
		  // 		      //		      cout << pixel_position << endl;
		  // 		      //		      cout << "distances are " << mc_distance << " " << grad_des_distance << endl;
		  // 		      //		      cout << "new " << tmp_position << field->value(tmp_position) << endl;
		  // 		      //		      cout << "old " << new_position << field->value(new_position) << endl;
		  // 		      mc_distance = grad_des_distance;
		  // 		      new_position = tmp_position;
		  // 		    }
		  //		}
		  else
		    {
		      field->computeScalarFieldParams(new_position, params);
                      cout << i << " " << j << " " << k << endl;
		      cout << "got bad return on project - pixel: " << pixel_position[0] << " " <<
                        pixel_position[1] << " " << pixel_position[2] << " - surface: " <<
                        new_position[0] << " " << new_position[1] << " " << new_position[2] <<
                        " value "  << params._F << endl;
		      //		      constraintProject(new_position, surface_position, field, tolerance_project, true);
		      //		      constraintProject(new_position, surface_position, field, tolerance_project, true);
		    }

		  vol_position.poke(i, j, k) = new_position;
		  if (return_error_code)
		    {
		      cout << "got bad return on minimize: " << new_position[0] << " " <<
                        new_position[1] << " " << new_position[2] << " "<< endl;
		      //		      cout << "return " << minimizeDistance(new_position, pixel_position, new_position, field, tolerance, tolerance_project, true) << endl;
		    }
		}
		//	      else
		//		vol_position.poke(i, j, k) = pixel_position;
		

// 	      if (field->value(pixel_position) > 0.0) 
// 		vol.poke(i, j, k) = mc_distance;
// 	      else 
// 		vol.poke(i, j, k) = -1.0*mc_distance;
	    }
	}
    }

  cout << "distance min " << vol.min() << " and max " << vol.max() << endl;
  im_file.write((vol.image()).becomeFlat(), "vol_dist.fts");

  vec<3> surface_position2, pixel_position2, difference_vector2;
  vec<3> normal_dx, normal_dy, normal_dz;
  matrix<3, 3> normal_gradient;
  double flux, distance, x, y, z;
  int num_bad_project = 0, num_bad_minimize = 0;

  vec<3> intersection_position, init_position;
  vec<3> gradient2;
  vec<3> vec_tmp;
  double correction;
  num_pts = 0;
  for (k = (z_lo + 1); k < (z_hi - 1); k++)
    {
      for (j = (y_lo + 1); j < (y_hi - 1); j++)
	for (i = (x_lo + 1); i < (x_hi - 1 ); i++)
//   for (k = 1; k < (d-1); k++)
//     {
//       //      cout << "about to do slice " <<  k << endl;
//       for (j = 1; j < (h - 1) ; j++)
// 	for (i = 1; i < (w - 1); i++)
	  {
	    pixel_position = vec<3>(spacing[0]*(double)i,
                                    spacing[1]*(double)j,
                                    spacing[2]*(double)k);
	    //	    if (field->value(pixel_position) > 0.0) 
	      {
		surface_position = vol_position.peek(i, j, k);
		difference_vector = pixel_position - surface_position;
		distance = difference_vector.length(1.0e-20);
		gradient = difference_vector/distance;
		intersection_position = findBoxIntersection(pixel_position, gradient);
		x = intersection_position[0]; y = intersection_position[1]; z = intersection_position[2]; 
		if (vol_position.checkBounds((float)(x/spacing[0]),
                                             (float)(y/spacing[1]),
                                             (float)(z/spacing[2])))
		  {
		    distance = (intersection_position - surface_position).length();
		    init_position = (octree.closestPoint(&Point(intersection_position))).position(); 
		    //		    new_position = vol_position.interp(x, y, z); 
		    //		    return_error_code = constraintProject(new_position, intersection_position, init_position, field, tolerance_project);
		    return_error_code = constraintProject(new_position, init_position, field, tolerance_project);
		    //		    if (!return_error_code)
		    if (return_error_code)
		      num_bad_project++;
		    return_error_code = minimizeDistance(new_position, intersection_position, new_position, field, tolerance, tolerance_project);
		    if (return_error_code)
		      num_bad_minimize++;
// 		    else
// 		      {
// 			cout << "*** got bad error code on constrain in flux calc" << endl;
// 			cout << pixel_position << endl;
			
// 			field->computeScalarFieldParams(pixel_position, params);
// 			cout << "value " << params._F << endl;
// 			//			cout << "value " << field->value(pixel_position) << endl;
// 			cout << "intersection position" << intersection_position << endl;
// 			field->computeScalarFieldParams(intersection_position, params);
// 			cout << "intersection value " << params._F << endl;
// 			//			cout << "intersection value " << field->value(intersection_position) << endl;
// 			cout << "new position" << new_position << endl;
// 			field->computeScalarFieldParams(new_position, params);
// 			cout << "foot point value " << params._F << endl;
// 			field->computeScalarFieldParams(init_position, params);
// 			cout << "initial foot point value " << params._F << endl;
// 			cout << "trace " << params._F << endl;
// 			//			constraintProject(new_position, intersection_position, init_position, field, tolerance_project,true);
// 			constraintProject(new_position, init_position, field, tolerance_project,true);
// 		      }
		    gradient2 = (intersection_position - new_position);
		    new_distance = (gradient2).length(1.0e-20);
		    gradient2 /= new_distance;
		    if (DotProduct(gradient2, gradient) < CURVATURE_TOLERANCE)
		      {
			vol.poke(i, j, k)  = 1.0;
			//   		    cout << "pixel is " << pixel_position << endl;
			//		    cout << "surface_position is " << surface_position << endl;
			//		    vec_tmp = (pixel_position - vec<3>((double)(w/2), (double)(h/2), (double)(d/2))).normalize();
			//   		    cout << "adjusted pixel is " << vec_tmp << endl;
			//   		    cout << "difference is  " << gradient << endl;
			//		    cout << "field gradient is " << field->gradient(surface_position) << endl;
			//   		    cout << "difference 2 is  " << gradient2 << endl;
			//		    cout << "field gradient 2 is " << field->gradient(new_position) << endl;

			//		    cout << "dot_product of diff 1 is " << DotProduct(vec_tmp, gradient) << endl;
			//		    cout << "dot_product of field 1 is " << DotProduct(vec_tmp , (field->gradient(surface_position)).normalize()) << endl;
		    
			//   		    cout << "intersection_position is " << intersection_position << endl;
			//		    cout << "interception pt is " << intersection_position << endl;
			correction = (distance - new_distance)/(1.0 - DotProduct(gradient2, gradient));
			//			if (correction < 2.0*(pixel_position - intersection_position).length()) 
			  {
			    num_pts++;
			    vec_tmp = intersection_position - gradient*correction;;
			    os_pts <<  vec_tmp[0] <<  " " <<  vec_tmp[1] << " "  <<  vec_tmp[2] << " ";
			    vec_tmp = (gradient + gradient2).normalize(1.0e-20);
			    vec_tmp = CrossProduct(CrossProduct(gradient, gradient2), vec_tmp);
			    vec_tmp = vec_tmp.normalize(1.0e-20);
			    os_pts <<  vec_tmp[0] <<  " " <<  vec_tmp[1] << " "  <<  vec_tmp[2] << " " << "0.0" << endl;
			  }
			//   		    cout << "surface position " << surface_position << endl;
			//   		    cout << "new position " << new_position << endl;
			//   		    cout << "distances " << distance << " " << new_distance << endl;
			// 	    else
			// 	      {
			// 		if (flux > 2.0)
			// 		  {
			// 		    minimizeDistance(new_position, pixel_position, surface_position, field, tolerance, tolerance_project, true);
			// 		    cout << "pixel_pos " << pixel_position; 
			// 		    cout << "surface_pos " << surface_position; 
			// 		    cout << "distance " << vol.peek(i, j, k) << endl;
			// 		    cout << value << endl;
			// 		  }
			// 		vol.poke(i, j, k) = 0.0;
			// 	      }
		      }
		    else
		      vol.poke(i, j, k) = 0.0;
		  }
		//		else 
		//		  vol.poke(i, j, k) = 0.0;

		//		vol.poke(i, j, k) = distance - new_distance;
	      }
	  }
    }
  cout << "num pts " << num_pts << endl;
  cout << "percentage bad project " << (float)num_bad_project/(float)((x_hi - x_lo)*(y_hi - y_lo)*(z_hi - z_lo)) << endl;;
  cout << "percentage bad minimize " << (float)num_bad_minimize/(float)((x_hi - x_lo)*(y_hi - y_lo)*(z_hi - z_lo)) << endl;;
  vol.setBorder(0.0, 1);

#ifdef this_one_does_down_wind_gradient
  double regularize = 1.0e-20;
  for (k = 1; k < (d-1); k++)
    {
      cout << "about to do slice " <<  k << endl;
      for (j = 1; j < (h - 1) ; j++)
	for (i = 1; i < (w - 1); i++)
	  {
	    surface_position = vol_position.peek(i, j, k);
	    pixel_position = vec<3>(spacing[0]*(double)i,
                                    spacing[1]*(double)j,
                                    spacing[2]*(double)k);
	    value = field->value(pixel_position);
	    difference_vector = pixel_position - surface_position;
	    distance = difference_vector.length(regularize);
	    gradient = difference_vector/distance;

	    //	      if ((difference_vector2.length() == 0.0)||(difference_vector.length() == 0.0))
	    //		cout << "got zero distance to surface value is " << value << endl;

	    //dx
	    if (gradient[0] > 0.0)
	      {
		surface_position = vol_position.peek(i+1, j, k);
		pixel_position = vec<3>(spacing[0]*(double)(i+1),
                                        spacing[1]*(double)j,
                                        spacing[2]*(double)k);
		difference_vector2 = pixel_position - surface_position;
		normal_dx = difference_vector2.normalize(regularize) - gradient;
	      }
	    else
	      {
		surface_position = vol_position.peek(i-1, j, k);
		pixel_position = vec<3>(spacing[0]*(double)(i-1),
                                        spacing[1]*(double)j,
                                        spacing[2]*(double)k);
		difference_vector2 = pixel_position - surface_position;
		normal_dx = difference_vector.normalize(regularize) - gradient;
	      }

	    //dy
	    if (gradient[1] > 0.0)
	      {
		surface_position = vol_position.peek(i, j+1, k);
		pixel_position = vec<3>(spacing[0]*(double)i,
                                        spacing[1]*(double)(j+1),
                                        spacing[2]*(double)k);
		difference_vector2 = pixel_position - surface_position;
		normal_dy = difference_vector2.normalize(regularize) - gradient;
	      }
	    else
	      {
		surface_position = vol_position.peek(i, j-1, k);
		pixel_position = vec<3>(spacing[0]*(double)i,
                                        spacing[1]*(double)(j-1),
                                        spacing[2]*(double)k);
		difference_vector2 = pixel_position - surface_position;
		normal_dy = difference_vector.normalize(regularize) - gradient;
	      }

	    //dz
	    if (gradient[2] > 0.0)
	      {
		surface_position = vol_position.peek(i, j, k+1);
		pixel_position = vec<3>(spacing[0]*(double)i,
                                        spacing[1]*(double)j,
                                        spacing[2]*(double)(k+1));
		difference_vector2 = pixel_position - surface_position;
		normal_dz = difference_vector2.normalize(regularize) - gradient;
	      }
	    else
	      {
		surface_position = vol_position.peek(i, j, k-1);
		pixel_position = vec<3>(spacing[0]*(double)i,
                                        spacing[1]*(double)j,
                                        spacing[2]*(double)(k-1));
		difference_vector2 = pixel_position - surface_position;
		normal_dz = difference_vector.normalize(regularize) - gradient;
	      }

	    normal_gradient.set(normal_dx, normal_dy, normal_dz);
	    flux = (normal_gradient*gradient).length();	
 	    // if (flux > 2.0)
//  	      {
//  		cout << "pixel_pos " << pixel_position; 
//  		cout << "surface_pos " << surface_position; 
//  		cout << value << endl;
//  	      }

// need to threshold distance, because it appears that too close to the surface is unreliable.
//	    if (distance > 0.2
//	    if (fabs(value) > 0.001)
	      vol.poke(i, j, k) = flux;
// 	    else
// 	      {
// 		if (flux > 2.0)
// 		  {
// 		    minimizeDistance(new_position, pixel_position, surface_position, field, tolerance, tolerance_project, true);
// 		    cout << "pixel_pos " << pixel_position; 
// 		    cout << "surface_pos " << surface_position; 
// 		    cout << "distance " << vol.peek(i, j, k) << endl;
// 		    cout << value << endl;
// 		  }
// 		vol.poke(i, j, k) = 0.0;
// 	      }
	  }
    }
  vol.setBorder(0.0, 1);
#endif

  //      volume_file.write_double(vol_out, "dist.vol");
  cout << "vol flux min " << vol.min() << " and max " << vol.max() << endl;
  im_file.write((vol.image()).becomeFlat(), outname_vol);

  // clean up fb
  //  if (bf != std::cout.rdbuf()) delete fb;

  cout << "done " << endl;
  
}


