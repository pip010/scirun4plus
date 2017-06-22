#include <morphology.h>
#include <volumefile.h>

#define MAX_REINIT_ITERATIONS (50)
// this is the value you need for the OLD reinit
//#define REINIT_TOLERANCE (1.0e-3)
// this is the value you need for the NEW reinit
#define REINIT_TOLERANCE (2.0e-3)
//#define REINIT_TOLERANCE (1.0e-2)
#define MINMAX (1.0f)

// discrete - BINARY morphology

/*
void
gaussDiffuse(VolumeScalar &vol, double sigma)
{
  double max_dt = 1.0/12.0, time = sigma*sigma/2.0, extra_time;
  int iterations = (int)floor(time/max_dt);
  extra_time = time-iterations*max_dt;
  double laplacian = 0.0;
  int w = vol.width(), h = vol.height(), d = vol.depth();
  VolumeScalar update(w, h, d);
  for (int l = 0; l < iterations+1; l++)
  {
    cout << "gaussDiffuse about to do iteration " << l << endl;
    for (int k = 0; k < d; k++)
      for (int j = 0; j < h; j++)
	      for (int i = 0; i < w; i++)
	      {
	         update.poke(i, j, k) = max_dt*(-6.0*vol.peek(i, j, k) +
					 vol.peek(VISmax(i - 1, 0), j, k) +
					 vol.peek(VISmin(i + 1, w - 1), j, k) +
					 vol.peek(i, VISmax(j - 1, 0), k) +
					 vol.peek(i, VISmin(j + 1, h - 1), k) +
					 vol.peek(i, j, VISmax(k - 1, 0)) +
					 vol.peek(i, j, VISmin(k + 1, d - 1)));
	     }
    if (l < iterations)
      vol += update;
    else
      vol += (float)(extra_time/max_dt)*update;
  }
}
*/

VolumeScalar
erode(const VolumeScalar& vol)
{
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

  int w = vol.width(), h = vol.height(), d = vol.depth();
  VolumeScalar ret = vol;
  int i, j, k, ii, jj, kk, l;

  for (k = 0; k < d; k++)
    for (j = 0; j < h; j++)
      for (i = 0; i < w; i++)
      {
	if (vol.peek(i, j, k) == 1.0f)
	{
	  for (l = 0; (l < 6); l++)
	  {
	    ii = i + N[l][0];
	    jj = j + N[l][1];
	    kk = k + N[l][2];
	    if ((vol.checkBounds(ii, jj, kk))&&
		(vol.peek(ii, jj, kk) == 0.0f))
	    {
	      ret.poke(i, j, k) = 0.0f;
	    }
	  }
	}
      }
  return(ret);
}

VolumeScalar
dilate(const VolumeScalar& vol)
{
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

  int w = vol.width(), h = vol.height(), d = vol.depth();
  VolumeScalar ret = vol;
  int i, j, k, ii, jj, kk, l;

  for (k = 0; k < d; k++)
    for (j = 0; j < h; j++)
      for (i = 0; i < w; i++)
      {
	if (vol.peek(i, j, k) == 0.0f)
	{
	  for (l = 0; (l < 6); l++)
	  {
	    ii = i + N[l][0];
	    jj = j + N[l][1];
	    kk = k + N[l][2];
	    if ((vol.checkBounds(ii, jj, kk))&&
		(vol.peek(ii, jj, kk) == 1.0f))
	    {
	      ret.poke(i, j, k) = 1.0f;
	    }
	  }
	}
      }
  return(ret);
}

VolumeScalar
open(const VolumeScalar& vol, int size)
{
  int i;
  VolumeScalar ret = vol;
  for (i = 0; i < size; i++)
    ret = erode(ret);
  for (i = 0; i < size; i++)
    ret = dilate(ret);
  return(ret);
}

VolumeScalar
close(const VolumeScalar& vol, int size)
{
  int i;
  VolumeScalar ret = vol;
  for (i = 0; i < size; i++)
    ret = dilate(ret);
  for (i = 0; i < size; i++)
    ret = erode(ret);
  return(ret);
}


// this does an advect by shifting the grey level.
// Expects that the function looks like signed distance for values < minmax
void
advect_new(VolumeScalar& vol, float distance, float target_grad, float minmax)
{
  std::cout << "advect-new..." << std::endl;
  // how much do you need to move the intensities, total...
  // how much can you move in each iteration...

  float d_update = VISmin<float>(target_grad*fabs(distance), minmax/2.0f);
  float d_distance = d_update/target_grad;
  int iterations = (int)ceil(fabs(distance)/d_distance);
  d_update = target_grad*distance/iterations;
  float init_update_max;
  float this_value;
  VISVolume<boolean> mask;

  int w, h, d;
  int num_updates;
  VISImageFile im_file;

  cout << "new advect update amount " << d_update
       << " num iterations " << iterations << endl;
  cout << "vol in min " << vol.min() << " and max " << vol.max() << endl;

  for (int l = 0; l < iterations; l++)
  {
    vol += d_update;
    if (d_update < 0.0)
      vol = vol.max(-minmax);
    else
      vol = vol.min(minmax);
    reinitDistance(vol, target_grad, minmax, REINIT_TOLERANCE,
		   MAX_REINIT_ITERATIONS);
  }
}



void
erode_grey(VISVolume<float> &volume, float t, float target_grad = -1.0)
{
  advect_new(volume, -t, target_grad, 1.0);
}

void
dilate_grey(VISVolume<float> &volume, float t, float target_grad = -1.0)
{
  advect_new(volume, t, target_grad, 1.0);
}


void
open_grey(VISVolume<float> &volume, float t, float target_grad = -1.0, bool writeIntermediateFiles = false)
{

  cout << "holyyyy sh!!!!!t" << endl;

  VISVolumeFile vol_file;
  VISImageFile im_file;
  //   int iterations = 0;
  //   float init_update_max = FLT_MAX;
  //   if (target_grad > 0.0)
  //     while ((init_update_max > 1.0e-3)&&(iterations++ < MAX_REINIT_ITERATIONS))
  //       {
  // 	init_update_max = reinitDistance(volume, target_grad);
  // 	cout << "open update max is " << init_update_max << endl;
  //       }
  cout << "about to erode in open" << endl;
  erode_grey(volume, t, target_grad);
  //  erode_grey(volume, t);
  if (writeIntermediateFiles)
  {
    //  vol_file.write_float(volume, "erode_open.vol");
    im_file.write((volume.image()).becomeFlat(), "vol_erode.fts");
  }
  cout << "about to dilate in open" << endl;
  dilate_grey(volume, t, target_grad);
  if (writeIntermediateFiles)
  {
    im_file.write((volume.image()).becomeFlat(), "vol_dilate.fts");
  }
  //  dilate_grey(volume, t);
  //  vol_file.write_float(volume, "dilate_open.vol");
  //  cout << "done open" << endl;
}

void
close_grey(VISVolume<float> &volume, float t, float target_grad = -1.0, bool writeIntermediateFiles = false)
{
  VISVolumeFile vol_file;
  dilate_grey(volume, t, target_grad);
  erode_grey(volume, t, target_grad);
  if (writeIntermediateFiles)
  {
    vol_file.write_float(volume, "close.vol");
  }
}
// this does an advect by shifting the grey level.
// Expects that the function looks like signed distance for values < minmax
void advect_new_vol(VolumeScalar& vol, VolumeScalar& vol_t, float distance, float target_grad, float minmax)
{
  std::cout << "advect-new-vol..." << std::endl;
  // how much do you need to move the intensities, total...
  // how much can you move in each iteration...

  float d_update = VISmin<float>(target_grad*fabs(distance), minmax/2.0f);
  float d_distance = d_update/target_grad;
  int iterations = (int)ceil(fabs(distance)/d_distance);
  d_update = target_grad*distance/iterations;
  float init_update_max;
  float this_value;
  VISVolume<boolean> mask;

  int w, h, d;
  int num_updates;
  VISImageFile im_file;

  cout << "new advect update amount " << d_update << " num iterations " << iterations << endl;
  cout << "vol in min " << vol.min() << " and max " << vol.max() << endl;

  for (int l = 0; l < iterations; l++)
  {
    vol += d_update;

    if (d_update < 0.0)
      vol = vol.max(-minmax);
    else
      vol = vol.min(minmax);

    reinitDistance(vol, vol_t, target_grad, minmax, REINIT_TOLERANCE,  MAX_REINIT_ITERATIONS);
  }
}

float
reinitDistance(VISVolume<float>& vol, VISVolume<float>& vol_t, float target_grad, float minmax,
	       float reinit_tolerance, int max_reinit_iterations)
{
  float init_update_max = FLT_MAX;
  int iterations = 0;
  VISVolIndexVISList mask_list;

  cout << "about to create mask " << endl;
  createMask(vol, mask_list, target_grad);

  while ((init_update_max > reinit_tolerance) &&
	 (iterations++ < max_reinit_iterations))
  {
    init_update_max = reinitDistance(vol, vol_t, target_grad, MINMAX, mask_list);
    cout << "reinit update max is " << init_update_max << endl;
  }
  return (init_update_max);
}

float reinitDistance(VISVolume<float>& vol, VISVolume<float>& vol_t ,float target_grad, float minmax, VISVolIndexVISList &mask_list)
{
  int w = vol.width(), h = vol.height(), d = vol.depth();
  VISVolume<float> update = vol.createToSize();
  float grad_forward, grad_back, grad_x, grad_y, grad_z, speed;
  float max_speed, max_update;
  float this_value, this_update;
  int i, j, k;
  VISVolIndex index;
  float dt;

  max_speed = 0.0f;
  max_update = 0.0f;
  mask_list.reset();
  float total = 0.0f;
  int number = 0;

  if (!mask_list.valid())
    cout << "reinitDistance got invalid mask" << endl;

  while (mask_list.valid())
  {
    index = mask_list.atCurrent();
    i = index.a(); j = index.b(); k = index.c();
    this_value = vol.peek(i, j, k);

    grad_forward = vol.peek(VISmin(i+1, w-1), j, k) - vol.peek(i, j, k);
    grad_back = vol.peek(i, j, k) - vol.peek(VISmax(i-1, 0), j, k);

    if (this_value > 0)
      grad_x = power(VISmin(grad_forward, 0.0f), 2) + power(VISmax(grad_back, 0.0f), 2);
    else
      grad_x = power(VISmax(grad_forward, 0.0f), 2) + power(VISmin(grad_back, 0.0f), 2);

    grad_forward = vol.peek(i, VISmin(j+1, h-1), k) - vol.peek(i, j, k);
    grad_back = vol.peek(i, j, k) - vol.peek(i, VISmax(j-1, 0), k);

    if (this_value > 0)
      grad_y = power(VISmin(grad_forward, 0.0f), 2) + power(VISmax(grad_back, 0.0f), 2);
    else
      grad_y = power(VISmax(grad_forward, 0.0f), 2) + power(VISmin(grad_back, 0.0f), 2);

    grad_forward = vol.peek(i, j, VISmin(k+1, d-1)) - vol.peek(i, j, k);
    grad_back = vol.peek(i, j, k) - vol.peek(i, j, VISmax(k-1, 0));

    if (this_value > 0)
      grad_z = power(VISmin(grad_forward, 0.0f), 2) + power(VISmax(grad_back, 0.0f), 2);
    else
      grad_z = power(VISmax(grad_forward, 0.0f), 2) + power(VISmin(grad_back, 0.0f), 2);


    this_update = this_value*(target_grad - sqrt(grad_x + grad_y + grad_z));

    if (this_update > 0.0f)
      this_update = VISmin(VISmin((this_value + this_update), minmax) - this_value, this_update);
    else
      this_update = VISmax(VISmax((this_value + this_update), -minmax) - this_value, this_update);

    update.poke(i, j, k) = this_update;
    max_update = VISmax((float)fabs(max_update), this_update);
    max_speed = VISmax((float)fabs(this_value), max_speed);
    total += fabs(this_update);
    number++;
    mask_list.stepForward();
  } // look over mask_list

  if (max_speed > 0.0f)
    dt = (1.0f/(6.0f*max_speed));
  else
    dt = 0.0f;
  dt = VISmin(dt, 0.5f/target_grad);
  cout << "reinite dt is " << dt << endl;

  mask_list.reset();
  while (mask_list.valid())
  {
    index = mask_list.atCurrent();
    i = index.a(); j = index.b(); k = index.c();
    vol.poke(i, j, k) += dt*update.peek(i, j, k);
    mask_list.stepForward();
  }

  cout << "reinit max update is " << max_update << endl;
  cout << "reinit avg update is " << total/(float)number << endl;

  //  cout << "reinitdistance min/max are " << vol.min() << " " << vol.max() << endl;
  //  return(max_update);
  // ****  return average update
  return(total/(float)number);
}

void erode_grey(VolumeScalar& volume, VolumeScalar& volt, float t , float target_grad)
{
  advect_new_vol(volume, volt, -t, target_grad, 1.0);
}
void dilate_grey(VolumeScalar& volume, VolumeScalar& volt, float t , float target_grad)
{
  advect_new_vol(volume, volt, t, target_grad, 1.0);
}
void open_grey(VolumeScalar& volume, VolumeScalar& volt, float t ,  float target_grad, bool writeIntermediateFiles)
{
    cout << "holyyyy sh!!!!!t   - new" << endl;

    VISVolumeFile vol_file;
    VISImageFile im_file;

    cout << "about to erode in open" << endl;

    erode_grey(volume,volt, t, target_grad);

    if (writeIntermediateFiles)
    {
      im_file.write((volume.image()).becomeFlat(), "vol_erode.fts");
    }

    cout << "about to dilate in open" << endl;

    dilate_grey(volume,volt, t, target_grad);

    if (writeIntermediateFiles)
    {
      im_file.write((volume.image()).becomeFlat(), "vol_dilate.fts");
    }

}
void close_grey(VolumeScalar& volume, VolumeScalar& volt, float t ,  float target_grad, bool writeIntermediateFiles)
{
  VISVolumeFile vol_file;
  dilate_grey(volume, volt, t, target_grad);
  erode_grey(volume, volt, t, target_grad);
  if (writeIntermediateFiles)
  {
    vol_file.write_float(volume, "close.vol");
  }
}




void
get_derivs (const float* n, float* data)
{
  data[DX] = 0.5f*(n[14] - n[12]);
  data[DY] = 0.5f*(n[16] - n[10]);
  data[DZ] = 0.5f*(n[22] - n[4]);
  data[DPX] = n[14] - n[13];
  data[DPY] = n[16] - n[13];
  data[DPZ] = n[22] - n[13];
  data[DMX] = n[13] - n[12];
  data[DMY] = n[13] - n[10];
  data[DMZ] = n[13] - n[4];
  data[DYPX] = 0.5f*(n[17] - n[11]);
  data[DZPX] = 0.5f*(n[23] - n[5]);
  data[DYMX] = 0.5f*(n[15] - n[9]);
  data[DZMX] = 0.5f*(n[21] - n[3]);
  data[DXPY] = 0.5f*(n[17] - n[15]);
  data[DZPY] = 0.5f*(n[25] - n[7]);
  data[DXMY] = 0.5f*(n[11] - n[9]);
  data[DZMY] = 0.5f*(n[19] - n[1]);
  data[DXPZ] = 0.5f*(n[23] - n[21]);
  data[DYPZ] = 0.5f*(n[25] - n[19]);
  data[DXMZ] = 0.5f*(n[5] - n[3]);
  data[DYMZ] = 0.5f*(n[7] - n[1]);
}

float
reinitDistance(VISImage<float>& image, float target_grad)
{
  int w = image.width(), h = image.height();
  VISImage<float> update = image.createToSize();
  float grad_forward, grad_back, grad_x, grad_y, speed;
  float max_speed, max_update;
  float this_value, this_update;
  int i, j;

  max_speed = 0.0f;
  max_update = 0.0f;
  for (j = 0; j < h; j++)
    for (i = 0; i < w; i++)
    {
      this_value = image.peek(i, j);
      //
      grad_forward = image.peek(VISmin(i+1, w-1), j) - image.peek(i, j);
      grad_back = image.peek(i, j) - image.peek(VISmax(i-1, 0), j);
      if (this_value > 0)
	grad_x = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
      else
	grad_x = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

      grad_forward = image.peek(i, VISmin(j+1, h-1)) - image.peek(i, j);
      grad_back = image.peek(i, j) - image.peek(i, VISmax(j-1, 0), 0);
      if (this_value > 0)
	grad_y = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
      else
	grad_y = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

      speed = this_value*(target_grad - sqrt(grad_y*grad_y + grad_x*grad_x));

      grad_forward = image.peek(VISmin(i+1, w-1), j) - image.peek(i, j);
      grad_back = image.peek(i, j) - image.peek(VISmax(i-1, 0), j);
      if (speed < 0.0f)
	grad_x = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
      else
	grad_x = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

      grad_forward = image.peek(i, VISmin(j+1, h-1)) - image.peek(i, j);
      grad_back = image.peek(i, j) - image.peek(i, VISmax(j-1, 0), 0);
      if (speed < 0)
	grad_y = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
      else
	grad_y = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

      this_update = speed*sqrt(grad_y*grad_y + grad_x*grad_x);
      update.poke(i, j) = this_update;

      max_update = VISmax((float)fabs(max_update), this_update);
      max_speed = VISmax((float)fabs(speed), max_speed);
    }

  if (max_speed > 0.0f)
    image += (1.0f/(4.0f*max_speed))*update;
  return(max_update);
}



float
reinitDistance(VISVolume<float>& vol, float target_grad, float minmax,
	       const VISVolume<boolean> &mask, float reinit_tolerance,
	       int max_reinit_iterations)
{
  float init_update_max = FLT_MAX;
  int iterations = 0;
  while ((init_update_max > reinit_tolerance) &&
	 (iterations++ < max_reinit_iterations))
  {
    init_update_max = reinitDistance(vol, target_grad, MINMAX, mask);
    //	init_update_max = reinitDistanceOld(vol, target_grad, mask);
    cout << "reinit update max is " << init_update_max << endl;
  }
  return (init_update_max);
}

float
reinitDistance(VISVolume<float>& vol, float target_grad, float minmax,
	       float reinit_tolerance, int max_reinit_iterations)
{
  float init_update_max = FLT_MAX;
  int iterations = 0;
  VISVolIndexVISList mask_list;

  cout << "about to create mask " << endl;
  createMask(vol, mask_list, target_grad);

  while ((init_update_max > reinit_tolerance) &&
	 (iterations++ < max_reinit_iterations))
  {
    init_update_max = reinitDistance(vol, target_grad, MINMAX, mask_list);
    //	init_update_max = reinitDistanceOld(vol, target_grad);
    cout << "reinit update max is " << init_update_max << endl;
  }
  return (init_update_max);
}

// creates a list of voxels that form a band around the zero set
// assumes a -1 to 1 spread
// sets all stuff outside that band to -1 or 1.
// changes
void
createMask(VISVolume<float>& vol, VISVolIndexVISList &list, float gradient)
{
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
  int w = vol.width(), h = vol.height(), d = vol.depth();
  int i, j, k;
  int ii, jj, kk, ll;
  float this_value, next_value;
  int num_layers = static_cast<int>(ceil(sqrt(3.0)/gradient) + 2);

  VISVolume<boolean> vol_mask(w, h, d);
  vol_mask = false;

  VISVolIndexVISList list_transition, list_new;
  list_new.clean(); list.clean(), list_transition;
  for (k = 0; k < (d - 1); k++)
    for (j = 0; j < (h - 1); j++)
      for (i = 0; i < (w - 1); i++)
      {
	for (ll = 0; ll < 6; ll+=2)
	{
	  ii = i + N[ll][0];
	  jj = j + N[ll][1];
	  kk = k+ N[ll][2];
	  this_value = vol.peek(i, j, k);
	  next_value = vol.peek(ii, jj, kk);
	  if ((this_value*next_value <= 0.0) &&
	      (!((next_value == 0.0f) && (this_value == 0.0f))))
	  {
	    if (!vol_mask.peek(i, j, k))
	    {
	      list_new.appendItem(VISVolIndex(i, j, k));
	      vol_mask.poke(i, j, k) = true;
	    }
	    if (!vol_mask.peek(ii, jj, kk))
	    {
	      list_new.appendItem(VISVolIndex(ii, jj, kk));
	      vol_mask.poke(ii, jj, kk) = true;
	    }
	  }
	}
      }

  VISVolIndex index;
  for (int iter = 0; iter < num_layers; iter++)
  {
    list_new.reset();
    while (list_new.valid())
    {
      list_transition.appendItem(list_new.atCurrent());
      list_new.removeCurrent();
    }
    list_transition.reset();
    while (list_transition.valid())
    {
      index = list_transition.atCurrent();
      list.appendItem(index);
      i = index.a();
      j = index.b();
      k = index.c();
      //	    printf("flood fill iteration %d x %d y %d z %d\n", n++, at_x,
      //		   at_y , at_z);
      for (ll = 0; ll < 6; ll++)
      {
	ii = i + N[ll][0];
	jj = j + N[ll][1];
	kk = k+ N[ll][2];
	if (vol.checkBounds(ii, jj, kk)&&(!vol_mask.peek(ii, jj, kk)))
	{
	  list_new.appendItem(VISVolIndex(ii, jj, kk));
	  vol_mask.poke(ii, jj, kk) = true;
	}
      }
      list_transition.removeCurrent();
    }
  }

  list_new.reset();
  while (list_new.valid())
  {
    list.appendItem(list_new.atCurrent());
    list_new.removeCurrent();
  }

  for (k = 0; k < d; k++)
    for (j = 0; j < h; j++)
      for (i = 0; i < w; i++)
      {
	if (!vol_mask.peek(i, j, k))
	  if (vol.peek(i, j, k) > 0.0f)
	    vol.poke(i, j, k) = 1.0f;
	  else
	    vol.poke(i, j, k) = -1.0f;
      }
}

float
reinitDistance(VISVolume<float>& vol, float target_grad, float minmax,
	       VISVolIndexVISList &mask_list)
{
  int w = vol.width(), h = vol.height(), d = vol.depth();
  VISVolume<float> update = vol.createToSize();
  float grad_forward, grad_back, grad_x, grad_y, grad_z, speed;
  float max_speed, max_update;
  float this_value, this_update;
  int i, j, k;
  VISVolIndex index;
  float dt;

  max_speed = 0.0f;
  max_update = 0.0f;
  mask_list.reset();
  float total = 0.0f;
  int number = 0;
  if (!mask_list.valid())
    cout << "reinitDistance got invalid mask" << endl;
  while (mask_list.valid())
  {
    index = mask_list.atCurrent();
    i = index.a(); j = index.b(); k = index.c();
    this_value = vol.peek(i, j, k);
    //
    grad_forward = vol.peek(VISmin(i+1, w-1), j, k) - vol.peek(i, j, k);
    grad_back = vol.peek(i, j, k) - vol.peek(VISmax(i-1, 0), j, k);
    if (this_value > 0)
      //	grad_x = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
      grad_x = power(VISmin(grad_forward, 0.0f), 2) + power(VISmax(grad_back, 0.0f), 2);
    else
      //	grad_x = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));
      grad_x = power(VISmax(grad_forward, 0.0f), 2) + power(VISmin(grad_back, 0.0f), 2);

    grad_forward = vol.peek(i, VISmin(j+1, h-1), k) - vol.peek(i, j, k);
    grad_back = vol.peek(i, j, k) - vol.peek(i, VISmax(j-1, 0), k);
    if (this_value > 0)
      //	grad_y = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
      grad_y = power(VISmin(grad_forward, 0.0f), 2) + power(VISmax(grad_back, 0.0f), 2);
    else
      //	grad_y = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));
      grad_y = power(VISmax(grad_forward, 0.0f), 2) + power(VISmin(grad_back, 0.0f), 2);

    grad_forward = vol.peek(i, j, VISmin(k+1, d-1)) - vol.peek(i, j, k);
    grad_back = vol.peek(i, j, k) - vol.peek(i, j, VISmax(k-1, 0));
    if (this_value > 0)
      //	grad_z = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
      grad_z = power(VISmin(grad_forward, 0.0f), 2) + power(VISmax(grad_back, 0.0f), 2);
    else
      //	grad_z = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));
      grad_z = power(VISmax(grad_forward, 0.0f), 2) + power(VISmin(grad_back, 0.0f), 2);

    //      this_update = this_value*(target_grad - sqrt(grad_y*grad_y + grad_x*grad_x + grad_z*grad_z));
    this_update = this_value*(target_grad - sqrt(grad_x + grad_y + grad_z));

    if (this_update > 0.0f)
      this_update = VISmin(VISmin((this_value + this_update), minmax) - this_value, this_update);
    else
      this_update = VISmax(VISmax((this_value + this_update), -minmax) - this_value, this_update);

    update.poke(i, j, k) = this_update;
    max_update = VISmax((float)fabs(max_update), this_update);
    max_speed = VISmax((float)fabs(this_value), max_speed);
    total += fabs(this_update);
    number++;
    mask_list.stepForward();
  } // look over mask_list

  if (max_speed > 0.0f)
    dt = (1.0f/(6.0f*max_speed));
  else
    dt = 0.0f;
  dt = VISmin(dt, 0.5f/target_grad);
  cout << "reinite dt is " << dt << endl;

  mask_list.reset();
  while (mask_list.valid())
  {
    index = mask_list.atCurrent();
    i = index.a(); j = index.b(); k = index.c();
    vol.poke(i, j, k) += dt*update.peek(i, j, k);
    mask_list.stepForward();
  }

  cout << "reinit max update is " << max_update << endl;
  cout << "reinit avg update is " << total/(float)number << endl;

  //  cout << "reinitdistance min/max are " << vol.min() << " " << vol.max() << endl;
  //  return(max_update);
  // ****  return average update
  return(total/(float)number);
} // reinitdistance

float
reinitDistance(VISVolume<float>& vol, float target_grad, float minmax,
	       const VISVolume<boolean> &mask)
{
  int w = vol.width(), h = vol.height(), d = vol.depth();
  VISVolume<float> update = vol.createToSize();
  float grad_forward, grad_back, grad_x, grad_y, grad_z, speed;
  float max_speed, max_update;
  float this_value, this_update;
  int i, j, k;
  float dt;
  boolean mask_invalid = !mask.isValid();

  max_speed = 0.0f;
  max_update = 0.0f;
  for (k = 0; k < d; k++)
    for (j = 0; j < h; j++)
      for (i = 0; i < w; i++)
	if ((mask_invalid)||(mask.peek(i, j, k)))
	{
	  this_value = vol.peek(i, j, k);
	  //
	  grad_forward = vol.peek(VISmin(i+1, w-1), j, k) - vol.peek(i, j, k);
	  grad_back = vol.peek(i, j, k) - vol.peek(VISmax(i-1, 0), j, k);
	  if (this_value > 0)
	    grad_x = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
	  else
	    grad_x = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

	  grad_forward = vol.peek(i, VISmin(j+1, h-1), k) - vol.peek(i, j, k);
	  grad_back = vol.peek(i, j, k) - vol.peek(i, VISmax(j-1, 0), k);
	  if (this_value > 0)
	    grad_y = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
	  else
	    grad_y = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

	  grad_forward = vol.peek(i, j, VISmin(k+1, d-1)) - vol.peek(i, j, k);
	  grad_back = vol.peek(i, j, k) - vol.peek(i, j, VISmax(k-1, 0));
	  if (this_value > 0)
	    grad_z = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
	  else
	    grad_z = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

	  this_update = this_value*(target_grad - sqrt(grad_y*grad_y + grad_x*grad_x + grad_z*grad_z));

	  if (this_update > 0.0f)
	    this_update = VISmin(VISmin((this_value + this_update), minmax) - this_value, this_update);
	  else
	    this_update = VISmax(VISmax((this_value + this_update), -minmax) - this_value, this_update);

	  update.poke(i, j, k) = this_update;
	  max_update = VISmax((float)fabs(max_update), this_update);
	  max_speed = VISmax((float)fabs(this_value), max_speed);
	} // if mask
	else
	  update.poke(i, j, k) = 0.0f;

  if (max_speed > 0.0f)
    dt = (1.0f/(8.0f*max_speed));
  else
    dt = 0.0f;
  dt = VISmin(dt, 0.5f/target_grad);
  cout << "reinite dt is " << dt << endl;
  vol += dt*update;

  cout << "reinit max update is " << max_update << endl;

  //  cout << "reinitdistance min/max are " << vol.min() << " " << vol.max() << endl;
  return(max_update);
} // reinitdistance






 inline VISMatrix calculate_curvature (const float *derivatives,
						float &curve_trace, float &curve_norm) {
  VISMatrix curve(3, 3);
  VISVector Nf(3), Nb(3), normal(3);
  float MIN_NORM = 1.0e-10;

  // do the x derivative
  Nf.poke(0) = derivatives[DPX];
  Nf.poke(1) = 0.5f*(derivatives[DYPX] + derivatives[DY]);
  Nf.poke(2) = 0.5f*(derivatives[DZPX] + derivatives[DZ]);
  Nf = Nf/(Nf.norm()+MIN_NORM);
  Nb.poke(0) = derivatives[DMX];
  Nb.poke(1) = 0.5f*(derivatives[DYMX] + derivatives[DY]);
  Nb.poke(2) = 0.5f*(derivatives[DZMX] + derivatives[DZ]);
  Nb = Nb/(Nb.norm()+MIN_NORM);
  curve.pokeROI(0, 0, Nf - Nb);

  // do the y derivative
  Nf.poke(0) = 0.5f*(derivatives[DXPY] + derivatives[DX]);
  Nf.poke(1) = derivatives[DPY];
  Nf.poke(2) = 0.5f*(derivatives[DZPY] + derivatives[DZ]);
  Nf = Nf/(Nf.norm()+MIN_NORM);
  Nb.poke(0) = 0.5f*(derivatives[DXMY] + derivatives[DX]);
  Nb.poke(1) = derivatives[DMY];
  Nb.poke(2) = 0.5f*(derivatives[DZMY] + derivatives[DZ]);
  Nb = Nb/(Nb.norm()+MIN_NORM);
  curve.pokeROI(0, 1, Nf - Nb);

  // do the z derivative
  Nf.poke(0) = 0.5f*(derivatives[DXPZ] + derivatives[DX]);
  Nf.poke(1) = 0.5f*(derivatives[DYPZ] + derivatives[DY]);
  Nf.poke(2) = derivatives[DPZ];
  Nf = Nf/(Nf.norm()+MIN_NORM);
  Nb.poke(0) = 0.5f*(derivatives[DXMZ] + derivatives[DX]);
  Nb.poke(1) = 0.5f*(derivatives[DYMZ] + derivatives[DY]);
  Nb.poke(2) = derivatives[DMZ];
  Nb = Nb/(Nb.norm()+MIN_NORM);
  curve.pokeROI(0, 2, Nf - Nb);

  curve_trace = curve.peek(0, 0) + curve.peek(1, 1) + curve.peek(2, 2);
  curve_norm = 0.0f;
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
    curve_norm += pow(curve.peek(i, j), 2);
  return(curve);
}

void get_neighborhood (int x, int y, int z, float* data, const VISVolume<float> &values) {
  int i, j, k, l = 0;

  for (k = -1; k <= 1; k++) for (j = -1; j <= 1; j++) for (i = -1; i <= 1; i++)
    data[l++] = values.peek(VISmax(VISmin(x + i, (int)values.width() - 1), 0),
			    VISmax(VISmin(y + j, (int)values.height() - 1), 0),
			    VISmax(VISmin(z + k, (int)values.depth() - 1), 0));
}

void antialias(VolumeScalar& vol, float epsilon)
{
  float neighborhood[NSIZE], derivs[NUM_DERIVS];
  int i, j, k, ii;
  int w = vol.width(), h = vol.height(), d = vol.depth();
  VolumeScalar update = vol.createToSize();
  float curve_mean, curve_norm;
  float grad, grad_p, grad_m;
  float value, max_change = FLT_MAX;
  float dt = 1.0/6.0;
  int iterations = 0;

    while ((iterations++ < 100)&&(max_change > 1.0e-6))
    {
      max_change = 0.0;
      for (k = 0; k < d; k++)
	for (j = 0; j < h; j++)
	  for (i = 0; i < w; i++)
	    {
	      get_neighborhood (i, j, k, neighborhood, vol);
	      get_derivs(neighborhood, derivs);
	      calculate_curvature(derivs, curve_mean, curve_norm);
	      grad = 0.0f;
	      if (curve_mean > 0.0f)
		for (ii = 0; ii < 3; ii++)
		  {
		    grad_p = VISmax(derivs[DPX + ii], 0.0f);
		    grad_m = VISmin(derivs[DMX + ii], 0.0f);
		    grad += grad_p*grad_p + grad_m*grad_m;
		  }
	      else if (curve_mean < 0.0f)
		{
		  for (ii = 0; ii < 3; ii++)
		    {
		      grad_p = VISmin(derivs[DPX + ii], 0.0f);
		      grad_m = VISmax(derivs[DMX + ii], 0.0f);
		      grad += grad_p*grad_p + grad_m*grad_m;
		    }
		}
	      update.poke(i, j, k) = grad*curve_mean;
	    }
      for (k = 0; k < d; k++)
	for (j = 0; j < h; j++)
	  for (i = 0; i < w; i++)
	    {
	      value = vol.peek(i, j, k);
	      if (value > 0.0)
		vol.poke(i, j, k) = VISmax(value + dt*update.peek(i, j, k), epsilon);
	      else
		vol.poke(i, j, k) = VISmin(value + dt*update.peek(i, j, k), -epsilon);
	      max_change = VISmax((float)fabs(vol.peek(i, j, k) - value), max_change);
	    }
    }
}



void
clamp_curvature(VolumeScalar& vol, float epsilon)
{
  float neighborhood[NSIZE], derivs[NUM_DERIVS];
  int ii;
  VISMatrix projection;
  int w = vol.width(), h = vol.height(), d = vol.depth();
  VolumeScalar update = vol.createToSize();
  float curve_trace = 0.0f, curve_norm = 0.0f, k1 = 0.0f, k2 = 0.0f;
  float grad = 0.0f, grad_p = 0.0f, grad_m = 0.0f;
  float value = 0.0f, max_update = FLT_MAX;
  float avg_update = 0.0f;
  float dt = 1.0f/8.0f;
  int iterations = 0;
  float curve_principle[2], flow = 0.0f;
  float tmp = 0.0f;
  float dxx = 0.0f, dyy = 0.0f, dzz = 0.0f, dxz = 0.0f, dyz = 0.0f, dxy = 0.0f;
  VISVector normal(3);
  float grad_mag = 0.0f;
  VISMatrix curv_adjust(3,3);

  float grad_forward, grad_back, grad_x, grad_y, grad_z, speed;
  float target_grad = 0.25;
  float max_speed, reinit_weight = 0.25;
  float this_value = 0.0f;
  float this_dt = 0.0f;

  while ((iterations++ < 200)&&(max_update > 1.0e-6))
  {
    max_speed = 0.0f;
    for (int j = 0; j < h; j++)
      for (int i = 0; i < w; i++)
	      for (int k = 0; k < d; k++)
	      {
	          get_neighborhood (i, j, k, neighborhood, vol);
	          get_derivs(neighborhood, derivs);

	          //redistancing

	          this_value = vol.peek(i, j, k);
	          //
	          grad_forward = derivs[DPX];
	          grad_back = derivs[DMX];
	          if (this_value > 0)
	            grad_x = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
	          else
	            grad_x = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

	          grad_forward = derivs[DPY];
	          grad_back = derivs[DMY];
	          if (this_value > 0)
	            grad_y = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
	          else
	            grad_y = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

	          grad_forward = derivs[DPZ];
	          grad_back = derivs[DMZ];
	          if (this_value > 0)
	            grad_z = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
	          else
	            grad_z = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));


	          speed = reinit_weight*this_value*(target_grad - sqrt(grad_y*grad_y + grad_x*grad_x + grad_z*grad_z));

	          grad_forward = derivs[DPX];
	          grad_back = derivs[DMX];
	          if (speed < 0.0f)
	            grad_x = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
	          else
	            grad_x = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

	          grad_forward = derivs[DPY];
	          grad_back = derivs[DMY];
	          if (speed < 0)
	            grad_y = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
	          else
	            grad_y = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

	          grad_forward = derivs[DPZ];
	          grad_back = derivs[DMZ];
	          if (speed < 0)
	            grad_z = VISMaxAbs(VISmin(grad_forward, 0.0f), VISmax(grad_back, 0.0f));
	          else
	            grad_z = VISMaxAbs(VISmax(grad_forward, 0.0f), VISmin(grad_back, 0.0f));

	          update.poke(i, j, k)  = speed*sqrt(grad_y*grad_y + grad_x*grad_x + grad_z*grad_z);
	          max_speed = VISmax((float)fabs(speed), max_speed);
	     }

    this_dt = 1.0f/(6.0f*max_speed);
    vol += this_dt*update;

    cout << "about to do iteration "  << iterations << endl;
    max_update = 0.0f;
    avg_update = 0.0f;
    for (int j = 0; j < h; j++)
      for (int i = 0; i < w; i++)
	      for (int k = 0; k < d; k++)
	      {
	        get_neighborhood (i, j, k, neighborhood, vol);
	        get_derivs(neighborhood, derivs);
	        for (ii = 0; ii < 3; ii++)
	          normal.poke(ii) = derivs[DX + ii];
	        grad_mag = normal.norm(1.0e-10);
	        normal /= grad_mag;
	        dxx = derivs[DPX] - derivs[DMX];
	        dyy = derivs[DPY] - derivs[DMY];
	        dzz = derivs[DPZ] - derivs[DMZ];
	        dxy = (derivs[DXPY] - derivs[DXMY])/2.0f;
	        dxz = (derivs[DXPZ] - derivs[DXMZ])/2.0f;
	        dyz = (derivs[DYPZ] - derivs[DYMZ])/2.0f;
	        // 	      curv_adjust.poke(0, 0) = dyy + dzz;
	        // 	      curv_adjust.poke(1, 1) = dxx + dzz;
	        // 	      curv_adjust.poke(2, 2) = dxx + dyy;
	        // 	      curv_adjust.poke(0, 1) = curv_adjust.poke(1, 0) = -dxy;
	        // 	      curv_adjust.poke(0, 2) = curv_adjust.poke(2, 0) = -dxz;
	        // 	      curv_adjust.poke(1, 2) = curv_adjust.poke(2, 1) = -dyz;

	        curv_adjust.poke(0, 0) = dxx;
	        curv_adjust.poke(1, 1) = dyy;
	        curv_adjust.poke(2, 2) = dzz;
	        curv_adjust.poke(0, 1) = curv_adjust.poke(1, 0) = dxy;
	        curv_adjust.poke(0, 2) = curv_adjust.poke(2, 0) = dxz;
	        curv_adjust.poke(1, 2) = curv_adjust.poke(2, 1) = dyz;

	        projection = VISIdentity(3) - normal.exterior(normal.t());

	        curv_adjust = (1.0/grad_mag)*projection*curv_adjust*projection;

	        curve_norm = curv_adjust.norm();
	        curve_trace = curv_adjust.trace();

	        tmp = -curve_trace*curve_trace + 2.0f*curve_norm;
	        if (tmp < 0.0)
	          cout << "got bad sqrt " << tmp << " " << curve_trace << " " << curve_norm << endl;
	        k1 = (curve_trace + sqrt(tmp))/2.0f;
	        k2 = (curve_trace - sqrt(tmp))/2.0f;

	        if (k1 > 0.0)
	          k1 = VISmax(k1 - epsilon, 0.0f);
	        else
	          k1 = VISmin(k1 + epsilon, 0.0f);

	        if (k2 > 0.0)
	          k2 = VISmax(k2 - epsilon, 0.0f);
	        else
	          k2 = VISmin(k2 + epsilon, 0.0f);


	        update.poke(i, j, k) = grad_mag*(k1+k2);

	        //curve_trace = (curv_adjust*normal).dot(normal);
	        //	      tmp = -curve_trace*curve_trace + 2.0f*curve_norm;
	        //	      if (tmp < 0.0)
	        //		cout << "got bad sqrt " << tmp << " " << curve_trace << " " << curve_norm << endl;
	        //	      curve_principle[0] = (curve_trace + sqrt(tmp))/2.0f;
	        //	      curve_principle[1] = (curve_trace - sqrt(tmp))/2.0f;
	        //	      for (ii = 0; ii < 2; ii++)
	        //		if (curve_principle[ii] > 0.0f)
	        //		  curve_principle[ii] = VISmax(curve_principle[ii] - epsilon, 0.0f);
	        //		else
	        //		  curve_principle[ii] = VISmin(curve_principle[ii] + epsilon, 0.0f);
	        //	      flow = curve_principle[1] + curve_principle[2];
	        max_update = VISmax<float>(fabs(update.peek(i, j, k)), (double)max_update);
	        avg_update += fabs(update.peek(i, j, k));

	      }
    vol += dt*update;

    cout << "max update " << max_update << endl;
    cout << "avg update " << avg_update/(w*h*d) << endl;

  }
}

  void tighten(VISVolume<float>& volume, VISVolume<float>& volume_radius, float radius)
  {
    float target_grad = 0.2f;
    float minmax = 1.0f;
    VISVolumeFile vol_file;
    VISImageFile im_file;

    //  gaussDiffuse(volume, 1.0f);
    //  gaussDiffuse(volume, 0.2f);

    reinitDistance(volume, target_grad, minmax, REINIT_TOLERANCE,
  		 MAX_REINIT_ITERATIONS);

    //  vol_file.write_float(volume, "init.vol");
    // writeScalarVolumeFile("init.nrrd" , volume);
    //  im_file.write((volume.image()).becomeFlat(), "init.fts");
    VISVolume<float> vol_open = volume, vol_close = volume;
    int w = volume.width(), h = volume.height(), d = volume.depth();
    int iterations;

    float neighborhood[NSIZE], derivs[NUM_DERIVS];

    open_grey(vol_open, volume_radius, radius, target_grad, false);

    // im_file.write((vol_open.image()).becomeFlat(), "vol_open.fts");
    // writeScalarVolumeFile("vol_open.nrrd", vol_open);
    // cout << "done im open" << endl;

    close_grey(vol_close, volume_radius, radius, target_grad, false);

    // im_file.write((vol_close.image()).becomeFlat(), "vol_close.fts");
    // writeScalarVolumeFile("vol_close.nrrd", vol_close);
    // cout << "done im close" << endl;

    VISVolume<float> update = volume.createToSize(), volume_new;
    float grad_mag;
    float dxx, dyy, dxy, dx, dy, dxz, dyz, dzz, dz;
    float value, max_update = FLT_MAX, avg_update;
    float dt = 1.0/4.0;
    int ii;
    int reinit_iterations;
    VISVolume<boolean> mask(w, h, d);

    VISMatrix ident = VISIdentity(3), projection, curvature(3, 3);
    VISVector normal(3);

    volume = volume;

    iterations = 0;
    while ((iterations++ < 50)&&(max_update > REINIT_TOLERANCE))
    {
      cout << "about to do iteration "  << iterations << endl;
      max_update = 0.0f;
      avg_update = 0.0f;
      for (int k = 0; k < d; k++)
        for (int j = 0; j < h; j++)
  	      for (int i = 0; i < w; i++)
  	      {
  	        get_neighborhood (i, j, k, neighborhood, volume);
  	        get_derivs(neighborhood, derivs);
  	        for (ii = 0; ii < 3; ii++)
  	          normal.poke(ii) = derivs[DX + ii];
  	        grad_mag = normal.norm(1.0e-10);
  	        normal /= grad_mag;
  	        dxx = derivs[DPX] - derivs[DMX];
  	        dyy = derivs[DPY] - derivs[DMY];
  	        dzz = derivs[DPZ] - derivs[DMZ];
  	        dxy = (derivs[DXPY] - derivs[DXMY])/2.0f;
  	        dxz = (derivs[DXPZ] - derivs[DXMZ])/2.0f;
  	        dyz = (derivs[DYPZ] - derivs[DYMZ])/2.0f;

  	        update.poke(i, j, k)  = (power(normal[1], 2) +
  				         power(normal[2], 2))*dxx
  	          + (power(normal[0], 2) + power(normal[2], 2))*dyy
  	          + (power(normal[0], 2) + power(normal[1], 2))*dzz
  	          - 2.0*(normal[0]*normal[1]*dxy + normal[0]*normal[2]*dxz +
  		         normal[1]*normal[2]*dyz);
  	      } // loop over image...

          volume_new = volume + dt*update;
          volume_new = volume_new.max(vol_open);
          volume_new = volume_new.min(vol_close);

          reinitDistance(volume_new, target_grad, MINMAX, REINIT_TOLERANCE, MAX_REINIT_ITERATIONS);

          update = (volume_new - volume).abs();
          max_update = update.max();
          avg_update = update.sum();
          cout << "max update " << max_update << endl;
          cout << "avg update " << avg_update/(w*h*d) << endl;
          volume = volume_new;

    } // loop over iterations...

    // fix it at the end to be a proper distance transform
    reinitDistance(volume_new, target_grad, MINMAX, REINIT_TOLERANCE,
  		 MAX_REINIT_ITERATIONS);
    volume = volume_new;
}

void
tighten(VISVolume<float>& volume, float radius)
{
  float target_grad = 0.2f;
  float minmax = 1.0f;
  VISVolumeFile vol_file;
  VISImageFile im_file;

  //  gaussDiffuse(volume, 1.0f);
  //  gaussDiffuse(volume, 0.2f);

  reinitDistance(volume, target_grad, minmax, REINIT_TOLERANCE,
		 MAX_REINIT_ITERATIONS);

  //  vol_file.write_float(volume, "init.vol");
  // writeScalarVolumeFile("init.nrrd" , volume);
  //  im_file.write((volume.image()).becomeFlat(), "init.fts");
  VISVolume<float> vol_open = volume, vol_close = volume;
  int w = volume.width(), h = volume.height(), d = volume.depth();
  int iterations;

  float neighborhood[NSIZE], derivs[NUM_DERIVS];

  open_grey(vol_open, radius, target_grad);

  // im_file.write((vol_open.image()).becomeFlat(), "vol_open.fts");
  // writeScalarVolumeFile("vol_open.nrrd", vol_open);
  // cout << "done im open" << endl;

  close_grey(vol_close, radius, target_grad);

  // im_file.write((vol_close.image()).becomeFlat(), "vol_close.fts");
  // writeScalarVolumeFile("vol_close.nrrd", vol_close);
  // cout << "done im close" << endl;

  VISVolume<float> update = volume.createToSize(), volume_new;
  float grad_mag;
  float dxx, dyy, dxy, dx, dy, dxz, dyz, dzz, dz;
  float value, max_update = FLT_MAX, avg_update;
  float dt = 1.0/4.0;
  int ii;
  int reinit_iterations;
  VISVolume<boolean> mask(w, h, d);

  VISMatrix ident = VISIdentity(3), projection, curvature(3, 3);
  VISVector normal(3);

  volume = volume;

  iterations = 0;
  while ((iterations++ < 50)&&(max_update > REINIT_TOLERANCE))
  {
    cout << "about to do iteration "  << iterations << endl;
    max_update = 0.0f;
    avg_update = 0.0f;
    for (int k = 0; k < d; k++)
      for (int j = 0; j < h; j++)
	      for (int i = 0; i < w; i++)
	      {
	        get_neighborhood (i, j, k, neighborhood, volume);
	        get_derivs(neighborhood, derivs);
	        for (ii = 0; ii < 3; ii++)
	          normal.poke(ii) = derivs[DX + ii];
	        grad_mag = normal.norm(1.0e-10);
	        normal /= grad_mag;
	        dxx = derivs[DPX] - derivs[DMX];
	        dyy = derivs[DPY] - derivs[DMY];
	        dzz = derivs[DPZ] - derivs[DMZ];
	        dxy = (derivs[DXPY] - derivs[DXMY])/2.0f;
	        dxz = (derivs[DXPZ] - derivs[DXMZ])/2.0f;
	        dyz = (derivs[DYPZ] - derivs[DYMZ])/2.0f;

	        update.poke(i, j, k)  = (power(normal[1], 2) +
				         power(normal[2], 2))*dxx
	          + (power(normal[0], 2) + power(normal[2], 2))*dyy
	          + (power(normal[0], 2) + power(normal[1], 2))*dzz
	          - 2.0*(normal[0]*normal[1]*dxy + normal[0]*normal[2]*dxz +
		         normal[1]*normal[2]*dyz);
	      } // loop over image...

        volume_new = volume + dt*update;
        volume_new = volume_new.max(vol_open);
        volume_new = volume_new.min(vol_close);
        reinitDistance(volume_new, target_grad, MINMAX, REINIT_TOLERANCE,
		       MAX_REINIT_ITERATIONS);

        update = (volume_new - volume).abs();
        max_update = update.max();
        avg_update = update.sum();
        cout << "max update " << max_update << endl;
        cout << "avg update " << avg_update/(w*h*d) << endl;
        volume = volume_new;

  } // loop over iterations...

  // fix it at the end to be a proper distance transform
  reinitDistance(volume_new, target_grad, MINMAX, REINIT_TOLERANCE,
		 MAX_REINIT_ITERATIONS);
  volume = volume_new;
}


void
clamp_curvature(VISImage<float>& image, float epsilon)
{
  int ii;
  int w = image.width(), h = image.height();
  VISImage<float> update = image.createToSize();
  float curve = 0.0f;
  float grad_mag_sq = 0.0f, grad_mag = 0.0f;
  float value = 0.0f, max_update = FLT_MAX, avg_update;;
  float dt = 1.0/8.0, this_dt;
  int iterations = 0;
  float tmp = 0.0f;
  float dxx = 0.0f, dyy = 0.0f, dxy = 0.0f, dx = 0.0f, dy = 0.0f;
  float target_grad = 0.25;

  float init_update_max = FLT_MAX;
  while (init_update_max > REINIT_TOLERANCE)
  {
    init_update_max = reinitDistance(image, target_grad);
    cout << "update max is " << init_update_max << endl;
  }

  while ((iterations++ < 1000)&&(max_update > 1.0e-4))
  {
    cout << "about to do iteration "  << iterations << endl;
    max_update = 0.0f;
    avg_update = 0.0f;
    for (int j = 0; j < h; j++)
      for (int i = 0; i < w; i++)
      {
	      dx = (image.peek(VISmin(i+1, w-1), j) -
	            image.peek(VISmax(i-1, 0), j))/2.0f;
	      dy = (image.peek(i, VISmin(j+1, h-1)) -
	            image.peek(i, VISmax(j-1, 0), 0))/2.0f;
	      dxx = image.peek(VISmin(i+1, w-1), j) +
	        image.peek(VISmax(i-1, 0), j) - 2.0*image.peek(i, j);
	      dyy = image.peek(i, VISmin(j+1, h-1)) +
	        image.peek(i, VISmax(j-1, 0), 0) - 2.0*image.peek(i, j);
	      dxy = 0.25*(image.peek(VISmin(i+1, w-1), VISmin(j+1, h-1)) +
		          image.peek(VISmax(i-1, 0), VISmax(j-1, 0))
		          - image.peek(VISmin(i+1, w-1), VISmax(j-1, 0)) -
		          image.peek(VISmax(i-1, 0), VISmin(j+1, h-1)));

	      grad_mag_sq = dx*dx + dy*dy + 1.0e-10;
	      grad_mag = sqrt(grad_mag_sq);

	      curve = (1.0/(grad_mag*grad_mag_sq))*
	        (dxx*dy*dy + dyy*dx*dx - 2.0f*dx*dy*dxy);

	      if (curve > 0.0)
	        curve = VISmax(curve - epsilon, 0.0f);
	      else
	        curve = VISmin(curve + epsilon, 0.0f);

	      update.poke(i, j) = grad_mag*curve;
	      //		update.poke(i, j) = 0.0;
	      //redistancing
	      max_update = VISmax<float>(fabs(update.peek(i, j)), (double)max_update);
	      avg_update += fabs(update.peek(i, j));
      }
    image += dt*update;

    init_update_max = FLT_MAX;
    while (init_update_max > 1.0e-4)
    {
      init_update_max = reinitDistance(image, target_grad);
      cout << "clamp curve update max is " << init_update_max << endl;
    }
    cout << "max update " << max_update << endl;
    cout << "avg update " << avg_update/(w*h) << endl;
    cout << "image min " << image.min() << " and max " << image. max() << endl;
  }
}
