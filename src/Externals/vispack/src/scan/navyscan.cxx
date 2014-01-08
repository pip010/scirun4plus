#include "scan/navyscan.h"
#include "image/imagefile.h"
#include <fstream>

#define AMP_ADJ (7.0f)
#define PHASE_ADJ (-0.1f)
//#define AMP_ADJ (0.0f)
//#define PHASE_ADJ (0.0f)

//*********************************************************************
//*********************************************************************
//********************** 3D NAVY    Scan ****************************
//*********************************************************************
//*********************************************************************

#define TRYKEYWORDREQ(var, pref, key) { \
sprintf(scan_keyword, "%s%s", pref, keywords[key]); \
      if (VPF::set(var, pf[scan_keyword][0]) != VPF::VALID) { \
      cout << "Navyscan::myload -- missing required keyword : " << scan_keyword << endl; \
      exit(-1); } \
}

#define TRYKEYWORD(var, pref, key) { \
      sprintf(scan_keyword, "%s%s", pref, keywords[key]); \
      VPF::set(var, pf[scan_keyword][0]); \
}



char* NavyScan::keywords[6] = 
    {
      "PHI0", "THETA0", "DELTA_THETA", "DELTA_PHI",
      "ALPHA", "BETA"
    };

NavyScan::NavyScan(const VISImage<float> &range_image):Scan(range_image)
{
//   setphi0(PHI_O);
//   settheta0(THETA_O);
//   setDeltaTheta(D_THETA); 
//   setDeltaPhi(D_PHI);
//   setalpha(ALPHA);
//   setbeta(BETA);
  setDefaultParams();

  VISImage<float> conf(range_image.width(), range_image.height());
  conf = 1.0f;
  setConfidence(conf);
}

NavyScan::NavyScan(const VISImage<float> &range_image, float normals_scale)
  :Scan(range_image, normals_scale)
{

//   setphi0(PHI_O);
//   settheta0(THETA_O);
//   setDeltaTheta(D_THETA); 
//   setDeltaPhi(D_PHI);
//   setDeltaPhi(D_PHI);
//   setalpha(ALPHA);
//   setbeta(BETA);

  setDefaultParams();

  VISImage<float> conf(range_image.width(), range_image.height());
  conf = 1.0f;
  setConfidence(conf);
}

void NavyScan::setParams(const VISVector &params)
{
  if (params.n() == NUMNAVYPARAMS)
    { 
      //      cout << " in set " << params.peek(NavScan::PHI0) << " " << 
      //	params.peek(NavScan::THETA0);
      setphi0(params.peek(PHI0)); 
      settheta0(params.peek(THETA0)); 
      setDeltaPhi(params.peek(DPHI)); 
      setDeltaTheta(params.peek(DTHETA)); 
      setalpha(params.peek(ALPHA));
      setbeta(params.peek(BETA));
      setRangeDelta(params.peek(RDELTA)); 
      setRangeOffset(params.peek(ROFFSET)); 
    }

  /* if (params.n() == 6)
    {
      setphi0(params.peek(0)); 
      settheta0(params.peek(1)); 
      setDeltaPhi(params.peek(2)); 
      setDeltaTheta(params.peek(3)); 
      setalpha(params.peek(4));
      setbeta(params.peek(5));
      }*/
}

VISVector NavyScan::getParams() const
{
  VISVector r(NUMNAVYPARAMS);
  r.poke(PHI0) = _phi_0; 
  r.poke(THETA0) = _theta_0; 
  r.poke(DPHI) = _delta_phi; 
  r.poke(DTHETA) = _delta_theta; 
  r.poke(ALPHA) = _alpha;
  r.poke(BETA) = _beta;
  r.poke(RDELTA) = _r_delta; 
  r.poke(ROFFSET) = _r_offset; 
  return(r);
}

VISVector NavyScan::getSurfaceNormal(int j,  int i) const
{
    VISVector r(4);

    float dRdi = (_r_delta*_range_map_orig_dy.peek(j, i));
    float dRdj = (_r_delta*_range_map_orig_dx.peek(j, i));
    float R = _r_offset + _r_delta*_range_map_orig.peek(j, i);

   float cos_t, sin_t, cos_p, sin_p;
   float cos2_t, sin2_t, cos2_p, sin2_p;
   float xi, yi, zi, xj, yj, zj;


    float phi = _delta_phi  *(j+1) + _phi_0,
        theta = _delta_theta*(i+1 + AMP_ADJ*sin(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
						/(float)(_range_map_orig).height())) 
      + _theta_0;
    float d_theta_di
      = _delta_theta*(1 + AMP_ADJ*(2.0f*M_PI/(float)(_range_map_orig).height())
		      *cos(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
			   /(float)(_range_map_orig).height()));
    cos_t = cos(theta);
    sin_t = sin(theta);
    cos_p = cos(phi);
    sin_p = sin(phi);

    xi = R*cos_p*cos_t*d_theta_di + dRdi*cos_p*sin_t;
    yi = dRdi*sin_p;
    zi = -R*cos_p*sin_t*d_theta_di + dRdi*cos_p*cos_t;

    xj = -R*sin_p*cos_t*_delta_phi + dRdj*cos_p*sin_t + _alpha;
    yj = R*cos_p*_delta_phi + dRdj*sin_p;;
    zj = -R*sin_p*cos_t*_delta_phi + dRdj*cos_p*cos_t + _beta;

    r.poke(0) = yi*zj - zi*yj;
    r.poke(1) = zi*xj - xi*zj;
    r.poke(2) = xi*yj - yi*xj;
    r.poke(3) = 0.0f;

//     cos2_t = cos_t*cos_t;
//     sin2_t = sin_t*sin_t;
//     cos2_p = cos_p*cos_p;
//     sin2_p = sin_p*sin_p;

//     r.poke(0) = -1*R*dRdi*cos_t*_delta_phi + _beta*dRdi*sin_p 
//       + R*R*cos2_p*sin_t*_delta_theta*_delta_phi +
// R*dRdj*sin_p*cos_p*sin_t*_delta_theta;
//     r.poke(1) = R*dRdj*cos2_p*_delta_phi -
// R*R*sin_p*cos_p*_delta_theta*_delta_phi
//       + R*cos_p*_delta_theta*(_beta*cos_t + _alpha*sin_t)
//       + dRdi*cos_p*(_beta*sin_t - _alpha*cos_t);
//     r.poke(2) = R*R*cos2_p*cos_t*_delta_theta*_delta_phi +
// R*dRdi*sin_t*_delta_phi
//       + R*dRdj*sin_p*cos_p*cos_t*_delta_theta - _alpha*dRdi*sin_p;
//     r.poke(3) = 0.0f;

// should this be positive or negative???
    r /= -r.peek(2);
    r /= r.norm();
    //r /= r.norm();
    return(r);
}

VISVector NavyScan::getSurfaceNormal(int j, int i, int scale) const
{

    VISVector r(4);
    //    float this_scale = _scale_factor[scale];

//     cout << "in get surf normal" << endl;

//     cout <<  (_range_maps_dx.peek(scale)).width() <<
    //" " (_range_maps_dx.peek(scale)).height() << endl;
//     cout << " j " << j << " i " << i << 

   float dRdi = (_r_delta*_range_maps_dy.peek(scale)).interp(j, i);
   float dRdj = (_r_delta*_range_maps_dx.peek(scale)).interp(j, i);
   float R = (_r_offset + _r_delta*_range_maps.peek(scale)).interp(j, i);

   float cos_t, sin_t, cos_p, sin_p;
   float cos2_t, sin2_t, cos2_p, sin2_p;
   float xi, yi, zi, xj, yj, zj;


    float phi = _delta_phi  *(j+1) + _phi_0,
         theta = _delta_theta*(i+1 + AMP_ADJ*sin(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
						 /(float)(_range_maps.peek(scale)).height()))
      + _theta_0;;
    float d_theta_di
      = _delta_theta*(1 + AMP_ADJ*(2.0f*M_PI/(float)(_range_maps.peek(scale)).height())
		      *cos(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
			   /(float)(_range_maps.peek(scale)).height()));
    cos_t = cos(theta);
    sin_t = sin(theta);
    cos_p = cos(phi);
    sin_p = sin(phi);

    xi = R*cos_p*cos_t*d_theta_di + dRdi*cos_p*sin_t;
    yi = dRdi*sin_p;
    zi = -R*cos_p*sin_t*d_theta_di + dRdi*cos_p*cos_t;

    xj = -R*sin_p*sin_t*_delta_phi + dRdj*cos_p*sin_t + _alpha;
    yj = R*cos_p*_delta_phi + dRdj*sin_p;;
    zj = -R*sin_p*cos_t*_delta_phi + dRdj*cos_p*cos_t + _beta;

    r.poke(0) = yi*zj - zi*yj;
    r.poke(1) = zi*xj - xi*zj;
    r.poke(2) = xi*yj - yi*xj;
    r.poke(3) = 0.0f;

//     cos2_t = cos_t*cos_t;
//     sin2_t = sin_t*sin_t;
//     cos2_p = cos_p*cos_p;
//     sin2_p = sin_p*sin_p;

//     r.poke(0) = -1*R*dRdi*cos_t*_delta_phi + _beta*dRdi*sin_p 
//       + R*R*cos2_p*sin_t*_delta_theta*_delta_phi +
// R*dRdj*sin_p*cos_p*sin_t*_delta_theta;
//     r.poke(1) = R*dRdj*cos2_p*_delta_phi -
// R*R*sin_p*cos_p*_delta_theta*_delta_phi
//       + R*cos_p*_delta_theta*(_beta*cos_t + _alpha*sin_t)
//       + dRdi*cos_p*(_beta*sin_t - _alpha*cos_t);
//     r.poke(2) = R*R*cos2_p*cos_t*_delta_theta*_delta_phi +
// R*dRdi*sin_t*_delta_phi
//       + R*dRdj*sin_p*cos_p*cos_t*_delta_theta - _alpha*dRdi*sin_p;
//     r.poke(3) = 0.0f;


// should this be positive or negative???
    r /= -r.peek(2);
    r /= r.norm();
    //r /= r.norm();

    return(r);
}

void NavyScan::readParams(VPF::ParameterFile pf, const char* prefix)
{

  char scan_keyword[80];

  TRYKEYWORDREQ(_phi_0, prefix, PHI0);
  TRYKEYWORDREQ(_theta_0, prefix, THETA0);
  TRYKEYWORDREQ(_delta_phi, prefix, DPHI);
  TRYKEYWORDREQ(_delta_theta, prefix, DTHETA);
  TRYKEYWORDREQ(_alpha, prefix, ALPHA);
  TRYKEYWORDREQ(_beta, prefix, BETA);
}

VISImage<float> NavyScan::transformRange(const VISImage<float> &im) const
{

  // for debugging 
  // REMOVE THIS!!!!
  //  return(im);

  int w = im.width(), h = im.height();
  VISImage<float> r = im, change;
  
  cout << "transformed range" << endl;

  //  for (j = 0; j < h; j++)
  //    for (i = 0; i < w; i++)
  //      r.poke(i, j) = navy_scan_ripple_adj[i%7]*im.peek(i, j);

// does Gaussian blurring using diffusion to handle boundary conditions
//


  float sigma = 2.0, k = 2.0f;
  float t = sigma*sigma/2.0f;
  int i = (int)ceil(::VISmax(4.0f, t/0.2f)), j, ii;
  float dt = t/(float)i;

//
//   for (j = 0; j < i; j++) {
//     change = r.dx(2);
//     for (ii = 0; ii < h; ii++) {
//       change.poke(0, ii) = r.peek(1, ii) - r.peek(0, ii);
//       change.poke(w - 1, ii) = r.peek(w - 2, ii) - r.peek(w - 1, ii);
//     }
//     r += dt*change;
//   }
//   return(r);
//

// do an ansitropic diffusion to correct problems in the scanner

    VISImage<float> f_x, f_y, 
	dx_half, dy_half, 
	dx, dy, 
	x_kernel, y_kernel;
    x_kernel = VISImage<float>(3, 1);
    x_kernel = 0.5f;
    x_kernel.at(0, 0) = 0.0f;

    for (j = 0; j < i; j++) {
      f_y = r.dy();
      dx_half = r.dxHalfForward();
      f_y = f_y.convolve(y_kernel);
      dx = dx_half*(((dx_half.power(2) + f_y.power(2))
		     /(-2.0f*k*k)).exp());
      change = dx.dxHalfBack();

      for (ii = 0; ii < h; ii++) {
	change.poke(0, ii) = dx.peek(0, ii);
	change.poke(w - 1, ii) = -1.0f*dx.peek(w - 1, ii);
      }
    
      r += dt*change;
    }
   return(r);
}


int NavyScan::getIJ(float &i, float &j, float X, float Y, float Z) const
{
  //This function uses the Newton method to
  //solve for i and j

  int max_iter = 200;
  int l;
  float theta, phi, sum1, sum2;
  VISMatrix J(2,2),F(2,1),dx(2,1),x(2,1);
  int cont = 1, iter = 0;
  float s,c;

  //  int ii, jj;
  float tol = 0.0001;

  j = (atan2(Y,sqrt(X*X+Z*Z)) - _phi_0) / _delta_phi;
  i = (atan2(X,Z) - _theta_0)/_delta_theta;
  VISMatrix J_inv(2,2);
   
  
  while(cont)
    {
      //      ii = (int) rint(i);
      //      jj = (int) rint(j);
      
      //      theta = _delta_theta*(i) + _theta_0;
          theta = _delta_theta*(i + AMP_ADJ*sin(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)(i - 1)
						/(float)(_range_map_orig).height()))
	    + _theta_0;
    float d_theta_di
      = _delta_theta*(1 + AMP_ADJ*(2.0f*M_PI/(float)(_range_map_orig).height())
		      *cos(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)(i - 1)
			   /(float)(_range_map_orig).height()));
      phi   = _delta_phi*(j) + _phi_0;
      s = sin(phi);
      c = cos(phi);

      F.poke(0,0) = (X-(Y/s)*c*sin(theta) - _alpha*(j)) * -1;
      F.poke(1,0) = (Z-(Y/s)*c*cos(theta) - _beta *(j)) * -1;


      sum1 = (fabs(F(0,0))  + fabs(F(1,0)))/2;
      if (sum1 < tol)
      break;
      J.poke(0,0) = (-1*(Y/sin(phi))*cos(phi)*cos(theta)*d_theta_di);
      J.poke(0,1) = ((Y/(sin(phi)*sin(phi)))*sin(theta)*_delta_phi - _alpha);
      J.poke(1,0) = ((Y/sin(phi))*cos(phi)*sin(theta)*d_theta_di);
      J.poke(1,1) = ((Y/(sin(phi)*sin(phi)))*cos(theta)*_delta_phi - _beta);
      float den;

      den = J(1,1)*J(0,0) - J(0,1)*J(1,0);

      if (den != 0.0f)
	{
	  J_inv.poke(0,0) =    J(1,1) / den;
	  J_inv.poke(0,1) = -1*J(0,1) / den;
	  J_inv.poke(1,0) = -1*J(1,0) / den;
	  J_inv.poke(1,1) =    J(0,0) / den;
	}
      else
	{
	  //	  cout << "WARNING:  getIJ did not converge - zero determinant" << endl;
	  i = j = -1;
	  return(-1);
	}

      dx = (J_inv*F);
      
      i += dx(0,0);
      j += dx(1,0);

      sum2 = (fabs(dx(0,0)) + fabs(dx(1,0)))/2;

      if (sum2 < tol)
	cont = 0;
      else if (iter == max_iter)
	cont = 0;
      
      iter++;
    }
  i -= 1;
  j -= 1;
  return(0);
  //  cout<<"get_ij :"<<i<<"   "<<j<<endl;
}  

VISVector NavyScan::myImageCoord(const VISVector &p) const
{
  float x = p.peek(0), y = p.peek(1), z = p.peek(2);
  float theta, phi;
  float r, c;

  if (z > 0)
    {
      getIJ(r,c,x,y,z);
      return(VISVector(c, r));
    }
  else
    {
      //      cout << "imageCoord: ERROR " << p << endl;
      return VISVector();
    }
}


VISVector NavyScan::imageCoord(const VISVector &p, int scale) const
{
    VISVector coord = imageCoord(p);    
    return(coord/_scale_factor[scale]);
}


float NavyScan::depth(const VISVector &p, int scale) const
{
  VISVector coord = imageCoord(p);
  float c, r;
  VISImage<float> this_map;
    
  if (coord.isValid()&&
      ((this_map = _range_maps.peek(scale))
       .checkBounds(c = (coord.peek(0)/
			 _scale_factor[scale]),
		    r = (coord.peek(1)/
			 _scale_factor[scale]))))
    {
//              cout << "point " << p << endl;
//              (this_map).print();
//              cout << "coord " << coord << endl;
//              cout << "c " << c << " r " << r << endl;
//              cout << "depth " << this_map.VISImage<float>::interp(c, r) << endl;
      return(_r_offset + _r_delta*this_map.VISImage<float>::interp(c, r));
    }
  else 		
    return(errorCode());
}

float NavyScan::depth(const VISVector &p) const
{
  // geometry of the range map is defined here.....
  VISVector coord = imageCoord(p);
  float c, r;

  if (coord.isValid()&&(_range_map_orig
			.checkBounds(c = coord.peek(0),r = coord.peek(1))))
    {
//                    cout << "point " << p << endl;
//              _range_map_orig.print();
//              cout << "coord " << coord << endl;
//              cout << "c " << c << " r " << r << endl;
//              cout << "depth " << _range_map_orig.VISImage<float>::interp(c, r) << endl;
    return(_r_offset + _r_delta*_range_map_orig.VISImage<float>::interp(c, r));
    }
  else 		
    return(errorCode());
}

float NavyScan::distance(const VISVector &p) const
{
  // geometry of the range map is defined here.....
  VISVector coord = imageCoord(p);
  float i, j, phi, theta;

  if (coord.isValid())    
    {
      i = coord(1);
      j = coord(0);
      phi   = _delta_phi*(j+1)   + _phi_0;
      theta = _delta_theta*(i+1 + AMP_ADJ*sin(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
					      /(float)(_range_map_orig).height()))
	+ _theta_0;
      //      theta = _delta_theta*(i+1) + _theta_0;
      return
	sqrt(pow((p.peek(0) - _alpha*(j+1)), 2)  
	     +  pow((p.peek(1)), 2)  
	     + pow((p.peek(2) - _beta*(j+1)), 2));
    }
  else 		
    return(errorCode());
}


VISVector NavyScan::grad_depth(const VISVector &p, int scale) const
{
  VISVector coord = myImageCoord(p);
  float c, r;
  float j, i;
  VISImage<float> this_map, im_dx, im_dy;
  float dRdi, dRdj;
  //  float x = p.peek(0), y = p.peek(1), z = p.peek(2);

  if (!(coord.isValid()&&
	((im_dx  = _range_maps_dx.peek(scale))
	 .checkBounds(c = ((j = coord.peek(0))/_scale_factor[scale]),
		      r = ((i = coord.peek(1))/_scale_factor[scale])))))
    return VISVector(0.0f, 0.0f, 0.0f, 0.0f);
  im_dy = _range_maps_dy.peek(scale);

  dRdj = _r_delta*im_dx.interp(c, r);
  dRdi = _r_delta*im_dy.interp(c, r);

  //  cout << "dc " << dc << " and dr " << dr << endl;
  float cos_t, sin_t, cos_p, sin_p;

  VISMatrix dT(3, 3);


  float phi = _delta_phi *(j+1) + _phi_0,
      theta = _delta_theta*(i+1 + AMP_ADJ*sin(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
					      /(float)_range_map_orig.height()))
    + _theta_0;

    float d_theta_di
      = _delta_theta*(1 + AMP_ADJ*(2.0f*M_PI/(float)(_range_map_orig).height())
		      *cos(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
			   /(float)(_range_map_orig).height()));

  //    theta = _delta_theta*(c+1) + _theta_0;
  // you must come up with expressions for derivatives of this!!!!!!!
  float R = sqrt(pow((p.peek(0) - _alpha*(j+1)), 2)  
	     +  pow((p.peek(1)), 2)  
 	     + pow((p.peek(2) - _beta*(j+1)), 2));

  cos_t = cos(theta);
  sin_t = sin(theta);
  cos_p = cos(phi);
  sin_p = sin(phi);

  dT.poke(0, 0) = R*cos_p*cos_t*d_theta_di;
  dT.poke(0, 1) = -R*sin_p*sin_t*_delta_phi + _alpha;
  dT.poke(0, 2) = cos_p*sin_t;

  //  dyd = ????;
  dT.poke(1, 0) = 0.0f;
  dT.poke(1, 1) = R*cos_p*_delta_phi;
  dT.poke(1, 2) = sin_p;

  dT.poke(2, 0) = -R*cos_p*sin_t*d_theta_di;
  dT.poke(2, 1) = -R*sin_p*cos_t*_delta_phi + _beta;
  dT.poke(2, 2) = cos_p*cos_t;

  dT = dT.inverseGJ();

  return VISVector
    (
     dRdi*dT.peek(0,0) + dRdj*dT.peek(1,0),
     dRdi*dT.peek(0,1) + dRdj*dT.peek(1,1),
     dRdi*dT.peek(0,2) + dRdj*dT.peek(1,2),
     0.0f
     );
}

VISVector NavyScan::grad_conf(const VISVector &p, int scale) const
{
  VISVector coord = myImageCoord(p);
  float c, r;
  float j, i;
  VISImage<float> this_map, im_dx, im_dy;
  float dCdi, dCdj;
  //  float x = p.peek(0), y = p.peek(1), z = p.peek(2);

  if (!(coord.isValid()&&
	((im_dx  = _conf_maps_dx.peek(scale))
	 .checkBounds(c = ((j = coord.peek(0))/_scale_factor[scale]),
		      r = ((i = coord.peek(1))/_scale_factor[scale])))))
    return VISVector(0.0f, 0.0f, 0.0f, 0.0f);
  im_dy = _conf_maps_dy.peek(scale);

  dCdj = im_dx.interp(c, r);
  dCdi = im_dy.interp(c, r);

  //  cout << "dc " << dc << " and dr " << dr << endl;
  float cos_t, sin_t, cos_p, sin_p;

  VISMatrix dT(3, 3);


  float phi = _delta_phi*(j+1) + _phi_0,
    theta = _delta_theta*(i+1 + AMP_ADJ*sin(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
					      /(float)(_range_map_orig).height()))
    + _theta_0;

  float d_theta_di
    = _delta_theta*(1 + AMP_ADJ*(2.0f*M_PI/(float)(_range_map_orig).height())
		    *cos(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
			 /(float)(_range_map_orig).height()));
  //    theta = _delta_theta*(c+1) + _theta_0;

  // you must come up with expressions for derivatives of this!!!!!!!
  float R = sqrt(pow((p.peek(0) - _alpha*(j+1)), 2)  
	     +  pow((p.peek(1)), 2)  
 	     + pow((p.peek(2) - _beta*(j+1)), 2));

  cos_t = cos(theta);
  sin_t = sin(theta);
  cos_p = cos(phi);
  sin_p = sin(phi);

  dT.poke(0, 0) = R*cos_p*cos_t*d_theta_di;
  dT.poke(0, 1) = -R*sin_p*sin_t*_delta_phi + _alpha;
  dT.poke(0, 2) = cos_p*sin_t;

  //  dyd = ????;
  dT.poke(1, 0) = 0.0f;
  dT.poke(1, 1) = R*cos_p*_delta_phi;
  dT.poke(1, 2) = sin_p;

  dT.poke(2, 0) = -R*cos_p*sin_t*d_theta_di;
  dT.poke(2, 1) = -R*sin_p*cos_t*_delta_phi + _beta;
  dT.poke(2, 2) = cos_p*cos_t;

  dT = dT.inverseGJ();

  return VISVector
    (
     dCdi*dT.peek(0,0) + dCdj*dT.peek(1,0),
     dCdi*dT.peek(0,1) + dCdj*dT.peek(1,1),
     dCdi*dT.peek(0,2) + dCdj*dT.peek(1,2),
     0.0f
     );
}




// VISVector NavyScan::grad_depth(const VISVector &p, int scale) const
// {
//   VISVector coord = myImageCoord(p);
//   float c, r;
//   VISImage<float> this_map, im_dx, im_dy;
//   float dRdi, dRdj;
//   //  float x = p.peek(0), y = p.peek(1), z = p.peek(2);

//   if (!(coord.isValid()&&
// 	((im_dx  = _range_maps_dx.peek(scale))
// 	 .checkBounds(c = (coord.peek(0)/_scale_factor[scale]),
// 		      r = (coord.peek(1)/_scale_factor[scale])))))
//     return VISVector(0.0f, 0.0f, 0.0f, 0.0f);
//   im_dy = _range_maps_dy.peek(scale);

//   dRdi = im_dx.interp(c, r);
//   dRdj = im_dy.interp(c, r);

//   //  cout << "dc " << dc << " and dr " << dr << endl;
//   float dxdi, didx, dydi, didy, dzdi, didz;
//   float dxdj, djdx, dydj, djdy, dzdj, djdz;
//   float cos_t, sin_t, cos_p, sin_p;
//   float R = (_range_maps.peek(scale)).interp(c, r);
//   float phi = _delta_phi  *(r+1) + _phi_0,
//         theta = _delta_theta*(c+1) + _theta_0;

//   cos_t = cos(theta);
//   sin_t = sin(theta);
//   cos_p = cos(phi);
//   sin_p = sin(phi);

//   dxdi = R*cos_p*cos_t*_delta_theta + dRdi*cos_p*sin_t;
//   dydi = dRdi*sin_p;
//   dzdi = -1*R*cos_p*sin_t*_delta_theta + dRdi*cos_p*cos_t;
//   dxdj = -1.0*R*sin_t*sin_p*_delta_phi + dRdj*sin_t*cos_p + _alpha;
//   dydj = R*cos_p*_delta_phi + dRdj*sin_p;
//   dzdj = -1*R*cos_t*sin_p*_delta_phi + dRdj*cos_t*cos_p + _beta;

//   if (dxdi != 0.0f)
//     didx = 1.0f/dxdi;
//   else 
//     didx = 0.0f;
//   if (dydi != 0.0f)
//     didy = 1.0f/dydi;
//   else
//     didy = 0.0f;
//   if (dzdi != 0.0f)
//     didz = 1.0f/dzdi;
//   else
//     didz = 0.0f;
//   if (dxdj != 0.0f)
//     djdx = 1.0f/dxdj;
//   else 
//     djdx = 0.0f;
//   if (dydj != 0.0f)
//     djdy = 1.0f/dydj;
//   else
//     djdy = 0.0f;
//   if (dzdj != 0.0f)
//     djdz = 1.0f/dzdj;
//   else
//     djdz = 0.0f;

// //   if ((dxdi == 0.0f)||
// //       (dxdj == 0.0f)||
// //       (dydi == 0.0f)||
// //       (dydj == 0.0f)||
// //       (dzdi == 0.0f)||
// //       (dzdj == 0.0f))
// //     {
// //       cout << coord << endl;
// //       cout << " " << cos_t << " " << 
// // 	sin_t << " " << 
// // 	cos_p << " " << 
// // 	sin_p << endl;
// //       cout << R << " " << dRdi << " " << dRdj << endl;
// //     }

//   return VISVector
//     (
//      dRdi*didx + dRdj*djdx, 
//      dRdi*didy + dRdj*djdy, 
//      dRdi*didz + dRdj*djdz, 
//      0.0f
//      );
// }


float NavyScan::confidence(const VISVector &p) const
{
  VISVector coord = imageCoord(p);
  float c, r;

  //  cout << "image coord" << coord << endl;
  //  cout << "point " << p << endl;
    
  if (coord.isValid()&&
	(_conf_map.checkBounds(c = coord.peek(0),r = coord.peek(1))))
	return(_conf_map.VISImage<float>::interp(c, r));
    else 		
	{
	  //	    cout << "got conf oob at " << c << " " << r << endl 
	  //		 << flush;
	    return(0.0f);
	}
}

float NavyScan::confidence(const VISVector &p, int scale) const
{

  return(confidence(p));
  
}

VISVector NavyScan::lineOfSight(const VISVector &p) const
{
  float 
    i,j,
    z = p[2],
    y = p[1],
    x = p[0];

  
  VISVector coord = myImageCoord(p);

  if (coord.isValid())
    {
      i = coord(1);
      j = coord(0);

      if (depth(p) == errorCode())
	{
	  //	  cout << "point " << p << " is an error" << endl ;
	  //	  cout << "coord " << coord  << endl ;
	  //
	  //return(VISVector(0.0f, 0.0f, 0.0f, 0.0f));
	}
  
	  float theta, phi;
  
	  phi   = _delta_phi*(j+1)   + _phi_0;
	  //	  theta = _delta_theta*(i+1) + _theta_0;
	  theta = _delta_theta*(i+1 + AMP_ADJ*sin(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
						  /(float)(_range_map_orig).height()))
	    + _theta_0;
	  VISVector ret;
      
	  ret = VISVector(cos(phi)*sin(theta), 
			  sin(phi),
			  cos(phi)*cos(theta),
			  0.0f);
	  // the convention is to have the line of sight pointing back toward the scanner
	  ret /= -1.0f*ret.norm();
	  return ret;
	}
      else
	{
	  //      cout << "lineofsight/imageCoord: ERROR " << p << endl;      
	  return VISVector();
	}
    }

VISVector NavyScan::get3DPoint(int j, int i) const
{
  float 
    phi = _delta_phi*(j+1) + _phi_0,
    theta = _delta_theta*(i+1 + AMP_ADJ*sin(M_PI*PHASE_ADJ + 2.0f*M_PI*(float)i
					    /(float)(_range_map_orig).height())) 
    + _theta_0,
  //    theta = _delta_theta*(i+1) + _theta_0,
    r = _r_offset + _r_delta*_range_map_orig.interp(j, i);

  float x,y,z;
  x = r*cos(phi)*sin(theta) + _alpha*(j+1);
  y = r*sin(phi);
  z = r*cos(phi)*cos(theta) + _beta *(j+1);
  return(VISVector(x, y, z, 1.0f)); 
}

VISVector NavyScan::get3DPoint(int i, int j, int scale) const
// should this assume i, i are scaled or not?
// for now, yes
{
    float this_scale = _scale_factor[scale];
    float phi, theta, r, x, y, z;
    // I don't think you should be scalling these like this --- RTW
    float ii = this_scale*i;
    float jj = this_scale*j;
    
    phi  = _delta_phi*(ii+1) + _phi_0;
    theta = _delta_theta*(jj+1  +  AMP_ADJ
			  *sin(M_PI*PHASE_ADJ 
			       + 2.0f*M_PI*(float)j/(float)(_range_maps.peek(scale)).height())) + _theta_0;
    r = _r_offset + _r_delta*(_range_maps.peek(scale)).peek(i, j);

    //    cout << "phi " << _delta_phi << " " << ii << " " << _phi_0/_delta_phi << endl;
    //    cout << "theta " << _delta_theta << " " << jj << " " << _theta_0/_delta_theta << endl;
    //    cout << "N phi " << phi << "N theta " << theta << "N r " << r << endl;

    
    x = r*cos(phi)*sin(theta) + _alpha*(ii+1);
    y = r*sin(phi);
    z = r*cos(phi)*cos(theta) + _beta *(ii+1);

    return(VISVector(x, y, z, 1.0f)); 
}

VISMatrix NavyScan::dx_dparams(int c, int r) const
{
  return VISMatrix();
}
