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
 *  GenerateStreamLinesWithPlacementHeuristic.h:  Create optimal stream lines
 *
 *  Written by:
 *   Frank B. Sachse
 *   CVRTI
 *   University of Utah
 *   February 2004, JULY 2004
 */

#include <Core/Thread/Thread.h>

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Core/Algorithms/Fields/SampleField/GeneratePointSamplesFromField.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>

namespace SCIRun {


class GenerateStreamLinesWithPlacementHeuristicData
{
  public:
    FieldHandle msl;
    FieldHandle cf;

    FieldHandle seed_fieldH;
    FieldHandle compare_fieldH;

    FieldHandle src_fieldH;
    FieldHandle weighting_fieldH;
    FieldHandle vfhandle;
    FieldHandle pts_fieldH;
    
    VField*  vf_field;
    VField*  src_field;
    VMesh*   src_mesh;
    VField*  weighting_field;
    VField*  seed_field;
    VField*  compare_field;
    VMesh*   compare_mesh;

    int numsl; // number of streamlines
    int numpts;  // number of trials
    double minper; // mininal rendering radius in percent of max. bounding box length
    double maxper; // maximal rendering radius in percent of max. bounding box length
    double ming; // minimum for clamping and normalizing of field magnitude
    double maxg;  // maximum for clamping and normalizing of field magnitude
    double maxlen; // max. bounding box length
    int numsamples; // number of samples for 1D-sampling of render function
    double stepsize; // step size for integration of streamlines
    int stepout; // increment for output of streamline nodes
    int maxsteps; // maximal number of steps
    int numsteps; // current number of steps
    double minmag; // minimal magnitude for stop of integration
    int direction; // direction of streamline integration
    int method; // numerical method of streamline integration
    double thresholddot; // threshold value for dot product of new and old direction for stop of integration

    double comparemean;
    double srcmean;

    GenerateStreamLinesWithPlacementHeuristicData()
    {
      seed_fieldH=0;
      compare_fieldH=0;
      src_fieldH=0;
      weighting_fieldH=0;
      pts_fieldH=0;
      vf_field=0;

      minper=0.;
      maxper=1.;
      maxlen=0;
      numsamples=1;

      direction=1;
      method=0;
      thresholddot=-.99;
    }

    //find the greatest dimension of the source field
    void boundingBoxMaxLength(VMesh* mesh) 
    {
      const BBox bbox = mesh->get_bounding_box();
      if (bbox.valid()) 
      {
        Vector size = bbox.diagonal();
        maxlen = std::max(size.x(), size.y());
        maxlen = std::max(maxlen, size.z());
      }
      else
      {
        maxlen=0;
      }
    }
};

class GenerateStreamLinesWithPlacementHeuristicGUI;

class GenerateStreamLinesWithPlacementHeuristicAlgo :
  public GenerateStreamLinesWithPlacementHeuristicData
{
  public:
    void execute(GenerateStreamLinesWithPlacementHeuristicData &, ProgressReporter * mod );

  private:
    double evaluateTrial();
    double evaluateRMSDif();
    double getMean(VField *field); 
    void renderPointSet(double factor);
    void renderPoint(const Point &seed, double factor);

    void FindNodes(Point x, double stepsize);
    bool interpolate(const Point &p, Vector &v);
    void addPoint(const int i, const Point &p);
    void FindEuler(Point x, double stepsize); 
    void FindRK2(Point x, double stepsize);
    void FindRK4(Point x, double stepsize);
    FieldHandle createStreamLine(Point seed);
    
};


void 
GenerateStreamLinesWithPlacementHeuristicAlgo::
execute(GenerateStreamLinesWithPlacementHeuristicData &MESL, ProgressReporter * mod)
{
  using std::cerr;
  using std::endl;
 cerr<<"check0"<<endl;
   *(GenerateStreamLinesWithPlacementHeuristicData *)this=MESL;  
   cerr<<"check0"<<endl;
  src_field = src_fieldH->vfield();
  src_mesh =  src_fieldH->vmesh();
  src_field->clear_all_values();
  vf_field=vfhandle->vfield();
  
 cerr<<"check0"<<endl; 
  compare_field = compare_fieldH->vfield();
  compare_mesh =  compare_field->vmesh();
  comparemean=getMean(compare_field);
  cerr<<"check0"<<endl;
  boundingBoxMaxLength(src_mesh);
  
  cerr<<"check1f"<<endl;

  weighting_field = weighting_fieldH->vfield();
  seed_field = seed_fieldH->vfield();

  FieldInformation fi("CurveMesh",1,"double");
  msl = CreateField(fi);
  VMesh* mslm= msl->vmesh();

  SCIRunAlgo::GeneratePointSamplesFromFieldAlgo SFalgo;
  unsigned int inumpts = (int)numpts;
  double rmsmin;

  SFalgo.set_int("num_seed_points",1);
  SFalgo.set_option("seed_method","uniuni");
  SFalgo.set_bool("clamp",false);
  SFalgo.set_progress_reporter(mod);

  //for each streamline, try placing it in a number of random points, and keep the best one
  for(int i=0; i<numsl; i++) 
  {
    mod->remark("streamline: "+to_string(i));
    rmsmin=HUGE;
    bool found=false;
    FieldHandle mslminH;

    for(unsigned int j=0 ; j<inumpts; ++j) 
    {
      //mod->remark("trial: "+to_string(j));
#ifdef _WIN32
      // standard Windows SDK only provides rand() and srand(..)
      SFalgo.set_int("rng_seed",(int)(pow(2.,31.)*rand()));
#else
      SFalgo.set_int("rng_seed",(int)(pow(2.,31.)*drand48()));
      //SFalgo.set_int("rng_seed",(int)(1));
#endif

      FieldHandle rph;

      SFalgo.run(seed_fieldH,rph);

      VMesh *rp_mesh = rph->vmesh();
      SCIRun::index_type mysize;
      mysize=rp_mesh->num_nodes();
      
      //cout<<"num nodes: "<<mysize<<endl;
      
      VMesh::Node::iterator rp_itr;
      rp_mesh->begin(rp_itr);
      Point seed;
      rp_mesh->get_point(seed, *rp_itr);
      //cout<<"seed: "<<seed<<endl;
      pts_fieldH = createStreamLine(seed);
      if (numsteps<2)
      {
        j--;
      }
      else 
      { // evaluate the new image
        double rms=evaluateTrial();
        if (rmsmin>rms)
        {
          rmsmin=rms;
          found=true;
          mslminH=pts_fieldH;
        }
      }

      pts_fieldH=0;
    }

    if (found) 
    {
      //          cerr << "GenerateStreamLinesWithPlacementHeuristicAlgoT::execute best assigned " << rmsmin << "\n";
      pts_fieldH=mslminH;
      renderPointSet(1);

      VMesh *mslminm=mslminH->vmesh();
      VMesh::Node::iterator n_itr, n_end;
      mslminm->begin(n_itr); mslminm->end(n_end);

      VMesh::Edge::iterator e_itr, e_end;
      mslminm->begin(e_itr);mslminm->end(e_end);
      VMesh::Node::array_type na(2);
      
      if (n_itr!=n_end) 
      {
        Point p;
        mslminm->get_center(p, *n_itr);
        VMesh::Node::index_type n1=mslm->add_node(p), n2;
        ++n_itr;
        while (n_itr!=n_end) 
        {
          mslminm->get_center(p, *n_itr);
          n2=mslm->add_node(p);
          if (e_itr!=e_end) 
          {
            VMesh::Node::array_type n;
            mslminm->get_nodes(n, *e_itr);
            if (n[0]==*n_itr-1 && n[1]==*n_itr) 
            {
              ++e_itr;
              na[0] = n1; na[1] = n2;
              mslm->add_elem(na);
            }
          }
          n1=n2;
          ++n_itr;
        }
      }
    }
    mod->update_progress(i+1, numsl);
  }
  msl->vfield()->resize_values();

  MESL=*(GenerateStreamLinesWithPlacementHeuristicData *)this;
}


double 
GenerateStreamLinesWithPlacementHeuristicAlgo::
evaluateTrial()
{
  renderPointSet(1);
  double rms=evaluateRMSDif();
  renderPointSet(-1);
  return rms;
}

double 
GenerateStreamLinesWithPlacementHeuristicAlgo::
evaluateRMSDif()
{
  VMesh::Node::iterator src_itr, src_end_itr;
  src_mesh->begin(src_itr);
  src_mesh->end(src_end_itr);
  VMesh::Node::iterator compare_itr, compare_end_itr;
  compare_mesh->begin(compare_itr);
  compare_mesh->end(compare_end_itr);

  srcmean=getMean(src_field);

  double rms=0., srcscale=srcmean, comparescale=comparemean;
  //    cerr << "srcmean: " << srcscale << "\t" << "comparemean: " << comparescale;
  if (srcscale==0) srcscale=1.;
  if (comparescale==0) comparescale=1.;
  while(src_itr != src_end_itr) 
  {
    double srcnorm, comparenorm;
    src_field->get_value(srcnorm,*src_itr);
    srcnorm=srcnorm/srcscale;
    compare_field->get_value(comparenorm,*compare_itr);
    comparenorm=comparenorm/comparescale;
    double tmp=srcnorm-comparenorm;
    rms+=tmp*tmp;

    ++src_itr;
    ++compare_itr;
  }
  //  cerr << "\t" << "rms: " << sqrt(rms) << endl;
  return sqrt(rms);
}


double 
GenerateStreamLinesWithPlacementHeuristicAlgo::
getMean(VField *field)
{
  VMesh::Node::iterator itr, end_itr;
  VMesh *mesh = field->vmesh();

  mesh->begin(itr);
  mesh->end(end_itr);

  double mean=0, dummy;
  size_type num_nodes = mesh->num_nodes();
  while(itr != end_itr) 
  {
    field->get_value(dummy, *itr);
    mean+=dummy;
    ++itr;
  }
  if (num_nodes == 0) return (0);
  return (mean/num_nodes);
}

void 
GenerateStreamLinesWithPlacementHeuristicAlgo::
renderPointSet(double factor)
{
  VMesh *pts_mesh = pts_fieldH->vmesh();

  VMesh::Node::iterator pts_itr, pts_end_itr;
  pts_mesh->begin(pts_itr);
  pts_mesh->end(pts_end_itr);

  for( ; pts_itr != pts_end_itr; ++pts_itr) 
  {
    Point p;
    pts_mesh->get_center(p, *pts_itr);
    renderPoint(p, factor);
  }
}

void 
GenerateStreamLinesWithPlacementHeuristicAlgo::
renderPoint(const Point &seed, double factor)
{
  Point p;
  double l;
  
  vf_field->interpolate(l, seed);
  if (l<ming) l=ming;
  else if (l>maxg) l=maxg;

  double maxp=(maxper/100.)*maxlen;
  double minp=(minper/100.)*maxlen;
  double lnorm=(ming==maxg ? 0.5 : (l-ming)/(maxg-ming));
  double radius=(numsamples==1 ? 0. : maxp-(maxp-minp)*lnorm);
  double dradius3=(radius ? 1./(3.*radius*radius) : 1.);

  double ffactor=factor*(radius ? (sqrt(2.*M_PI)/radius) : 1.);
  double inc=(numsamples==1 ? 0. : 2.*radius/(numsamples-1));

  double pz=seed.z();
  double py=seed.y();
  double px=seed.x();

  double az=pz-radius;

  VMesh::ElemInterpolate ei;
  std::vector<double> values;

  for(int i=0; i<numsamples; i++, az+=inc) 
  {
    p.z(az);
    double rz2=pz-az;
    rz2*=rz2;
    double ay=py-radius;

    for(int i=0; i<numsamples; i++, ay+=inc) 
    {
      p.y(ay);
      double ry2=py-ay;
      ry2*=ry2;
      double ax=px-radius;

      for(int i=0; i<numsamples; i++, ax+=inc) 
      {
        p.x(ax);
        double rx2=px-ax;
        rx2*=rx2;
        double val=ffactor*exp(-(rx2+ry2+rz2)*dradius3);
        
        src_mesh->get_interpolate_weights(p,ei,1);
        src_field->get_values(values,ei.node_index);
        for (size_t j=0; j<ei.weights.size();j++) values[j] += val*ei.weights[j];
        src_field->set_values(values,ei.node_index);
      }
    }
  }
}

void 
GenerateStreamLinesWithPlacementHeuristicAlgo::
FindNodes(Point x, double stepsize)
{
  switch(method) {
  case 0:
    FindEuler(x, stepsize);
    
    break;

  case 1:
    FindRK2(x, stepsize);
    break;

  case 2:
    FindRK4(x, stepsize);
    break;

  default:
    ASSERT(0);
    break;
  }
}

//! interpolate using the generic linear interpolator
bool 
GenerateStreamLinesWithPlacementHeuristicAlgo::
interpolate(const Point &p, Vector &v)
{
   VMesh::index_type idx;
   idx=0;
  if (vf_field->interpolate(v, p)) 
  {
    double a=v.safe_normalize();
    //cerr<<"safe norm "<<a<<endl;
    // cerr<<"minmag "<<minmag<<endl;
    return a > minmag;
  }

  return (false);
}

//! add point to streamline
void 
GenerateStreamLinesWithPlacementHeuristicAlgo::
addPoint(const int i, const Point &p)
{
  if (!(i%stepout)) 
  {
    numsteps++;
    VMesh::Node::index_type n1 = cf->vmesh()->add_node(p);
    VMesh::Node::array_type na(2);
    na[0] = static_cast<VMesh::Node::index_type>(n1-1); 
    na[1] = n1;
    if (i && n1) cf->vmesh()->add_elem(na);
  }
}

void 
GenerateStreamLinesWithPlacementHeuristicAlgo::
FindEuler(Point x, double stepsize)
{
  Vector v0, v0old;
  int i;
  for (i=0; i < maxsteps; i++)
  {
    if (!interpolate(x, v0)) 
    {
      //cerr<<"Didn't work"<<endl;
    break;
    }
    addPoint(i, x);

    if (i) if (Dot(v0, v0old)<thresholddot)   
    {
      break;
    }
    v0old=v0;
    x += stepsize*v0 ;
  }
}


void 
GenerateStreamLinesWithPlacementHeuristicAlgo::
FindRK2(Point x, double stepsize)
{
  Vector v0, v0old;
  int i;

  for (i=0; i < maxsteps; i ++) 
  {
    if (!interpolate(x, v0)) break;
    addPoint(i, x);

    if (!interpolate(x + v0*stepsize*0.5, v0)) break;

    if (i) if (Dot(v0, v0old)<thresholddot) 
    {
      /*        cerr << "Stopped Dot\n"; */
      break;
    }
    v0old=v0;
    x += stepsize*v0;
  }
  /*       cerr << "Stopped after " << i << "\n"; */
}


void 
GenerateStreamLinesWithPlacementHeuristicAlgo::
FindRK4(Point x, double stepsize)
{
  Vector f[4], v0, v0old;
  int i;

  for (i = 0; i < maxsteps; i++) 
  {
    if (!interpolate(x, f[0])) break;
    addPoint(i, x);

    if (!interpolate(x + f[0] * stepsize*0.5, f[1])) break;
    if (!interpolate(x + f[1] * stepsize*0.5, f[2])) break;
    if (!interpolate(x + f[2] * stepsize, f[3])) break;

    v0= f[0] + 2.0 *(f[1] + f[2]) + f[3];
    v0.safe_normalize();

    if (i) if (Dot(v0, v0old)<thresholddot) 
    {
      break;
    }
    v0old=v0;
    x += stepsize*v0;
  }
}

FieldHandle 
GenerateStreamLinesWithPlacementHeuristicAlgo::
createStreamLine(Point seed)
{
  numsteps=0;
  FieldInformation fi("CurveMesh",1,"double");
  cf = CreateField(fi);

  // Find the negative streamlines.
  if( direction <= 1 ) FindNodes(seed, -stepsize);

  // Append the positive streamlines.
  if( direction >= 1 ) FindNodes(seed, stepsize);

  return cf;
}


class GenerateStreamLinesWithPlacementHeuristicGUI : virtual public GenerateStreamLinesWithPlacementHeuristicData
{
public:
  GuiInt numsl_;
  GuiInt numpts_;
  GuiDouble minper_;
  GuiDouble maxper_;
  GuiDouble ming_;
  GuiDouble maxg_;
  GuiInt numsamples_;
  GuiInt method_;
  GuiDouble stepsize_;
  GuiInt stepout_;
  GuiInt maxsteps_;
  GuiDouble minmag_;
  GuiInt direction_;

  GenerateStreamLinesWithPlacementHeuristicGUI(GuiContext* ctx):
    numsl_(ctx->subVar("numsl")),
    numpts_(ctx->subVar("numpts")),
    minper_(ctx->subVar("minper")),
    maxper_(ctx->subVar("maxper")),
    ming_(ctx->subVar("ming")),
    maxg_(ctx->subVar("maxg")),
    numsamples_(ctx->subVar("numsamples")),
    method_(ctx->subVar("method")),
    stepsize_(ctx->subVar("stepsize")),
    stepout_(ctx->subVar("stepout")),
    maxsteps_(ctx->subVar("maxsteps")),
    minmag_(ctx->subVar("minmag")),
    direction_(ctx->subVar("direction"))
  {
  };

  void LoadGuiVariables()
  {
    numsl=numsl_.get();
    numpts=numpts_.get();
    minper=minper_.get();
    maxper=maxper_.get();
    ming=ming_.get();
    maxg=maxg_.get();
    numsamples=numsamples_.get();
    stepsize = stepsize_.get();
    stepout = stepout_.get();
    maxsteps = maxsteps_.get();
    minmag = minmag_.get();
    direction = direction_.get();
    method = method_.get();
  }
};



class GenerateStreamLinesWithPlacementHeuristic
  : public GenerateStreamLinesWithPlacementHeuristicGUI,
    public Module
{
  public:
    GenerateStreamLinesWithPlacementHeuristic(GuiContext* ctx);
    virtual ~GenerateStreamLinesWithPlacementHeuristic() {}
    virtual void execute();

};


DECLARE_MAKER(GenerateStreamLinesWithPlacementHeuristic)

GenerateStreamLinesWithPlacementHeuristic::
GenerateStreamLinesWithPlacementHeuristic(GuiContext* ctx) :
  GenerateStreamLinesWithPlacementHeuristicGUI(ctx),
  Module("GenerateStreamLinesWithPlacementHeuristic", ctx, Source, "Visualization", "SCIRun")
{
}

void
GenerateStreamLinesWithPlacementHeuristic::execute()
{
  
  
  get_input_handle("Flow",vfhandle,true);
  get_input_handle("Source",src_fieldH,true);
  get_input_handle("Weighting",weighting_fieldH,true);
  get_input_handle("Compare",compare_fieldH,true);
  get_input_handle("Seed points",seed_fieldH,false);

  // Check that the flow field input is a vector field.
  if (!(vfhandle->vfield()->is_vector())) 
  {
    error("Flow is not a Vector field.");
    return;
  }
  
  // Check that the source field input is a scalar field.
  if (!(src_fieldH->vfield()->is_scalar())) 
  {
    error("Source is not a Scalar field.");
    return;
  }

  src_fieldH.detach();

  if (!(weighting_fieldH->vfield()->is_scalar())) 
  {
    error("Weighting is not a scalar field.");
    return;
  }

  if (!(compare_fieldH->vfield()->is_scalar())) 
  {
    error("Compare is not a Scalar field.");
    return;
  }

  // Inform module that execution started
  update_state(Executing);
  
  if (!( seed_fieldH.get_rep())) seed_fieldH=src_fieldH;

  GenerateStreamLinesWithPlacementHeuristicGUI::LoadGuiVariables();

  GenerateStreamLinesWithPlacementHeuristicAlgo algo;
  algo.execute(*(GenerateStreamLinesWithPlacementHeuristicData *)this, (ProgressReporter *)this);

  send_output_handle("Streamlines",msl,true);
  send_output_handle("Render",src_fieldH,true);
}
  
} // End namespace SCIRun
