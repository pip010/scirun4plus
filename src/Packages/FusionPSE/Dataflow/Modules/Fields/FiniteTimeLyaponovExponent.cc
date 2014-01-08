/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
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
 *  FiniteTimeLyaponovExponent.cc:
 *
 *  Written by:
 *   Allen R. Sanderson
 *   SCI Institute
 *   University of Utah
 *   September 2005
 *
 */

#include <Core/Thread/Mutex.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Core/Containers/Handle.h>

#include <Core/Algorithms/Fields/StreamLines/StreamLineIntegrators.h>

namespace FusionPSE {

using namespace SCIRun;
using namespace SCIRunAlgo;



class FiniteTimeLyaponovExponent : public Module
{
public:
  FiniteTimeLyaponovExponent(GuiContext* ctx);
  virtual ~FiniteTimeLyaponovExponent();

  virtual void execute();

protected:
  GuiDouble gui_stepsize_;
  GuiInt    gui_numsteps_;
  GuiDouble gui_delta_;
  GuiInt    gui_dim_;
  GuiInt    gui_plane_;
};

class FiniteTimeLyaponovExponentAlgo
{
  public:
    Mutex lock;
    
    FiniteTimeLyaponovExponentAlgo() :
      lock("FiniteTimeLyaponovExponentAlgo") {}
  
    //! virtual interface. 
    void execute(FieldHandle& vector_src,
		 FieldHandle& point_src,
		 FieldHandle& ftle_dst,
		 double stepSize,
		 int numSteps,
 		 double delta,
		 int dim,
		 int plane,
		 ProgressReporter* pr );

  virtual ~FiniteTimeLyaponovExponentAlgo() {}

  short SIGN( double val ) { return (val < 0 ? -1 : 1 ); };
};


DECLARE_MAKER(FiniteTimeLyaponovExponent)

FiniteTimeLyaponovExponent::FiniteTimeLyaponovExponent(GuiContext* context)
  : Module("FiniteTimeLyaponovExponent", context, Source, "Fields", "FusionPSE"),
    gui_stepsize_(context->subVar("stepsize"), 0.01),
    gui_numsteps_(context->subVar("numsteps"), 250),
    gui_delta_(context->subVar("delta"), 0.01),
    gui_dim_(context->subVar("dim"), 2),
    gui_plane_(context->subVar("plane"), 1)
{
}

FiniteTimeLyaponovExponent::~FiniteTimeLyaponovExponent()
{
}

void
FiniteTimeLyaponovExponent::execute()
{
  // The vector field input is required.
  FieldHandle vector_field_input_handle;

  get_input_handle( "Input Vector", vector_field_input_handle, true );
  FieldInformation vector_fi(vector_field_input_handle);
  
  if (!(vector_fi.is_vector()))
  {
    error("Only available for vector data.");
    return;
  }

  // The seed field input is required.
  FieldHandle seed_field_input_handle;

  get_input_handle( "Input Slice", seed_field_input_handle, true );
  FieldInformation seed_fi(seed_field_input_handle);
  
//   if (!(seed_fi.is_quadsurfmesh()) && !(seed_fi.is_structquadsurfmesh()))
//   { 
//     error("Only available for quad surfaces.");
//     return;
//   }


  std::cerr << "FiniteTimeLyaponovExponent executing " << std::endl;

  // If no data or a changed recalcute.
  if( inputs_changed_ ||

      !oport_cached("Output FTLE") ||

      gui_stepsize_.changed( true ) ||
      gui_numsteps_.changed( true ) ||
      gui_delta_.changed( true ) ||
      gui_dim_.changed( true ) ||
      gui_plane_.changed( true ) )
  {
    update_state( Executing );

    std::vector< std::pair< unsigned int, unsigned int > > topology;

    FieldHandle ftle_field_output_handle;

    FiniteTimeLyaponovExponentAlgo algo;

    algo.execute(vector_field_input_handle,
		 seed_field_input_handle,
		 ftle_field_output_handle,
		 gui_stepsize_.get(),
		 gui_numsteps_.get(),
		 gui_delta_.get(),
		 gui_dim_.get(),
		 gui_plane_.get(),
		 this);

    // Send the data downstream
    send_output_handle( "Output FTLE", ftle_field_output_handle );

    std::cerr << "FiniteTimeLyaponovExponent done " << std::endl;
  }
}


void
FiniteTimeLyaponovExponentAlgo::
execute(FieldHandle& vfield_h,
	FieldHandle& sfield_h,
	FieldHandle& ofield_h,
	double stepSize,
	int numSteps,
	double delta,
	int dim,
	int plane,
	ProgressReporter* pr)
{
  FieldInformation fo(sfield_h);
  fo.make_double();
  
  ofield_h = CreateField(fo,sfield_h->mesh());
  ofield_h->copy_properties(sfield_h.get_rep());
  ofield_h->resize_fdata();

  vfield_h->vmesh()->synchronize(Mesh::EPSILON_E |
				 Mesh::ELEM_LOCATE_E |
				 Mesh::EDGES_E |
				 Mesh::FACES_E);

  VField *ofield = ofield_h->vfield();
  VMesh  *omesh =  ofield->vmesh();

  unsigned int numnodes = omesh->num_nodes();

  VMesh::Node::index_type n1 = 0;
  VMesh::Node::iterator onodeItr, onodeEnd;

  omesh->begin( onodeItr );
  omesh->end( onodeEnd );

  if(onodeItr == onodeEnd) return;

  Vector vecList[3][2];
  vecList[0][0] = Vector(-1, 0, 0);
  vecList[0][1] = Vector( 1, 0, 0);
  vecList[1][0] = Vector( 0,-1, 0);
  vecList[1][1] = Vector( 0, 1, 0);
  vecList[2][0] = Vector( 0, 0,-1);
  vecList[2][1] = Vector( 0, 0, 1);

  int offset;

  if( dim == 2 ) 
    offset = 2;
  else // if( dim == 3 ) 
    offset = 1;

  Vector ptList[3][2];

  double grad[3][3], D[3][3];

  StreamLineIntegrators SI;
  SI.nodes_.reserve(numSteps);              // storage for points
  SI.max_steps_    = numSteps;              // max number of steps
  SI.vfield_       = vfield_h->vfield();    // the vector field  

  while (onodeItr != onodeEnd) 
  {
    // Get the next point in mesh.
    Point basePt;
    omesh->get_center(basePt, *onodeItr);

    double ftle = 0;

    // Search in the negative direction to find joining structures and
    // search in the posative direction to find diverging structures
    for( int dir = -1; dir <= 1; dir+=2 )
    {
      SI.step_size_ = (double) dir * stepSize; // step size

      // Get an inital streamline to see how many integration steps
      // are possible.
      SI.seed_ = basePt;
      SI.max_steps_ = numSteps;
      SI.nodes_.clear();
      SI.integrate( 0 ); // Adams Bashforth

      // If less than the desired number of integration steps then use
      // 10% than that for the neighbor integration.
      if( SI.nodes_.size() < (unsigned int) SI.max_steps_ )
      {
	SI.max_steps_ = (unsigned int) ((double) SI.nodes_.size() * 0.9);
      }
    
      if( n1 == numnodes/2 )
      {
	cerr << "basePt  " << basePt << endl;
      }

      bool valid = true;
    
      // Advect four or six points surrounding the base point.
      for( unsigned int j=0; j<3; j+=offset )
      {
	// Skip the direction not in the plane
	if( dim == 2 && (int) j == plane )
	  continue;

	for( unsigned int i=0; i<2; ++i )
        {
	  SI.seed_ = basePt + (delta * vecList[j][i]);
	  
	  SI.nodes_.clear();
	  SI.integrate( 0 ); // Adams Bashforth
	  
	  // Each streamline must complete the same number of
	  // integration steps if not then skip it.
	  if( SI.nodes_.size() == (unsigned int) SI.max_steps_ )
	  {
	    ptList[j][i] = (Vector) SI.nodes_[SI.nodes_.size()-1];
	    
	    if( n1 == numnodes/2 )
	    {
	      cerr << j << "  " << i << "  " << ptList[j][i] << endl;
	    }
	  }
	  else
	  {
	    cerr << "Not enough points from the neighbor streamline "
		 << SI.max_steps_ << "  " << SI.nodes_.size() << endl;
	    
	    i = j = 2;
	    valid = false;
	  }
	}
      }

      if( dim == 2 && valid )
      {
	if( plane == 2 ) // X-Y plane
	{
	  // Calculate the gradient using central differences.
	  grad[0][0] = (ptList[0][1].x() - ptList[0][0].x()) / (2.0*delta);
	  grad[0][1] = (ptList[1][1].x() - ptList[2][0].x()) / (2.0*delta);
	  
	  grad[1][0] = (ptList[0][1].y() - ptList[0][0].y()) / (2.0*delta);
	  grad[1][1] = (ptList[1][1].y() - ptList[1][0].y()) / (2.0*delta);
	}

	else if( plane == 1 ) // X-Z plane.
	{
	  // Calculate the gradient using central differences.
	  grad[0][0] = (ptList[0][1].x() - ptList[0][0].x()) / (2.0*delta);
	  grad[0][1] = (ptList[2][1].x() - ptList[2][0].x()) / (2.0*delta);
	  
	  grad[1][0] = (ptList[0][1].z() - ptList[0][0].z()) / (2.0*delta);
	  grad[1][1] = (ptList[2][1].z() - ptList[2][0].z()) / (2.0*delta);
	}

	else // if( plane == 0 ) // Y-Z plane.
	{
	  // Calculate the gradient using central differences.
	  grad[0][0] = (ptList[1][1].y() - ptList[1][0].y()) / (2.0*delta);
	  grad[0][1] = (ptList[2][1].y() - ptList[2][0].y()) / (2.0*delta);
	  
	  grad[1][0] = (ptList[1][1].z() - ptList[1][0].z()) / (2.0*delta);
	  grad[1][1] = (ptList[2][1].z() - ptList[2][0].z()) / (2.0*delta);
	}

	// Calculate the finite time Lyapunov exponent.
	// sigma = ln( sqrt( max(eigen( gradtransposed * grad ) ) ) ) / numSteps

	// Multiply the transposed gradient and gradient together.
	D[0][0] = grad[0][0] * grad[0][0] + grad[0][1] * grad[0][1];
	D[0][1] = grad[0][0] * grad[1][0] + grad[0][1] * grad[1][1];

	D[1][0] = grad[1][0] * grad[0][0] + grad[1][1] * grad[0][1];
	D[1][1] = grad[1][0] * grad[1][0] + grad[1][1] * grad[1][1];

	double a = D[0][0]; 	double b = D[0][1];
	double c = D[1][0]; 	double d = D[1][1];

	// Get the eigen values for the matrix. For a 2x2 matrix we
	// just need to solve for the roots of a quadratic eqation.
	double B = (a + d);   // trace
	double C = a*d - b*c; // determinant

	double l0 = (B + sqrt( B*B - 4*C ) ) / 2.0;
	double l1 = (B - sqrt( B*B - 4*C ) ) / 2.0;

	// Get the natural log of the maximum eigen value
	double le = log( sqrt( max(l0,l1) ) ) / (double) SI.max_steps_;

	if( le < 0 ) 
	  le = 0;

	// Get the greatest Lyapunov exponent between the two
	// directions. Negative direction means joining so positive
	// Lyapunov exponent positive dir means diverging so negative
	// Lyapunov exponent.
	if( fabs( ftle ) < le )
	  ftle = -1.0 * dir * le;

	if( n1 == numnodes/2 )
        {
	  cerr << l0 << "  " << l1 << "    " << le << endl;
	}
      }

      else if( dim == 3 && valid )
      {
	// Calculate the gradient using central differences.
	grad[0][0] = (ptList[0][1].x() - ptList[0][0].x()) / (2.0*delta);
	grad[0][1] = (ptList[1][1].x() - ptList[1][0].x()) / (2.0*delta);
	grad[0][2] = (ptList[2][1].x() - ptList[2][0].x()) / (2.0*delta);

	grad[1][0] = (ptList[0][1].y() - ptList[0][0].y()) / (2.0*delta);
	grad[1][1] = (ptList[1][1].y() - ptList[1][0].y()) / (2.0*delta);
	grad[1][2] = (ptList[2][1].y() - ptList[2][0].y()) / (2.0*delta);

	grad[2][0] = (ptList[0][1].z() - ptList[0][0].z()) / (2.0*delta);
	grad[2][1] = (ptList[1][1].z() - ptList[1][0].z()) / (2.0*delta);
	grad[2][2] = (ptList[2][1].z() - ptList[2][0].z()) / (2.0*delta);
      
	// Calculate the finite time Lyapunov exponent.
	// sigma = ln( sqrt( max(eigen( gradtransposed * grad ) ) ) ) / numSteps

	// Multiply the transposed gradient and gradient together.
	D[0][0] = grad[0][0] * grad[0][0] + grad[0][1] * grad[0][1] + grad[0][2] * grad[0][2];
	D[0][1] = grad[0][0] * grad[1][0] + grad[0][1] * grad[1][1] + grad[0][2] * grad[1][2];
	D[0][2] = grad[0][0] * grad[2][0] + grad[0][1] * grad[2][1] + grad[0][2] * grad[2][2];

	D[1][0] = grad[1][0] * grad[0][0] + grad[1][1] * grad[0][1] + grad[1][2] * grad[0][2];
	D[1][1] = grad[1][0] * grad[1][0] + grad[1][1] * grad[1][1] + grad[1][2] * grad[1][2];
	D[1][2] = grad[1][0] * grad[2][0] + grad[1][1] * grad[2][1] + grad[1][2] * grad[2][2];

	D[2][0] = grad[2][0] * grad[0][0] + grad[2][1] * grad[0][1] + grad[2][2] * grad[0][2];
	D[2][1] = grad[2][0] * grad[1][0] + grad[2][1] * grad[1][1] + grad[2][2] * grad[1][2];
	D[2][2] = grad[2][0] * grad[2][0] + grad[2][1] * grad[2][1] + grad[2][2] * grad[2][2];

	double a = D[0][0]; double b = D[0][1]; double c = D[0][2];
	double d = D[1][1]; double e = D[1][1]; double f = D[1][1];
	double g = D[1][1]; double h = D[1][1]; double i = D[1][1];

	// Get the eigen values for the matrix. For a 3x3 matrix we
	// just need to solve for the roots of a cubic eqation.
	// A - trace
	// B - trace(matrix^2) - trace ^ 2;
	// C - determinant;

	double A = -(a + e + i);
	double B = -(d*b + g*c + f*h - a*e - a*i - e*i);
	double C = -(a*e*i - a*f*h - d*b*i + d*c*h + g*b*f - g*c*e);

	// Roots of a cubic equation.
	double Q = (A*A - 3.0*B) / 9.0;
	double R = (2.0*A*A*A - 9.0*A*B + 27.0*C) / 54.0;

	double theta = acos( R / sqrt (Q*Q*Q) );

	double l0 = -2.0 * sqrt( Q ) * cos (theta / 3) - A / 3.0;
	double l1 = -2.0 * sqrt( Q ) * cos ((theta+2.0*M_PI) / 3) - A / 3.0;
	double l2 = -2.0 * sqrt( Q ) * cos ((theta-2.0*M_PI) / 3) - A / 3.0;

	// Get the natural log of the maximum eigen value
	double le = log( sqrt( max(max(l0,l1),l2) ) ) / (double) SI.max_steps_;

	if( le < 0 ) 
	  le = 0;

	// Get the greatest Lyapunov exponent between the two
	// directions.  Negative direction means joining so positive
	// Lyapunov exponent positive dir means diverging so negative
	// Lyapunov exponent.
	if( fabs( ftle ) < le )
	  ftle = -1.0 * dir * le;

	if( n1 == numnodes/2 )
        {
	  string tmpstr;
	
	  cerr << tmpstr << endl
	       << l0 << "  " << l1 << "  " << l2 << "    " << le << endl;
	}
      }

      ofield->set_value( ftle, n1 );
    }

    ++onodeItr;
    ++n1;

    pr->update_progress(n1,numnodes);
  }


}

} // End namespace FusionPSE
