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
 *  SimulateDipoleInSphereAnalytically: Calculation of potential on 
 *                  conducting sphere due to the dipole sources
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   October 1994
 *
 *  Modified by:
 *   Samsonov Alexei
 *   Department of Computer Science
 *   University of Utah
 *   March 2001
 *
 */


#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/DenseMatrix.h>

#include <Core/Containers/Array1.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/BBox.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>

#include <iostream>
#include <stdio.h>
#include <math.h>

namespace BioPSE {
using namespace SCIRun;

class SimulateDipoleInSphereAnalytically : public Module {

  private:
    //! Private Methods
    // -- fills in the surface with potentials for single sphere uniform model
    void fillOneSphere(DenseMatrix&, VField*, VField*);

  public:
    
    SimulateDipoleInSphereAnalytically(GuiContext *context);
    virtual ~SimulateDipoleInSphereAnalytically();
    virtual void execute();
};

DECLARE_MAKER(SimulateDipoleInSphereAnalytically)


SimulateDipoleInSphereAnalytically::SimulateDipoleInSphereAnalytically(GuiContext *context)
  : Module("SimulateDipoleInSphereAnalytically", context, Filter, "Forward", "BioPSE")
{
}


SimulateDipoleInSphereAnalytically::~SimulateDipoleInSphereAnalytically()
{
}


void
SimulateDipoleInSphereAnalytically::execute()
{
  update_state(NeedData);
  
  FieldHandle field_handle;
  get_input_handle("Sphere", field_handle, true);
 
  FieldHandle dip_handle;
  get_input_handle("Dipole Sources", dip_handle, false);
 
  FieldInformation fi(field_handle);
  FieldInformation dip_fi(dip_handle);
  
  
  if (fi.is_trisurfmesh() && fi.is_lineardata())
  {
    VMesh* hMesh = field_handle->vmesh();

    fi.make_double();
    FieldHandle pfield = CreateField(fi,field_handle->mesh());
    VField* hNewSurf = pfield->vfield();

    fi.make_vector();
    FieldHandle mfield = CreateField(fi,field_handle->mesh());
    VField* hBSurf =  mfield->vfield();

    if (dip_fi.is_pointcloudmesh() && dip_fi.is_vector())
    {  
      VField* pDips = dip_handle->vfield();
      VMesh* hMesh = dip_handle->vmesh();
      
      VMesh::Node::iterator ii;
      VMesh::Node::iterator ii_end;
      
      Point p;
      Vector qdip;
      std::vector<Vector> dips;
      std::vector<Point>  pos;
      hMesh->begin(ii);
      hMesh->end(ii_end);
      for (; ii != ii_end; ++ii) 
      {
        pDips->get_value(qdip, *ii);
        dips.push_back(qdip);
        hMesh->get_point(p, *ii);
        pos.push_back(p);
      }
      
      DenseMatrix dip_mtrx((int)pos.size(), 6);
      unsigned int i;
      msg_stream_ << "Dipoles: " << endl;
      for (i=0; i<pos.size(); ++i)
      {
        qdip = dips[i];
        p = pos[i];
        dip_mtrx[i][0] = p.x(); dip_mtrx[i][1] = p.y();  dip_mtrx[i][2] = p.z();
        dip_mtrx[i][3] = qdip.x(); dip_mtrx[i][4] = qdip.y();  dip_mtrx[i][5] = qdip.z();
        msg_stream_ << "Pos: " << p << ", moment: " << qdip << endl;
      }
      
      update_state(JustStarted);
      fillOneSphere(dip_mtrx, hNewSurf, hBSurf);

      send_output_handle("SphereWithPots", pfield);
      send_output_handle("SphereWithMagneticField", mfield);
    }
    else 
    {
      warning("No dipole info found in the mesh supplied or supplied field is not of type PointCloudField<Vector>.");
    }
  }
  else 
  {
    error("The supplied field is not of type TriSurfField<double>.");
    return;
  }
}


void
SimulateDipoleInSphereAnalytically::fillOneSphere(DenseMatrix& dips, VField* hSurf, 
                                                  VField* hBSurf) 
{  
  VMesh* hMesh = hSurf->vmesh();

//  std::vector<double>& data = hSurf->fdata();
//  std::vector<Vector>& bdata = hBSurf->fdata();
  
  VMesh::Node::size_type nsize; hMesh->size(nsize);

  BBox bbox = hMesh->get_bounding_box();
  
  if (!bbox.valid())
  {
    error("No valid input mesh.");
    return;
  }

  double R = 0.5*bbox.longest_edge();
  
  double gamma=1;
  double E[3];
  msg_stream_ << "Radius of the sphere is " << R << std::endl;
  Point p;

  VMesh::Node::iterator niter; hMesh->begin(niter);
  VMesh::Node::iterator niter_end; hMesh->end(niter_end);

  // -- for every point
  while (niter != niter_end) 
  {
    hSurf->set_value(0.0,*niter);
    hBSurf->set_value(Vector(0,0,0),*niter);
    hMesh->get_point(p, *niter);
      
    // -- for every dipole
    int id;
    for (id = 0; id < dips.nrows(); ++id)
    {
      double V = 0.0;
      E[0] = p.x();
      E[1] = p.y();
      E[2] = p.z();

      double rho = sqrt( pow((E[0] - dips[id][0]),2) + pow((E[1] - dips[id][1]),2) + pow((E[2] - dips[id][2]),2));
      double S = E[0]*dips[id][0] + E[1]*dips[id][1] + E[2]*dips[id][2];
	
      for(int k=0;k<3;k++) 
      {
        double F[3];
        F[k] = (1/(4*M_PI*gamma*rho)) * (2*(E[k]-dips[id][k])/pow(rho,2) +
                 (1/pow(R,2)) * (E[k] + (E[k]*S/R - R*dips[id][k])/(rho+R-S/R)));
        V += F[k]*dips[id][k+3];
      }
	
      double data;
      hSurf->get_value(data,*niter);
      hSurf->set_value(data+V,*niter);

      // magnetic
      Vector r_vec = p.asVector()*1.01;
      Vector r0_vec = Vector(dips[id][0], dips[id][1], dips[id][2]);
      double r_mag = r_vec.length();
      Vector a_vec = r_vec-r0_vec;
      double a_mag = a_vec.length();
      Vector Q_vec(dips[id][3], dips[id][4], dips[id][5]);
      double F_mag = a_mag*(r_mag*a_mag + r_mag*r_mag - Dot(r0_vec, r_vec));
      Vector gradF_vec((a_mag*a_mag/r_mag + Dot(a_vec,r_vec)/a_mag + 2*a_mag +
                        2*r_mag)*r_vec - (a_mag + 2*r_mag + Dot(a_vec,r_vec)/a_mag)*r0_vec);
      
      Vector vec;
      hBSurf->get_value(vec,*niter);
      vec += (F_mag*(Cross(Q_vec,r0_vec)) - 
              Dot(Cross(Q_vec,r0_vec),r_vec)*gradF_vec)/(4*M_PI*F_mag*F_mag);
      hBSurf->set_value(vec,*niter);
    }
    ++niter;
  }
}

} // End namespace BioPSE
