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
 *  CreateImage.cc:  Make an ImageField that fits the source field.
 *
 *  Written by:
 *   Michael Callahan
 *   Department of Computer Science
 *   University of Utah
 *   March 2001
 *
 */

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Core/Math/MiscMath.h> // for M_PI
#include <iostream>

namespace SCIRun {

class CreateImage : public Module
{
public:
  CreateImage(GuiContext* ctx);
  virtual ~CreateImage() {}

  virtual void execute();

private:

  GuiInt size_x_;
  GuiInt size_y_;
  GuiInt size_z_;
  GuiInt z_value_;
  GuiInt auto_size_;
  GuiInt axis_;
  GuiDouble padpercent_;
  GuiDouble position_;
  GuiString data_at_;
  GuiString update_type_;
  GuiPoint custom_origin_;
  GuiVector custom_normal_;

  enum DataTypeEnum { SCALAR, VECTOR, TENSOR };
};


DECLARE_MAKER(CreateImage)
  
CreateImage::CreateImage(GuiContext* ctx) : 
  Module("CreateImage", ctx, Filter, "NewField", "SCIRun"),
  size_x_(get_ctx()->subVar("sizex"), 20),
  size_y_(get_ctx()->subVar("sizey"), 20),
  size_z_(get_ctx()->subVar("sizez"), 2),
  z_value_(get_ctx()->subVar("z_value"), 0),
  auto_size_(get_ctx()->subVar("auto_size"), 0),
  axis_(get_ctx()->subVar("axis"), 0),
  padpercent_(get_ctx()->subVar("padpercent"), 0),
  position_(get_ctx()->subVar("pos"), 0),
  data_at_(get_ctx()->subVar("data-at"), "Nodes"),
  update_type_(get_ctx()->subVar("update_type"), "On Release"),
  custom_origin_(get_ctx()->subVar("corigin"), Point(0.0, 0.0, 0.0)),
  custom_normal_(get_ctx()->subVar("cnormal"), Vector(1.0, 1.0, 1.0))
{
}

void
CreateImage::execute()
{
  update_state(NeedData);
  const int axis = Min(2, Max(0, axis_.get()));

  Transform trans;
  trans.load_identity();
	
	//checking for input matrices
	MatrixHandle size_matrix, ov_matrix;
	const bool size_matrix_p =get_input_handle("Image Size",size_matrix,false);
	const bool ov_matrix_p =get_input_handle("Custom Center and Normal",ov_matrix,false);
	FieldHandle ifieldhandle;
	const bool ifield_p=get_input_handle("Input Field", ifieldhandle, false);
	
	if (size_matrix_p)
	{
		if (size_matrix->get_data_size() == 1)
		{
			double* data = size_matrix->get_data_pointer();
			unsigned int size = static_cast<unsigned int>(data[0]);
			size_x_.set(size);
			size_y_.set(size);
			//get_ctx()->reset();
			
		}
		else if (size_matrix->get_data_size() == 2)
		{
			double* data = size_matrix->get_data_pointer();
			unsigned int size1 = static_cast<unsigned int>(data[0]);
			unsigned int size2 = static_cast<unsigned int>(data[1]);
			size_x_.set(size1);
			size_y_.set(size2);
			//get_ctx()->reset();
		}
		else 
		{
			error("Image Size matrix must have only 1 or 2 elements");
		}
	}
  
	if (ov_matrix_p)
	{
		if (ov_matrix->nrows() != 2 || ov_matrix->ncols() != 3)
		{
			error("Custom Center and Normal matrix must be of size 2x3.  The Center is the first row and the normal is the second");
		}
		
		double* data=ov_matrix->get_data_pointer();
		custom_origin_.set(Point(data[0],data[1],data[2]));
		Vector vect=Vector(data[3],data[4],data[5]);
		vect.safe_normalize();
		custom_normal_.set(vect);
		//get_ctx()->reset();
	}	
	
	double angle = 0;
  Vector axis_vector(0.0, 0.0, 1.0);
  switch (axis)
  {
  case 0:
    angle = M_PI * -0.5; 
    axis_vector = Vector(0.0, 1.0, 0.0);
    break;

  case 1:
    angle = M_PI * 0.5; 
    axis_vector = Vector(1.0, 0.0, 0.0);
    break;

  case 2:
    angle = 0.0;
    axis_vector = Vector(0.0, 0.0, 1.0);
    break;

  default:
    break;
  }
  trans.pre_rotate(angle, axis_vector);

  if (axis_.get() == 3)
  {
    Vector tmp_normal(-custom_normal_.get());
    Vector fakey(Cross(Vector(0.0, 0.0, 1.0), tmp_normal));
    if (fakey.length2() < 1.0e-6)
    {
      fakey = Cross(Vector(1.0, 0.0, 0.0), tmp_normal);
    }
    Vector fakex(Cross(tmp_normal, fakey));
    tmp_normal.safe_normalize();
    fakex.safe_normalize();
    fakey.safe_normalize();
		double dg=1.0;
		if (ifield_p)
		{
			BBox box = ifieldhandle->vmesh()->get_bounding_box();
			Vector diag(box.diagonal());
			dg=diag.maxComponent();
			trans.pre_scale(Vector(dg,dg,dg));
		}
		Transform trans2;
    trans2.load_identity();
		trans2.load_basis(Point(0, 0, 0), fakex, fakey, tmp_normal);
		trans2.invert();
		trans.change_basis(trans2);
    const Vector &origin(custom_origin_.get().asVector());
		trans.pre_translate(origin);
  }

  
  DataTypeEnum datatype;
  unsigned int sizex, sizey, sizez;
  if (!ifield_p)
  {
    update_state(Executing);
    datatype = SCALAR;  
    // Create blank mesh.
    sizex = Max(2, size_x_.get());
    sizey = Max(2, size_y_.get());
  }
  else
  {
    update_state(Executing);
    datatype = SCALAR;
    FieldInformation fi(ifieldhandle);
    if (fi.is_tensor())
    {
      datatype = TENSOR;
    }
    else if (fi.is_vector())
    {
      datatype = VECTOR;
    }
  
    int basis_order = 1;
    if( auto_size_.get() ) 
    {   // Guess at the size of the sample plane.
      // Currently we have only a simple algorithm for LatVolFields.

      if (fi.is_latvolmesh())
      {
        VMesh *lvm = ifieldhandle->vmesh();
        basis_order = ifieldhandle->vfield()->basis_order();
        
        switch( axis ) 
        {
          case 0:
            sizex = Max(2, (int)lvm->get_nj());
            size_x_.set( sizex );
            sizey = Max(2, (int)lvm->get_nk());
            size_y_.set( sizey );
            sizez = Max(2, (int)lvm->get_ni());
            if( basis_order == 0 )
            {
              size_z_.set( sizez - 1 );
            } 
            else 
            {
              size_z_.set( sizez );
            }
            TCLInterface::execute(get_id()+" edit_scale");
            break;
          case 1: 
            sizex =  Max(2, (int)lvm->get_ni());
            size_x_.set( sizex );
            sizey =  Max(2, (int)lvm->get_nk());
            size_y_.set( sizey );
            sizez = Max(2, (int)lvm->get_nj());
            if( basis_order == 0 )
            {
              size_z_.set( sizez - 1 );
            } 
            else 
            {
              size_z_.set( sizez );
            }
            TCLInterface::execute(get_id()+" edit_scale");
            break;
          case 2:
            sizex =  Max(2, (int)lvm->get_ni());
            size_x_.set( sizex );
            sizey =  Max(2, (int)lvm->get_nj());
            size_y_.set( sizey );
            sizez =  Max(2, (int)lvm->get_nk());
            if( basis_order == 0 )
            {
              size_z_.set( sizez - 1 );
            } 
            else 
            {
              size_z_.set( sizez );
            }
            TCLInterface::execute(get_id()+" edit_scale");
            break;
          default:
            warning("Custom axis, resize manually.");
            sizex = Max(2, size_x_.get());
            sizey = Max(2, size_y_.get());
            break;
        }
      } 
      else 
      {
        warning("No autosize algorithm for this field type, resize manually.");
        sizex = Max(2, size_x_.get());
        sizey = Max(2, size_y_.get());
        auto_size_.set(0);
        TCLInterface::execute(get_id()+" edit_scale");
      }
    } 
    else 
    {
      // Create blank mesh.
      sizex = Max(2, size_x_.get());
      sizey = Max(2, size_y_.get());
    }

    

    if (axis_.get() != 3)
    {
			// Compute Transform.
			BBox box = ifieldhandle->vmesh()->get_bounding_box();
			Vector diag(box.diagonal());
			trans.pre_scale(diag);
			
      Point loc(box.center());
      position_.reset();
      double dist;
      if ( !auto_size_.get() ) 
      {
        dist = position_.get()/2.0;
      } 
      else 
      {
        if( basis_order == 0 ) 
        {
          dist = double( z_value_.get() )/size_z_.get() + 0.5/size_z_.get();
					position_.set( ( dist) * 2.0 );
        } 
        else 
        {
          dist = double( z_value_.get() )/size_z_.get();
					position_.set( ( dist) * 2.0 );
        }
      }
      switch (axis)
      {
      case 0:
        loc.x(loc.x() + diag.x() * dist);
        break;

      case 1:
        loc.y(loc.y() + diag.y() * dist);
        break;

      case 2:
        loc.z(loc.z() + diag.z() * dist);
        break;
      
      default:
        break;
      }
      trans.pre_translate(Vector(loc));
    }
  }

  Point minb(-0.5, -0.5, 0.0);
  Point maxb(0.5, 0.5, 0.0);
  Vector diag((maxb.asVector() - minb.asVector()) * (padpercent_.get()/100.0));
  minb -= diag;
  maxb += diag;
	
  int basis_order;
  if (data_at_.get() == "Nodes") basis_order = 1;
  else if (data_at_.get() == "Faces") basis_order = 0;
  else if (data_at_.get() == "None") basis_order = -1;
  else 
  {
    error("Unsupported data_at location " + data_at_.get() + ".");
    return;
  }
  
  FieldInformation ifi("ImageMesh",basis_order,"double");
  
  if (datatype == VECTOR) ifi.make_vector();
  else if (datatype == TENSOR) ifi.make_tensor();  
  
  MeshHandle imagemesh = CreateMesh(ifi,sizex, sizey, minb, maxb);
  FieldHandle ofh = CreateField(ifi,imagemesh);
  
  // Transform field.
  ofh->vmesh()->transform(trans);

  send_output_handle("Output Sample Field", ofh);
}


} // End namespace SCIRun

