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


#include <Core/Algorithms/Fields/ROI/MeshROI.h>

//! For mapping matrices
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Geometry/Vector.h>
//#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/FieldInformation.h>

//! STL classes needed
#include <sci_hash_map.h>
#include <algorithm>
#include <set>
#include <cmath>
#include <limits>

//! BOOST classes needed
#include <boost/lexical_cast.hpp>



namespace SCIRunAlgo {

	// Distance function

	double distance(Point from, Point to)
	{
		sqrt( pow(from[0]-to[0],2) + pow(from[1]-to[1],2) + pow(from[2]-to[2],2) );
	}
	
	void trace_neigb(VMesh* mesh, VMesh::Elem::index_type elem_idx, size_t top, std::set<VMesh::Elem::index_type>& traced_elems)
	{
		VMesh::Elem::array_type neigb;
		
		// bound of our topologival distance range
		if(top)
		{

			if(traced_elems.insert(elem_idx).second)
			{
			  //std::cout << "non-dublicate" << std::endl;
			}
			else
			{
			  //std::cout << "dublicate" << std::endl;
			}

			mesh->get_neighbors(neigb,elem_idx);


			for (size_t p = 0; p < neigb.size(); p++)
			{
				trace_neigb(mesh,neigb[p],top-1,traced_elems);
			}
		  
		}
	}

	// General access function
	bool RoiMeshAlgo::run(FieldHandle& input, MatrixHandle& ref, FieldHandle& output)
	{
		algo_start("ROIMesh");

		if (input.get_rep() == 0)
		{
		error("No input field was given.");
		algo_end(); return (false);
		}

		std::string select; 
		get_option("select",select);

		int topo_dist;
		get_int("distance",topo_dist);

		if (select == "")
		{
			// just passing through
			output = input;
			algo_end(); 
			return (true); 
		}
		
		if (select == "line")
		{
			error("NA");
			algo_end();
			return false;
		}


		Matrix<double>* matref = ref.get_rep();
		VField* infield   = input->vfield();
		VMesh*  inmesh    = input->vmesh();
		//VMesh*  outmesh = output->vmesh();
		//VField* rfield  = output->vfield();
		  
		VMesh::Elem::array_type elems;
		VMesh::Elem::iterator it, eit;
		VMesh::Elem::size_type sz;

		VMesh::Node::array_type nodearray;
		//VMesh::Node::array_type onodes(4);

		VMesh::size_type num_nodes = inmesh->num_nodes();
		VMesh::size_type num_elems = inmesh->num_elems();


		status( "--DEBUG: topodist : " + boost::lexical_cast<std::string>(topo_dist) );
		status( "--DEBUG: selected ROI : " + select );
		status( "--DEBUG: elements : " + boost::lexical_cast<std::string>(num_elems) );

		std::set<VMesh::Elem::index_type> traced_elems;
		std::vector<char> filter_values(num_elems,0);


		////////////////////////////////////////////////////////////////////
		inmesh->synchronize(Mesh::ALL_ELEMENTS_E);
		////////////////////////////////////////////////////////////////////


		inmesh->begin(it);
		inmesh->end(eit);
		inmesh->size(sz);
		//VMesh::index_type cnt = 0, c = 0;

		double min_l2_dist = std::numeric_limits<double_t>::max();

		VMesh::Elem::index_type elem_rori = 0;

		

		if(select == "point")
		{
			Point from(matref->get(0,0),matref->get(0,1),matref->get(0,2));

			while (it != eit)
			{
				inmesh->get_nodes(nodearray, *it);

				size_t nsize = nodearray.size();

				Point to;

				inmesh->get_center(to, *it);

				double d = distance(from,to);
				//std::cout << from << to << d << std::endl;

				if(d < min_l2_dist)
				{
					min_l2_dist = d;
					elem_rori = *it;	  
				}
			  

			  ++it;
			  
			  
			  //cnt++; 
			  //if (cnt==1000) 
			  //{ 
			  //  cnt=0; c+=1000; 
			   // algo->update_progress(c,sz); 
			   // if (algo->check_abort()) break;
			  //}
			}
		}
		else if(select == "line")
		{
			while (it != eit)
			{
				inmesh->get_nodes(nodearray, *it);

				size_t nsize = nodearray.size();

				
				
				Point to;

				inmesh->get_center(to, *it);
				
				Point from1(matref->get(0,0),matref->get(0,1),matref->get(0,2));
				Point from2(matref->get(1,0),matref->get(1,1),matref->get(1,2));				

				double d = distance(from1,to);
				
				Point p0, p1, p2;
				inmesh->get_point(p0,nodearray[0]);
				inmesh->get_point(p1,nodearray[1]);
				inmesh->get_point(p2,nodearray[2]);
				
				Vector v1(p1 - p0);
				Vector v2(p2 - p0);
				Vector vn = Cross(v1,v2);
				
				double t = - Dot(from1 - p0,vn) / Dot(from2, vn);
				
				//intersection poin on the plane of the triangle polygon
				Vector x = from1 + t*from2;
								
				if(d < 2*min_l2_dist)
				{
					//ray intersection test
					if( Dot(Cross(p1 - p0, x - p0), vn) >= 0.0 &&
						Dot(Cross(p2 - p1, x - p1), vn) >= 0.0 &&
						Dot(Cross(p0 - p2, x - p2), vn) >= 0.0)
					{
						elem_rori = *it;
						std::cout << " --DEBUG inside triangle \n";
					}

					if(d < min_l2_dist)
					{
						min_l2_dist = d;
						//elem_rori = *it;	  
					}
			  
				}
				
			  ++it;
			  
			  
			  //cnt++; 
			  //if (cnt==1000) 
			  //{ 
			  //  cnt=0; c+=1000; 
			   // algo->update_progress(c,sz); 
			   // if (algo->check_abort()) break;
			  //}
			}
		}
		else
		{
			error(" Unknown selected methood! Only point and line are supported.");
			return false;
		}

		// find the topolicaly nearest elements
		trace_neigb(inmesh,elem_rori,topo_dist,traced_elems);
		
		status( " --DEBUG number of ROI elements :  " + boost::lexical_cast<std::string>(traced_elems.size()) );

		// ctreate filter mesh data
		inmesh->begin(it);
		inmesh->end(eit);
		size_t iii = 0;

		while (it != eit)
		{
		  const std::set<VMesh::Elem::index_type>::iterator loc = traced_elems.find(*it);
		  
		  filter_values[iii] = loc == traced_elems.end() ? 0 : 1;
		  
		  ++iii;
		  ++it;
		}
		 

		////////////////////////////////////////////////////////////////////
		inmesh->unsynchronize(Mesh::ALL_ELEMENTS_E);
		///////////////////////////////////////////////////////////////////


		FieldInformation fis(input);
		fis.set_data_type("char");
		fis.set_basis_type(0);//constatn data on elements
		FieldHandle  sfield = CreateField(fis,input->mesh()->clone());

		// THROWS!
		//sfield->vmesh()->copy_elems(inmesh);
		sfield->vfield()->set_values(filter_values);

		if (!(sfield.get_rep()))
		{
		  error("Could not allocate output field.");
		  return (false);
		}

		MatrixHandle mapping;

		clipping_algo_.set_option("method","element");

		if(!(clipping_algo_.run(input,sfield,output,mapping))) 
		{
			return false;
		}
  

	  algo_end(); 
	  return true;
	}                           
		 
} // namespace    
