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

#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Mesh.h>

#include <Packages/CardioWaveInterface/Core/Model/BuildMembraneTable.h>

#include <sci_hash_map.h>
#include <algorithm>

namespace CardioWaveInterface {

using namespace SCIRun;

bool BuildMembraneTableAlgo::BuildMembraneTable(ProgressReporter *pr, FieldHandle ElementType, FieldHandle MembraneModel, MatrixHandle CompToGeom, MatrixHandle NodeLink, MatrixHandle ElemLink, MembraneTable& membranetable, MatrixHandle& MappingMatrix, MatrixHandle& MappingMatrix2, Matrix::index_type &offset)
{
  if (ElementType.get_rep() == 0)
  {
    pr->error("BuildMembraneTable: No element type field");
    return (false);
  }

  if (MembraneModel.get_rep() == 0)
  {
    pr->error("BuildMembraneTable: No membrane model field");
    return (false);
  }

  // no precompiled version available, so compile one

  FieldInformation fi(ElementType);
  FieldInformation fi2(MembraneModel);
  
  if (!(fi.is_constantdata()))
  {
    pr->error("BuildMembraneTable: The element type field needs to have one data value assigned to each element");
    return (false);
  }
  
  if (!fi.is_unstructuredmesh())
  {
    pr->error("BuildMembraneTable: This function is not defined for structured meshes");
    return (false);
  }  

  if (!((fi.is_volume()&&fi2.is_surface())||(fi.is_surface()&&fi2.is_curve())))
  {
    pr->error("BuildMembraneTable: The element type field and the MembraneModelace needs to be a volumea and a surface OR a surface and a curve");
    return (false);
  }  
  
  if (!fi2.is_unstructuredmesh())
  {
    pr->error("BuildMembraneTable: This function is not defined for structured meshes");
    return (false);
  }  

  // Setup dynamic files

  VField* elementtypefield = ElementType->vfield();
  if (elementtypefield == 0)
  { 
    pr->error("BuildMembraneTable: Could not obtain element class description field");
    return (false);
  }

  // 2) Check whether we can obtain the mesh of the elementtype field
  VMesh* elementtypemesh = ElementType->vmesh();
  if (elementtypemesh == 0)
  {
    pr->error("BuildMembraneTable: No mesh associated with element class description field");
    return (false);
  }

  // 3) Check whether the membranefield is valid
  VField* membranefield = MembraneModel->vfield();
  if (membranefield == 0)
  { 
    pr->error("BuildMembraneTable: Could not obtain model description field");
    return (false);
  }

  // 4) Same for the mesh
  VMesh* membranemesh = MembraneModel->vmesh();
  if (membranemesh == 0)
  {
    pr->error("BuildMembraneTable: No mesh associated with model description field");
    return (false);
  }

  //Get the sizes of both fields
  
  if (elementtypemesh->dimensionality() == 3) elementtypemesh->synchronize(Mesh::NODE_NEIGHBORS_E|Mesh::FACES_E);
  if (elementtypemesh->dimensionality() == 2) elementtypemesh->synchronize(Mesh::NODE_NEIGHBORS_E|Mesh::EDGES_E);
  if (elementtypemesh->dimensionality() == 1) elementtypemesh->synchronize(Mesh::NODE_NEIGHBORS_E|Mesh::NODES_E);
  
  VMesh::Node::size_type numelementnodes;
  VMesh::Node::size_type nummembranenodes;
  VMesh::DElem::size_type numelementdelems;

  elementtypemesh->size(numelementdelems);
  elementtypemesh->size(numelementnodes);
  membranemesh->size(nummembranenodes);
  
  // If there are no nodes there is no point in persuing this effort...
  if ((numelementnodes == 0)||(nummembranenodes == 0))
  {
    pr->error("BuildMembraneTable: One of the input meshes has no elements");
    return (false);  
  }

  // CompToGeom: Geometric nodes to computational nodes
  // ElemLink: Border elements to Opposing border elements
  // In this case only where we have two different domains in the outerboundary link
  // these need to be membranes as well.

  // Check whether CompToGeom was defined on the input, if not
  // maybe it is a property, else we ignore it...
  // it is an optional parameter
  if (CompToGeom.get_rep() == 0)
  {
    elementtypefield->get_property("CompToGeom",CompToGeom);
  }
  
  // Same for ElemLink
  if (ElemLink.get_rep() == 0)
  {
    elementtypefield->get_property("ElemLink",ElemLink);
  }

  if (NodeLink.get_rep() == 0)
  {
    elementtypefield->get_property("NodeLink",NodeLink);
  }  

  // VVV ALGORITHM STARTS HERE VVV

  // Make sure we have everything the mesh needs, compute all internal structures
  // needed
  

#ifdef HAVE_HASH_MAP
  typedef hash_multimap<int,VMesh::Node::index_type> nodemap_type;
#else
  typedef multimap<int,VMesh::Node::index_type> nodemap_type;
#endif

  int key;
  
  nodemap_type nodemap;
  Point point, point2, point3;

  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  double zmin = 0.0;
  double zmax = 0.0;
    
  VMesh::Node::iterator it, it_end;
  elementtypemesh->begin(it);
  elementtypemesh->end(it_end);
  
  // Compute a bounding box...
	
  // Initiate the parameters in a single loop
  if (it != it_end)
  {
    elementtypemesh->get_center(point,*it);
    xmin = point.x(); xmax = point.x();
    ymin = point.y(); ymax = point.y();
    zmin = point.z(); zmax = point.z();
    ++it;
  }
    
   // now compute the full bounding box.
  while (it != it_end)
  {
    elementtypemesh->get_center(point,*it);
    if (point.x() < xmin) xmin = point.x();
    if (point.x() > xmax) xmax = point.x();
    if (point.y() < ymin) ymin = point.y();
    if (point.y() > ymax) ymax = point.y();
    if (point.z() < zmin) zmin = point.z();
    if (point.z() > zmax) zmax = point.z();
    ++it;
  }

  // TO DO:
  // Need to check whether bounding box is already available in all meshes..
  // It is a small overhead but may help reduce computational times

  // Define multipliers for the grid we are putting on top
  // Unfortunately the internal mesh grid is not accessible for every mesh type
  // hence we have to redo it. (Design error in SCIRun Mesh classes!)
  double xmul = 0.0; if ((xmax-xmin) > 0.0 ) xmul = 250/(xmax-xmin);
  double ymul = 0.0; if ((ymax-ymin) > 0.0 ) ymul = 250/(ymax-ymin);
  double zmul = 0.0; if ((zmax-zmin) > 0.0 ) zmul = 250/(zmax-zmin);

  xmin -= (xmax-xmin)*0.01;
  ymin -= (ymax-ymin)*0.01;
  zmin -= (zmax-zmin)*0.01;
                  
  elementtypemesh->begin(it);
  elementtypemesh->end(it_end);
  
  
  // Compute a key for each node in the mesh, so we can quickly find nodes when
  // we are looking for them. This is memory overhead, but otherwise this procedure
  // is extremely slow and we cannot use the internal locate function as we allow for
  // multiple nodes to be at the same position, hence the accelarator functions
  // in SCIRun are useless here, hence we need to do our own accelation.
  while (it != it_end)
  {
    elementtypemesh->get_center(point,*it);
    
    key = static_cast<int>((point.x()-xmin)*xmul);
    key += (static_cast<int>((point.y()-ymin)*ymul))<<8;
    key += (static_cast<int>((point.z()-zmin)*zmul))<<16;  
  
    nodemap.insert(nodemap_type::value_type(key,*it));
    ++it;
  }

  // Process ElemLink: this matrix tells which faces are connected by a 
  // membrane to faces on the opposite side of the mesh

  // Assume we do not have it
  bool iselemlink = false;
  Matrix::index_type* elemlinkrr = 0;
  Matrix::index_type* elemlinkcc = 0;

  
  if (ElemLink.get_rep())
  {
    // We have a ElemLink Matrix
    
    // Do a sanity check, if not return a proper error
    if ((numelementdelems != ElemLink->nrows())&&(numelementdelems != ElemLink->ncols()))
    {
      pr->error("BuildMembraneTable: The ElemLink property is not of the right dimensions");
      return (false);        
    }
    
    // Get the SparseMatrix, if not we will not do this operation
    SparseRowMatrix *spr = dynamic_cast<SparseRowMatrix *>(ElemLink.get_rep());
    if (spr)
    {
      elemlinkrr = spr->get_rows();
      elemlinkcc = spr->get_cols();
      iselemlink = true;
    }
    else
    {
      // Inform the user that they tried something ill-fated
      pr->error("BuildMembraneTable: The ElemLink matrix is not a sparse matrix");
      return (false);       
    }
    
    for (size_t q=0; q<numelementdelems; q++)
    {
      if (elemlinkrr[q+1]-elemlinkrr[q] > 1)
      {
        // Inform the user that they tried something ill-fated
        pr->error("BuildMembraneTable: Each derivative element can only be connected to one other");
        return (false);             
      }
    }
  }  


  // Process NodeLink: this matrix tells which nodes are connected by a 
  // membrane to faces on the opposite side of the mesh

  // Assume we do not have it
  bool isnodelink = false;
  Matrix::index_type* nodelinkrr = 0;
  Matrix::index_type* nodelinkcc = 0;
  
	
  if (NodeLink.get_rep())
  {
    // We have a NodeLink Matrix
    
    // Do a sanity check, if not return a proper error
    if ((numelementnodes != NodeLink->nrows())&&(numelementnodes != NodeLink->ncols()))
    {
      pr->error("BuildMembraneTable: The NodeLink property is not of the right dimensions");
      return (false);        
    }
    
    // Get the SparseMatrix, if not we will not do this operation
    SparseRowMatrix *spr = dynamic_cast<SparseRowMatrix *>(NodeLink.get_rep());
    if (spr)
    {
      nodelinkrr = spr->rows;
      nodelinkcc = spr->columns;
      isnodelink = true;
    }
    else
    {
      // Inform the user that they tried something ill-fated
      pr->error("BuildMembraneTable: The NodeLink matrix is not a sparse matrix");
      return (false);       
    }
  }  

  if (NodeLink.get_rep() && !(ElemLink.get_rep()))
  {
    pr->error("BuildMembraneTable: This function needs both ElemLink and NodeLink, not only a NodeLink");
    return (false);         
  }

  if (!(NodeLink.get_rep()) && ElemLink.get_rep())
  {
    pr->error("BuildMembraneTable: This function needs both ElemLink and NodeLink, not only a ElemLink");
    return (false);         
  }

  // The GeomToComp: Matrix for conversion of node numbers from geometrical domain
  // to computational domain.

  bool isgeomtocomp = false;
  Matrix::index_type* geomrr = 0;
  Matrix::index_type* geomcc = 0;
  Matrix::size_type numgeomtocomp = 0;
  
  if (CompToGeom.get_rep())
  {
    if ((numelementnodes != CompToGeom->nrows()))
    {
      pr->error("BuildMembraneTable: The CompToGeom matrix property is not of the right dimensions");
      return (false);        
    }
    SparseRowMatrix *spr = dynamic_cast<SparseRowMatrix *>(CompToGeom.get_rep());
    if (spr)
    {
      geomrr = spr->rows;
      geomcc = spr->columns;
      
      isgeomtocomp = true;
    }
    else
    {
      // Inform the user that they tried something ill-fated
      pr->error("BuildMembraneTable: The CompToGeom matrix is not a sparse matrix");
      return (false);       
    }    
    
    // Check whether it is a one-on-one conversion. Each geometrical node has to
    // map to one computational node (the opposite does not need to be true)
    for (Matrix::index_type r=0; r<numelementnodes+1; r++)
    {
      if (geomrr[r] != r)
      {
        pr->error("BuildMembraneTable: The CompToGeom matrix property maps a geometric position to multiple computational nodes. This is not allowed");
        return (false);      
      }
    }
  }

  VMesh::Elem::iterator eit, eit_end;
  VMesh::Node::array_type nodes;  

  membranemesh->begin(eit);
  membranemesh->end(eit_end);  

  membranemesh->get_nodes(nodes,*eit); 
  
  Matrix::size_type nodespersurf =nodes.size();
  std::vector<Point> points(nodespersurf);
  std::vector<double> surfaces(nodespersurf);

  bool foundsurf1;
  bool foundsurf2;

  std::pair<nodemap_type::iterator,nodemap_type::iterator> lit;
  VMesh::Node::index_type idx;
  VMesh::Elem::array_type elems;
  VMesh::Elem::array_type elems2;
  VMesh::DElem::array_type delems;
  VMesh::DElem::array_type delems2;
  VMesh::DElem::array_type delemlist;
  VMesh::Elem::array_type elemlist;
  VMesh::Node::array_type enodes;
  VMesh::Node::array_type enodes2;
  double modeltype;
  double classvalue1,classvalue2;
  double   surface;
  
  Matrix::index_type k = 0;
  VMesh::Elem::size_type numelems;
  membranemesh->size(numelems);
  // Create a table that lists how two nodes in the split model are connected
  // and to which membrane node they correspond and which surface connects them
  // The table is build of node per surface components and later the individual
  // parts are combined. 
  MembraneTable membranetablelisttemp; 
  membranetablelisttemp.resize(nodespersurf*numelems);
	
  // Informing where we are
  std::ostringstream oss;
  oss << nodespersurf*numelems;
  pr->remark("BuildMembraneTable: Building a table with "+oss.str()+" entries");
  
  // Loop over all elements in the membrane model
  while (eit != eit_end)
  {
    // We want to find two surfaces in the elementtype model for each element in
    // membrane model. This function is similar to locate(), but will find multiple
    // entries. Since membranes are infinitely thin, meshes touch eachother. And
    // unfortunately locate() was not designed for that purpose. Thus here is an
    // algorithm that avoids the problems with locate()
    
    // Assume we have not found the surfaces yet
    foundsurf1 = false;
    foundsurf2 = false;
    
    // Get the node numbers of the membrane element
    membranemesh->get_nodes(nodes,*eit); 
    
    surface = (membranemesh->get_size(*eit))/(nodes.size());
    
    // As we need to compare physical locations, start with getting the node locations
    for (size_t p =0; p < nodes.size(); p++) 
    {
      membranemesh->get_center(point2,nodes[p]);
      points[p] = point2;
      surfaces[p] = surface; // <<<< NEED TO ALTER THIS
    }

    // points contains locations of membrane patch
    // surface contains surface weight factorss

    // Calculate the key for the first point to find the corresponding node
    // in the hash_map
    key = static_cast<int>((points[0].x()-xmin)*xmul);
    key += (static_cast<int>((points[0].y()-ymin)*ymul))<<8;
    key += (static_cast<int>((points[0].z()-zmin)*zmul))<<16;
    
    // Find all nodes that can be used in the elementtype mesh
    // This to make sure we do not need to do an exhaustive search through the
    // whole mesh.
    lit = nodemap.equal_range(key);

    // loop through all possibilities
    while ((lit.first != lit.second)&&(foundsurf2 == false))
    {
      idx = (*(lit.first)).second;
      elementtypemesh->get_center(point2,idx);

      // check 1: target point needs to be same as point we are looking for
      // SearchGrid is of finite dimensions, thus two different points can have same key
      if (point2 == points[0])
      {
        // get all elements it is connected to:
        elementtypemesh->get_elems(elems,idx);
        delemlist.clear();
        elemlist.clear();
        
        // Get all unique faces that connect 
        for (size_t p = 0; p < elems.size(); p++)
        {
          // for each connected element get the faces
          elementtypemesh->get_delems(delems,elems[p]);
          for (size_t r = 0; r < delems.size();  r++)
          {
             // Check if it was not already in the list
             size_t s;
             for (s=0; s<delemlist.size(); s++) if (delems[r] == delemlist[s]) break;
             if (s < delemlist.size()) continue;
             elementtypemesh->get_nodes(enodes,delems[r]);
             // Only add it if it contains the target node
             size_t u;
             for (u=0;u<enodes.size();u++) if (enodes[u] == idx) break;
             if (u < enodes.size()) { delemlist.push_back(delems[r]); elemlist.push_back(elems[p]); }
          }
        }

        // check whether  the faces correspond to the face we are looking for
        bool isequal = false;
        for (size_t r = 0; (r < delemlist.size())&&(!isequal);  r++)
        {
          // Get the nodes of the face for testing whether it contains the other locations
          // we are looking for
          elementtypemesh->get_nodes(enodes,delemlist[r]);
      
          isequal = true;
          for (size_t q = 0; q < enodes.size(); q++)
          {
            elementtypemesh->get_center(point2,enodes[q]);
            size_t t=0;
            for (; t< points.size(); t++) { if (point2 == points[t]) break; } 
            // if we did not find the node in this face it was definitely not
            // the one we were looking for.
            if (t == points.size()) 
            {
              isequal = false;
            }
          }

          // in case we found a surface
          if (isequal)
          {
            
            // is it the first one?
            if (foundsurf1 == false)
            {
              // Get the element type do decide the orientation of the membrane model
              elementtypefield->value(classvalue1,elemlist[r]);

              // Add the newly found surface to the list
              for (size_t q = 0; q < enodes.size(); q++)
              {           
                membranetablelisttemp[k+q].node1 = enodes[q];
              }
              foundsurf1 = true;

              // check whether we can figure out whether the link matrices can
              // provide us with the other surface
              if (iselemlink&&isnodelink)
              {    
                // if so elemlink so point to another face
                // Check whether there is an entry in elemlink
                if (elemlinkrr[delemlist[r]] < elemlinkrr[delemlist[r]+1])
                {
                  // get the index of this other face
                  VMesh::DElem::index_type fidx = elemlinkcc[elemlinkrr[delemlist[r]]];
                  
                  // Now the ugly part, we need to access the element type of that face and there is no
                  // function from face to element.
                  // So we need to use the nodes 
                  
                  elementtypemesh->get_nodes(enodes2,fidx);
                  elementtypemesh->get_elems(elems2,enodes2[0]);
                  
                  bool foundelem = false;
                  for (size_t q = 0; (q < elems2.size())&&(!foundelem); q++)
                  {
                    elementtypemesh->get_delems(delems2,elems2[q]);
                    for (size_t u = 0; u < delems2.size(); u++)
                    {
                      if (delems2[u] == fidx) { elementtypefield->value(classvalue2,elems2[q]); foundelem = true; break; }
                    }
                  }
                  
                  if (foundelem = false)
                  {
                    pr->error("BuildMembraneTable: Could not find elem connected to face, internal error in mesh");
                    return (false);
                  }

                  // if they are in the same domain we should need to link them  
                  //if (classvalue1 != classvalue2)
                  //{  
                    // Now find the how the nodes connect over the boundary,
                    // loop through the nodelink matrix and figure out which link
                    // finds the nodes of the faces. This way we can establish the
                    // order as faces may show a different ordering of the nodes
                    for (size_t q = 0; q < enodes.size(); q++)
                    {           
                      bool foundnode = false;
                      for (size_t u=nodelinkrr[enodes[q]]; (u < nodelinkrr[enodes[q]+1])&&(!foundnode); u++)
                      {
                        for (size_t t=0; t < enodes2.size(); t++)
                        {
                          if (enodes2[t] == nodelinkcc[u]) 
                          { 
                            // add the node numbers in the proper order
                            membranetablelisttemp[k+q].node2 = enodes2[t];
                            foundnode = true; 
                            break;
                          }
                        }
                      } 
                    }
                    
                    // Now find how it connects to the model surface:
                    // Here we rely on physical point locations as that is the
                    // only way, we did not keep track of node numbers when splitting
                    // out the surfaces. As membranes are subdivided a couple of 
                    // times, keeping track is a GUI nightmare, hence we reconstruct
                    // it here.
                    										
                    for (size_t q = 0; q < enodes.size(); q++)
                    {  
                      elementtypemesh->get_center(point2,static_cast<VMesh::Node::index_type>(membranetablelisttemp[k+q].node1));
                      size_t s=0;
                      for (; s < points.size(); s++)
                      {
                        // TO DO: at tolerance here
                        if (points[s] == point2) break;
                      }
                      if (s == points.size())
                      {
                        pr->error("BuildMembraneTable: Membrane face and ElementType face do not line up");
                        return (false);
                      }
                      
                      // Surfaces for the connection are also defined based on the 
                      // membrane surface, hence we use this one here as well.
                      membranetablelisttemp[k+q].node0 = nodes[s];
                      membranetablelisttemp[k+q].surface = surfaces[s];                    
                     
                      // In case the membrane should be defined the other way around:
                      // Order matters: for instance a lot of models have an ICS and
                      // ECS side. Similarly gap junctions have different sides
                      // depending on cell type: e.g. fibroblast and myocyte                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                      if (classvalue1 > classvalue2)
                      {
                        unsigned int temp = membranetablelisttemp[k+q].node2;
                        membranetablelisttemp[k+q].node2 = membranetablelisttemp[k+q].node1;
                        membranetablelisttemp[k+q].node1 = temp;
                      }
                      else if (classvalue1 == classvalue2)
                      {                        
                        if (membranetablelisttemp[k+q].node2 < membranetablelisttemp[k+q].node1)
                        {
                          unsigned int temp;
                          temp = membranetablelisttemp[k+q].node1;
                          membranetablelisttemp[k+q].node1 = membranetablelisttemp[k+q].node2;
                          membranetablelisttemp[k+q].node2 = temp;                        
                        }
                      }
                    }
                    
                    // finish by saying that we do not need to process second surface
                    foundsurf2 = true; 
                    k += enodes.size();
                  //}
                }
              }
            }
            else if (foundsurf2 == false)
            {
              // we already have one surface but the other one is still missing
              // Get the element value of this surface
              elementtypefield->value(classvalue2,elemlist[r]);
              
              // find out how surfaces are connected node wise
              for (size_t q = 0; q < enodes.size(); q++)
              { 
                // get the physical location to match it with the other surface
                // and the model surface
                
                // get physical location of first face nodes
                elementtypemesh->get_center(point2,static_cast<VMesh::Node::index_type>(membranetablelisttemp[k+q].node1));
                
                // compare with currently found face
                size_t t=0;
                for (; t < enodes.size(); t++)
                {
                  elementtypemesh->get_center(point3,enodes[t]);
                  if (point3 == point2) break;
                }
                // t is now relative index on how different the nodes of the face line up

                // compare with membrane model itself
                size_t s=0;
                for (; s < points.size(); s++)
                {
                  if (points[s] == point2) break;
                }
                // s is now relative index
                            
                // depending on order we swap node2 and node1                                                                                                                                                                                  
                if (classvalue1 < classvalue2)
                {
                  membranetablelisttemp[k+q].node2 = enodes[t];
                }
                else if (classvalue1 == classvalue2)
                {
                  membranetablelisttemp[k+q].node2 = enodes[t];
                  if (membranetablelisttemp[k+q].node2 < membranetablelisttemp[k+q].node1)
                  {
                    unsigned int temp;
                    temp = membranetablelisttemp[k+q].node1;
                    membranetablelisttemp[k+q].node1 = membranetablelisttemp[k+q].node2;
                    membranetablelisttemp[k+q].node2 = temp;
                  }
                }
                else
                {
                  membranetablelisttemp[k+q].node2 = membranetablelisttemp[k+q].node1;
                  membranetablelisttemp[k+q].node1 = enodes[t];
                }
                membranetablelisttemp[k+q].node0 = nodes[s];
                membranetablelisttemp[k+q].surface = surfaces[s];
              }
              
              foundsurf2 = true;
              k += enodes.size();
            }
          }
        }
      }
      ++(lit.first);
    }
      
    if ((foundsurf1 == false)||(foundsurf2 == false))
    {
      pr->error("BuildMembraneTable: Not every surface/curve in the membrane model can be found in the element type mesh");
      return (false); 
   }
    ++eit;
  }
  
	
  // we have the complete list
  // now merge it together:
  // we have double entries as we searched per face, adjoining faces will share
  // connections. Remove these to speed up computations in simulation phase


  Matrix::index_type *mrr = new Matrix::index_type[nummembranenodes+1];
  Matrix::index_type *mcc = new Matrix::index_type[2*nummembranenodes];
  double* mvv = new double[2*nummembranenodes];

  Matrix::index_type *mrr2 = new Matrix::index_type[nummembranenodes+1];
  Matrix::index_type *mcc2 = new Matrix::index_type[nummembranenodes];
  double* mvv2 = new double[nummembranenodes];


  if ((mrr == 0)||(mcc == 0)||(mvv == 0)||(mrr2 == 0)||(mcc2 == 0)||(mvv2 == 0))
  {
    delete[] mrr;
    delete[] mcc;
    delete[] mvv;
    delete[] mrr2;
    delete[] mcc2;
    delete[] mvv2;
    
    pr->error("BuildMembraneTable: Could not generate mapping matrix");
    return (false);
  }

  for (index_type r=0; r < nummembranenodes; r++) mrr[r] = r;
  
  if (isgeomtocomp)
  {
    for (size_t q=0; q < membranetablelisttemp.size(); q++)
    {                  
      if (membranetablelisttemp[q].node1 == membranetablelisttemp[q].node2)
      {
        std::cerr << "node is itself\n";
      }
      membranetablelisttemp[q].node1 = geomcc[membranetablelisttemp[q].node1];
      membranetablelisttemp[q].node2 = geomcc[membranetablelisttemp[q].node2];          
    }
  }
  
  std::sort(membranetablelisttemp.begin(),membranetablelisttemp.end());
  
  membranetable.clear();

  // Informing where we are
  std::ostringstream oss2;
  oss2 << membranetablelisttemp.size();
  pr->remark("BuildMembraneTable: Membrane table has "+oss2.str()+" entries");
  
  // Find unique entries and count them and move surface areas to add surfaces
  // of duplicate entries
  if (membranetablelisttemp.size() > 0)
  {
    size_t tablesize = 1;
    size_t q = 0;
    for (size_t p=1; p < membranetablelisttemp.size(); p++)
    {
      if (membranetablelisttemp[p] == membranetablelisttemp[q])
      {
        membranetablelisttemp[q].surface += membranetablelisttemp[p].surface;
        membranetablelisttemp[p].surface = 0.0;
        if (membranetablelisttemp[q].node0 > membranetablelisttemp[p].node0)
        {
          mrr[membranetablelisttemp[q].node0] = membranetablelisttemp[p].node0;
          membranetablelisttemp[q].node0 = membranetablelisttemp[p].node0;
        }
        else
        {
          mrr[membranetablelisttemp[p].node0] = membranetablelisttemp[q].node0;
        }
      }
      else if (membranetablelisttemp[p].surface != 0.0)
      {
        q = p;
        tablesize++;
      }
      else
      {
        std:cerr << "Detected synapse with surface area equal to zero\n";
      }
    } 

    for (index_type r=0; r<nummembranenodes; r++)
    {
      index_type p = r;
      while (mrr[p] != p) p = mrr[p];
      mrr[r] = p;      
    }

    // Build the final list with only unique entries
    q = 0;
    membranetable.resize(tablesize);
    for (size_t p=0; p < membranetablelisttemp.size(); p++)
    {
      if (membranetablelisttemp[p].surface != 0.0)
      {
        membranetable[q] = membranetablelisttemp[p];
        membranetable[q].node0 = mrr[membranetable[q].node0];
        membranetable[q].snode = offset+q; 
        q++;
      }    
    }
    offset += q;
  }
  else
  {
    pr->error("BuildMembraneTable: The Membrane geometry does not correspond to any of the element faces/edges of the element type field");
    return (true);    
  }

  for (int r=0; r < nummembranenodes; r++)
  {
    mcc[2*r] = 0;
    mcc[2*r+1] = 0;
    mvv[2*r] = 0.0;
    mvv[2*r+1] = 0.0;
    mcc2[r] = 0;
    mvv2[r] = 0.0;
	}
	
  size_t nummembraneentries = membranetable.size();
  for (index_type r=0; r < nummembraneentries; r++)
  {
    index_type q = membranetable[r].node0; 

		if (membranetable[r].node2>membranetable[r].node1)
		{
			mcc[2*q] = membranetable[r].node1;
			mcc[2*q+1] = membranetable[r].node2;
			mvv[2*q] = -1.0;
			mvv[2*q+1] = 1.0;		
		}
		else
		{
			mcc[2*q] = membranetable[r].node2;
			mcc[2*q+1] = membranetable[r].node1;
			mvv[2*q] = 1.0;
			mvv[2*q+1] = -1.0;
		}
    
    mcc2[q] = membranetable[r].snode;
    mvv2[q] = 1.0;
	}

  for (index_type r=0; r < nummembranenodes; r++)
  {
    if (mrr[r] < r)
    {
      index_type q = 2*mrr[r];
      mcc[2*r] =   mcc[q];
      mcc[2*r+1] = mcc[q+1];
      mvv[2*r] =   mvv[q];
      mvv[2*r+1] = mvv[q+1];
      q = mrr[r];
      mcc2[r] = mcc2[q];
      mvv2[r] = mvv2[q];
    }
    mrr[r] = 2*r;
    mrr2[r] = r;
  }
 
  for (index_type r=0; r < nummembranenodes; r++)
  {
    if (mvv[2*r] == 0.0) std::cerr << "Membrane node " << r << "does not have any potentials assigned to it\n";
  }
 
  mrr[nummembranenodes] = 2*nummembranenodes; // An extra entry goes on the end of rr.
  mrr2[nummembranenodes] = nummembranenodes; // An extra entry goes on the end of rr.
  

  MappingMatrix = new SparseRowMatrix(nummembranenodes, offset, mrr, mcc, 2*nummembranenodes, mvv);
  if (MappingMatrix.get_rep() == 0)
  {
    pr->error("LinkToCompGrid: Could build geometry to computational mesh mapping matrix");
    return (false);
  }

  MappingMatrix2 = new SparseRowMatrix(nummembranenodes, offset, mrr2, mcc2, nummembranenodes, mvv2);
  if (MappingMatrix2.get_rep() == 0)
  {
    pr->error("LinkToCompGrid: Could build geometry to computational mesh mapping matrix");
    return (false);
  }
  
  // Success:
  return (true);
}


} // End namespace ModelCreation
