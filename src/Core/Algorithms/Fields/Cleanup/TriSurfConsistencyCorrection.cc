/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2013 Scientific Computing and Imaging Institute,
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

#include <algorithm>
#include <iostream>
#include <sstream>
#include <set>
#include <sstream>
#include <queue>
#include <vector>

#include <Core/Datatypes/VField.h>

#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/VMesh.h>

#include <Core/Datatypes/ColumnMatrix.h>

#include <Core/Algorithms/Fields/Cleanup/TriSurfConsistencyCorrection.h>
#include <Core/Datatypes/FieldInformation.h>

namespace SCIRunAlgo {

using namespace SCIRun;

bool 
TriSurfConsistencyCorrectionAlgo::run(FieldHandle inputField,
                                  FieldHandle& outputField, MatrixHandle& invertedElementsListMatrix)
{
  algo_start("TriSurfConsistencyCorrection");

  // Handle: the function get_rep() returns the pointer contained in the handle
  if (inputField.get_rep() == 0)
  {
    // If we encounter a null pointer we return an error message and return to
    // the program to deal with this error. 
    error("No input field");
    algo_end(); return (false);
  }
  
  // Step 1: determine the type of the input fields and determine what type the
  // output field should be.
  
  FieldInformation fi(inputField);
  FieldInformation fo(inputField);

  // Here we test whether the class is part of any of these newly defined 
  // non-linear classes. If so we return an error.
  if (fi.is_nonlinear())
  {
    error("This algorithm has not yet been defined for non-linear elements yet");
    algo_end(); return (false);
  }
  
  // This one
  if (! fi.is_trisurfmesh() ) 
  {
    // Notify the user that no action is done  
    error("This algorithm only works on a TriSurfMesh");
    algo_end(); return (false);
  }
  
  VField* inputVField = inputField->vfield();
  //VMesh* inputVMesh = inputField->vmesh();
  
  fo.make_trisurfmesh();
  outputField = CreateField( fo, inputField->mesh() );
  outputField.detach();
  outputField->mesh_detach();

  if (outputField.get_rep() == 0)
  {
    error("Could not allocate output field");
    algo_end(); return (false);
  }

  if (! fo.is_trisurfmesh() )
  {
    // Notify the user that no action is done  
    error("Output is not a TriSurfMesh");
    algo_end(); return (false);
  }

  VMesh *outputVMesh = outputField->vmesh();  
  const VMesh::size_type n = outputVMesh->num_nodes();
  const VMesh::size_type m = outputVMesh->num_elems();
  
  size_type i, j, k;
  std::vector<int> invertedElements;
  
  std::vector< std::set<int> >elemsOfVert;
  std::vector< std::set<int> >elemNeighbors;
  std::vector< std::vector<int> > elem;
  std::set<int> dummy;
  std::vector<int> tempVec;

  for (i = 0; i < 3; i++)
    tempVec.push_back(0);
  
  // initialize vector of vectors
  for (i = 0; i < m; i++)
    elem.push_back(tempVec);
  
  VMesh::Face::iterator meshFaceIter;
  VMesh::Face::iterator meshFaceEnd;
  
  VMesh::Node::array_type nodesFromFace(3);
  
  outputVMesh->end(meshFaceEnd);
  
  for (outputVMesh->begin(meshFaceIter); meshFaceIter != meshFaceEnd; ++meshFaceIter)
  {
    // get nodes from mesh element
    VMesh::Face::index_type elemID = *meshFaceIter;
    outputVMesh->get_nodes(nodesFromFace, elemID);
 
    for (i = 0; i < 3; i++)
      elem[elemID][i] = nodesFromFace[i];
  }
  
  //m=elem.size(); n=nodes.size();
  for (i = 0; i < m; i++)
    elemNeighbors.push_back(dummy);

  for (i = 0; i < n; i++)
    elemsOfVert.push_back(dummy);
  
  // get elements
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < 3; j++)
      elemsOfVert[elem[i][j]].insert(i);
  }

  // get edges
  // create graph
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < 3; j++)
    {
      k = (j+1)%3;
      std::set<int>::iterator it;
      for (it = elemsOfVert[elem[i][j]].begin(); it != elemsOfVert[elem[i][j]].end(); it++)
      {
        if (elemsOfVert[elem[i][k]].find(*it) != elemsOfVert[elem[i][k]].end())
        {
          if (*it != i)
          {
            elemNeighbors[i].insert(*it);
            elemNeighbors[i].insert(*it);
          }
        }
      }
    }
  }

  // traverse graph
  std::vector<bool> elemTraversed;
  for (i = 0; i < m; i++)
    elemTraversed.push_back(false);

  std::queue<int> bfs;
  bfs.push(0);
  elemTraversed[0] = true;
  std::set<std::pair<int,int> > edges;
  
  while (! bfs.empty() )
  {
    //insert elems not traversed
    std::set<int>::iterator it;
    i = bfs.front();
    for (it = elemNeighbors[i].begin(); it != elemNeighbors[i].end(); it++)
    {
      if (! elemTraversed[*it])
      {
        bfs.push(*it);
        elemTraversed[*it] = true;
      }
    }
    // reordered elem i if desired
    bool flag = false;
    for (j = 0; j < 3 && flag == false; j++)
    {
      k = (j+1)%3;
      if (edges.find(std::make_pair(elem[i][j], elem[i][k])) != edges.end())
      {
        invertedElements.push_back(i);
        // invert ordering
        flag = true;
        int temp = elem[i][j];
        elem[i][j]=elem[i][k];
        elem[i][k]=temp;
        // invert in the output field
        VMesh::Face::index_type elemID = i;
        outputVMesh->get_nodes(nodesFromFace, elemID);
        nodesFromFace[j] = elem[i][j];
        nodesFromFace[k] = elem[i][k];
      }
    }
    for (j = 0; j < 3; j++)
    {
      k=(j+1)%3;
      edges.insert(std::make_pair(elem[i][j], elem[i][k]));
    }
    
    bfs.pop();
  }

  VField* outputVField = outputField->vfield();
  outputVField->resize_values();
  std::vector<VMesh::index_type> order;

  //! Copy field data (non-linear not supported, check made upstream)
  int basis_order = inputVField->basis_order();
  if (basis_order == 0)
  {
    VField::size_type size = order.size();
    for(VField::index_type idx = 0; idx < size; idx++)
      outputVField->copy_value(inputVField, order[idx], idx);
  }
  else if (basis_order == 1)
  {
    outputVField->copy_values(inputVField);
  }
  
  //! Copy properties of the property manager
  outputField->copy_properties(inputField.get_rep());

  if (invertedElements.size() > 0)
  {
    std::ostringstream oss;
    oss << invertedElements.size() << " elements were found to be incorrectly oriented in input TriSurf mesh.";
    warning(oss.str());
    
    bool output_inverted_element_list = get_bool("output_inverted_element_list");
    if (output_inverted_element_list)
    {
      invertedElementsListMatrix = new ColumnMatrix(invertedElements.size());    
      if (invertedElementsListMatrix.get_rep() == 0)
      {
        error("Could not allocate inverted element matrix");
        algo_end(); return (false);  
      }
      std::copy(invertedElements.begin(), invertedElements.end(), invertedElementsListMatrix->begin());
    }
  }
  else
  {
    remark("Input TriSurf mesh is consistent.");
  }

  // Success:
  algo_end(); return (true);
}

} // End namespace SCIRunAlgo
