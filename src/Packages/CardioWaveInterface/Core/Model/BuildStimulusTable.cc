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

#include <Packages/CardioWaveInterface/Core/Model/BuildStimulusTable.h>

#include <sci_hash_map.h>
#include <algorithm>

namespace CardioWaveInterface {

using namespace SCIRun;

bool 
BuildStimulusTableAlgo::
BuildStimulusTable(ProgressReporter *pr,  FieldHandle ElementType, 
                   FieldHandle Stimulus, MatrixHandle CompToGeom, 
                   double domaintypemin, double domaintypemax, 
                   bool selectbynode, StimulusTable& stimulustable)
{
  if (ElementType.get_rep() == 0)
  {
    pr->error("BuildStimulusTable: No domain nodetype field");
    return (false);
  }

  if (Stimulus.get_rep() == 0)
  {
    pr->error("BuildStimulusTable: No Stimulus model field");
    return (false);
  }

  // no precompiled version available, so compile one

  FieldInformation fi(ElementType);
  FieldInformation fi2(Stimulus);
  
  if (!(fi.is_constantdata()))
  {
    pr->error("BuildStimulusTable: The ElementType field needs to have one data value assigned to each element");
    return (false);
  }
  
  if (!fi.is_unstructuredmesh())
  {
    pr->error("BuildStimulusTable: This function is not defined for structured meshes");
    return (false);
  }  

  if (!(fi.is_volume()||fi.is_surface()))
  {
    pr->error("BuildStimulusTable: The domain nodetype field needs to be a volume or surface");
    return (false);
  }  

  if (fi2.is_pointcloud()) 
  {
    VField *elementtypefield = ElementType->vfield();
    VField *stimfield = Stimulus->vfield();

    if (stimfield == 0)
    {
      pr->error("BuildStimulusTable: The Stimulus field is empty");
      return (false);
    }

    if (elementtypefield == 0)
    {
      pr->error("BuildStimulusTable: The domaintype field is empty");
      return (false);
    }

    VMesh* elementtypemesh = elementtypefield->vmesh();
    if (elementtypemesh == 0)
    {
      pr->error("BuildStimulusTable: The domaintype mesh is empty");
      return (false);
    }

    VMesh* stimmesh = stimfield->vmesh();
    if (stimmesh == 0)
    {
      pr->error("BuildStimulusTable: The stimmesh is empty");
      return (false);
    }

    VMesh::Node::size_type nnodes;
    bool isgeomtocomp = false;
    Matrix::index_type  *geomtocomprr = 0;
    Matrix::index_type  *geomtocompcc = 0;

    if (CompToGeom.get_rep())
    {
      elementtypemesh->size(nnodes);
      if (nnodes != CompToGeom->nrows())
      {
        pr->error("BuildStimulusTable: The number of rows in CompToGeom is not of the right size");
        return (false);    
      }
      
      SparseRowMatrix* mat = dynamic_cast<SparseRowMatrix *>(CompToGeom.get_rep());
      if (mat == 0)
      {
        pr->error("BuildStimulusTable: CompToGeom is not a sparse matrix");
        return (false);    
      }
    
      geomtocomprr = mat->rows;
      geomtocompcc = mat->columns;
      
      for (int r=0; r<nnodes+1;r++)
      {
        if (geomtocomprr[r] != r)
        {
          pr->error("BuildStimulusTable: CompToGeom is mapping a geometrical node to multiple computational nodes");
          return (false);     
        }
      }
      isgeomtocomp = true;
    }
    VMesh::Elem::array_type elems;
    VMesh::Node::iterator it, it_end;
    VMesh::Node::size_type sz;
    VMesh::Node::index_type ci;
    double val, dval;
    Point point;

    elementtypemesh->size(sz);
    stimmesh->begin(it);
    stimmesh->end(it_end);  

    std::vector<bool> indomain(sz,false);

    double dvalmin = domaintypemin;
    double dvalmax = domaintypemax;
    
    elementtypemesh->synchronize(Mesh::LOCATE_E|Mesh::EDGES_E);

    while (it != it_end)
    {    
      stimmesh->get_center(point,*it);
      if (elementtypemesh->locate(ci,point))
      {
        elementtypemesh->get_elems(elems,ci);
        if (elems.size() > 0)
        {
          elementtypefield->get_value(val,elems[0]);
          if (val >= dvalmin && val <= dvalmax)
          {
            stimulusparam_type stimitem;
            stimitem.node = static_cast<unsigned int>(ci);
            stimitem.weight = 1.0;     
            stimulustable.push_back(stimitem);
          }
        }
      }
      ++it;
    }

    if (isgeomtocomp)
    {
      for (size_t p=0; p < stimulustable.size(); p++) 
      {
        stimulustable[p].node = geomtocompcc[stimulustable[p].node];
      }
    }

    std::sort(stimulustable.begin(),stimulustable.end());

    size_t k = 0;
    for (size_t p=1; p < stimulustable.size(); p++) 
    {
      if (stimulustable[p].node == stimulustable[k].node)
      {
        stimulustable[k].weight += stimulustable[p].weight;
        stimulustable[p].weight = 0.0;
      }
      else
      {
        k++;
        stimulustable[k] = stimulustable[p];
      }
    }
       
    if (stimulustable.size() > 0)       
    {
      stimulustable.resize(k+1);
    }
                       
    // Success:
    return (true);
  }
  else
  {
    if (!selectbynode)
    {
      if (fi2.is_curve())
      {
        VField *elementtypefield = ElementType->vfield();
        VField *stimfield = Stimulus->vfield();

        if (stimfield == 0)
        {
          pr->error("BuildStimulusTable: The Stimulus field is empty");
          return (false);
        }

        if (elementtypefield == 0)
        {
          pr->error("BuildStimulusTable: The domaintype field is empty");
          return (false);
        }

        VMesh* elementtypemesh = elementtypefield->vmesh();
        if (elementtypemesh == 0)
        {
          pr->error("BuildStimulusTable: The domaintype mesh is empty");
          return (false);
        }

        VMesh* stimmesh = stimfield->vmesh();
        if (stimmesh == 0)
        {
          pr->error("BuildStimulusTable: The stimmesh is empty");
          return (false);
        }

        VMesh::Node::size_type nnodes;
        bool isgeomtocomp = false;
        Matrix::index_type  *geomtocomprr = 0;
        Matrix::index_type  *geomtocompcc = 0;

        if (CompToGeom.get_rep())
        {
          elementtypemesh->size(nnodes);
          if (nnodes != CompToGeom->nrows())
          {
            pr->error("BuildStimulusTable: The number of rows in CompToGeom is not of the right size");
            return (false);    
          }
          
          SparseRowMatrix* mat = dynamic_cast<SparseRowMatrix *>(CompToGeom.get_rep());
          if (mat == 0)
          {
            pr->error("BuildStimulusTable: CompToGeom is not a sparse matrix");
            return (false);    
          }
        
          geomtocomprr = mat->rows;
          geomtocompcc = mat->columns;
          
          for (int r=0; r<nnodes+1;r++)
          {
            if (geomtocomprr[r] != r)
            {
              pr->error("BuildStimulusTable: CompToGeom is mapping a geometrical node to multiple computational nodes");
              return (false);     
            }
          }
          isgeomtocomp = true;
        }

        VMesh::Elem::array_type elems;
        VMesh::Node::iterator it, it_end;
        VMesh::Node::size_type sz;
        VMesh::Elem::index_type ci;
        double val, dval;
        Point point;

        elementtypemesh->size(sz);
        elementtypemesh->begin(it);
        elementtypemesh->end(it_end);  

        std::vector<bool> indomain(sz);

        double dvalmin = domaintypemin;
        double dvalmax = domaintypemax;
        
        stimmesh->synchronize(Mesh::LOCATE_E|Mesh::EDGES_E);

        while (it != it_end)
        {
          indomain[(static_cast<unsigned int>(*it))] = false;
          
          elementtypemesh->get_center(point,*it);
          if (stimmesh->locate(ci,point))
          {
            elementtypemesh->get_elems(elems,*it);
            if (elems.size() > 0)
            {
              elementtypefield->get_value(val,elems[0]);
              if (val >= dvalmin && val <= dvalmax)
              {
                indomain[static_cast<unsigned int>(*it)] = true;        
              }
            }
          }
          ++it;
        }

        VMesh::Edge::iterator cit, cit_end;
        VMesh::Node::array_type nodes;  
          
        elementtypemesh->begin(cit);
        elementtypemesh->end(cit_end);  
       
        stimulustable.clear();
       
        while (cit != cit_end)
        {
          elementtypemesh->get_nodes(nodes,*cit);
          size_t p = 0;
          for (; p < nodes.size(); p++)
          {
            if (indomain[static_cast<unsigned int>(nodes[p])] == false) break;
          }
          if (p == nodes.size())
          {
            elementtypemesh->get_center(point,*cit);
            if (stimmesh->locate(ci,point))
            {    
              for (p = 0; p < nodes.size(); p++)
              {
                stimulusparam_type stimitem;
                stimitem.node = static_cast<unsigned int>(nodes[p]);
                stimitem.weight = elementtypemesh->get_size(*cit)/elementtypemesh->num_nodes_per_edge();
                stimulustable.push_back(stimitem);
              }
            }
          }
          ++cit;
        }

        if (isgeomtocomp)
        {
          for (size_t p=1; p < stimulustable.size(); p++) 
          {
            stimulustable[p].node = geomtocompcc[stimulustable[p].node];
          }
        }

        std::sort(stimulustable.begin(),stimulustable.end());

        size_t k = 0;
        for (size_t p=0; p < stimulustable.size(); p++) 
        {
          if (stimulustable[p].node == stimulustable[k].node)
          {
            stimulustable[k].weight += stimulustable[p].weight;
            stimulustable[p].weight = 0.0;
          }
          else
          {
            k++;
            stimulustable[k] = stimulustable[p];
          }
        }
           
        if (stimulustable.size() > 0)       
        {
          stimulustable.resize(k+1);
        }

                       
        // Success:
        return (true);
      
      }
      else if (fi2.is_surface())
      {

        VField *elementtypefield = ElementType->vfield();
        VField *stimfield = Stimulus->vfield();

        if (stimfield == 0)
        {
          pr->error("BuildStimulusTable: The Stimulus field is empty");
          return (false);
        }

        if (elementtypefield == 0)
        {
          pr->error("BuildStimulusTable: The domaintype field is empty");
          return (false);
        }

        VMesh* elementtypemesh = ElementType->vmesh();
        if (elementtypemesh == 0)
        {
          pr->error("BuildStimulusTable: The domaintype mesh is empty");
          return (false);
        }

        VMesh* stimmesh = Stimulus->vmesh();
        if (stimmesh == 0)
        {
          pr->error("BuildStimulusTable: The stimmesh is empty");
          return (false);
        }


        VMesh::Node::size_type nnodes;
        bool isgeomtocomp = false;
        Matrix::index_type  *geomtocomprr = 0;
        Matrix::index_type  *geomtocompcc = 0;

        if (CompToGeom.get_rep())
        {
          elementtypemesh->size(nnodes);
          if (nnodes != CompToGeom->nrows())
          {
            pr->error("BuildStimulusTable: The number of rows in CompToGeom is not of the right size");
            return (false);    
          }
          
          SparseRowMatrix* mat = dynamic_cast<SparseRowMatrix *>(CompToGeom.get_rep());
          if (mat == 0)
          {
            pr->error("BuildStimulusTable: CompToGeom is not a sparse matrix");
            return (false);    
          }
        
          geomtocomprr = mat->rows;
          geomtocompcc = mat->columns;
          
          for (int r=0; r<nnodes+1;r++)
          {
            if (geomtocomprr[r] != r)
            {
              pr->error("BuildStimulusTable: CompToGeom is mapping a geometrical node to multiple computational nodes");
              return (false);     
            }
          }
          isgeomtocomp = true;
        }

        VMesh::Node::iterator it, it_end;
        VMesh::Node::size_type sz;
        VMesh::Elem::array_type elems;
        double val, dval;
        VMesh::Elem::index_type ci;
        Point point;

        elementtypemesh->size(sz);
        elementtypemesh->begin(it);
        elementtypemesh->end(it_end);  

        std::vector<bool> indomain(sz);

        double dvalmin = domaintypemin;
        double dvalmax = domaintypemax;

        stimmesh->synchronize(Mesh::LOCATE_E|Mesh::FACES_E);

        while (it != it_end)
        {
          indomain[(static_cast<unsigned int>(*it))] = false;

          elementtypemesh->get_center(point,*it);
          if (stimmesh->locate(ci,point))
          {
            elementtypemesh->get_elems(elems,*it);
            if (elems.size() > 0)
            {
              elementtypefield->get_value(val,elems[0]);
              if (val >= dvalmin && val <= dvalmax)
              {
                indomain[(static_cast<unsigned int>(*it))] = true;        
              }
            }
          }
          ++it;
        }

        VMesh::Face::iterator cit, cit_end;
        VMesh::Node::array_type nodes;  
          
        elementtypemesh->begin(cit);
        elementtypemesh->end(cit_end);  
       
        stimulustable.clear();
       
        while (cit != cit_end)
        {
          elementtypemesh->get_nodes(nodes,*cit);
          size_t p = 0;
          for (; p < nodes.size(); p++)
          {
            if (indomain[static_cast<unsigned int>(nodes[p])] == false) break;
          }
          if (p == nodes.size())
          {
            elementtypemesh->get_center(point,*cit);
            if (stimmesh->locate(ci,point))
            {
              for (p = 0; p < nodes.size(); p++)
              {
                stimulusparam_type stimitem;
                stimitem.node = static_cast<unsigned int>(nodes[p]);
                stimitem.weight = elementtypemesh->get_size(*cit)/elementtypemesh->num_nodes_per_face();
                stimulustable.push_back(stimitem);
              }
            }
          }
          ++cit;
        }

        if (isgeomtocomp)
        {
          for (size_t p=0; p < stimulustable.size(); p++) 
          {
            stimulustable[p].node = geomtocompcc[stimulustable[p].node];
          }
        }

        std::sort(stimulustable.begin(),stimulustable.end());

        size_t k = 0;
        for (size_t p=1; p < stimulustable.size(); p++) 
        {
          if (stimulustable[p].node == stimulustable[k].node)
          {
            stimulustable[k].weight += stimulustable[p].weight;
            stimulustable[p].weight = 0.0;
          }
          else
          {
            k++;
            stimulustable[k] = stimulustable[p];
          }
        }
       
        if (stimulustable.size() > 0)       
        {
          stimulustable.resize(k+1);
        }
                       
        // Success:
        return (true);

      }
      else if (fi2.is_volume())
      {
        // Check whether we have all pointers
        VField* elementtypefield = ElementType->vfield();
        VField* stimfield = Stimulus->vfield();

        if (stimfield == 0)
        {
          pr->error("BuildStimulusTable: The Stimulus field is empty");
          return (false);
        }

        if (elementtypefield == 0)
        {
          pr->error("BuildStimulusTable: The domaintype field is empty");
          return (false);
        }

        VMesh* elementtypemesh = ElementType->vmesh();
        if (elementtypemesh == 0)
        {
          pr->error("BuildStimulusTable: The domaintype mesh is empty");
          return (false);
        }

        VMesh* stimmesh = Stimulus->vmesh();
        if (stimmesh == 0)
        {
          pr->error("BuildStimulusTable: The stimmesh is empty");
          return (false);
        }

        // We have pointers and handles

        // Setup algorithm to renumber nodes as the domain has linked nodes somewhere
        VMesh::Node::size_type nnodes;
        bool isgeomtocomp = false;
        Matrix::index_type  *geomtocomprr = 0;
        Matrix::index_type  *geomtocompcc = 0;

        // If it is there we setup the tables
        if (CompToGeom.get_rep())
        {
          // Sanity check..
          elementtypemesh->size(nnodes);
          if (nnodes != CompToGeom->nrows())
          {
            pr->error("BuildStimulusTable: The number of rows in CompToGeom is not of the right size");
            return (false);    
          }
          
          SparseRowMatrix* mat = dynamic_cast<SparseRowMatrix *>(CompToGeom.get_rep());
          if (mat == 0)
          {
            pr->error("BuildStimulusTable: CompToGeom is not a sparse matrix");
            return (false);    
          }
        
          geomtocomprr = mat->rows;
          geomtocompcc = mat->columns;
          
          // WE only support 1-on-1 mappings here
          for (Matrix::index_type r=0; r<nnodes+1;r++)
          {
            if (geomtocomprr[r] != r)
            {
              pr->error("BuildStimulusTable: CompToGeom is mapping a geometrical node to multiple computational nodes");
              return (false);     
            }
          }
          isgeomtocomp = true;
        }


        VMesh::Node::iterator it, it_end;
        VMesh::Node::size_type sz;
        double val, dval;
        VMesh::Elem::array_type elems;
        VMesh::Elem::index_type ci;
        
        Point point;

        elementtypemesh->size(sz);
        elementtypemesh->begin(it);
        elementtypemesh->end(it_end);  
        
        std::vector<bool> indomain(sz);
        double dvalmin = domaintypemin;
        double dvalmax = domaintypemax;
        
        stimmesh->synchronize(Mesh::LOCATE_E|Mesh::CELLS_E);

        // First make a table of which nodes are actually in the domain.

        while (it != it_end)
        {
          indomain[(static_cast<unsigned int>(*it))] = false;
          
          elementtypemesh->get_center(point,*it);
          if (stimmesh->locate(ci,point))
          {
            elementtypemesh->get_elems(elems,*it);
            if (elems.size() > 0)
            {
              elementtypefield->get_value(val,elems[0]);
              if (val >= dvalmin && val <= dvalmax)
              {
                indomain[(static_cast<unsigned int>(*it))] = true;        
              }
            }
          }
          ++it;
        }

        VMesh::Cell::iterator cit, cit_end;
        VMesh::Node::array_type nodes;  
          
        elementtypemesh->begin(cit);
        elementtypemesh->end(cit_end);  
       
        stimulustable.clear();

        // now iterate over each element

        while (cit != cit_end)
        {
          elementtypemesh->get_nodes(nodes,*cit);
          size_t p = 0;
          for (; p < nodes.size(); p++)
          {
            if (indomain[static_cast<unsigned int>(nodes[p])] == false) break;
          }
          if (p == nodes.size())
          {
            elementtypemesh->get_center(point,*cit);
            if (stimmesh->locate(ci,point))
            {
            
              for (p = 0; p < nodes.size(); p++)
              {
                stimulusparam_type stimitem;
                stimitem.node = static_cast<unsigned int>(nodes[p]);
                stimitem.weight = elementtypemesh->get_size(*cit)/elementtypemesh->num_nodes_per_elem();
                stimulustable.push_back(stimitem);
              }
            }
          }
          ++cit;    
        }
        
        
        if (isgeomtocomp)
        {
          for (size_t p=0; p < stimulustable.size(); p++) 
          {
            stimulustable[p].node = geomtocompcc[stimulustable[p].node];
          }
        }

        if (stimulustable.size() > 0)       
        {  
          std::sort(stimulustable.begin(),stimulustable.end());

          size_t k = 0;
          for (size_t p=1; p < stimulustable.size(); p++) 
          {
            if (stimulustable[p].node == stimulustable[k].node)
            {
            stimulustable[k].weight += stimulustable[p].weight;
            stimulustable[p].weight = 0.0;
            }
            else
            {
              k++;
              stimulustable[k] = stimulustable[p];
            }
          }
             
          stimulustable.resize(k+1);
        }

        // Success:
        return (true);      
      }
      else
      {
        return (false);
      }
    }
    else
    {
      VField* elementtypefield = ElementType->vfield();
      VField* stimfield = Stimulus->vfield();

      if (stimfield == 0)
      {
        pr->error("BuildStimulusTable: The Stimulus field is empty");
        return (false);
      }

      if (elementtypefield == 0)
      {
        pr->error("BuildStimulusTable: The domaintype field is empty");
        return (false);
      }

      VMesh* elementtypemesh = elementtypefield->vmesh();
      if (elementtypemesh == 0)
      {
        pr->error("BuildStimulusTable: The domaintype mesh is empty");
        return (false);
      }

      VMesh* stimmesh = stimfield->vmesh();
      if (stimmesh == 0)
      {
        pr->error("BuildStimulusTable: The stimmesh is empty");
        return (false);
      }

      VMesh::Node::size_type nnodes;
      bool isgeomtocomp = false;
      Matrix::index_type  *geomtocomprr = 0;
      Matrix::index_type  *geomtocompcc = 0;

      if (CompToGeom.get_rep())
      {
        elementtypemesh->size(nnodes);
        if (nnodes != CompToGeom->nrows())
        {
          pr->error("BuildStimulusTable: The number of rows in CompToGeom is not of the right size");
          return (false);    
        }
        
        SparseRowMatrix* mat = dynamic_cast<SparseRowMatrix *>(CompToGeom.get_rep());
        if (mat == 0)
        {
          pr->error("BuildStimulusTable: CompToGeom is not a sparse matrix");
          return (false);    
        }
      
        geomtocomprr = mat->rows;
        geomtocompcc = mat->columns;
        
        for (Matrix::index_type r=0; r<nnodes+1;r++)
        {
          if (geomtocomprr[r] != r)
          {
            pr->error("BuildStimulusTable: CompToGeom is mapping a geometrical node to multiple computational nodes");
            return (false);     
          }
        }
        isgeomtocomp = true;
      }

      VMesh::Elem::array_type elems;
      VMesh::Node::iterator it, it_end;
      VMesh::Node::size_type sz;
      VMesh::Elem::index_type ci;
      double val, dval;
      Point point;

      elementtypemesh->size(sz);
      elementtypemesh->begin(it);
      elementtypemesh->end(it_end);  

      std::vector<bool> indomain(sz);

      double dvalmin = domaintypemin;
      double dvalmax = domaintypemax;
      
      stimmesh->synchronize(Mesh::LOCATE_E|Mesh::EDGES_E);

      while (it != it_end)
      {
        indomain[(static_cast<unsigned int>(*it))] = false;
        
        elementtypemesh->get_center(point,*it);
        if (stimmesh->locate(ci,point))
        {
          elementtypemesh->get_elems(elems,*it);
          if (elems.size() > 0)
          {
            elementtypefield->get_value(val,elems[0]);
            if (val >= dvalmin && val <= dvalmax)
            {
              indomain[static_cast<unsigned int>(*it)] = true;        
            }
          }
        }
        ++it;
      }

      VMesh::Elem::iterator cit, cit_end;
      VMesh::Node::array_type nodes;  
        
      elementtypemesh->begin(cit);
      elementtypemesh->end(cit_end);  
     
      stimulustable.clear();
     
      while (cit != cit_end)
      {
        elementtypemesh->get_nodes(nodes,*cit);
        for (size_t p=0; p<nodes.size();p++)
        {
          if (indomain[static_cast<unsigned int>(nodes[p])] == true) 
          {
            stimulusparam_type stimitem;
            stimitem.node = static_cast<unsigned int>(nodes[p]);
            stimitem.weight = elementtypemesh->get_size(*cit)/elementtypemesh->num_nodes_per_elem();
            stimulustable.push_back(stimitem);
          }
        }
        ++cit;
      }

      if (isgeomtocomp)
      {
        for (size_t p=0; p < stimulustable.size(); p++) 
        {
          stimulustable[p].node = geomtocompcc[stimulustable[p].node];
        }
      }

      std::sort(stimulustable.begin(),stimulustable.end());

      size_t k = 0;
      for (size_t p=1; p < stimulustable.size(); p++) 
      {
        if (stimulustable[p].node == stimulustable[k].node)
        {
          stimulustable[k].weight += stimulustable[p].weight;
          stimulustable[p].weight = 0.0;
        }
        else
        {
          k++;
          stimulustable[k] = stimulustable[p];
        }
      }
         
      if (stimulustable.size() > 0)       
      {
        stimulustable.resize(k+1);
      }
                         
      // Success:
      return (true);    
    }
  }

  return (false);
}


} // End namespace ModelCreation
