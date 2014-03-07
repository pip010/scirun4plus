//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : GenerateMPMData.cc
//    Author : Martin Cole
//    Date   : Wed Jun  2 16:08:19 2004

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Persistent/Pstreams.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sys/stat.h>
#include <map>
#include <sstream>

using namespace SCIRun;

std::vector<std::ofstream*> *files = 0;
std::vector<std::ofstream*> *aux_files = 0;
std::vector<unsigned> fcount;
std::vector<unsigned> num_files;
std::vector<std::string> prename;

const unsigned int max_lines = 25000000;
double part_per_mm = 4.;

#ifdef _WIN32
#define round(x) (int) (x+.5)
#endif

void
write_MPM_fibdir(FieldHandle& field_h, const std::string &outdir)
{
  // this outputs only one material
  VField *ifield = field_h->vfield();
  VMesh* mesh = field_h->vmesh();

  files = new std::vector<std::ofstream*>;
  aux_files = new std::vector<std::ofstream*>;
  fcount.push_back(0);
  num_files.push_back(1);
  std::stringstream fname(std::stringstream::in | std::stringstream::out);
  std::string nm;
  fname << outdir << "/" << "heart-diastole";
  fname >> nm;
  prename.push_back(nm);

  fname.clear();
  fname << nm << "-"; 
  fname.fill('0');
  fname.width(4);
  fname << std::right << 1;
  std::string nm1;
  fname >> nm1;
  std::string aux_nm = nm1 + ".elems";
  files->push_back(new std::ofstream(nm1.c_str(), std::ios_base::out));
  aux_files->push_back(new std::ofstream(aux_nm.c_str(), std::ios_base::out));

  mesh->synchronize(Mesh::NODES_E | Mesh::EDGES_E | Mesh::FACES_E | Mesh::CELLS_E); 
	VMesh::Elem::iterator iter, end;
  mesh->begin( iter );
  mesh->end( end );
  
  while (iter != end) 
  {
    Point o, s, t, u;
    VMesh::coords_type coord(3, 0.0);
    mesh->interpolate(o, coord, *iter);
    coord[0] = 1.0;
    mesh->interpolate(s, coord, *iter);
    coord[0] = 0.0;
    coord[1] = 1.0;
    mesh->interpolate(t, coord, *iter);
    coord[1] = 0.0;
    coord[2] = 1.0;
    mesh->interpolate(u, coord, *iter);
    BBox elem_bbox;
    elem_bbox.extend(o);
    elem_bbox.extend(s);
    elem_bbox.extend(t);
    elem_bbox.extend(u);


    // assume length in mm's
    double dx = (s - o).length();
    double dy = (t - o).length();
    double dz = (u - o).length();

    double vol = dx * dy * dz; 

    const int div_per_x = (int)round(part_per_mm * dx);
    const int div_per_y = (int)round(part_per_mm * dy);
    const int div_per_z = (int)round(part_per_mm * dz);

    const double part_vol = vol / (div_per_x * div_per_y * div_per_z);
    
    dx = 1. / (double)div_per_x;
    dy = 1. / (double)div_per_y;
    dz = 1. / (double)div_per_z;

    for (double  k = dz * 0.5; k < 1.0; k += dz) 
    {
      for (double j = dy * 0.5; j < 1.0; j += dy) 
      {
        for (double i = dx * 0.5; i < 1.0; i += dx) {
          
          coord[0] = i;
          coord[1] = j;
          coord[2] = k;
          //cerr << "at : [" << i << "," << j << "," << k << "]" << endl;
          Point c;
          mesh->interpolate(c, coord, *iter);
          //sanity check..
          // 	  if (! elem_bbox.inside(c)) {
          // 	    cerr << "Error: interpolated to a point outside the element!" 
          // 		 << endl;
          // 	  }

          Vector v;
          ifield->interpolate(v, coord, *iter);

          if (fcount[0] > max_lines) {
            std::string nm;
            std::stringstream fname(std::stringstream::in | std::stringstream::out);
            fname << prename[0] << "-";
            fname.fill('0');
            fname.width(4);
            fname << std::ios::right << ++num_files[0];
            fname >> nm;
            std::string aux_nm = nm + ".elems";
            delete (*files)[0];
            (*files)[0] = new std::ofstream(nm.c_str(), std::ios_base::out);
            (*aux_files)[0] = new std::ofstream(aux_nm.c_str(), std::ios_base::out);
            fcount[0] = 0;
          }

          std::ofstream* str = (*files)[0];
          std::ofstream* aux_str = (*aux_files)[0];
          (*str) << std::setprecision (9) 
           << c.x() << " " << c.y() << " " << c.z() << " " 
           << part_vol << " "
           << v.x() << " " << v.y() << " " << v.z() << std::endl;

          // aux file just writes the elem number that the particle belongs to.
          (*aux_str) << *iter << std::endl;
          fcount[0]++;
        }
      }
    }
    std::cout << "." << *iter << ".";
    std::cout.flush();
    ++iter;
  }
  std::cout << std::endl;
}

void
write_MPM(FieldHandle& field_h, const std::string &outdir)
{
  VField *ifield = field_h->vfield();
  VMesh* mesh =    field_h->vmesh();
  
  std::vector<std::pair<std::string, Tensor> > field_tensors;
  files = new std::vector<std::ofstream*>;
  
  if (ifield->get_property("conductivity_table", field_tensors)) 
  {
    std::vector<std::pair<std::string, Tensor> >::iterator titer = field_tensors.begin();
    while (titer != field_tensors.end()) 
    {
      fcount.push_back(0);
      num_files.push_back(1);
      std::stringstream fname(std::stringstream::in | std::stringstream::out);
      std::string nm;
      fname << outdir << "/" << titer++->first;
      fname >> nm;
      prename.push_back(nm);

      fname.clear();
      fname << nm << "-"; 
      fname.fill('0');
      fname.width(4);
      fname << std::ios::right << 1;
      std::string nm1;
      fname >> nm1;
      files->push_back(new std::ofstream(nm1.c_str(), std::ios::out));

    }    
  } 
  else 
  {
    // no properties..
    fcount.push_back(0);
    num_files.push_back(1);
    std::stringstream fname(std::stringstream::in | std::stringstream::out);
    std::string nm;
    fname << outdir << "/unknown_material";
    fname >> nm;
    prename.push_back(nm);

    fname.clear();
    fname << nm << "-"; 
    fname.fill('0');
    fname.width(4);
    fname << std::ios::right << 1;
    std::string nm1;
    fname >> nm1;
    files->push_back(new std::ofstream(nm1.c_str(), std::ios::out));
  }

  const double max_vol_s = 1.0L;
  BBox bbox = field_h->vmesh()->get_bounding_box();
  Point minb, maxb;
  minb = bbox.min();
  maxb = bbox.max();

  int sizex = (int)(round(maxb.x() - minb.x()) / max_vol_s);
  int sizey = (int)(round(maxb.y() - minb.y()) / max_vol_s);
  int sizez = (int)(round(maxb.z() - minb.z()) / max_vol_s);
  
  FieldInformation fi("LatVolMesh",1,"double");
  MeshHandle meshH = CreateMesh(fi,sizex, sizey, sizez, minb, maxb);

  VMesh *cntrl = meshH->vmesh();
  
  double vol = cntrl->get_size(VMesh::Elem::index_type(0));
  
  VMesh::Cell::iterator iter, end;
  cntrl->begin( iter );
  cntrl->end( end );
  int count = 0;
  
  while (iter != end) 
  {
    Point c;
    cntrl->get_center(c, *iter);
    ++iter; count++;
    VMesh::Elem::index_type ci;
    if (mesh->locate(ci, c)) 
    {
      int file_index = 0;
      if (ifield->basis_order() > 0) 
      {
        ifield->get_value(file_index,ci);
      }
      if (fcount[file_index] > max_lines) 
      {
        std::string nm;
        std::stringstream fname(std::stringstream::in | std::stringstream::out);
        fname << prename[file_index] << "-";
        fname.fill('0');
        fname.width(4);
        fname << std::ios::right << ++num_files[file_index];
        fname >> nm;
        delete (*files)[file_index];
        (*files)[file_index] = new std::ofstream(nm.c_str(), std::ios::out);
        fcount[file_index] = 0;
      }

      std::ofstream* str = (*files)[file_index];
      (*str) << std::setprecision (9) <<c.x() << " " << c.y() << " " << c.z() 
	     << " " << vol << std::endl;
      fcount[file_index]++;
    }
    
    if (count % ((int)(sizex*sizey*sizez*0.01)) == 0) { std::cout << "."; }
  }    
  std::cout << std::endl;
}





int
main(int argc, char **argv)
{
  if (argc != 4)
  {
    std::cout << "Usage:  GenerateMPMData <input-field> <dest-dir> <num particles per millimeter>" << std::endl;
    exit(0);
  }

  struct stat buff;
  if (stat(argv[1], &buff) == -1)
  {
    std::cout << "File " << argv[1] << " not found\n";
    exit(99);
  }

  if (stat(argv[2], &buff) == -1)
  {
    std::cout << "Directory  " << argv[2] << " not found\n";
    exit(100);
  }
  std::string ppm(argv[3]);
  part_per_mm = atof(ppm.c_str());


  std::string outdir(argv[2]);

  FieldHandle field_handle;
  
  PiostreamPtr in_stream = auto_istream(argv[1]);
  if (!in_stream)
  {
    std::cout << "Error reading file " << argv[1] << ".\n";
    exit(101);
  }

  Pio(*in_stream, field_handle);

  write_MPM_fibdir(field_handle, outdir);

  return 0;
}    
