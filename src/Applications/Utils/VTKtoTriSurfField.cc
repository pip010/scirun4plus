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

//    File   : VTKtoTriSurfField.cc
//    Author : Martin Cole
//    Date   : Fri May  7 10:23:05 2004


// Warning!: this converter is only partially implemented. It supports only 
//           a subset of the vtk format. Please extend it as you need more

#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/Field.h>
#include <Core/Persistent/Pstreams.h>

#include <Core/Init/init.h>
#include <Core/Util/FileUtils.h>

#include <sci_deprecated.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

using namespace SCIRun;


#define check_error(str) \
  if (str.fail()) { \
    std::cerr << "fail state at line: " << __LINE__ << " " << std::endl; \
    return 66; }

bool bin_out;
bool swap_endian;

void setDefaults() 
{
  bin_out=false;
  swap_endian=false;
}

int parseArgs(int argc, char *argv[]) 
{
  int currArg = 3;
  while (currArg < argc) 
  {
    if (!strcmp(argv[currArg], "-bin_out")) 
    {
      bin_out=true;
      currArg++;
    } 
    else if (!strcmp(argv[currArg], "-swap_endian")) 
    {
      swap_endian=true;
      currArg++;
    } 
    else 
    {
      std::cerr << "Error - unrecognized argument: "<<argv[currArg]<<"\n";
      return 0;
    }
  }
  return 1;
}

void printUsageInfo(char *progName) 
{
  std::cerr << std::endl 
       << "Usage: " << progName 
       << " VTKFile TriSurfField [-bin_out] [-swap_endian]" 
       << std::endl << std::endl
       << "\t This program will read in a .vtk file and output a TriSurfField"
       << std::endl
       << "\t -bin_out specifies binare vs. ASCII output. ASCII by default."
       << std::endl 
       << "\t -swap_endian swaps the data from big to little endian or" 
       << std::endl
       << "\t vice versa. false by default."
       << std::endl << std::endl;

  std::cerr << "Warning!: this converter is only partially implemented. It supports "
       << "only a subset of the vtk format. Please extend it as you need more."
       << std::endl;
}


void swap_endianess_4(unsigned *dw)
{
  if (!swap_endian) return;
  register unsigned tmp;
  tmp =  (*dw & 0x000000FF);
  tmp = ((*dw & 0x0000FF00) >> 0x08) | (tmp << 0x08);
  tmp = ((*dw & 0x00FF0000) >> 0x10) | (tmp << 0x08);
  tmp = ((*dw & 0xFF000000) >> 0x18) | (tmp << 0x08);
  *dw = tmp;
}

int
read_n_points(VMesh *tsm, int n, std::ifstream &str) 
{
  float arr[3];
  check_error(str);

  for(int i = 0; i < n; i++) 
  {
    str.read((char*)arr, sizeof(float) * 3);
    check_error(str);
    swap_endianess_4((unsigned*)&arr[0]);
    swap_endianess_4((unsigned*)&arr[1]);
    swap_endianess_4((unsigned*)&arr[2]);
    tsm->add_point(Point(arr[0], arr[1], arr[2]));
    //std::cout << arr[0] << ", " << arr[1] << ", " << arr[2] << std::endl;
  }
  return 0;
}


int
read_n_polys(VMesh *tsm, int n, std::ifstream &str) 
{
  unsigned int arr[3];
  check_error(str);
  VMesh::Node::array_type nodes(3);
  for(int i = 0; i < n; i++) 
  {
    unsigned int val;
    str.read((char*)&val, sizeof(float));
    check_error(str);
    swap_endianess_4(&val);

    if (val != 3) 
    {
      std::cout << "ERROR: can only handle triangle polys atm..." << std::endl;
      exit(1);
    }

    str.read((char*)arr, sizeof(unsigned int) * 3);
    check_error(str);

    swap_endianess_4((unsigned int*)&arr[0]);
    swap_endianess_4((unsigned int*)&arr[1]);
    swap_endianess_4((unsigned int*)&arr[2]);
    
    nodes[0] = static_cast<VMesh::Node::index_type>(arr[0]);
    nodes[1] = static_cast<VMesh::Node::index_type>(arr[1]);
    nodes[2] = static_cast<VMesh::Node::index_type>(arr[2]);
    tsm->add_elem(nodes);
  }
  return 0;
}

int
read_n_strips(VMesh *tsm, int n, std::ifstream &str) 
{
  unsigned int arr[3];
  check_error(str);
  std::cerr << "reading " << n << " strips." << std::endl;
  for(int i = 0; i < n; i++) 
  {
    unsigned int val;
    str.read((char*)&val, sizeof(unsigned int));
    check_error(str);
    swap_endianess_4(&val);

    std::cout << val << " points in strip." << std::endl; 
    std::vector<int> strip;
    strip.resize(val);
    // read in all the indicies
    str.read((char*)&strip[0], sizeof(unsigned int) * val);
    check_error(str);

    bool ccw = true;

    VMesh::Node::array_type nodes(3);

    for (int j = 2; j < static_cast<int>(val); ++j)
    {
      arr[0] = strip[j-2];
      arr[1] = strip[j-1];
      arr[2] = strip[j];
      
      swap_endianess_4((unsigned*)&arr[0]);
      swap_endianess_4((unsigned*)&arr[1]);
      swap_endianess_4((unsigned*)&arr[2]);
    

      if (ccw) 
      {
        nodes[0] = static_cast<VMesh::Node::index_type>(arr[0]);
        nodes[1] = static_cast<VMesh::Node::index_type>(arr[1]);
        nodes[2] = static_cast<VMesh::Node::index_type>(arr[2]);
      } 
      else 
      {
        nodes[0] = static_cast<VMesh::Node::index_type>(arr[0]);
        nodes[1] = static_cast<VMesh::Node::index_type>(arr[2]);
        nodes[2] = static_cast<VMesh::Node::index_type>(arr[1]);
      }

      tsm->add_elem(nodes);
      ccw = !ccw;
    }
  }
  return 0;
}


int
read_scalar_lookup(VField *fld, VField::size_type n, std::ifstream &str) 
{
  fld->resize_values();
  //val_t last = 0;
  int vset = 0;
  check_error(str);

  for(VField::index_type i = 0; i < n; i++) 
  {
    float val;
    str.read((char*)&val, sizeof(float));
    check_error(str);
    swap_endianess_4((unsigned int*)&val);
    vset++;
    fld->set_value(val, i);    
  }
  return 0;
}


int
main(int argc, char **argv) 
{

  int status = 0;
  if (argc < 3 || argc > 5) 
  {
    printUsageInfo(argv[0]);
    return 2;
  }

  SCIRunInit();
  setDefaults();

  MeshHandle tsm = CreateMesh(TRISURFMESH_E);
  //exit(5);
  char *in = argv[1];
  char *out = argv[2];
  if (!parseArgs(argc, argv)) 
  {
    printUsageInfo(argv[0]);
    status = 1;
    return 1;
  }
  std::ifstream vtk(in, std::ios::binary);
  check_error(vtk);

  char id[256], header[256], format[256];

  vtk.getline(id, 256);
  vtk.getline(header, 256);

  vtk >> format;
  std::cout << format << std::endl;  

  std::string dataset;
  vtk >> dataset;

  if (dataset != "DATASET") 
  {
    std::cerr << "ERROR: expected DATASET keyword." << std::endl;
    status = 5;
    return 5;
  } 
  else 
  {    
    vtk >> dataset;
    if (dataset != "POLYDATA") 
    {
      std::cerr << "ERROR: can only handle POLYDATA type vtk files at this time. "
        << "got: " << dataset << std::endl;
      return 5;
    }
  }
  check_error(vtk);
  
  std::string attrib;
  vtk >> attrib;

  int n;
  vtk >> n;

  std::string type;
  vtk >> type;
  //std::cout << "attrib is : " << attrib << " " << n << " type is " << type << std::endl;
  check_error(vtk);
  vtk.get(); // eat a newline

  read_n_points(tsm->vmesh(), n, vtk);
  
  check_error(vtk);

  std::string poly;
  vtk >> poly;

  if (poly == "POLYGONS") 
  {
    vtk >> n;
    int sz;
    vtk >> sz;
    std::cout << poly << " " << n << " " << sz << std::endl;
    
    vtk.get(); // eat a newline
    read_n_polys(tsm->vmesh(), n, vtk);
    check_error(vtk);
    
  } 
  else if (poly == "TRIANGLE_STRIPS") 
  {
    vtk >> n;
    int sz;
    vtk >> sz;
    std::cout << poly << " " << n << " " << sz << std::endl;
    
    vtk.get(); // eat a newline
    read_n_strips(tsm->vmesh(), n, vtk);
    check_error(vtk);
    
  } 
  else 
  {
    std::cerr << "ERROR: can only handle POLYGONS or TRIANGLE_STRIPS type polygonal data " 
      << " files at this time. got: " << poly << std::endl;
    return 1;
  }
  
  std::string dat;
  std::string d_attrib;
  vtk >> dat; 
  vtk >> n;

  check_error(vtk)
  std::cout << dat << " " << n << std::endl;

  if (dat == "CELL_DATA") 
  {
    vtk >> d_attrib;
    check_error(vtk);

    if (d_attrib == "POINT_DATA") 
    {	
      dat = d_attrib;
      vtk >> n;
    }
    std::cout << dat << " " << n << std::endl;
  }
  
  FieldHandle ts_handle;
  std::string data, name;
  vtk >> data;
  check_error(vtk);
  std::cout << data << std::endl;
  vtk >> name;
  check_error(vtk);
  std::cout << name << std::endl;
  vtk >> type;
  check_error(vtk);
  std::cout << type << std::endl;

  if (dat == "CELL_DATA") 
  {
    if (type != "float") 
    {
      std::cerr << "supporting float only atm..." << std::endl;
      return 1;
    }
    
    FieldInformation fi(TRISURFMESH_E,CONSTANTDATA_E,FLOAT_E);
    ts_handle = CreateField(fi,tsm);
    
    ts_handle->vfield()->resize_values();
    std::string table, tname;
    vtk >> table >> tname;
    vtk.get(); // eat a newline
    read_scalar_lookup(ts_handle->vfield(), n, vtk);
  } 
  else 
  {
    // node centered data ...
    if (type != "float") 
    {
      std::cerr << "supporting float only atm..." << std::endl;
      std::cerr << "got " << type << std::endl;
      return 1;
    }
    if (data == "SCALARS") 
    {
      FieldInformation fi(TRISURFMESH_E,LINEARDATA_E,FLOAT_E);
      ts_handle = CreateField(fi,tsm);
      
      ts_handle->vfield()->resize_values();
      std::string table, tname;
      vtk >> table >> tname;
      //std::cout << table << " " << tname << std::endl;
      vtk.get(); // eat a newline
      read_scalar_lookup(ts_handle->vfield(), n, vtk);
    } 
    else
    {
      FieldInformation fi(TRISURFMESH_E,NODATA_E,FLOAT_E);
      ts_handle = CreateField(fi,tsm);
      
      ts_handle->vfield()->resize_values();
    }
  }

  while (!vtk.eof()) 
  {
    vtk.get();
  }

  if (bin_out) 
  {
    BinaryPiostream out_stream(out, Piostream::Write);
    Pio(out_stream, ts_handle);
  } 
  else 
  {
    TextPiostream out_stream(out, Piostream::Write);
    Pio(out_stream, ts_handle);
  }

  return status;  
}    
