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

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/Field.h>
#include <Core/Geometry/Transform.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Point.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <Core/Algorithms/Fields/TransformMesh/TransformMeshWithTransform.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace SCIRun;

std::string input_field_fname;
std::string output_field_fname;
std::string input_matrix_fname;

bool has_input_field = false;
bool has_output_field = false;
bool has_input_matrix = false;

bool parseArgs(int argc, char *argv[]) 
{
  int curr=1;
  while (curr<argc-1) 
  {
    if (argv[curr] == std::string("-input")) 
    {
      input_field_fname = argv[curr+1];
      has_input_field = true;
      curr+=2;
    } 
    else if (argv[curr] == std::string("-output")) 
    {
      output_field_fname = argv[curr+1];
      has_output_field = true;
      curr+=2;
    } 
    else if (argv[curr] == std::string("-transform")) 
    {
      input_matrix_fname = argv[curr+1];
      has_input_matrix = true;
      curr+=2;
    } 
    else 
    {
      std::cerr << "Error -- invalid argument: "<<argv[curr]<<"\n";
      return false;
    }
  }
  return true;
}


void
printUsage(char *appName) 
{
  std::cerr << "Usage: "<<appName <<" -input field_file -output field_file -transform transform_file\n";
  return;
}


int
main(int argc, char *argv[]) 
{
  if (!parseArgs(argc, argv)) 
  {
    printUsage(argv[0]);
    exit(1);
  }

  if (!has_input_field || !has_output_field || !has_input_matrix)
  {
    printUsage(argv[0]);
    exit(1);
  }

  std::cout << "Reorienting the output field to match the original nrrd." << std::endl;
               
  SCIRunAlgo::TransformMeshWithTransformAlgo algo;
  
  FieldHandle mesh;
  PiostreamPtr mesh_stream = auto_istream(input_field_fname);
  if (!mesh_stream) 
  {
    std::cerr << "Error: couldn't open field file "<<input_field_fname<< std::endl;
    exit(1);
  }

  Pio(*mesh_stream, mesh);
  if (!mesh.get_rep()) 
  {
    std::cerr << "Error reading field from file " << input_field_fname
              << ".  Exiting..." << std::endl;
    exit(1);
  }

  MatrixHandle tf;
  PiostreamPtr matrix_stream = auto_istream(input_matrix_fname);
  if (!matrix_stream) 
  {
    std::cerr << "Error: couldn't open transform matrix file "<<input_matrix_fname << std::endl;
    exit(1);
  }

  Pio(*matrix_stream, tf);
  if (!tf.get_rep()) 
  {
    std::cerr << "Error reading matrix transform from file " << input_matrix_fname << ".  Exiting..." << std::endl;
    exit(1);
  }  
  
  FieldHandle trans_mesh;
  if(!(algo.run(mesh,tf,trans_mesh)))
  {
    std::cerr << "Transformation algorithm failed" << std::endl;
  }

  ProgressReporter* pr = new ProgressReporter;
  PiostreamPtr stream = auto_ostream(output_field_fname, "Binary", pr);
  
  if (stream->error())
  {
    std::cerr << "Could not open file for writing " << output_field_fname << std::endl;
    return (1);
  }
  else
  {
    // Write the file
    Pio(*stream, trans_mesh);
  } 

  // Yeah it worked !!
  return (0);  
}    

