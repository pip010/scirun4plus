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
//    File   : InterfaceWithTetGen.cc
//    Author : David Weinstein
//    Date   : Mon Jan  5 16:22:37 MST 2009

#include <Core/Algorithms/Fields/MedialAxis/MedialAxis2.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Core/Persistent/Pstreams.h>
#include <Core/Util/StringUtil.h>

#include <iostream>
#include <fstream>

#include <stdlib.h>

using std::cerr;
using std::ifstream;
using std::endl;

using namespace SCIRun;

SCIRunAlgo::MedialAxis2Algo algo;
std::string surface_fname;
std::string output_fname;

bool parseArgs(int argc, char *argv[]) 
{
  int curr=1;
  while (curr<argc-1) 
  {
    if (argv[curr] == std::string("-surface")) 
    {
      surface_fname = std::string(argv[curr+1]);
      curr+=2;
    } 
    else if (argv[curr] == std::string("-output")) 
    {
      output_fname = std::string(argv[curr+1]);
      curr+=2;
    } 
    else if (argv[curr] == std::string("-levels")) 
    {
      int levels;
      from_string(std::string(argv[curr+1]),levels);
      algo.set_int("refinement_levels",levels);
      curr+=2;
    } 
    else if (argv[curr] == std::string("-refinement_levels")) 
    {
      int levels;
      from_string(std::string(argv[curr+1]),levels);
      algo.set_int("refinement_levels",levels);
      curr+=2;
    } 
    else if (argv[curr] == std::string("-start_num_bins")) 
    {
      int start_num_bins;
      from_string(std::string(argv[curr+1]),start_num_bins);
      algo.set_int("start_num_bins",start_num_bins);
      curr+=2;
    } 
    else if (argv[curr] == std::string("-normal_minimum_angle")) 
    {
      double normal_minimum_angle;
      from_string(std::string(argv[curr+1]),normal_minimum_angle);
      algo.set_int("normal_minimum_angle",normal_minimum_angle);
      curr+=2;
    } 
    else if (argv[curr] == std::string("-axis_minimum_angle")) 
    {
      double axis_minimum_angle;
      from_string(std::string(argv[curr+1]),axis_minimum_angle);
      algo.set_int("axis_minimum_angle",axis_minimum_angle);
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
  std::cerr << "Usage: "<<appName<<" -surface surface_filename -output points_filename"
  << " [-levels (4)]" << " [-start_num_bins (100)]"
  << " [-normal_minimum_angle (50)]" << " [-axis_minimum_angle (50)]"
  << "\n";
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
  
  if (!surface_fname.length()) 
  {
    std::cerr << "Error -- no '-surface' argument found on command line\n";
    printUsage(argv[0]);
    exit(1);
  }
  
  if (!output_fname.length()) 
  {
    std::cerr << "Error -- no '-output' argument found on command line\n";
    printUsage(argv[0]);
    exit(1);
  }
  
  FieldHandle surface;
  PiostreamPtr surface_stream = auto_istream(surface_fname);
  if (!surface_stream) 
  {
    std::cerr << "Error: couldn't open surface file "<<surface_fname<<std::endl;
    exit(1);
  }
  
  Pio(*surface_stream, surface);
  
  if (!surface.get_rep()) 
  {
    std::cerr << "Error reading field from file " << surface_fname
      << ".  Exiting..." << std::endl;
    exit(1);
  }
  
  FieldHandle output;
  if(!(algo.run(surface, output)))
  {
    std::cerr << "algorithm failed. Exiting...." << std::endl;
    exit(1);
  }

  BinaryPiostream output_stream(output_fname, Piostream::Write);
  Pio(output_stream, output);

  return 0;  
}    
