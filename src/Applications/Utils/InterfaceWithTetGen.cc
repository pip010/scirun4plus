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

#include <Core/Algorithms/Fields/TetGen/InterfaceWithTetGen.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Core/Persistent/Pstreams.h>

#include <iostream>
#include <fstream>

#include <stdlib.h>

using std::cerr;
using std::ifstream;
using std::endl;

using namespace SCIRun;

std::string main_surface_fname;
std::vector<std::string> region_fnames;
std::string points_fname;
std::string cmmd_ln;
std::string region_attribs_fname;
std::string tetvol_fname;

bool parseArgs(int argc, char *argv[]) {
  int curr=1;
  while (curr<argc-1) {
    if (argv[curr] == std::string("-cmmd_line")) {
      cmmd_ln = argv[curr+1];
      curr+=2;
    } else if (argv[curr] == std::string("-main")) {
      main_surface_fname = argv[curr+1];
      curr+=2;
    } else if (argv[curr] == std::string("-region_attribs")) {
      region_attribs_fname = argv[curr+1];
      curr+=2;
    } else if (argv[curr] == std::string("-points")) {
      points_fname = argv[curr+1];
      curr+=2;
    } else if (argv[curr] == std::string("-out")) {
      tetvol_fname = argv[curr+1];
      curr+=2;
    } else if (argv[curr] == std::string("-regions")) {
      int num_regions;
      if (sscanf(argv[curr+1], "%d", &num_regions) != 1) {
	cerr << "Error -- couldn't parse "<<argv[curr+2]<<" as a number (num regions)\n";
	return false;
      }
      curr+=2;
      if (curr+num_regions > argc) {
	cerr << "Error -- listed "<<num_regions<<" regions, but only "<<argc-curr<<" arguments left on command line\n";
	return false;
      }
      for (int i=0; i<num_regions; i++, curr++) {
	region_fnames.push_back(std::string(argv[curr]));
      }
    } else {
      cerr << "Error -- invalid argument: "<<argv[curr]<<"\n";
      return false;
    }
  }
  return true;
}

void
printUsage(char *appName) {
  cerr << "Usage: "<<appName<<" -cmmd_line cmmd_line_string -main main_surface_filename -out tetvol_filename [-region_attribs region_attribs_filename] [-points points_filename] [-regions num_regions region1_filename region2_filename ...]\n";
  return;
}

int
main(int argc, char *argv[]) {
  if (!parseArgs(argc, argv)) {
    printUsage(argv[0]);
    exit(1);
  }
  if (!cmmd_ln.length()) {
    cerr << "Error -- no '-cmmd_line' argument found on command line\n";
    printUsage(argv[0]);
    exit(1);
  }
  if (!main_surface_fname.length()) {
    cerr << "Error -- no '-main' argument found on command line\n";
    printUsage(argv[0]);
    exit(1);
  }
  if (!tetvol_fname.length()) {
    cerr << "Error -- no '-out' argument found on command line\n";
    printUsage(argv[0]);
    exit(1);
  }
  FieldHandle main_surface;
  PiostreamPtr main_surface_stream = auto_istream(main_surface_fname);
  if (!main_surface_stream) {
    cerr << "Error: couldn't open main_surface file "<<main_surface_fname<<endl;
    exit(1);
  }
  Pio(*main_surface_stream, main_surface);
  if (!main_surface.get_rep()) {
    cerr << "Error reading field from file " << main_surface_fname
	 << ".  Exiting..." << endl;
    exit(1);
  }
  std::vector<FieldHandle> regions;
  for (unsigned int i=0; i<region_fnames.size(); i++) {
    PiostreamPtr region_stream = auto_istream(region_fnames[i]);
    if (!region_stream) {
      cerr << "Error: couldn't open region file "<<region_fnames[i]<<endl;
      exit(1);
    }
    FieldHandle temp;
    Pio(*region_stream, temp);
    if (!temp.get_rep()) {
      cerr << "Error reading field from file " << region_fnames[i]
	   << ".  Exiting..." << endl;
      exit(1);
    }
    regions.push_back(temp);
  }
  FieldHandle points;
  if (points_fname.length()) {
    PiostreamPtr points_stream = auto_istream(points_fname);
    if (!points_stream) {
      cerr << "Error: couldn't open points_fname file "<<points_fname<<endl;
      exit(1);
    }
    Pio(*points_stream, points);
    if (!points.get_rep()) {
      cerr << "Error reading field from file " << points_fname
	   << ".  Exiting..." << endl;
      exit(1);
    }
  }
  FieldHandle region_attribs;
  if (region_attribs_fname.length()) {
    PiostreamPtr region_attribs_stream = auto_istream(region_attribs_fname);
    if (!region_attribs_stream) {
      cerr << "Error: couldn't open region_attribs file "<<region_attribs_fname<<endl;
      exit(1);
    }
    Pio(*region_attribs_stream, region_attribs);
    if (!region_attribs.get_rep()) {
      cerr << "Error reading field from file " << region_attribs_fname
	   << ".  Exiting..." << endl;
      exit(1);
    }
  }
  SCIRunAlgo::InterfaceWithTetGenAlgo algo;
  FieldHandle tetvol_out;
  algo.set_string("cmmd_ln", cmmd_ln);
  algo.run(main_surface, regions, points, region_attribs, tetvol_out);
  BinaryPiostream tetvol_stream(tetvol_fname, Piostream::Write);
  Pio(tetvol_stream, tetvol_out);
  return 0;  
}    
