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
 *  VistaToNrrd.cc: create a nrrd header for a vista .v file
 *
 *  Written by:
 *   Mark Hartner
 *   SCI Institute
 *   University of Utah
 *   December 2004
 *
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <sci_deprecated.h>



int main(int argc, char **argv){

  int nrrdSkip = 0;
  char *in = argv[1];
  std::ifstream vistaFileStream(in, std::ios::binary);
  char *out = argv[2];
  std::ofstream nrrdHeader(out, std::ios::binary);
  std::string z_dim("0");
  std::string y_dim("0");
  std::string x_dim("0");
  
  if (argc < 3) {
    std::cerr << "Usage: "<<argv[0]<<" inputFile.v outputFile.nhdr \n";
    return 2;
  }

	
  if (! vistaFileStream){
    std::cout<<"could not input vista file "<<in<<std::endl;
    return 1;
  }
  if (! nrrdHeader){
    std::cout<<"could not create nrrd header file "<<out<<std::endl;
    return 1;
  }
	
  int lineCount = 0;
  int foundFormFeed = 0;	
  while (! vistaFileStream.eof() && !foundFormFeed){
    std::string current_line;
    getline(vistaFileStream,current_line);
    lineCount++;

    if (current_line[0] == '\f'){
      nrrdSkip = lineCount;
      foundFormFeed = 1;
    }

    {
      //look for the Z factor line because we need that for NRRD spacing
      std::string z_marker("nframes: ");
      std::string::size_type z_factor_pos = current_line.find(z_marker);
      if ( z_factor_pos != std::string::npos ){
	z_dim = current_line.substr(z_factor_pos + z_marker.size());
      }
    }
    
    //look for the Y factor line because we need that for NRRD spacing
    {
      std::string y_marker("nrows: ");
      std::string::size_type y_factor_pos = current_line.find(y_marker);
      if ( y_factor_pos != std::string::npos ){
	y_dim = current_line.substr(y_factor_pos + y_marker.size());
      }
    }

    //look for the X factor line because we need that for NRRD spacing
    {
      std::string x_marker("ncolumns: ");
      std::string::size_type x_factor_pos = current_line.find(x_marker);
      if ( x_factor_pos != std::string::npos ){
	x_dim = current_line.substr(x_factor_pos + x_marker.size());
      }
    }
    


  }
  vistaFileStream.close();

  nrrdHeader<<"NRRD0001"<<std::endl
	    <<"type: uchar"<<std::endl
	    <<"dimension: 3"<<std::endl
	    <<"spacings: 1 1 1"<<std::endl
	    <<"sizes: "<<x_dim<<" "<<y_dim<<" "<<z_dim<<std::endl
	    <<"data file: "<<in<<std::endl
	    <<"endian: big"<<std::endl
	    <<"encoding: raw"<<std::endl
	    <<"line skip: "<<nrrdSkip<<std::endl;

  nrrdHeader.close();
  return 0;
}
   
