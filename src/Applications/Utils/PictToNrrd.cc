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
 *  pict2nrrd.cc: create a nrrd header for a pict file
 *  this convrter is called by the Nrrd Reader Module
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


const int FIRST_HEADER_BYTESIZE = 76;
const int SECOND_HEADER_MAX_BYTESIZE = 480;

int main(int argc, char **argv){

  int x=0, y=0, z=0;
  char *in = argv[1];
  std::ifstream picFileStream(in, std::ios::binary);
  char *out = argv[2];
  std::ofstream nrrdHeader(out, std::ios::binary);

  if (argc < 3) {
    std::cerr << "Usage: "<<argv[0]<<" inputFile.pic outputFile.nhdr \n";
    return 2;
  }

	
  if (! picFileStream){
    std::cout<<"could not input pic file "<<in<<std::endl;
    return 1;
  }
  if (! nrrdHeader){
    std::cout<<"could not create nrrd header file "<<out<<std::endl;
    return 1;
  }
	
  unsigned char tmp_uchar1;
  unsigned char tmp_uchar2;

  //read data from the first header
  picFileStream>>tmp_uchar1;
  picFileStream>>tmp_uchar2;
  x = tmp_uchar1 + (tmp_uchar2 * 256);
  picFileStream>>tmp_uchar1;
  picFileStream>>tmp_uchar2;
  y = tmp_uchar1 + (tmp_uchar2 * 256);
  picFileStream>>tmp_uchar1;
  picFileStream>>tmp_uchar2;
  z = tmp_uchar1 + (tmp_uchar2 * 256);
  
  //seek to the 2nd header
  int dataSize = x*y*z;
  int startOf2ndHeader = FIRST_HEADER_BYTESIZE +  dataSize;
  picFileStream.seekg(startOf2ndHeader);
  int endOf2ndHeader = startOf2ndHeader + SECOND_HEADER_MAX_BYTESIZE;
  //there is actuall more data beyond the end of the 2nd header, but
  //it looks like data on the settings for the aquisition machine.

  //parse the 2nd header
  std::string z_factor;
  std::string comments;
  char out_tmp;
  int end_of_header = 0;
  picFileStream>>out_tmp;
  while (!picFileStream.eof() &&
	 !end_of_header &&
	 (picFileStream.tellg() < endOf2ndHeader)
	 ){

    if ((int)out_tmp != -1){
      end_of_header = 1;
      std::cout<<"Found end of header!"<<std::endl;
    }
    else{
    
      //the first 16 char of the line are garbage
      int i = 0;
      while ( (i < 15) && (!picFileStream.eof()) ){
	picFileStream>>out_tmp;
	i++;
      }

      //now comes the text
      std::string current_line;
      int j=0;
      picFileStream>>out_tmp;
      while ( (((int)out_tmp) != 0) &&  (!picFileStream.eof()) ){
	current_line.push_back(out_tmp);
	j++;
	picFileStream>>out_tmp;
      }

      //save off the text to use as a comments in our NRRD header file
      comments += '#';
      comments.append(current_line);
      comments += '\n';

      //look for the Z factor line because we need that for NRRD spacing
      std::string::size_type z_factor_pos = current_line.find("Z_CORRECT_FACTOR=");
      if ( z_factor_pos != std::string::npos ){
	z_factor.assign(current_line,25,32);
      }

      //eat trailing NULLS (each line is a fixed 96? bytes long, padded with nulls)
      while ((out_tmp == 0) && (!picFileStream.eof())){
	picFileStream>>out_tmp;
      }
      
    }
  }

  //there is actuall more data beyond the end of the 2nd header, but
  //it looks like data on the settings for the aquisition machine.

  //just close the stream, we have everything for the nrrd header
  picFileStream.close();

  
  nrrdHeader<<"NRRD0001"<<std::endl
	    <<"type: uchar"<<std::endl
	    <<"dimension: 3"<<std::endl
	    <<"sizes: "<<x<<" "<<y<<" "<<z<<std::endl
	    <<"spacings: 1 1 "<<z_factor<<std::endl
	    <<"data file: "<<in<<std::endl
	    <<"endian: big"<<std::endl
	    <<"encoding: raw"<<std::endl
            <<"byte skip: "<<FIRST_HEADER_BYTESIZE<<std::endl;
  nrrdHeader<<"# Information contained in the pict header:"<<std::endl;
  nrrdHeader<<comments;
  nrrdHeader.close();

  return 0;
}

