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
#include <Core/Datatypes/NrrdData.h>
#include <Core/Geometry/Transform.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Point.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace SCIRun;

std::string input_nrrd_fname;
std::string output_nrrd_fname;
std::string output_matrix_fname;

bool has_input_nrrd = false;
bool has_output_nrrd = false;
bool has_output_matrix = false;

bool parseArgs(int argc, char *argv[]) 
{
  int curr=1;
  while (curr<argc-1) 
  {
    if (argv[curr] == std::string("-input")) 
    {
      input_nrrd_fname = argv[curr+1];
      has_input_nrrd = true;
      curr+=2;
    } 
    else if (argv[curr] == std::string("-output")) 
    {
      output_nrrd_fname = argv[curr+1];
      has_output_nrrd = true;
      curr+=2;
    } 
    else if (argv[curr] == std::string("-transform")) 
    {
      output_matrix_fname = argv[curr+1];
      has_output_matrix = true;
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
  std::cerr << "Usage: "<<appName <<" -input nrrd_file -output nrrd_file -transform transform_file\n";
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

  if (!has_input_nrrd || !has_output_nrrd || !has_output_matrix)
  {
    printUsage(argv[0]);
    exit(1);
  }

  std::cout << "Unorienting the Nrrd file so it will be axis aligned and the" <<
               " origin of the coordinate system is at (0,0,0)" << std::endl;
               

  // Step 1 -- Read the nrrd

  Nrrd* nrrd = nrrdNew();
  if ( nrrdLoad(nrrd, airStrdup(input_nrrd_fname.c_str()), 0) )
  {
    char *err = biffGetDone(NRRD);
    
    std::cerr << "Could not read Nrrd file" << std::endl;
    std::cerr << "Nrrd Error: " << err << std::endl;

    free(err);
    return (1);
  }
  
  unsigned int dim = nrrd->dim;
  if (dim != 3)
  {

    std::cerr << "Nrrd file needs to contain a 3D nrrd" << std::endl;
    std::cerr << "The input nrrd has a dimension of : " << dim << std::endl;
    nrrdNuke(nrrd);
    // abort
    return (1);
  }
  
 // Step 2 -- Determine the transformation matrix
 
  Transform tf;
  
  // This is the old way of defining a coordinate system in a NRRD
  std::vector<int>    size(3);
  std::vector<double> min(3);
  std::vector<double> max(3);

  // This is the new way of defining a coordinate system in a NRRD
  std::vector<Vector> SpaceDir(3);
  Vector Origin;
		
    
  // We assume the data is on the NODES
  // TODO: We need to fix this  
  
  for (size_t p=0; p<3; p++) 
  {
    // We get the size, this one is always defined
    size[p] = nrrd->axis[p].size;
    
    // min is not always defined. In case it is not it will be NaN
    if (airExists(nrrd->axis[p].min)) 
    {
      min[p] = nrrd->axis[p].min;
    }
    else
    {
      min[p] = 0.0;
    }
    
    // max is not always defined. In case it is not it will be NaN
    if (airExists(nrrd->axis[p].max)) 
    {
      max[p] = nrrd->axis[p].max;
    }
    else
    {
      if (airExists(nrrd->axis[p].spacing)) 
      {
        max[p] = nrrd->axis[p].spacing*(size[p]-1) + min[p];
      }
      else
      {
        max[p] = static_cast<double>(size[p]-1) + min[p];        
      }
    }
  }  
  
  // After this min and max for the three coordinates should be set
  
  // SpaceDim takes priority over min, max and spacing.....
  // Read the Nrrd documentation...
  
  if (nrrd->spaceDim > 0)
  {
    int sd = nrrd->spaceDim;
    if (sd != 3)
    {
      std::cerr << "Encountered an invalid nrrd file, the space dimension should" <<
       " match the nrrd dimension" << std::endl;
      
      // We fail, no need to clean up memory
      exit(-1);
    }

    for (size_t p=0; p< 3; p++)   
    {
      min[p] = 0.0;
      max[p] = static_cast<double>(size[p]-1);
      // Deal with voxel centered data
      if (nrrd->axis[p].center != nrrdCenterNode)
      {
        double cor = (max[p]-min[p])/(2*(size[p]-1));
        min[p] += cor;
        max[p] += cor;
        nrrd->axis[p].center = nrrdCenterNode;
      }
    } 
        
    // New way of defining origin

    if (airExists(nrrd->spaceOrigin[0])) Origin.x(nrrd->spaceOrigin[0]); else Origin.x(0.0);
    if (airExists(nrrd->spaceOrigin[1])) Origin.y(nrrd->spaceOrigin[1]); else Origin.y(0.0);
    if (airExists(nrrd->spaceOrigin[2])) Origin.z(nrrd->spaceOrigin[2]); else Origin.z(0.0);

    for (size_t p=0; p < 3;p++)
    {
      if (airExists(nrrd->axis[p].spaceDirection[0])) SpaceDir[p].x(nrrd->axis[p].spaceDirection[0]); else SpaceDir[p].x(0.0);
      if (airExists(nrrd->axis[p].spaceDirection[1])) SpaceDir[p].y(nrrd->axis[p].spaceDirection[1]); else SpaceDir[p].y(0.0);
      if (airExists(nrrd->axis[p].spaceDirection[2])) SpaceDir[p].z(nrrd->axis[p].spaceDirection[2]); else SpaceDir[p].z(0.0);
    }
  }
  else
  {
    Origin.x(min[0]);
    Origin.y(min[1]);
    Origin.z(min[2]);
    
    SpaceDir[0] = Vector((max[0]-min[0])/static_cast<double>(size[0]-1),0.0,0.0);
    SpaceDir[1] = Vector(0.0,(max[1]-min[1])/static_cast<double>(size[1]-1),0.0);
    SpaceDir[2] = Vector(0.0,0.0,(max[2]-min[2])/static_cast<double>(size[2]-1));
  }
  
  
  // Making sure that the coordinate system has a spacing close to 1.0 for the smallest
  // spacing. Due to arbitrary defined epsilons in the numeric code: things would go wrong
  // if spacing is small or really large. Compensate for that here, by redoing the scale issue
  // and setting the final transform to get us back to the original space dimensions.
  
  double minspace = SpaceDir[0].length();
  if (SpaceDir[1].length() < minspace) minspace = SpaceDir[1].length();
  if (SpaceDir[2].length() < minspace) minspace = SpaceDir[2].length();
  
  SpaceDir[0] = SpaceDir[0]*(1.0/minspace);
  SpaceDir[1] = SpaceDir[1]*(1.0/minspace);
  SpaceDir[2] = SpaceDir[2]*(1.0/minspace);
  
  // This defines the transform to the new coordinate system for the final field 
  tf.load_basis(Point(Origin),SpaceDir[0],SpaceDir[1],SpaceDir[2]);
  
// Correct for spacing
  Transform  spacing;
  spacing.pre_scale(Vector(minspace/(SpaceDir[0].length()),minspace/(SpaceDir[1].length()),minspace/(SpaceDir[2].length())));	
  tf.post_trans (spacing);
  
  
  // Clean up the nrrd file
  
  nrrd->spaceDim = 3;
  for (int p=0; p < 3; p++)
  {
    nrrd->axis[p].min     = static_cast<double>(AIR_NAN);
    nrrd->axis[p].max     = static_cast<double>(AIR_NAN);
    nrrd->axis[p].spacing = static_cast<double>(AIR_NAN);
    nrrd->axis[p].center = nrrdCenterNode;
  
    nrrd->axis[0].spaceDirection[p] = 0.0;
    nrrd->axis[1].spaceDirection[p] = 0.0;
    nrrd->axis[2].spaceDirection[p] = 0.0;
    
    nrrd->spaceOrigin[p] = 0.0;
  }

  nrrd->axis[0].spaceDirection[0] = SpaceDir[0].length();
  nrrd->axis[1].spaceDirection[1] = SpaceDir[1].length();
  nrrd->axis[2].spaceDirection[2] = SpaceDir[2].length();



  NrrdIoState *nio = nrrdIoStateNew();
  // set encoding to be gzip
  nio->encoding = nrrdEncodingArray[nrrdEncodingTypeGzip];
  // set format to be nrrd
  nio->format = nrrdFormatArray[nrrdFormatTypeNRRD];
  // set endian to be endian of machine
  nio->endian = AIR_ENDIAN;
  nio->zlibLevel = 6;
       
  if (nrrdSave( airStrdup(output_nrrd_fname.c_str()), nrrd,nio)) 
  {
    // Ugly...
    char *err = biffGetDone(NRRD);
    std::cerr << "Could not save nrrd '" << output_nrrd_fname << "' because teem crashed for the following reason: " << err << std::endl;
    free(err);
    
    return (1);
  }

  MatrixHandle matrix = new DenseMatrix(tf);
  
  ProgressReporter* pr = new ProgressReporter;
  PiostreamPtr stream = auto_ostream(output_matrix_fname, "Binary", pr);
  
  if (stream->error())
  {
    std::cerr << "Could not open file for writing " << output_matrix_fname << std::endl;
    return (1);
  }
  else
  {
    // Write the file
    Pio(*stream, matrix);
  } 
  
  nrrdNuke(nrrd);

  // Yeah it worked !!
  return (0);  
}    

