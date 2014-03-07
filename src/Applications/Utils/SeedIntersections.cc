/*
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

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <Core/Geometry/Point.h>
#include <Core/Datatypes/NrrdData.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Algorithms/Converter/NrrdToField.h>
#include <Core/Geometry/Transform.h>
#include <Core/Algorithms/Fields/TransformMesh/TransformMeshWithTransform.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>


using namespace SCIRun;

typedef std::map<std::string, std::vector<Point> > seed_map_t;

void
do_seedQ(std::vector<Point> &pnts, size_t xi, size_t yi, size_t zi, Transform T)
{
  const float f = .75;
  for (float z = -0.5; z <= 0.5; z+=f) 
	{
    float zs = z + (float)zi;
    for (float y = -0.5; y <= 0.5; y+=f) 
		{
      float ys = y + (float)yi;
      for (float x = -0.5; x <= 0.5; x+=f)
			{
        float xs = x + (float)xi;
        //ostringstream str;
        //str << xs << " " << ys << " " << zs << " 0.0 0.0 0.0 0.0";
        pnts.push_back(T.project(Point(xs, ys, zs)));
        //pnts.push_back((Point(xs, ys, zs)));
      }
    }
  }
}

void
do_seed(std::vector<Point> &pnts, size_t xi, size_t yi, size_t zi, Transform T)
{
  pnts.push_back(T.project(Point(xi, yi, zi)));
  //pnts.push_back((Point(xi, yi, zi)));
}

void
seed_dbls(std::vector<NrrdDataHandle> &nrrds, seed_map_t &smap, Transform T)
{
  // Expect all nrrds to have the same size
  NrrdDataHandle nh0 = nrrds[0];
  unsigned int x_sz = nh0->nrrd_->axis[0].size;
  unsigned int y_sz = nh0->nrrd_->axis[1].size;
  unsigned int z_sz = nh0->nrrd_->axis[2].size;
  unsigned int n = nrrds.size();
  const int skip = 1;
  for (size_t i = 0; i < n; ++i) 
	{
    float *idat = (float*)nrrds[i]->nrrd_->data;
    for (size_t j = i+1; j < n; ++j) 
		{
      float *jdat = (float*)nrrds[j]->nrrd_->data;
      std::ostringstream nm;
      nm << "d" << "_" << i << "_" << j;
      smap[nm.str()] = std::vector<Point>();

      for (size_t s = 0; s < z_sz; s+=skip) 
			{
        for (size_t r = 0; r < y_sz; r+=skip)
				{
          for (size_t p = 0; p < x_sz; p+=skip) 
					{
            int s1 = s; if (s1 > 0) s1 = s - rand() % skip;
            int r1 = r; if (r1 > 0) r1 = r - rand() % skip;
            int p1 = p; if (p1 > 0) p1 = p - rand() % skip;
            size_t off = s1 * x_sz * y_sz + r1 * x_sz + p1;
            if (idat[off] > 0 && jdat[off] > 0) 
            {
              // possible double junction here
              do_seed(smap[nm.str()], p1, r1, s1, T);
            }
          }
        }
      }
    }
  }
}

void
seed_trips(std::vector<NrrdDataHandle> &nrrds, seed_map_t &smap, Transform T)
{
  // Expect all nrrds to have the same size
  NrrdDataHandle nh0 = nrrds[0];
  unsigned int x_sz = nh0->nrrd_->axis[0].size;
  unsigned int y_sz = nh0->nrrd_->axis[1].size;
  unsigned int z_sz = nh0->nrrd_->axis[2].size;
  unsigned int n = nrrds.size();
  const int skip = 1;
  for (size_t i = 0; i < n; ++i)
	{
    float *idat = (float*)nrrds[i]->nrrd_->data;
    for (size_t j = i+1; j < n; ++j) 
		{
      float *jdat = (float*)nrrds[j]->nrrd_->data;
      for (size_t k = j+1; k < n; ++k) 
			{
        float *kdat = (float*)nrrds[k]->nrrd_->data;
        std::ostringstream nm;
        nm << "t" << "_" << i << "_" << j << "_" << k;
        smap[nm.str()] = std::vector<Point>();
       
        for (size_t s = 0; s < z_sz; s+=skip) 
				{
          for (size_t r = 0; r < y_sz; r+=skip) 
					{
            for (size_t p = 0; p < x_sz; p+=skip) 
						{
              int s1 = s; if (s1 > 0) s1 = s - rand() % skip;
              int r1 = r; if (r1 > 0) r1 = r - rand() % skip;
              int p1 = p; if (p1 > 0) p1 = p - rand() % skip;
              size_t off = s1 * x_sz * y_sz + r1 * x_sz + p1;
              if (idat[off] > 0 && jdat[off] > 0 && 
                  kdat[off] > 0) 
              {
                // possible triple junction here
                do_seed(smap[nm.str()], p1, r1, s1, T);
              }
            }
          }
        }
      }
    }
  }
}


void
seed_quads(std::vector<NrrdDataHandle> &nrrds, seed_map_t &smap, Transform T)
{
  // Expect all nrrds to have the same size
  NrrdDataHandle nh0 = nrrds[0];
  unsigned int x_sz = nh0->nrrd_->axis[0].size;
  unsigned int y_sz = nh0->nrrd_->axis[1].size;
  unsigned int z_sz = nh0->nrrd_->axis[2].size;
  unsigned int n = nrrds.size();
  for (size_t i = 0; i < n; ++i)
	{
    float *idat = (float*)nrrds[i]->nrrd_->data;
    for (size_t j = i+1; j < n; ++j) 
		{
      float *jdat = (float*)nrrds[j]->nrrd_->data;
      for (size_t k = j+1; k < n; ++k) 
			{
        float *kdat = (float*)nrrds[k]->nrrd_->data;
        for (size_t l = k+1; l < n; ++l) 
				{
          float *ldat = (float*)nrrds[l]->nrrd_->data;
          std::ostringstream nm;
          nm << "q" << "_" << i << "_" << j << "_" << k << "_" << l;
          smap[nm.str()] = std::vector<Point>();
          
          for (size_t s = 0; s < z_sz; ++s) 
					{
            for (size_t r = 0; r < y_sz; ++r) 
						{
              for (size_t p = 0; p < x_sz; ++p) 
							{
                size_t off = s * x_sz * y_sz + r * x_sz + p;
                if (idat[off] > 0 && jdat[off] > 0 && 
                    kdat[off] > 0 && ldat[off] > 0) 
								{
                  // possible quad junction here
                  do_seedQ(smap[nm.str()], p, r, s, T);
                }
              }
            }
          }
        }
      }
    }
  }
}

int main(int argc, char *argv[])
{
  if (argc < 3)
	{
    std::cerr << "Usage: "<<argv[0]<<" ptclPath transform nrrdFile1 [nrrdFile2] [nrrdFile3] [...]\n";
    exit(1);
  }
  std::string ptcl_path(argv[1]);
  
  std::string input_matrix_fname;
  input_matrix_fname = argv[2];
  
  std::vector<NrrdDataHandle> nrrds; // get the names of these from the command line
  for (int i=3; i<argc; i++)
	{
    NrrdDataHandle nrrdH = new NrrdData;
    if (nrrdLoad(nrrdH->nrrd_, airStrdup(argv[i]), 0))
		{
      char *err = biffGetDone(NRRD);      
      std::cerr << "Error reading nrrd " << argv[i] << ". Error: "<< err << std::endl;
      free(err);
      biffDone(NRRD);
      exit(2);
    }
    if (nrrdH->nrrd_->type != nrrdTypeFloat) {
      std::cerr << "Error -- nrrds have to be of type 'float'\n";
      exit(3);
    }
    if (nrrdH->nrrd_->dim != 3) {
      std::cerr << "Error -- nrrds have to be 3D\n";
      exit(4);
    }
    if (nrrds.size() > 0) 
		{
      if (nrrdH->nrrd_->axis[0].size != nrrds[0]->nrrd_->axis[0].size ||
					nrrdH->nrrd_->axis[1].size != nrrds[0]->nrrd_->axis[1].size ||
					nrrdH->nrrd_->axis[2].size != nrrds[0]->nrrd_->axis[2].size) 
			{ 
				std::cerr << "Error -- nrrds all have to have the same axis sizes.\n";
				exit(5);
      }
    }
				
    nrrds.push_back(nrrdH);
  }
	
	
	NrrdDataHandle nrrdH=nrrds[1];
	
	std::vector<Vector> SpaceDir(3);
	for (size_t p=0; p < 3;p++)
	{
		if (airExists(nrrdH->nrrd_->axis[p].spaceDirection[0])) SpaceDir[p].x(nrrdH->nrrd_->axis[p].spaceDirection[0]); else SpaceDir[p].x(0.0);
		if (airExists(nrrdH->nrrd_->axis[p].spaceDirection[1])) SpaceDir[p].y(nrrdH->nrrd_->axis[p].spaceDirection[1]); else SpaceDir[p].y(0.0);
		if (airExists(nrrdH->nrrd_->axis[p].spaceDirection[2])) SpaceDir[p].z(nrrdH->nrrd_->axis[p].spaceDirection[2]); else SpaceDir[p].z(0.0);
	}

	if (nrrds.size() == 0) 
	{
    std::cerr << "Error: must have at least one nrrd input.";
    exit(6);
  }
	


	SCIRun::Vector spacing;
	
	spacing[0]=(SpaceDir[0].length());
	spacing[1]=(SpaceDir[1].length());
	spacing[2]=(SpaceDir[2].length());
	
	
	Transform tftrans, trans;

	tftrans.load_identity();
	tftrans.pre_scale(spacing);
	
	
	
  NrrdDataHandle nh0 = nrrds[0];
  seed_map_t seeds;
  seed_quads(nrrds, seeds,tftrans);
  seed_trips(nrrds, seeds,tftrans);
  seed_dbls(nrrds, seeds,tftrans);
	
	

  seed_map_t::iterator iter = seeds.begin();
  while (iter != seeds.end())
	{
    const std::string &key = (*iter).first;
    //cerr << key << endl;
    std::string fname(ptcl_path + "/" + key + "_seed.ptcl");
    std::vector<Point> &points = (*iter).second;
    if (points.size())
		{
      std::ofstream pts(fname.c_str());
      std::vector<Point>::iterator piter = points.begin();
      pts << points.size() << std::endl;
      while (piter != points.end()) 
			{
				Point &p = *piter++;
				pts << p.x() << " " << p.y() << " " << p.z() << " 0.0 0.0 0.0 0.0" << std::endl;
      }
      pts.close();
    }
    ++iter;
  }
  return 0;
}
