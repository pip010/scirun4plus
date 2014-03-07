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
//    File   : ComputeTightenedLabels.cc
//    Author : Bill Martin
//    Date   : Wed Dec 3 2009

#include <Core/Containers/Array3.h>
#include <Core/Datatypes/NrrdData.h>

#include <iostream>
#include <fstream>
#include <set>

#include <stdlib.h>

using std::cerr;
using std::ifstream;
using std::endl;
using std::set;

using namespace SCIRun;

class Label
{
public:
	Label(unsigned char id, float val): id_(id), value_(val) {}

	unsigned char id_;
	float value_;
};

// sorted so that dominant label is first.
bool operator<(const Label& l1, const Label& l2)
{
	return ((l1.value_ > l2.value_) ||
					((l1.value_ == l2.value_) && (l1.id_ < l2.id_)));
}

class LabelMap
{
public:
	LabelMap(size_t nx=0, size_t ny=0, size_t nz=0)
	{
		labels_.resize(nx,ny,nz);
	}
	void resize(size_t nx, size_t ny, size_t nz) { labels_.resize(nx,ny,nz); }
	std::set<Label> &operator()(size_t i, size_t j, size_t k)
	{
		return labels_(i,j,k);
	}
	size_t size(unsigned int axis)
	{
		if (axis == 0)
		{
			return labels_.dim1();
		} else if (axis == 1)
		{
			return labels_.dim2();
		} else if (axis == 2)
		{
			return labels_.dim3();
		}

		return 0;
	}
private:
	Array3< std::set<Label> > labels_;
};

int
main(int argc, char *argv[]) {
  if (argc < 2) {
    cerr << "Usage: "<<argv[0]<<" nrrdFileNames\n";
    exit(1);
  }

	LabelMap labels;
	char curr_label=0;
	
	// this is the output dominant label map.
	Nrrd* labelNrrd = nrrdNew();
	
	for (int m=1; m<argc; m++)
	{
		curr_label = m-1; // this is not technically ok -- should be
						  // the label assigned for the particular material
						  // not relative to the input order
		// open each material in turn
		NrrdDataHandle nrrdH = new NrrdData;
		if (nrrdLoad(nrrdH->nrrd_, airStrdup(argv[m]), 0)) {

			char *err = biffGetDone(NRRD);
			cerr << "Error reading nrrd " << argv[m] << ". Error: "<< err << endl;
			free(err);
			biffDone(NRRD);
			exit(2);
		}
		// for each voxel in this material, add it to the sorted list
		// of material contributions

		if (nrrdH->nrrd_->type != nrrdTypeFloat)
		{
			cerr << "Bad type " << argv[m] << "." << endl;
			exit(2);
		}

		// first time through, set the size of the volume to the
		// size of the first material
		if (m==1)
		{
				labels.resize(nrrdH->nrrd_->axis[0].size,
											nrrdH->nrrd_->axis[1].size,
											nrrdH->nrrd_->axis[2].size);
			
			// copy the header information from the input NRRDs to the output
			// label map
			labelNrrd->dim = 3;
			
			labelNrrd->type = nrrdTypeChar;
			labelNrrd->space = nrrdH->nrrd_->space;
			labelNrrd->spaceDim = nrrdH->nrrd_->spaceDim;
			
			for (int j=0; j<NRRD_SPACE_DIM_MAX; j++)
			{
				labelNrrd->spaceOrigin[j] = nrrdH->nrrd_->spaceOrigin[j];
			}
			
			for (int i=0; i<NRRD_DIM_MAX; i++)
			{
				labelNrrd->axis[i].size = nrrdH->nrrd_->axis[i].size;
				labelNrrd->axis[i].spacing = nrrdH->nrrd_->axis[i].spacing;
				labelNrrd->axis[i].min = nrrdH->nrrd_->axis[i].min;
				labelNrrd->axis[i].max = nrrdH->nrrd_->axis[i].max;
				labelNrrd->axis[i].center = nrrdH->nrrd_->axis[i].center;
				
				for (int j=0; j<NRRD_SPACE_DIM_MAX; j++)
				{
					labelNrrd->axis[i].spaceDirection[j] = nrrdH->nrrd_->axis[i].spaceDirection[j];
				}			
			}
			
			
		} else {
			if (nrrdH->nrrd_->axis[0].size != labels.size(0) ||
					nrrdH->nrrd_->axis[1].size != labels.size(1) ||
					nrrdH->nrrd_->axis[2].size != labels.size(2))
			{
				cerr << "Volumes must all be the same size!!" << endl;
				exit(2);
			}
		}

		printf("before label creation....\n");

		float* data = static_cast<float*>(nrrdH->nrrd_->data);
		for (size_t i=0; i<nrrdH->nrrd_->axis[0].size; i++)
		{
			for (size_t j=0; j<nrrdH->nrrd_->axis[1].size; j++)
			{
				for (size_t k=0; k<nrrdH->nrrd_->axis[2].size; k++)
				{
					// only store a label if its value is > -1.
					if (*data > -1.0f) {
						labels(i,j,k).insert(Label(curr_label,*data));
					}
					data++; // move to the next data value
				}
			}
		}
  }

	std::string base_dir_str = argv[1];

	std::string::size_type idx = base_dir_str.find_last_of ( "/\\" );
	if (idx == std::string::npos)
	{
		base_dir_str = "";
	} else {
		base_dir_str =  base_dir_str.substr(0,idx);
	}

	printf("before debugging\n");
	std::string labelmap_str = base_dir_str + "labelmap.txt";
	// BEGIN: DEBUG
	FILE* fp = fopen(labelmap_str.c_str(),"w");
	if (!fp) printf("Dummy, forgot to initialize this :)");
	for (size_t i=0; i<labels.size(0); i++)
	{
		for (size_t j=0; j<labels.size(1); j++)
		{
			for (size_t k=0; k<labels.size(2); k++)
			{
				fprintf(fp,"(%zu %zu %zu): ",i,j,k);
				std::set<Label>::iterator label_iter = labels(i,j,k).begin();
				while (label_iter !=  labels(i,j,k).end())
				{
					fprintf(fp,"(%d, %f) ",label_iter->id_,label_iter->value_);
					label_iter++;
				}
				fprintf(fp,"\n");
			}
		}
	}
	fclose(fp);

	printf("after debugging\n");

	// END: DEBUG

	NrrdIoState *nio = nrrdIoStateNew();
	nio->encoding = nrrdEncodingGzip;


	char *data;


	labelNrrd->data = data = new char[labels.size(0)*labels.size(1)*labels.size(2)];

	std::string label_file_str =  base_dir_str + "dominant_labelmap.nrrd";

	printf("before creating label map\n");

	for (size_t i=0; i<labels.size(0); i++)
	{
		for (size_t j=0; j<labels.size(1); j++)
		{
			for (size_t k=0; k<labels.size(2); k++)
			{
				std::set<Label>::iterator label_iter = labels(i,j,k).begin();
				// just write out the label who contribution is maximal at this point.
				if (label_iter != labels(i,j,k).end()) {
					*data = label_iter->id_;
				} else {
					*data = 0;
				}
				data++;
			}
		}
	}

	printf("before writing dominant label map\n");

	printf("label: %s\n",label_file_str.c_str());

	if (nrrdSave(label_file_str.c_str(),labelNrrd, nio)) {
		char *err = biffGetDone(NRRD);
		cerr << "Error writing nrrd " << label_file_str << ". Error: "<< err << endl;
		free(err);
		biffDone(NRRD);
		exit(2);
	}

	nrrdNuke(labelNrrd);

	printf("after dominant labelmap\n");

	// Now go through every material and adjust its tightened value
	// to account for the presence of other labels
	curr_label=0;
	for (int m=1; m<argc; m++)
	{
		curr_label = m-1; // technically this should be the label assigned to the particular material
		// open each material in turn
		NrrdDataHandle nrrdH = new NrrdData;
		if (nrrdLoad(nrrdH->nrrd_, airStrdup(argv[m]), 0)) {

			char *err = biffGetDone(NRRD);
			cerr << "Error reading nrrd " << argv[m] << ". Error: "<< err << endl;
			free(err);
			biffDone(NRRD);
			exit(2);
		}
		float* data = static_cast<float*>(nrrdH->nrrd_->data);
		for (size_t i=0; i<nrrdH->nrrd_->axis[0].size; i++)
		{
			for (size_t j=0; j<nrrdH->nrrd_->axis[1].size; j++)
			{
				for (size_t k=0; k<nrrdH->nrrd_->axis[2].size; k++)
				{
					std::set<Label>::iterator label_iter = labels(i,j,k).begin();

					// if there's a max and a second max, both not equal to -1
					if (labels(i,j,k).size() >= 2)
					{
						// if the current label is the max label, subtract the
						// next largest label
						if (curr_label == label_iter->id_)
						{
							*data -= (++label_iter)->value_;
						} else // otherwise, subtract the largest label from us
						{
							*data -= label_iter->value_;
						}
					} else if (labels(i,j,k).size() == 1) // only one value neq -1
					{
						if (m == label_iter->id_)
						{
							*data -= -1.0f;
						} else // otherwise, subtract the largest label from us
						{
							*data -= label_iter->value_;
						}
					} else { // this should never happen -- all voxels should have a label
						*data -= -1.0f;
					}
					data++; // move to the next data value
				}
			}
		}

		// now write the values back to a file
		std::string out_file_str = argv[m];
		std::string::size_type idx = out_file_str.find_last_of ( '.' );
		out_file_str =  out_file_str.substr(0,idx) + "-corrected" +
			out_file_str.substr(idx);

		if (nrrdSave(out_file_str.c_str(),nrrdH->nrrd_, nio)) {
			char *err = biffGetDone(NRRD);
			cerr << "Error writing nrrd " << out_file_str << ". Error: "<< err << endl;
			free(err);
			biffDone(NRRD);
			exit(2);
		}
	}
  return 0;
}
