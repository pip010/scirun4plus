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
//    File   : LabelTets.cc
//    Author : Jeroen Stinstra and Bill Martin
//    Date   : Wed Dec 2 2009

#include <Core/Algorithms/Converter/ConverterAlgo.h>
#include <Core/Algorithms/Fields/ClipMesh/ClipMeshBySelection.h>
#include <Core/Algorithms/Fields/DomainFields/MatchDomainLabels.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/NrrdData.h>
#include <Core/Persistent/Pstreams.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <iostream>
#include <fstream>

#include <stdlib.h>

using std::cerr;
using std::ifstream;
using std::endl;

using namespace SCIRun;

Point circumcirclecenter( VMesh::points_type& pts )
{
//  Vector ba = pts[1]-pts[0];
//  Vector ca = pts[2]-pts[0];
//  Vector da = pts[3]-pts[0];
//
//  return pts[0] +
//    ( da.length2()*Cross( ba, ca ) +
//      ca.length2()*Cross( da, ba ) +
//      ba.length2()*Cross( ca, da ) ) /
//    ( 2.0*Dot( ba, Cross(ca,da) ) );
	return (.25*(pts[0]+pts[1]+pts[2]+pts[3])).asPoint();
}

int
main(int argc, char *argv[]) {
  if (argc < 4) {
    cerr << "Usage: "<<argv[0]<<" nrrdLabelFile tetField labeledTetField\n";
    exit(1);
  }
  NrrdDataHandle nrrdH = new NrrdData;
  if (nrrdLoad(nrrdH->nrrd_, airStrdup(argv[1]), 0)) {
    char *err = biffGetDone(NRRD);
    cerr << "Error reading nrrd " << argv[1] << ". Error: "<< err << endl;
    free(err);
    biffDone(NRRD);
    exit(2);
  }

  ProgressReporter pr;
  SCIRunAlgo::ConverterAlgo algo(&pr);
  FieldHandle labelFieldH;
  // ADD OPTION FOR "Element"
  if (!algo.NrrdToField(nrrdH, labelFieldH, "Element", "Auto", "Do Not Correct")) {
    cerr << "Error converting nrrd to field.\n";
    exit(3);
  }

  FieldHandle tetFieldH;
  PiostreamPtr tetFieldStream = auto_istream(argv[2]);
  if (!tetFieldStream)
  {
    std::cerr << "Couldn't open file " << argv[2] << ".  Exiting..." << std::endl;
    exit(4);
  }
  Pio(*tetFieldStream, tetFieldH);
  if (!tetFieldH.get_rep())
  {
    std::cerr << "Error reading field from file " << argv[2] << ".  Exiting..."
              << std::endl;
    exit(5);
  }

  FieldHandle labelTetFieldH;
  SCIRunAlgo::MatchDomainLabelsAlgo label_algo;
  label_algo.run(tetFieldH, labelFieldH, labelTetFieldH);

//  FieldInformation fi(tetFieldH);
//  fi.set_basis_type(0);
//  fi.set_data_type("char");
//
//  FieldHandle labelTetFieldH=CreateField(fi, tetFieldH->mesh());
//  FieldHandle maskTetFieldH=CreateField(fi, tetFieldH->mesh());
//
//  VMesh *imesh = tetFieldH->vmesh();
//  VMesh *lmesh = labelFieldH->vmesh();
//  VField *lfield = labelFieldH->vfield();
//  VField *mfield = maskTetFieldH->vfield();
//  VField *ofield = labelTetFieldH->vfield();
//
//  VMesh::Node::array_type nodes;
//  VMesh::points_type points;
//
//  char remove_label = 0;
//
//  for (VMesh::Elem::index_type i=0; i<imesh->num_elems(); i++)
//  {
//	  imesh->get_nodes(nodes,i);
//	  imesh->get_centers(points,nodes);
//	  Point ccc = circumcirclecenter(points);
//	  char val;
//	  VMesh::Elem::index_type elem;
//	  lmesh->locate(elem,ccc);
//	  lfield->get_value(val,elem);
//	  ofield->set_value(val,i);
//
//	  if (val == remove_label)
//		  mfield->set_value(static_cast<char>(0.0),i);
//	  else
//		  mfield->set_value(static_cast<char>(1.0),i);
//  }

//  FieldHandle clippedLabelTetFieldH;
//  SCIRunAlgo::ClipMeshBySelectionAlgo calgo;
//  calgo.set_option("method","element");
//  bool success = calgo.run(labelTetFieldH,maskTetFieldH,clippedLabelTetFieldH);
//  if (!success)
//  {
//	  std::cerr << "Clip By Selection Failed.  Exiting..."
//	              << std::endl;
//	  exit(6);
//  }

  BinaryPiostream labelTetFieldStream(argv[3], Piostream::Write);
  Pio(labelTetFieldStream, labelTetFieldH);
  return 0;
}
