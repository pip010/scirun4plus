/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2013 Scientific Computing and Imaging Institute,
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

#include <gtest/gtest.h>

#include <boost/filesystem.hpp>

#include <iostream>
#include <string>

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Persistent/Pstreams.h>
#include <Core/Util/ProgressReporter.h>

#include <Core/Algorithms/Fields/Cleanup/TriSurfConsistencyCorrection.h>

using namespace ::testing;
using namespace SCIRun;

TEST(TriSurfConsistencyCorrectionTest, NullInputTest)
{
  FieldHandle input, output;
  MatrixHandle outputInvertedElementList;

  SCIRunAlgo::TriSurfConsistencyCorrectionAlgo algo_;

  // algorithm shouldn't succeed if input is null
  bool retval = algo_.run(input, output, outputInvertedElementList);
  EXPECT_EQ(retval, false);
}

TEST(TriSurfConsistencyCorrectionTest, NullInputTestWithSetParameter)
{
  FieldHandle input, output;
  MatrixHandle outputInvertedElementList;
  
  SCIRunAlgo::TriSurfConsistencyCorrectionAlgo algo_;
  algo_.set_bool("output_inverted_element_list", 0);
  
  // algorithm shouldn't succeed if input is null
  bool retval = algo_.run(input, output, outputInvertedElementList);
  EXPECT_EQ(retval, false);
}

TEST(TriSurfConsistencyCorrectionTest, InconsistentMeshTest)
{
  FieldHandle input, output;
  MatrixHandle outputInvertedElementList;

  // TODO: move field reading code to helper class or method...
  boost::filesystem::path inconsistentTriSurf = boost::filesystem::path(UNIT_TEST_DATA_ROOT) / "Fields" / "inconsistent-mesh.ts.fld";
  PiostreamPtr fileStream = auto_istream(inconsistentTriSurf.string());
  ASSERT_TRUE(fileStream != 0);
  
  Pio(*fileStream, input);
  ASSERT_TRUE( input.get_rep() != 0 );

  VMesh *inputVMesh = input->vmesh();
  ASSERT_TRUE( inputVMesh->is_trisurfmesh() );
  
  SCIRunAlgo::TriSurfConsistencyCorrectionAlgo algo_;
  algo_.set_bool("output_inverted_element_list", 1);

  bool retval = algo_.run(input, output, outputInvertedElementList);
  EXPECT_EQ(retval, true);
  
  std::string str;
  outputInvertedElementList->print(str);
  std::cerr << str;

  // compare meshes
  VMesh *outputVMesh = output->vmesh();
  ASSERT_TRUE( outputVMesh->is_trisurfmesh() );

  VMesh::dimension_type inputDims, outputDims;
  inputVMesh->get_dimensions(inputDims);
  outputVMesh->get_dimensions(outputDims);

  ASSERT_TRUE( inputVMesh->synchronize(Mesh::EDGES_E) );
  ASSERT_TRUE( outputVMesh->synchronize(Mesh::EDGES_E) );

  // TODO: replace with Google Mock container comparison...
  for(size_t i = 0; i < inputDims.size() && i < outputDims.size(); ++i)
  {
    EXPECT_EQ(inputDims[i], outputDims[i]);
  }

  EXPECT_EQ(inputVMesh->num_nodes(), outputVMesh->num_nodes());
  EXPECT_EQ(inputVMesh->num_edges(), outputVMesh->num_edges());
  EXPECT_EQ(inputVMesh->num_faces(), outputVMesh->num_faces());
}