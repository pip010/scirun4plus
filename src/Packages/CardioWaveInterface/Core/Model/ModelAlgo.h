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


#ifndef PACKAGES_CARDIOWAVE_CORE_MODEL_MODELALGO_H
#define PACKAGES_CARDIOWAVE_CORE_MODEL_MODELALGO_H 1

#include <Core/Algorithms/Util/AlgoLibrary.h>

#include <Core/Datatypes/Bundle.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/DenseMatrix.h>

#include <Dataflow/Network/Module.h>
#include <Packages/CardioWaveInterface/Core/Model/BuildMembraneTable.h>
#include <Packages/CardioWaveInterface/Core/Model/BuildStimulusTable.h>


#include <string>
#include <sstream>


namespace CardioWaveInterface {

using namespace SCIRun;

class ModelAlgo : public AlgoLibrary {

  public:
    ModelAlgo(ProgressReporter* pr); // normal case

    bool DMDOptimizeByCoordinate(std::string optim_coord, FieldHandle ElementType, std::vector<FieldHandle>& Membranes, std::vector<MembraneTable>& MembraneTables, MatrixHandle CompToGeom, MatrixHandle& Reorder);
    bool DMDBuildMembraneTable(FieldHandle ElementType, FieldHandle MembraneModel, MatrixHandle CompToGeom, MatrixHandle NodeLink, MatrixHandle ElemLink, MembraneTable& Table, MatrixHandle& MappingMatrix, MatrixHandle& MappingMatrix2, Matrix::index_type &offset);
    bool DMDBuildMembraneMatrix(std::vector<MembraneTable>& membranetable, std::vector<double>& nodetypes, std::vector<double>& cmvalues, size_type num_volumenodes, size_type num_synnodes, MatrixHandle& NodeType, MatrixHandle& Volume, MatrixHandle& Capacitance, MatrixHandle& MembaneMatrix);
    bool DMDMembraneTableToMatrix(MembraneTable MemTable, MatrixHandle& MemMatrix);

    bool DMDBuildStimulusTable(FieldHandle ElementType, FieldHandle StimulusModel, MatrixHandle CompToGeom, double domainmin, double domainmax, StimulusTable& Table);
    bool DMDBuildStimulusTableByElement(FieldHandle ElementType, FieldHandle StimulusModel, MatrixHandle CompToGeom, double domainmin, double domainmax, StimulusTable& Table);
    bool DMDStimulusTableToMatrix(StimulusTable StimTable, MatrixHandle& StimulusMatrix);

    bool DMDBuildReferenceTable(FieldHandle ElementType, FieldHandle ReferenceModel, MatrixHandle CompToGeom, double domainmin, double domainmax, ReferenceTable& Table);
    bool DMDBuildReferenceTableByElement(FieldHandle ElementType, FieldHandle ReferenceModel, MatrixHandle CompToGeom, double domainmin, double domainmax, ReferenceTable& Table);
    bool DMDReferenceTableToMatrix(ReferenceTable StimTable, MatrixHandle& ReferenceMatrix);

    bool DMDBuildSimulation(BundleHandle SimulationBundle, StringHandle FileName, BundleHandle& VisualizationBundle, StringHandle& Script);
};

}

#endif
