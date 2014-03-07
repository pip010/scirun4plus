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

#include <Packages/CardioWaveInterface/Core/Model/ModelAlgo.h>
#include <Packages/CardioWaveInterface/Core/Model/BuildMembraneTable.h>
#include <Packages/CardioWaveInterface/Core/Model/BuildStimulusTable.h>
#include <Core/Algorithms/DataIO/DataIOAlgo.h>
#include <Core/Algorithms/Converter/ConverterAlgo.h>
#include <Core/Algorithms/Regression/RegressionAlgo.h>
#include <Core/Datatypes/MatrixOperations.h>

#include <Core/Algorithms/Math/MathAlgo.h>

#include <Core/Algorithms/Math/SortMatrix/SortMatrix.h>
#include <Core/Algorithms/Math/ReorderMatrix/ReorderMatrix.h>
#include <Core/Algorithms/Math/SelectMatrix/SelectMatrix.h>
#include <Core/Algorithms/Math/FindMatrix/FindMatrix.h>
#include <Core/Algorithms/Math/AppendMatrix/AppendMatrix.h>
#include <Core/Algorithms/Math/MappingMatrix/ConvertMappingMatrixIntoMappingOrder.h>
#include <Core/Algorithms/Math/MappingMatrix/ConvertMappingOrderIntoMappingMatrix.h>
#include <Core/Algorithms/FiniteElements/Mapping/BuildFEGridMapping.h>
#include <Core/Algorithms/FiniteElements/BuildMatrix/BuildFEMatrix.h>
#include <Core/Algorithms/Fields/FindNodes/FindClosestNode.h>
#include <Core/Algorithms/Fields/FindNodes/FindClosestNodeByValue.h>
#include <Core/Algorithms/Fields/MeshData/GetMeshNodes.h>
#include <Core/Algorithms/Fields/FieldData/GetFieldData.h>
#include <Core/Algorithms/Fields/Mapping/MapFieldDataFromElemToNode.h>

#include <fstream>

namespace CardioWaveInterface {

using namespace SCIRun;
using namespace SCIRunAlgo;


ModelAlgo::ModelAlgo(ProgressReporter* pr) :
  AlgoLibrary(pr)
{
}

bool
ModelAlgo::DMDBuildMembraneTable( FieldHandle     ElementType,
                                  FieldHandle     MembraneModel,
                                  MatrixHandle    CompToGeom,
                                  MatrixHandle    NodeLink,
                                  MatrixHandle    ElemLink,
                                  MembraneTable & MembraneTable,
                                  MatrixHandle &  MappingMatrix,
                                  MatrixHandle &  MappingMatrix2,
                                  Matrix::size_type & offset)
{
  BuildMembraneTableAlgo algo;
  return(algo.BuildMembraneTable(pr_,ElementType,MembraneModel,CompToGeom,NodeLink, ElemLink,MembraneTable,MappingMatrix,MappingMatrix2,offset));
}


bool ModelAlgo::DMDOptimizeByCoordinate(std::string optim_coord, FieldHandle ElementType, std::vector<FieldHandle>& Membranes, std::vector<MembraneTable>& MembraneTables, MatrixHandle CompToGeom, MatrixHandle& Reorder)
{
  remark("DMDBuildSimulation: Optimize by coordinate");
  DataIOAlgo dalgo(pr_);
 
  SCIRunAlgo::ConvertMappingMatrixIntoMappingOrderAlgo mapping_algo_;
  mapping_algo_.set_progress_reporter(pr_);
  SCIRunAlgo::ConvertMappingOrderIntoMappingMatrixAlgo imapping_algo_;
  imapping_algo_.set_progress_reporter(pr_);  
  SCIRunAlgo::SelectMatrixAlgo select_algo_;
  select_algo_.set_progress_reporter(pr_);
  SCIRunAlgo::AppendMatrixAlgo append_algo_;
  append_algo_.set_progress_reporter(pr_);
  SCIRunAlgo::SortMatrixAlgo sort_algo_;
  sort_algo_.set_progress_reporter(pr_); 
  SCIRunAlgo::GetMeshNodesAlgo nodes_algo_;
  nodes_algo_.set_progress_reporter(pr_); 
   
  MatrixHandle Coords;
  
  if(!(nodes_algo_.run(ElementType,Coords))) return (false);

  if (CompToGeom.get_rep())
  {
    std::vector<Matrix::index_type> order;
    if(!(mapping_algo_.run(CompToGeom,order))) return (false);
    select_algo_.set_option("method","select_rows");
    if(!(select_algo_.run(Coords,Coords,order))) return (false);
  }
  
  for (size_t p=0; p<Membranes.size();p++)
  {
    MatrixHandle MCoords;
    if(!(nodes_algo_.run(Membranes[p],MCoords))) return (false);
    
    std::vector<Matrix::index_type> nodesused;
    size_t memtable_size = MembraneTables[p].size();
    nodesused.reserve(memtable_size);
    for (size_t q=0; q<memtable_size; q++)
    {
      nodesused.push_back(MembraneTables[p][q].node0);
    }
    
    select_algo_.set_option("method","select_rows");
    if(!(select_algo_.run(MCoords,MCoords,nodesused))) return (false);
    
    append_algo_.set_option("method","append_rows");
    if(!(append_algo_.run(Coords,Coords,MCoords))) return (false);
  }

  MatrixHandle S;
  /*
  std::vector<Matrix::index_type> sel(3);
  if (optim_coord == "x") { sel[0] = 0; sel[1] = 1; sel[2] = 2; }
  if (optim_coord == "y") { sel[0] = 1; sel[1] = 2; sel[2] = 0; }
  if (optim_coord == "z") { sel[0] = 2; sel[1] = 0; sel[2] = 1; }
*/
  std::vector<Matrix::index_type> sel(1); sel[0] = 2; 
  
  select_algo_.set_option("method","select_columns");
  if(!(select_algo_.run(Coords,Coords,sel))) return (false);

  sort_algo_.set_option("method","sort_rows");
  std::vector<Matrix::index_type> order;
  if(!(sort_algo_.run(Coords,order))) return (false);  

  if(!(imapping_algo_.run(order,Reorder))) return (false);
  
  remark("DMDBuildSimulation: Optimize by coordinate done");
  
  return (true);
}


bool ModelAlgo::DMDBuildMembraneMatrix(std::vector<MembraneTable>& membranetable, std::vector<double>& nodetypes, std::vector<double>& cmvalues, Matrix::index_type num_volumenodes, Matrix::index_type num_synnodes, MatrixHandle& NodeType, MatrixHandle& Volume, MatrixHandle& Capacitance,MatrixHandle& MembraneMatrix)
{
  SCIRunAlgo::MathAlgo mathalgo(pr_);
  Matrix::index_type num_totalnodes = num_volumenodes + num_synnodes;
  
  // Build a vector for the surface areas
  Volume = dynamic_cast<Matrix*>(new DenseMatrix(num_totalnodes,1));
  if (Volume.get_rep() == 0)
  {
    error("DMDBuildSimulation: Could not allocate Volume Vector Matrix");
    return (false);
  }
  double* volumeptr = Volume->get_data_pointer();


  Capacitance = dynamic_cast<Matrix*>(new DenseMatrix(num_totalnodes,1));
  if ( Capacitance.get_rep() == 0)
  {
    error("DMDBuildSimulation: Could not allocate Capacitance Vector Matrix");
    return (false);
  }
  double* capacitanceptr = Capacitance->get_data_pointer();
  
  // Build a vector for the node types
  NodeType = dynamic_cast<Matrix*>(new DenseMatrix(num_totalnodes,1));
  if (NodeType.get_rep() == 0)
  {
    error("DMDBuildSimulation: Could not allocate Volume Vector Matrix");
    return (false);
  }
  double* nodetypeptr = NodeType->get_data_pointer();
  
  for (Matrix::index_type p=0; p<num_totalnodes;p++) nodetypeptr[p] = 0.0;
  for (Matrix::index_type p=0; p<num_totalnodes;p++) volumeptr[p] = 0.0;
  for (Matrix::index_type p=0; p<num_totalnodes;p++) capacitanceptr[p] = 0.0;
  
  // Build the Membrane connections

  Matrix::size_type synnum = num_volumenodes;
  SparseElementVector sev(num_synnodes*7);
  Matrix::index_type k = 0;
  for (size_t p=0; p<membranetable.size();p++)
  {
    for (size_t q=0; q< membranetable[p].size(); q++)
    { 
      if (membranetable[p][q].snode != synnum)
        std::cerr << "There is an issue in accounting for the row numbers\n";
      sev[k].row = synnum;
      sev[k].col = synnum;
      sev[k].val = 1.0;
      volumeptr[synnum] = membranetable[p][q].surface;
      capacitanceptr[synnum] = cmvalues[p];
      nodetypeptr[synnum] = nodetypes[p];
      k++;  
      sev[k].row = membranetable[p][q].node1;
      sev[k].col = membranetable[p][q].node2;
      sev[k].val = 1000.0;
      k++;    
      sev[k].row = membranetable[p][q].node2;
      sev[k].col = membranetable[p][q].node1;
      sev[k].val = 1000.0;                                                                 
      k++; 
      sev[k].row = membranetable[p][q].node2;
      sev[k].col = synnum;
      sev[k].val = 1000.0;
      k++; 
      sev[k].row = synnum;
      sev[k].col = membranetable[p][q].node2;
      sev[k].val = 1000.0;
      k++; 
      sev[k].row = membranetable[p][q].node1;
      sev[k].col = synnum;
      sev[k].val = 1000.0;
      k++; 
      sev[k].row = synnum;
      sev[k].col = membranetable[p][q].node1;
      sev[k].val = 1000.0;
      k++; 
      synnum++;
   }    
  }

  if(!(mathalgo.CreateSparseMatrix(sev, MembraneMatrix, num_totalnodes, num_totalnodes)))
  {
    error("DMDBuildSimulation: Could not build synapse sparse matrix");
    return (false);
  }    

  return (true);
}

bool ModelAlgo::DMDBuildStimulusTable(FieldHandle ElementType, FieldHandle StimulusModel, MatrixHandle CompToGeom,  double stimulusdomainmin, double stimulusdomainmax, StimulusTable& StimulusTable)
{
  BuildStimulusTableAlgo algo;
  return(algo.BuildStimulusTable(pr_,ElementType,StimulusModel,CompToGeom,stimulusdomainmin,stimulusdomainmax,true,StimulusTable));
}

bool ModelAlgo::DMDBuildStimulusTableByElement(FieldHandle ElementType, FieldHandle StimulusModel, MatrixHandle CompToGeom,  double stimulusdomainmin, double stimulusdomainmax, StimulusTable& StimulusTable)
{
  BuildStimulusTableAlgo algo;
  return(algo.BuildStimulusTable(pr_,ElementType,StimulusModel,CompToGeom,stimulusdomainmin,stimulusdomainmax,false,StimulusTable));
}

bool ModelAlgo::DMDBuildReferenceTableByElement(FieldHandle ElementType, FieldHandle ReferenceModel, MatrixHandle CompToGeom,  double referencedomainmin, double referencedomainmax, ReferenceTable& ReferenceTable)
{
  BuildStimulusTableAlgo algo;
  return(algo.BuildStimulusTable(pr_,ElementType,ReferenceModel,CompToGeom,referencedomainmin,referencedomainmax,false,ReferenceTable));
}

bool ModelAlgo::DMDBuildReferenceTable(FieldHandle ElementType, FieldHandle ReferenceModel, MatrixHandle CompToGeom,  double referencedomainmin, double referencedomainmax, ReferenceTable& ReferenceTable)
{
  BuildStimulusTableAlgo algo;
  return(algo.BuildStimulusTable(pr_,ElementType,ReferenceModel,CompToGeom,referencedomainmin,referencedomainmax,true,ReferenceTable));
}

bool ModelAlgo::DMDReferenceTableToMatrix(ReferenceTable ReferenceTable,MatrixHandle& M)
{
  M = dynamic_cast<Matrix *>(new DenseMatrix(ReferenceTable.size(),2));
  if (M.get_rep() == 0) return(false);
  double *dataptr = M->get_data_pointer();
  
  Matrix::index_type p =0;
  for (size_t k=0; k <ReferenceTable.size(); k++)
  {
    dataptr[p++] = static_cast<double>(ReferenceTable[k].node);
    dataptr[p++] = static_cast<double>(ReferenceTable[k].weight);
  }
  
  return (true);
}


bool ModelAlgo::DMDStimulusTableToMatrix(StimulusTable StimulusTable,MatrixHandle& M)
{
  M = dynamic_cast<Matrix *>(new DenseMatrix(StimulusTable.size(),2));
  if (M.get_rep() == 0) return(false);
  double *dataptr = M->get_data_pointer();
  
  Matrix::index_type p =0;
  for (size_t k=0; k <StimulusTable.size(); k++)
  {
    dataptr[p++] = static_cast<double>(StimulusTable[k].node);
    dataptr[p++] = static_cast<double>(StimulusTable[k].weight);
  }
  
  return (true);
}


bool ModelAlgo::DMDMembraneTableToMatrix(MembraneTable MembraneTable,MatrixHandle& M)
{
  M = dynamic_cast<Matrix *>(new DenseMatrix(MembraneTable.size(),4));
  if (M.get_rep() == 0) return(false);
  double *dataptr = M->get_data_pointer();
  
  Matrix::index_type p =0;
  for (size_t k=0; k <MembraneTable.size(); k++)
  {
    dataptr[p++] = static_cast<double>(MembraneTable[k].snode);
    dataptr[p++] = static_cast<double>(MembraneTable[k].node1);
    dataptr[p++] = static_cast<double>(MembraneTable[k].node2);
    dataptr[p++] = static_cast<double>(MembraneTable[k].surface);
  }
  
  return (true);
}

bool ModelAlgo::DMDBuildSimulation(BundleHandle SimulationBundle, StringHandle FileName, BundleHandle& VisualizationBundle, StringHandle& Script)
{
  // Define all the dynamic algorithms.
  // Forward the ProgressReporter so everything can forward an error

  SCIRunAlgo::MathAlgo        mathalgo(pr_);
  SCIRunAlgo::ConverterAlgo   converteralgo(pr_);
  SCIRunAlgo::DataIOAlgo      dataioalgo(pr_);
  SCIRunAlgo::RegressionAlgo  regalgo(pr_); // Need this for field comparisons
  
  // Step 0: Some sanity checks
  
  if (SimulationBundle.get_rep() == 0)
  {
    error("DMDBuildSimulation: SimulationBundle is empty");
    return (false);
  }

  // Check whether the filename was defined
  if (FileName.get_rep() == 0)
  {
    error("DMDBuildSimulation: FileName is empty");
    return (false);  
  }

  bool debug           = true;
  bool visbundle       = true;
  bool optimizesystemx = false;
  bool optimizesystemy = false;
  bool optimizesystemz = true;
  
  if (SimulationBundle->is_property("enable_debug")) SimulationBundle->get_property("enable_debug",debug);
  if (SimulationBundle->is_property("build_visualization_bundle")) SimulationBundle->get_property("build_visualization_bundle",visbundle);
  if (SimulationBundle->is_property("optimize_systemx")) SimulationBundle->get_property("optimize_systemx",optimizesystemx);
  if (SimulationBundle->is_property("optimize_systemy")) SimulationBundle->get_property("optimize_systemy",optimizesystemy);
  if (SimulationBundle->is_property("optimize_systemz")) SimulationBundle->get_property("optimize_systemz",optimizesystemz);
  
  // Display options in reporter
  std::ostringstream oss;
  oss << "debug="<<debug<< ", build_vis_bundle="<< visbundle << ", optimizex="<< optimizesystemx<<", optimizesystemy="<<optimizesystemy<<", optimizesystemz="<<optimizesystemz<<"\n";
  remark(oss.str());
  
  // Step 1: Define the filenames
  
  std::string filename_simbundle;   // Store the full configuration to desk so we can retrieve it for debugging or database purposes
  std::string filename_visbundle;   // Visualization Bundle, as we need to renumber the system we need a back projection all information is contained in this bundle
  std::string filename_sysmatrix;   // The system matrix
  std::string filename_mapping;     // The system matrix
  std::string filename_nodetype;    // A vector with the nodetype for each node in the system
  std::string filename_potential0;  // A vector with the domaintype for each node in the system
  std::string filename_surface;     // Surface factors for the cell membrane
  std::string filename_capacitance; // Capacitance factors for the cell membrane
  std::string filename_in;          // Parameter file
  std::string filename_script;      // The script to build simulation
  std::string filename_membrane;    // Membrane description (geometry)
  std::string filename_stimulus;    // Stimulus nodes in a file
  std::string filename_reference;   // Reference nodes in a file   
  std::string filename_electrode;   // Reference nodes in a file   
  std::string filename_output;      // name of output file   
  
  // Try to derive the filename
  // Remove the .sim.bdl extension of the simulation setup bundle file
  // First remove .bdl and then .sim 
  filename_simbundle  = FileName->get();
  std::string filenamebase = filename_simbundle;
  if ((filename_simbundle.size() > 4)&&(filename_simbundle.substr(filename_simbundle.size()-4)==std::string(".bdl")))
  {
    filenamebase = filename_simbundle.substr(0,filename_simbundle.size()-4);
  }
  if ((filenamebase.size() > 4)&&(filenamebase.substr(filenamebase.size()-4)==std::string(".sim")))
  {
    filenamebase = filenamebase.substr(0,filenamebase.size()-4);
  }
  
  // We now have the filenamebase and now define each of the new filenames
 
  filename_visbundle = filenamebase + ".vis.bdl";
  filename_sysmatrix = filenamebase + ".fem.spr";
  filename_mapping   = filenamebase + ".map.spr";
  filename_nodetype  = filenamebase + ".nt.bvec";
  filename_potential0= filenamebase + ".vm0.vec";
  filename_surface   = filenamebase + ".area.vec";
  filename_capacitance= filenamebase + ".cap.vec";
  filename_in        = filenamebase + ".in";
  filename_script    = filenamebase + ".script.sh";
  filename_membrane  = filenamebase + ".mem.txt";
  filename_stimulus  = filenamebase + ".stim.txt";
  filename_reference = filenamebase + ".ref.txt";
  filename_electrode = filenamebase + ".elec.txt";
  filename_output    = filenamebase + ".out";
   
  // Step 2: Save the model to file, we need to have an archived version
  // This will allow us to go back to the original data
  


  BundleHandle Domain = SimulationBundle->getBundle("Domain");
  if (Domain.get_rep() == 0)
  {
    error("DMDBuildSimulation: Domain is not defined");
    return (false);        
  }
  
	try
  {
    std::ofstream infile;
    infile.open(filename_in.c_str());

    StringHandle Parameters = SimulationBundle->getString("Parameters");
    if (Parameters.get_rep() == 0)
    {
      error("DMDBuildSimulation: Could not find Parameters String");
      return (false);     
    }
    
    infile << Parameters->get();
    infile << "\n";
    std::string rel_filename; 
    size_t pos;
    
    rel_filename = filename_nodetype; pos = 0; while (pos!=string::npos) { pos = rel_filename.find('/'); rel_filename = rel_filename.substr(pos+1); }
    infile << "nodefile=" << rel_filename << "\n";    
    if (Domain->isField("InitialPotential"))
    {
      rel_filename = filename_potential0; pos = 0; while (pos!=string::npos) { pos = rel_filename.find('/'); rel_filename = rel_filename.substr(pos+1); }
      infile << "vm0file=" << rel_filename<< "\n";
    }
    rel_filename = filename_membrane; pos = 0; while (pos!=string::npos) { pos = rel_filename.find('/'); rel_filename = rel_filename.substr(pos+1); }
    infile << "synapsefile=" << rel_filename << "\n";
    rel_filename = filename_sysmatrix; pos = 0; while (pos!=string::npos) { pos = rel_filename.find('/'); rel_filename = rel_filename.substr(pos+1); }
    infile << "grid_int=" << rel_filename << "\n";
    rel_filename = filename_surface; pos = 0; while (pos!=string::npos) { pos = rel_filename.find('/'); rel_filename = rel_filename.substr(pos+1); }
    infile << "grid_area=" << rel_filename << "\n";
    rel_filename = filename_reference; pos = 0; while (pos!=string::npos) { pos = rel_filename.find('/'); rel_filename = rel_filename.substr(pos+1); }
    infile << "reffile="<< rel_filename << "\n";
    rel_filename = filename_output; pos = 0; while (pos!=string::npos) { pos = rel_filename.find('/'); rel_filename = rel_filename.substr(pos+1); }
    infile << "outputfile=" << rel_filename << "\n";
    rel_filename = filename_stimulus; pos = 0; while (pos!=string::npos) { pos = rel_filename.find('/'); rel_filename = rel_filename.substr(pos+1); }
    infile << "stimfile=" << rel_filename << "\n";
    rel_filename = filename_electrode; pos = 0; while (pos!=string::npos) { pos = rel_filename.find('/'); rel_filename = rel_filename.substr(pos+1); }
    infile << "electrodefile="<< rel_filename << "\n";
    rel_filename = filename_capacitance; pos = 0; while (pos!=string::npos) { pos = rel_filename.find('/'); rel_filename = rel_filename.substr(pos+1); }
    infile << "cmfile=" << rel_filename << "\n";

    if (debug) infile << "debug=4" << "\n"; else infile << "debug=0" << "\n";
    
    remark("Created parameter file");
    
  }
  catch (...)
  {
    error("DMDBuildSimulation: Could not write input file");
    return (false);  
  }   

  try
  {
    std::ofstream scriptfile;
    scriptfile.open(filename_script.c_str());

    StringHandle SourceFile = SimulationBundle->getString("SourceFile");
    if (SourceFile.get_rep() == 0)
    {
      error("DMDBuildSimulation: Could not find SourceFile String");
      return (false);     
    }
    
    scriptfile << "perl nw_make.pl " << SourceFile->get() << "\n";    
    
    Script = new String("perl nw_make.pl "+ string(SourceFile->get()));

    remark("Created script file");
  }
  catch (...)
  {
    error("DMDBuildSimulation: Could not write script");
    return (false);  
  } 


  // Step 3: Dissemble bundle and get pointers to model organized

  // Get the pointers to the ElementType and Conductivity fields
  // This should not be a big memory overhead as we are just reorganizing data
  
  // Preferably ElementType and Conductivity are the same field with shorts on the
  // data to make it more memory efficient.
  

  
  FieldHandle InitialPotential = Domain->getField("InitialPotential");
  FieldHandle ElementType      = Domain->getField("ElementType");
  FieldHandle Conductivity     = Domain->getField("Conductivity");
  MatrixHandle ConductivityTable = Domain->getMatrix("ConductivityTable");
  MatrixHandle ElemLink = Domain->getMatrix("ElemLink");
  MatrixHandle NodeLink = Domain->getMatrix("NodeLink");

  
  // Check the validity of the matrices and fields

  if (ElementType.get_rep() == 0)
  {
    error("DMDBuildSimulation: No field called 'ElementType' was defined");
    return (false);
  }

  if (ElemLink.get_rep() == 0)
  {
    if (ElementType->is_property("ElemLink")) ElementType->get_property("ElemLink",ElemLink);
    if (ElementType->is_property("NodeLink")) ElementType->get_property("NodeLink",NodeLink);
  }

  if (Conductivity.get_rep() == 0)
  {
    error("DMDBuildSimulation: No field called 'Conductivity' was defined");
    return (false);
  }

  // Get all the information from the sub bundles
  // in the model

  std::vector<FieldHandle> Membranes; // Collect all the Membrane Geometries in one vector
  std::vector<double>      nodetypes; // Get the corresponding nodetype as well
  std::vector<double>      cmvalues; // Get the capacitance values

  std::vector<FieldHandle> References;      // Collect all the Geometries of the References in here
  std::vector<double>      referencevalues;
  std::vector<double>      referencedomainmin;
  std::vector<double>      referencedomainmax;
  std::vector<bool>        usereferencevalues;
  std::vector<bool>        referenceuseelements;

  std::vector<FieldHandle> Electrodes;
  std::vector<double>      electrodedomainmin;
  std::vector<double>      electrodedomainmax;
  std::vector<int>         electrodemembrane;

  std::vector<FieldHandle> Stimulus;
  std::vector<double>      stimulusdomainmin;
  std::vector<double>      stimulusdomainmax;
  std::vector<double>      stimuluscurrent;
  std::vector<double>      stimulusstart;
  std::vector<double>      stimulusend;
  std::vector<double>      stimulusduration;
  std::vector<bool>        stimulusiscurrentdensity;
  std::vector<bool>        stimulususeelements;

  size_t num_membranes = 0;
  size_t num_stimulus = 0;
  size_t num_references = 0;
  size_t num_electrodes = 0;

  // Loop through all sub bundles in the parent bundle to get all the individual
  // components of the model. Next to the volume each piece of membrane has its
  // own geometry. Similarly every electrode (stimulus and reference has its own geometry)
  
  for (size_t p=0; p < static_cast<size_t>(SimulationBundle->numBundles()); p++)
  {
    // Do we have a membrane definition
    std::string bundlename = SimulationBundle->getBundleName(p);
    if (bundlename.substr(0,9) == "Membrane_")
    {
      BundleHandle MembraneBundle = SimulationBundle->getBundle(bundlename);
      FieldHandle Geometry = MembraneBundle->getField("Geometry");
      MatrixHandle Capacitance = MembraneBundle->getMatrix("Capacitance");
      
      std::istringstream iss(bundlename.substr(9));
      double num;
      iss >> num;
      
      if (Geometry.get_rep())
      {
        // Move the information in the vectors
        // We use vectors here as we are only shifting pointers
        // It is just book keeping
        size_type numnodes, numelems;
        numnodes = Geometry->vmesh()->num_nodes();
        numelems = Geometry->vmesh()->num_elems();
        
        if (numnodes != 0 && numelems != 0)
        {
          Membranes.push_back(Geometry);
          nodetypes.push_back(num);
          num_membranes++;
        }
      }
      else
      {
        warning("DMDBuildSimulation: One of the membrane fields misses a Geometry or a NodeType");
      }
      
      if (Capacitance.get_rep())
      {
        cmvalues.push_back(Capacitance->get(0,0));  
      }
      else
      {
        cmvalues.push_back(1.0);
      }
    }

    if (bundlename.substr(0,10) == "Reference_")
    {
      BundleHandle ReferenceBundle =  SimulationBundle->getBundle(bundlename);
      FieldHandle Geometry = ReferenceBundle->getField("Geometry");
      MatrixHandle RefValue = ReferenceBundle->getMatrix("RefValue");
      MatrixHandle DomainMin = ReferenceBundle->getMatrix("DomainMin");      
      MatrixHandle DomainMax = ReferenceBundle->getMatrix("DomainMax");      
      MatrixHandle UseElements = ReferenceBundle->getMatrix("UseElements");      

      if (Geometry.get_rep()&&DomainMin.get_rep())
      {
        size_type numnodes, numelems;
        numnodes = Geometry->vmesh()->num_nodes();
        numelems = Geometry->vmesh()->num_elems();
        
        if (numnodes != 0)
        {

          References.push_back(Geometry);
          // If no reference value is given we assume one
          // If one is given we just need to convert it.
          if (RefValue.get_rep())
          {
            double refvalue = 0.0;
            converteralgo.MatrixToDouble(RefValue,refvalue);
            referencevalues.push_back(refvalue);
            usereferencevalues.push_back(true);
          }
          else
          {
            referencevalues.push_back(0.0);
            usereferencevalues.push_back(false);
          }
          
          if (UseElements.get_rep())
          {
            double ue;
            converteralgo.MatrixToDouble(UseElements,ue);
            if (ue) referenceuseelements.push_back(true); else referenceuseelements.push_back(false);
          }
          else referenceuseelements.push_back(false);
          
          double domain;
          converteralgo.MatrixToDouble(DomainMin,domain);
          referencedomainmin.push_back(domain);        
          converteralgo.MatrixToDouble(DomainMax,domain);
          referencedomainmax.push_back(domain);        
          num_references++;
        }
      }
      else
      {
        warning("DMDBuildSimulation: One of the reference fields misses a Geometry or a Domain");
      }
    }

    if (bundlename.substr(0,9) == "Stimulus_")
    {
      BundleHandle StimulusBundle = SimulationBundle->getBundle(bundlename);
      FieldHandle  Geometry = StimulusBundle->getField("Geometry");
      MatrixHandle DomainMin = StimulusBundle->getMatrix("DomainMin");
      MatrixHandle DomainMax = StimulusBundle->getMatrix("DomainMax");
      MatrixHandle Current = StimulusBundle->getMatrix("Current");
      MatrixHandle Start = StimulusBundle->getMatrix("Start");
      MatrixHandle End = StimulusBundle->getMatrix("End");
      MatrixHandle Duration = StimulusBundle->getMatrix("Duration");
      MatrixHandle CurrentDensity = StimulusBundle->getMatrix("CurrentDensity");
      MatrixHandle UseElements = StimulusBundle->getMatrix("UseElements");
      
      // Check whether we have all the information needed
      if ((Geometry.get_rep())&&(DomainMin.get_rep())&&(Start.get_rep())&&(End.get_rep()))
      {
        size_type numnodes, numelems;
        numnodes = Geometry->vmesh()->num_nodes();
        numelems = Geometry->vmesh()->num_elems();
        
        if (numnodes != 0)
        {
          double domainmin, domainmax;
          double start;
          double end;
          double current;
          double duration = -1.0;

          Stimulus.push_back(Geometry);
          converteralgo.MatrixToDouble(DomainMin,domainmin);
          converteralgo.MatrixToDouble(DomainMax,domainmax);
          converteralgo.MatrixToDouble(Start,start);
          converteralgo.MatrixToDouble(End,end);
          if (Duration.get_rep())
          {
            converteralgo.MatrixToDouble(Duration,duration);
          }
          stimulusdomainmin.push_back(domainmin);
          stimulusdomainmax.push_back(domainmax);
          stimulusstart.push_back(start);     
          stimulusend.push_back(end);
          stimulusduration.push_back(duration);
          num_stimulus++;

          if (UseElements.get_rep())
          {
            double ue;
            converteralgo.MatrixToDouble(UseElements,ue);
            if (ue) stimulususeelements.push_back(true); else stimulususeelements.push_back(false);
          }
          else stimulususeelements.push_back(false);
           
          if (Current.get_rep())
          {
            converteralgo.MatrixToDouble(Current,current);
            stimuluscurrent.push_back(current);
            stimulusiscurrentdensity.push_back(false);
          }
          else if (CurrentDensity.get_rep())
          {
            converteralgo.MatrixToDouble(CurrentDensity,current);
            stimuluscurrent.push_back(current);
            stimulusiscurrentdensity.push_back(true);
          }
          else
          {
            stimuluscurrent.push_back(0.0);
            stimulusiscurrentdensity.push_back(false);
          }
        }
      }
      else
      {
        warning("DMDBuildSDimulation: One of the stimulus fields misses a Geometry or a stimulation Domain, or a Starting time or a Ending time");
      }
    }
  }

  for (size_t p=0; p < static_cast<size_t>(SimulationBundle->numBundles()); p++)
  {
    // Do we have a membrane definition
    std::string bundlename = SimulationBundle->getBundleName(p);

    if (bundlename.substr(0,10) == "Electrode_")
    {
      BundleHandle ElectrodeBundle =  SimulationBundle->getBundle(bundlename);
      FieldHandle Geometry = ElectrodeBundle->getField("Electrodes");
      FieldHandle Membrane = ElectrodeBundle->getField("Membrane");
      MatrixHandle DomainMin = ElectrodeBundle->getMatrix("DomainMin");    
      MatrixHandle DomainMax = ElectrodeBundle->getMatrix("DomainMax");    
      if (Geometry.get_rep() && DomainMin.get_rep())
      {
        double domainmin, domainmax;
        converteralgo.MatrixToDouble(DomainMin,domainmin);
        converteralgo.MatrixToDouble(DomainMax,domainmax);
        electrodedomainmin.push_back(domainmin);
        electrodedomainmax.push_back(domainmax);
        electrodemembrane.push_back(-1);
        Electrodes.push_back(Geometry);
        num_electrodes++;
      }
      else if (Geometry.get_rep() && Membrane.get_rep())
      {
        Matrix::index_type membranenum = -1;
        for (size_t q=0; q< Membranes.size(); q++)
        {
          if (regalgo.CompareFields(Membrane,Membranes[q])) { membranenum = q; break; }
        }
        
        if (membranenum >= 0)
        {
          electrodedomainmin.push_back(0.0);
          electrodedomainmax.push_back(0.0);
          electrodemembrane.push_back(membranenum);
          Electrodes.push_back(Geometry);
          num_electrodes++;
        }
        else
        {
          warning("DMDBuildSimulation: Membrane Electrode is defined for a membrane not part of the model");        
        }
      }
      else
      {
        warning("DMDBuildSimulation: One of the electrode fields misses a Membrane or a Domain");
      }      
    }    
  }

  {
    std::ostringstream oss;
    oss << "Number of membranes="<<num_membranes<<" ,number of stimulus="<<num_stimulus<<" ,number of references="<<num_references << " ,number of electrodes="<< num_electrodes;
    remark(oss.str());
  }
  
  
  remark("Verified Simulation Bundle");


  // Step 4: Build the finite element system

  // If we are linking surfaces we have an additional memory overhead of the linkage
  // matrices and renumbering matrices.
  
  MatrixHandle GeomToComp;
  MatrixHandle CompToGeom;
  
  // Optional linkage of boundaries
  if (NodeLink.get_rep())
  {
    BuildFEGridMappingAlgo build_mapping_;
    build_mapping_.set_progress_reporter(pr_);
    MatrixHandle Dummy1, Dummy2;
    
    build_mapping_.set_bool("build_potential_geomtogrid",false);
    build_mapping_.set_bool("build_current_gridtogeom",false);
    if(!(build_mapping_.run(ElementType,NodeLink,Dummy1,CompToGeom,GeomToComp,Dummy2)))
    {
      error("DMDBuildSimulation: Could not build computational grid to geometrical mesh linkage matrices");
      return (false);  
    }
        
    Domain->setMatrix("GeomToComp",GeomToComp);    
    Domain->setMatrix("CompToGeom",CompToGeom);    
        
    remark("Created Computational grid");
  }

  if (!(dataioalgo.WriteBundle(std::string(filename_simbundle),SimulationBundle)))
  {
    error("DMDBuildSimulation: Could not save simulation bundle to disk");
    return (false);    
  }
  
  MatrixHandle GeomToComp2;
  
  if (GeomToComp.get_rep())
  {
    GeomToComp2 = GeomToComp;
    GeomToComp2.detach();
    SparseRowMatrix* spr = GeomToComp2->sparse();
  
    Matrix::size_type nrows = spr->nrows();
    Matrix::size_type nocls = spr->ncols();
    Matrix::index_type *rr = spr->rows;
    Matrix::index_type *cc = spr->columns;
    double *vv = spr->a;
       
    for (Matrix::index_type r=0; r< nrows; r++)
    {
      Matrix::index_type start = rr[r];
      Matrix::index_type end = rr[r+1];
      
      if (end > start)
      {
        double factor = 1.0/static_cast<double>(end-start);
        for (int c=start; c<end;c++)
        {
          vv[c] = factor*vv[c];
        }
      }
    }
    remark("Created GeomToComp2");
  }
  
  // We should have ordered all the geometric information now
  
  num_membranes = Membranes.size();
  num_stimulus  = Stimulus.size();
  num_references = References.size();
  num_electrodes = Electrodes.size();
  
  // To be able to map data back onto the geometry we will be needed mapping matrices
  // Hence define a lot of properties. Most of it we need for the model or the
  // visualization bundle
  std::vector<MembraneTable>     membranetable(num_membranes);
  std::vector<MatrixHandle>      membranemapping(num_membranes);
  std::vector<MatrixHandle>      membranemapping2(num_membranes);
  std::vector<size_type>  membranenumnodes(num_membranes);
  std::vector<size_type>  membranenumelems(num_membranes);
  
  // A synapse node contains the parameters of the membrane in CW they are separate nodes
  // They are not contained in the domain, they are just computational.
  Matrix::size_type num_synnodes = 0;

  
  // Build the tables that tell how the different components can be linked into the
  // system.
  // We have a MembraneTable to specify how the membranes are linked
  // We have a StimulusTable and a ReferenceTable

  // Get the information from the domain field
  size_type num_domainnodes, num_domainelems;
  num_domainnodes = ElementType->vmesh()->num_nodes();
  num_domainelems = ElementType->vmesh()->num_elems();

  // Get the offset at which we have to insert the synapse nodes
  index_type offset = num_domainnodes;
  if (CompToGeom.get_rep())
  {
    offset = CompToGeom->ncols();
  }
  
  // Fill out the membrane tables
  for (size_t p=0; p<membranetable.size();p++)
  {
    // Get the size
    // This function will dynamically compile if needed
    membranenumnodes[p] = Membranes[p]->vmesh()->num_nodes();
    membranenumelems[p] = Membranes[p]->vmesh()->num_elems();

    // This is the expensive function, it figures out how the membrane panels
    // are linked into the volumetric model. As we need not save this information
    // we are reconstructing it here. The reason for not keeping it is that this
    // way the user can clip and delete panels without having to remember each
    // change, as they all change the node numbering.
    if (!(DMDBuildMembraneTable(ElementType,Membranes[p],CompToGeom,NodeLink,ElemLink,membranetable[p],membranemapping[p],membranemapping2[p],offset)))
    {
      error("DMDBuildSimulation: Could build membrane model");
      return(false);
    }

    // Figure out how many synapse nodes are needed.
    num_synnodes += membranetable[p].size();
  }
  
  remark("Created Membrane tables");

  
  // Resize the mapping matrices to have the proper sizes of the full system.
  for (size_t p=0; p<membranetable.size();p++)
  {
    if (!(mathalgo.ResizeMatrix(membranemapping[p],membranemapping[p],membranemapping[p]->nrows(),offset)))
    {
      error("DMDBuildSimulation: Could not resize mapping matrix");
      return (false);        
    }  
    if (!(mathalgo.ResizeMatrix(membranemapping2[p],membranemapping2[p],membranemapping2[p]->nrows(),offset)))
    {
      error("DMDBuildSimulation: Could not resize mapping matrix");
      return (false);        
    }  
  }
  
  remark("Located nodes that form the membranes in volume and computed synapse nodes");

  
  std::vector<StimulusTable>     stimulustable(num_stimulus);
  
  for (size_t p=0; p<stimulustable.size();p++)
  {  
    if (stimulususeelements[p])
    {
      if(!(DMDBuildStimulusTableByElement(ElementType, Stimulus[p], CompToGeom,stimulusdomainmin[p], stimulusdomainmax[p],stimulustable[p])))
      {
        error("DMDBuildSimulation: Could build stimulus model");
        return(false);   
      }    
    }
    else
    {
      if(!(DMDBuildStimulusTable(ElementType, Stimulus[p], CompToGeom,stimulusdomainmin[p], stimulusdomainmax[p], stimulustable[p])))
      {
        error("DMDBuildSimulation: Could build stimulus model");
        return(false);   
      }
    }
  }
  
  remark("Located nodes that form the stimulus nodes");  
  
  std::vector<ReferenceTable>    referencetable(num_references);
  
  for (size_t p=0; p<referencetable.size();p++)
  { 
    if (referenceuseelements[p])
    { 
      if(!(DMDBuildReferenceTableByElement(ElementType, References[p], CompToGeom,referencedomainmin[p], referencedomainmax[p], referencetable[p])))
      {
        error("DMDBuildSimulation: Could build reference model");
        return(false);   
      }
    }
    else
    {
      if(!(DMDBuildReferenceTable(ElementType, References[p], CompToGeom,referencedomainmin[p],referencedomainmax[p], referencetable[p])))
      {
        error("DMDBuildSimulation: Could build reference model");
        return(false);   
      }
    }
  }
  
  remark("Located nodes that form the reference nodes");  

  // Step 5: Build the real Stiffness matrix
  
  // Build the Finite Element Model here
  // We enter optional fields here, if they are empty handles the routine will skip
  // them.
  MatrixHandle fematrix;
  MatrixHandle synmatrix;
  
  SCIRunAlgo::BuildFEMatrixAlgo femalgo;
  femalgo.set_progress_reporter(pr_);
  
  if(!(femalgo.run(Conductivity,ConductivityTable,fematrix)))
  {
    error("DMDBuildSimulation: Could not build FE Matrix");
    return(false);
  }
  
  fematrix = GeomToComp*fematrix*CompToGeom;
  
  remark("Created the stiffness matrix");
  
  // Figure out how many nodes we actually have
  size_type num_volumenodes = fematrix->nrows();
  size_type num_totalnodes = num_volumenodes + num_synnodes;

  // Resize the matrix so we can use it to store the synapse nodes as well
  if(!(mathalgo.ResizeMatrix(fematrix,fematrix,num_totalnodes,num_totalnodes)))
  {
    error("DMDBuildSimulation: Could not resize FE matrix");
    return (false);
  }

  remark("Resized the stiffness matrix");
  
  MatrixHandle VolumeVec;
  MatrixHandle CapacitanceVec;
  MatrixHandle NodeType;
  MatrixHandle Potential0;

  SCIRunAlgo::MapFieldDataFromElemToNodeAlgo map_algo;
  map_algo.set_progress_reporter(pr_);
  map_algo.set_option("method","max");
  
  if(!(map_algo.run(ElementType,ElementType)))
  {
    error("DMDBuildSimulation: Could not move elementtype to nodes");
    return (false);   
  }

  if (InitialPotential.get_rep())
  {
    // Get
    SCIRunAlgo::GetFieldDataAlgo get_algo;
    get_algo.set_progress_reporter(pr_);
    
    if(!(get_algo.run(InitialPotential,Potential0)))
    {
      error("DMDBuildSimulation: Could not extract InitialPotential");
      return (false);   
    }

    if (GeomToComp2.get_rep())
    {
      Potential0 = GeomToComp2*Potential0;
    }
    
    
    if(!(mathalgo.ResizeMatrix(Potential0,Potential0,num_totalnodes,1)))
    {
      error("DMDBuildSimulation: Could not resize DomainType matrix");
      return (false);   
    }
  }

  
  if (!(DMDBuildMembraneMatrix(membranetable,nodetypes,cmvalues,num_volumenodes,num_synnodes,NodeType,VolumeVec,CapacitanceVec,synmatrix)))
  {
    error("DMDBuildSimulation: Could not build Synapse matrix");
    return (false); 
  }

  // Somehow CardioWave uses the negative matrix
  MatrixHandle sysmatrix = synmatrix - fematrix;
  if (sysmatrix.get_rep() == 0)
  {
    error("DMDBuildSimulation: Could not build system sparse matrix");
    return (false);  
  }
  
  synmatrix = 0;
  fematrix = 0;

  // What do we have now:
  // We have sysmatrix a combined matrix with synapses and the volume fe matrix
  // We have the nodetype
  // We have the volvec
  // We have the membranetable
  // We have the stimulustable
  // We have the referencetable

  remark("Created the full linear system");  

  // Step 5:
  // Now optimize the system:
  MatrixHandle mapping;
  MatrixHandle tmapping;
    
  //optimizesystemz = false;
  //optimizesystemy = false;
  //optimizesystemx = false;  
  
  if (optimizesystemz)
  {
    if(!(DMDOptimizeByCoordinate("z",ElementType,Membranes,membranetable,GeomToComp,mapping)))
    {
      error("DMDBuildSimulation: Matrix reordering failed");
      return (false);    
    }
    tmapping = mapping->transpose();
    sysmatrix = mapping*sysmatrix*tmapping;
  }
  else if (optimizesystemy)
  {
    if(!(DMDOptimizeByCoordinate("y",ElementType,Membranes,membranetable,GeomToComp,mapping)))
    {
      error("DMDBuildSimulation: Matrix reordering failed");
      return (false);    
    }
    tmapping = mapping->transpose();
    sysmatrix = mapping*sysmatrix*tmapping;
  }
  else if (optimizesystemx)
  {
    if(!(DMDOptimizeByCoordinate("x",ElementType,Membranes,membranetable,GeomToComp,mapping)))
    {
      error("DMDBuildSimulation: Matrix reordering failed");
      return (false);    
    }
    tmapping = mapping->transpose();
    sysmatrix = mapping*sysmatrix*tmapping;
  }  
  else
  {
    if(!(mathalgo.IdentityMatrix(num_totalnodes,mapping)))
    {
      error("DMDBuildSimulation: Matrix reordering failed");
      return (false);    
    }
    tmapping = mapping->transpose();
  }
  
  remark("Calculated reordering of system");  


  // For visualization bundle
  MatrixHandle ElementMapping, InverseElementMapping;
  
  if (CompToGeom.get_rep())
  {
    if(!(mathalgo.ResizeMatrix(CompToGeom,ElementMapping,CompToGeom->nrows(),num_totalnodes)))
    {
      error("DMDBuildSimulation: Could not resize ElementMapping matrix");
      return (false);   
    }  

    if(!(mathalgo.ResizeMatrix(GeomToComp,InverseElementMapping,num_totalnodes,GeomToComp2->ncols())))
    {
      error("DMDBuildSimulation: Could not resize InverseElementMapping matrix");
      return (false);   
    }  

  }
  else
  {
    if(!(mathalgo.IdentityMatrix(num_totalnodes,ElementMapping)))
    {
      error("DMDBuildSimulation: Could not build ElementMapping Matrix");
      return (false);    
    }
     
    if(!(mathalgo.ResizeMatrix(ElementMapping,ElementMapping,num_volumenodes,num_totalnodes)))
    {
      error("DMDBuildSimulation: Could not resize ElementMapping matrix");
      return (false);   
    }  

    if(!(mathalgo.IdentityMatrix(num_totalnodes,InverseElementMapping)))
    {
      error("DMDBuildSimulation: Could not build InverseElementMapping Matrix");
      return (false);    
    }
     
    if(!(mathalgo.ResizeMatrix(InverseElementMapping,InverseElementMapping,num_totalnodes,num_volumenodes)))
    {
      error("DMDBuildSimulation: Could not resize ElementMapping matrix");
      return (false);   
    }  


  }
  
  ElementMapping = ElementMapping*tmapping;
  InverseElementMapping = mapping*InverseElementMapping;

  // All the tables we just build should have been converted to the new numbering
  GeomToComp = 0;
  CompToGeom = 0;


  // Reorder domain properties
  NodeType = mapping*NodeType;
  if (Potential0.get_rep()) Potential0 = mapping*Potential0;
  CapacitanceVec = mapping*CapacitanceVec;
  VolumeVec = mapping*VolumeVec;

  if (!(dataioalgo.WriteMatrix(filename_mapping,mapping,"CardioWave Sparse Matrix")))
  {
    error("DMDBuildSimulation: Could not write mapping");  
    return (false);
  }
  
  if (!(dataioalgo.WriteMatrix(filename_sysmatrix,sysmatrix,"CardioWave Sparse Matrix")))
  {
    error("DMDBuildSimulation: Could not write system matrix");  
    return (false);
  }

  if (!(dataioalgo.WriteMatrix(filename_surface,VolumeVec,"CardioWave FP Vector")))
  {
    error("DMDBuildSimulation: Could not write volume vector");  
    return (false);
  }

  if (!(dataioalgo.WriteMatrix(filename_capacitance,CapacitanceVec,"CardioWave FP Vector")))
  {
    error("DMDBuildSimulation: Could not write volume vector");  
    return (false);
  }

  if (!(dataioalgo.WriteMatrix(filename_nodetype,NodeType,"CardioWave Byte Vector")))
  {
    error("DMDBuildSimulation: Could not write nodetype vector");  
    return (false);
  }

  if (Potential0.get_rep())
  {
    if (!(dataioalgo.WriteMatrix(filename_potential0,Potential0,"CardioWave FP Vector")))
    {
      error("DMDBuildSimulation: Could not write initial potential vector");  
      return (false);
    }
  }
  
  if (!(dataioalgo.WriteMatrix(filename_surface+".mat",VolumeVec,"Matlab Matrix")))
  {
    error("DMDBuildSimulation: Could not write volume vector");  
    return (false);
  }

  if (!(dataioalgo.WriteMatrix(filename_nodetype+".mat",NodeType,"Matlab Matrix")))
  {
    error("DMDBuildSimulation: Could not write nodetype vector");  
    return (false);
  }

  if (!(dataioalgo.WriteMatrix(filename_nodetype+".mapping.mat",mapping,"Matlab Matrix")))
  {
    error("DMDBuildSimulation: Could not write nodetype vector");  
    return (false);
  }
  
  if (Potential0.get_rep())
  {
    if (!(dataioalgo.WriteMatrix(filename_potential0+".mat",Potential0,"Matlab Matrix")))
    {
      error("DMDBuildSimulation: Could not write initial potential vector");  
      return (false);
    }
  }  
  
  
  remark("Created domain files");  
 
  // clean memory
  sysmatrix = 0;
  NodeType = 0;
  VolumeVec = 0;
  
  if (tmapping.get_rep() == 0)
  {
    error("DMDBuildSimulation: Could not build mapping matrix");
    return (false);      
  }

  SparseRowMatrix* spr = dynamic_cast<SparseRowMatrix *>(tmapping.get_rep());
  if (spr == 0)
  {
    error("DMDBuildSimulation: Could not build mapping matrix");
    return (false);     
  }

  Matrix::index_type* renumber = spr->columns;
  
  // Build membrane table
  for (size_t p=0; p<membranetable.size();p++)
  {
    for (size_t q=0; q< membranetable[p].size(); q++)
    { 
      membranetable[p][q].snode = renumber[membranetable[p][q].snode];
      membranetable[p][q].node1 = renumber[membranetable[p][q].node1];
      membranetable[p][q].node2 = renumber[membranetable[p][q].node2];
    }
  }
  
  
  for (size_t p=0; p<stimulustable.size(); p++)
  {
    for (size_t q=0; q< stimulustable[p].size(); q++)
    {
      stimulustable[p][q].node = renumber[stimulustable[p][q].node];
    }
  }

  for (size_t p=0; p<referencetable.size(); p++)
  {
    for (size_t q=0; q< referencetable[p].size(); q++)
    {
      referencetable[p][q].node = renumber[referencetable[p][q].node];
    }
  }

  remark("Renumbered Membrane table, Stimulus Table, and Reference Table");   
  
  // Build membrane table file
  
  try
  {
    std::ofstream memfile;
    memfile.open(filename_membrane.c_str());
    for (size_t p=0; p<membranetable.size();p++)
    {
      for (size_t q=0; q< membranetable[p].size(); q++)
      { 
        memfile << membranetable[p][q].snode << " " << membranetable[p][q].node2 << " " << membranetable[p][q].node1 << "\n";
      }
    }
  }
  catch (...)
  {
    error("DMDBuildSimulation: Could not write membrane connection file");
    return (false);  
  }


  try
  {
    std::ofstream stimfile;
    stimfile.open(filename_stimulus.c_str());
    for (size_t p=0; p<stimulustable.size();p++)
    {
      for (size_t q=0; q< stimulustable[p].size(); q++)
      { 
        if (stimulusduration[p] == -1.0)
        { 
          stimfile << stimulustable[p][q].node <<  " " << stimuluscurrent[p]*stimulustable[p][q].weight << ", " << stimulusstart[p] << ", " << stimulusend[p] << ", I\n";
        }
        else
        {
          stimfile << stimulustable[p][q].node <<  " " << stimuluscurrent[p]*stimulustable[p][q].weight << ", " << stimulusstart[p] << ", " << stimulusend[p] << ", " << stimulusduration[p] << ", I\n";        
        }
      }
    }
  }
  catch (...)
  {
    error("DMDBuildSimulation: Could not write stimulus file");
    return (false);  
  }

  try
  {
    std::ofstream reffile;
    reffile.open(filename_reference.c_str());
    
    for (size_t p=0; p<referencetable.size();p++)
    {
      for (size_t q=0; q< referencetable[p].size(); q++)
      { 
        reffile << referencetable[p][q].node <<  " " << referencevalues[p] << ", I\n";
      }
    }
  }
  catch (...)
  {
    error("DMDBuildSimulation: Could not write reference file");
    return (false);  
  }

  remark("Wrote stimulus, reference and membrane table"); 

  for (size_t p=0; p < num_membranes; p++)
  {
    membranemapping[p] = membranemapping[p]*tmapping;
    membranemapping2[p] = membranemapping2[p]*tmapping;
  }


  remark("Renumbered membrane mapping matrices"); 


  VisualizationBundle = new Bundle();
  BundleHandle VolumeField = new Bundle();
  BundleHandle ElementTypeField = new Bundle();
  BundleHandle ConductivityField = new Bundle();


  VolumeField->setMatrix("Mapping",ElementMapping);
  VolumeField->setMatrix("InverseMapping",InverseElementMapping);
  VisualizationBundle->setBundle("Tissue",VolumeField);
  VisualizationBundle->setBundle("ElementType",ElementTypeField);
  VisualizationBundle->setBundle("Conductivity",ConductivityField);
  ElementTypeField->setField("Field",ElementType);
  ConductivityField->setField("Field",Conductivity);
   
  for (size_t p=0; p <num_membranes; p++)
  {
    std::ostringstream oss;
    oss << "Membrane_" << p;
    
    BundleHandle MembraneBundle = new Bundle;
    if (MembraneBundle.get_rep() == 0)
    { 
      error("DMDBuildSimulation: Could not allocate new output bundle");
      return (false);
    }
    MembraneBundle->setField("Field",Membranes[p]);
    MembraneBundle->setMatrix("Mapping",membranemapping[p]);
    MembraneBundle->setMatrix("MappingFromMemNode",membranemapping2[p]);
    VisualizationBundle->setBundle(oss.str(),MembraneBundle);
  }

  MatrixHandle Mapping;
  std::vector<Matrix::index_type> columnselection;
  std::vector<std::vector<Matrix::index_type> > indices(num_electrodes);   

  SCIRunAlgo::SelectMatrixAlgo select_algo_;
  select_algo_.set_progress_reporter(pr_);
  
  SCIRunAlgo::AppendMatrixAlgo append_algo_;
  append_algo_.set_progress_reporter(pr_);
  
  SCIRunAlgo::FindClosestNodeAlgo fnode_algo_;
  fnode_algo_.set_progress_reporter(pr_);
 
  SCIRunAlgo::FindClosestNodeByValueAlgo fnodev_algo_;
  fnodev_algo_.set_progress_reporter(pr_); 
  for (size_t p=0; p < num_electrodes; p++)
  {
    MatrixHandle LocalMapping;
    std::vector<Matrix::index_type> rowselection;
    std::ostringstream oss;
    oss << "electrode " << p;
  
    if (electrodemembrane[p] >= 0)
    {    
      if(!(fnode_algo_.run(Membranes[electrodemembrane[p]],Electrodes[p],rowselection)))
      {
        error("DMDBuildSimulation: Could not find electrodes within specified membrane for "+oss.str());
        return (false);
      }

      select_algo_.set_option("method","select_rows");
      if(!(select_algo_.run(membranemapping[electrodemembrane[p]],LocalMapping,rowselection)))
      {
        error("DMDBuildSimulation: Could not select rows from mapping matrix for electrodes of "+oss.str());
        return (false);
      }

      append_algo_.set_option("method","append_rows");
      if(!(append_algo_.run(Mapping,Mapping,LocalMapping,indices[p])))
      {
        error("DMDBuildSimulation: Could not append rows to mapping matrix for electrodes of "+oss.str());
        return (false);      
      }
    }
    else
    {
      fnodev_algo_.set_scalar("min",electrodedomainmin[p]);
      fnodev_algo_.set_scalar("max",electrodedomainmax[p]);
      if(!(fnodev_algo_.run(Membranes[electrodemembrane[p]],Electrodes[p],rowselection)))
      {
        error("DMDBuildSimulation: Could not find electrodes within specified domain for "+oss.str());
        return (false);      
      }

      select_algo_.set_option("method","select_rows");
      if(!(select_algo_.run(ElementMapping,LocalMapping,rowselection)))
      {
        error("DMDBuildSimulation: Could not select rows from mapping matrix for electrodes of "+oss.str());
        return (false);
      }      

      append_algo_.set_option("method","append_rows");
      if(!(append_algo_.run(Mapping,Mapping,LocalMapping,indices[p])))
      {
        error("DMDBuildSimulation: Could not append rows to mapping matrix for electrodes of "+oss.str());
        return (false);      
      } 
			
    }
  }

  SCIRunAlgo::FindMatrixAlgo find_algo_;
  find_algo_.set_progress_reporter(pr_);
  
  if (num_electrodes >0)
  {
    find_algo_.set_option("method","find_nonzerocolumns");
    if (!(find_algo_.run(Mapping,columnselection)))
    {
      error("DMDBuildSimulation: Could not find non zero rows in mapping matrix");
      return (false);      
    }
    
    select_algo_.set_option("method","select_columns");
    if (!(select_algo_.run(Mapping,Mapping,columnselection)))
    {
      error("DMDBuildSimulation: Could not select non zero rows in mapping matrix");
      return (false);      
    }

    if (!(dataioalgo.WriteMatrix("MappingMatrix.mat",Mapping)))
    {
      error("DMDBuildSimulation: Could not write mapping matrix");  
      return (false);
    }

  }

  FieldInformation fi_et(ElementType); fi_et.make_lineardata();
  ElementType = CreateField(fi_et,ElementType->mesh());
  if (!(ElementType.get_rep()))
  {
    error("DMDBuildSimulation: Could not clear field");
    return (false);            
  }
  VolumeField->setField("Field",ElementType);


  for (size_t p=0; p <num_electrodes; p++)
  {
    std::ostringstream oss;
    oss << "Electrode_" << p;  
    
    MatrixHandle EMapping;
    select_algo_.set_option("method","select_rows");
    select_algo_.run(Mapping,EMapping,indices[p]);
  
    BundleHandle ElectrodeBundle = new Bundle;
    ElectrodeBundle->setField("Field",Electrodes[p]);
    ElectrodeBundle->setMatrix("Mapping",EMapping);
    VisualizationBundle->setBundle(oss.str(),ElectrodeBundle);
  }  

  if (!(dataioalgo.WriteBundle(filename_visbundle,VisualizationBundle)))
  {
    error("DMDBuildSimulation: Could not write visualization information");
    return (false);     
  }

  remark("Created Visualization Bundle"); 
  
  try
  {
    std::ofstream electrodefile;
    electrodefile.open(filename_electrode.c_str());
    for (size_t p=0; p<columnselection.size();p++)
    {
      electrodefile << columnselection[p] << ", I\n";
    }
  }
  catch (...)
  {
    error("DMDBuildSimulation: Could not write electrodes file");
    return (false);  
  }

  remark("Wrote electrode table");   
  return (true);
}
 
  
   
     

} // end namespace SCIRun
