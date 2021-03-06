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
 *  InterfaceWithNeuroFEMForward.cc:
 *
 *  Written by:
 *   Michael Callahan
 *   Department of Computer Science
 *   University of Utah
 *   April 2005
 *
 */

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>

#include <Core/ImportExport/Field/FieldIEPlugin.h>
#include <Core/ImportExport/Matrix/MatrixIEPlugin.h>
#include <Core/ImportExport/ExecConverter.h>

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <Core/Util/StringUtil.h>
#include <Core/Util/Dir.h> // for MKDIR

#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>

#ifdef _WIN32
#define mode_t int
#include <io.h>
#endif

namespace SCIRun {
 
class InterfaceWithNeuroFEMForward : public Module {

  public:
    GuiString ipm_pathTCL_;
    GuiInt associativityTCL_;
    GuiDouble eps_matTCL_;
    
    InterfaceWithNeuroFEMForward(GuiContext* ctx);
    virtual ~InterfaceWithNeuroFEMForward() {}

    virtual void execute();

  private:
    bool write_par_file(const std::string& filename);
    bool write_pebbles_file(const std::string& filename);
};


DECLARE_MAKER(InterfaceWithNeuroFEMForward)

InterfaceWithNeuroFEMForward::InterfaceWithNeuroFEMForward(GuiContext* ctx)
  : Module("InterfaceWithNeuroFEMForward", ctx, Filter, "NeuroFEM", "BioPSE"),
    ipm_pathTCL_(get_ctx()->subVar("ipm_pathTCL")),
    associativityTCL_(get_ctx()->subVar("associativityTCL")),
    eps_matTCL_(get_ctx()->subVar("eps_matTCL"))
{
}


static bool
write_geo_tet(ProgressReporter *pr, FieldHandle field, const char *filename)
{
  FieldInformation fi(field);

  if (!fi.is_tetvolmesh())
  {
    pr->error("Field does not contain a TetVolMesh.");
    return false;
  }

  VMesh* mesh = field->vmesh();
  
  VMesh::Node::iterator nitr, neitr;
  VMesh::Node::size_type nsize; 
  mesh->begin(nitr);
  mesh->end(neitr);
  mesh->size(nsize);

  VMesh::Elem::iterator citr, ceitr;
  VMesh::Elem::size_type csize; 
  mesh->begin(citr);
  mesh->end(ceitr);
  mesh->size(csize); 

  FILE *f = fopen(filename, "w");
  if (f == NULL)
  {
    pr->error(std::string("Unable to open file '") + filename + "' for writing.");
    return false;
  }

  // Write header.
  fprintf(f, "BOI - GEOMETRIEFILE\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "BOI - STEUERKARTE\n");
  fprintf(f, "  ANZAHL DER KNOTEN             : %d\n",(unsigned)(nsize));
  fprintf(f, "  ANZAHL DER ELEMENTE           : %d\n",(unsigned)(csize));
  fprintf(f, "  GEOMETR. STRUKTUR - DIMENSION :      3\n");
  fprintf(f, "EOI - STEUERKARTE\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "===================================================================\n");    
  fprintf(f, "BOI - KOORDINATENKARTE\n");  

  // write it out
  int odd=1;
  while (nitr != neitr)
  {
    Point p;
    mesh->get_point(p, *nitr);

    if(odd==1) 
    {
      fprintf(f, "   %11.5e %11.5e %11.5e ", p.x(), p.y(), p.z());
      odd=0;
    }
    else
    {
      fprintf(f, "   %11.5e %11.5e %11.5e\n", p.x(), p.y(), p.z());
      odd=1;
    }
    // Write p;

    ++nitr;
  }
  // write middle part
  fprintf(f, "\nEOI - KOORDINATENKARTE\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "===================================================================\n"); 
  fprintf(f, "BOI - ELEMENTKNOTENKARTE\n");

  while (citr != ceitr)
  {
    VMesh::Node::array_type nodes; // 0 based
    mesh->get_nodes(nodes, *citr);

    // Write nodes - Please notice the indexing order differs.
    fprintf(f,  "  303: %6d%6d%6d%6d\n",
            static_cast<int>(nodes[1]+1), static_cast<int>(nodes[0]+1), 
            static_cast<int>(nodes[2]+1), static_cast<int>(nodes[3]+1));

    ++citr;
  }

  // Write tail, bottom part.
  fprintf(f, "EOI - ELEMENTKNOTENKARTE\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "EOI - GEOMETRIEFILE\n");
  
  fclose(f);
  
  return true;
}


static bool
write_geo_hex(ProgressReporter *pr, FieldHandle field, const char *filename)
{
  FieldInformation fi(field);

  if (!fi.is_hexvolmesh())
  {
    pr->error("Field does not contain a HexVolMesh.");
    return false;
  }

  VMesh* mesh = field->vmesh();

  // TODO: WRITE HEX MESH
  VMesh::Node::iterator nitr, neitr;
  VMesh::Node::size_type nsize; 
  mesh->begin(nitr);
  mesh->end(neitr);
  mesh->size(nsize);

  VMesh::Elem::iterator citr, ceitr;
  VMesh::Elem::size_type csize; 
  mesh->begin(citr);
  mesh->end(ceitr);
  mesh->size(csize); 

  FILE *f = fopen(filename, "w");
  if (f == NULL)
  {
    pr->error(std::string("Unable to open file '") + filename + "' for writing.");
    return false;
  }

  // Write header.
  fprintf(f, "BOI - GEOMETRIEFILE\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "BOI - STEUERKARTE\n");
  fprintf(f, "  ANZAHL DER KNOTEN             : %d\n",(unsigned)(nsize));
  fprintf(f, "  ANZAHL DER ELEMENTE           : %d\n",(unsigned)(csize));
  fprintf(f, "  GEOMETR. STRUKTUR - DIMENSION :      3\n");
  fprintf(f, "EOI - STEUERKARTE\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "===================================================================\n");    
  fprintf(f, "BOI - KOORDINATENKARTE\n");  

  // write it out
  int odd=1;
  while (nitr != neitr)
  {
    Point p;
    mesh->get_point(p, *nitr);

    if(odd==1) {
      fprintf(f, "   %11.5e %11.5e %11.5e ", p.x(), p.y(), p.z());
      odd=0;
    }
    else{
      fprintf(f, "   %11.5e %11.5e %11.5e\n", p.x(), p.y(), p.z());
      odd=1;
    }
    // Write p;

    ++nitr;
  }
  // write middle part
  fprintf(f, "\nEOI - KOORDINATENKARTE\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "===================================================================\n"); 
  fprintf(f, "BOI - ELEMENTKNOTENKARTE\n");

  while (citr != ceitr)
  {
    VMesh::Node::array_type nodes; // 0 based
    mesh->get_nodes(nodes, *citr);

    fprintf(f,  "  323: %6d%6d%6d%6d%6d%6d%6d%6d\n",
            static_cast<int>(nodes[0]+1), static_cast<int>(nodes[1]+1), 
            static_cast<int>(nodes[2]+1), static_cast<int>(nodes[3]+1), 
            static_cast<int>(nodes[4]+1), static_cast<int>(nodes[5]+1), 
            static_cast<int>(nodes[6]+1), static_cast<int>(nodes[7]+1));

    ++citr;
  }

  // Write tail, bottom part.
  fprintf(f, "EOI - ELEMENTKNOTENKARTE\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "===================================================================\n");
  fprintf(f, "EOI - GEOMETRIEFILE\n");
  
  fclose(f);

  return true;
}


static bool
write_knw_tet(ProgressReporter *pr, FieldHandle field, const char *filename)
{
  FieldInformation fi(field);

  if (!fi.is_tetvolmesh())
  {
    pr->error("Field does not contain a TetVolMesh.");
    return false;
  }
  
  VField *tfield = field->vfield();
  VMesh *mesh = field->vmesh();
  
  FILE *f = fopen(filename, "w");
  if (f == NULL)
  {
    pr->error(std::string("Unable to open file '") + filename + "' for writing.");
    return false;
  }

  // Write header.
  fprintf(f, "BOI - TENSORVALUEFILE\n");
  fprintf(f, "========================================================\n");
  fprintf(f, "========================================================\n");
  fprintf(f, "BOI - TENSOR\n");  

  // Write it out.
  VMesh::Elem::iterator itr, eitr;
  mesh->begin(itr);
  mesh->end(eitr);
  
  while (itr != eitr)
  {
    Tensor t;
    tfield->value(t, *itr);

    fprintf(f, "%10d  %11f  %11f  %11f\n           %11f  %11f  %11f\n",
            static_cast<int>(1+*itr),  // fortran index, +1
            t.mat_[0][0], t.mat_[1][1], t.mat_[2][2],
            t.mat_[0][1], t.mat_[0][2], t.mat_[1][2]);

    ++itr;
  }

  // Write tail, bottom part.
  fprintf(f, "EOI - TENSOR\n");
  fprintf(f, "========================================================\n");
  fprintf(f, "========================================================\n");
  fprintf(f, "EOI - TENSORVALUEFILE\n");

  fclose(f);

  return true;
}


static bool
write_knw_hex(ProgressReporter *pr, FieldHandle field, const char *filename)
{

  FieldInformation fi(field);

  if (!fi.is_hexvolmesh())
  {
    pr->error("Field does not contain a HexVolMesh.");
    return false;
  }
  
  VField *tfield = field->vfield();
  VMesh *mesh = field->vmesh();

  FILE *f = fopen(filename, "w");
  if (f == NULL)
  {
    pr->error(std::string("Unable to open file '") + filename + "' for writing.");
    return false;
  }

  // Write header.
  fprintf(f, "BOI - TENSORVALUEFILE\n");
  fprintf(f, "========================================================\n");
  fprintf(f, "========================================================\n");
  fprintf(f, "BOI - TENSOR\n");  

  // Write it out.
  VMesh::Elem::iterator itr, eitr;
  mesh->begin(itr);
  mesh->end(eitr);
  
  while (itr != eitr)
  {
    Tensor t;
    tfield->value(t, *itr);

    fprintf(f, "%10d  %11f  %11f  %11f\n           %11f  %11f  %11f\n",
            static_cast<int>(1+*itr),  // fortran index, +1
            t.mat_[0][0], t.mat_[1][1], t.mat_[2][2],
            t.mat_[0][1], t.mat_[0][2], t.mat_[1][2]);

    ++itr;
  }

  // Write tail, bottom part.
  fprintf(f, "EOI - TENSOR\n");
  fprintf(f, "========================================================\n");
  fprintf(f, "========================================================\n");
  fprintf(f, "EOI - TENSORVALUEFILE\n");

  fclose(f);

  return true;
}


bool
write_geo_file(ProgressReporter *pr, FieldHandle field, const char *filename)
{
  FieldInformation fi (field);
  if (fi.is_tetvolmesh())
  {
    return write_geo_tet(pr, field, filename);
  }
  else if (fi.is_hexvolmesh())
  {
    return write_geo_hex(pr, field, filename);
  }
  else
  {
    pr->error("In write_geo_file, expected a TVMesh or HVMesh.");
    return false;
  }
}


bool
write_knw_file(ProgressReporter *pr, FieldHandle field, const char *filename)
{
  FieldInformation fi(field);
  
  if (fi.is_tetvolmesh())
  {
    return write_knw_tet(pr, field, filename);
  }
  else if (fi.is_hexvolmesh())
  {
    return write_knw_hex(pr, field, filename);
  }
  else
  {
    pr->error("In write_knw_file, expected a TVMesh or HVMesh.");
    return false;
  }
}



bool
write_elc_file(ProgressReporter *pr, FieldHandle fld, const char *filename)
{
  FieldInformation fi(fld);
  
  if (!fi.is_pointcloudmesh());
  {
    pr->error("Field does not contain a PCMesh.");
    return false;
  }

  VMesh *mesh = fld->vmesh();


  VMesh::Node::iterator niter; 
  VMesh::Node::iterator niter_end; 
  VMesh::Node::size_type nsize; 
  mesh->begin(niter);
  mesh->end(niter_end);
  mesh->size(nsize);
  
  FILE *f = fopen(filename, "w");
  if (f == NULL)
  {
    pr->error(std::string("Unable to open file '") + filename + "' for writing.");
    return false;
  }

  //write header
  fprintf(f, "# %d Electrodes\n", (unsigned)(nsize));
  fprintf(f, "ReferenceLabel	avg\n");
  fprintf(f, "NumberPositions=\t%d\n",(unsigned)(nsize));
  fprintf(f, "UnitPosition	mm\n");
  fprintf(f, "Positions\n");

  //write it out
  while(niter != niter_end) 
  {
    Point p;
    mesh->get_center(p, *niter);
    fprintf(f, "%9.3f %9.3f %9.3f\n", p.x(), p.y(), p.z());
    ++niter;
  }     

  fprintf(f,"Labels\n");
  /*
  for (unsigned int c=0; c<nsize; c++){
    fprintf(f,"E%d ",c+1);
  }
  */
  // for testing model
  /*  fprintf(f,"FPz	Cz	Oz	Pz	Fz	T10	T7	C3	T8	C4	F3	F4	Fp1	Fp2	F7	F8	O1	O2	P3	P4	P7	P8	FT7	FT8	FC3	FC4	CP5	CP6	AFz	AF3	AF4	AF7	AF8	F5	F6	F9	F10	FC5	FC6	FT9	FT10	C5	C6	CPz	CP3	TP7	TP9	CP4	TP8	TP10	P5	P9	P6	P10	POz	PO3	PO7	PO4	PO8	FCz	F1	F2	FC1	FC2	C1	C2	CP1	CP2	P1	P2	T9	\n");
   */
  
  // Only for epilepsy model
  fprintf(f, "CZ	A1	A2	A3	A4	A5	A6	A7	A8	A9	A10	A11	A12\n");
  fprintf(f, "A13	A14	A15	A16	A17	A18	A19	A20	A21	A22	A23	A24	B1\n");
  fprintf(f, "B2	B3	B4	B5	B6	B7	B8	B9	B10	B11	B12	B13	B14\n");
  fprintf(f, "B15	B16	B17	B18	B19	B20	B21	B22	B23	B24	B25	B26	B27\n");
  fprintf(f, "B28	B29	B30	B31	B32	B35	B36	B37	B38	B39	B40	B41	B42\n");
  fprintf(f, " B43	B44	B45	B46	B47	B48	STR1	STR2	STR3	STR4	STR5	STR6	STR7	STR8\n");
 
  fclose(f);
  
  return true;
}


static bool
write_dip_file(ProgressReporter *pr, MatrixHandle mh1, MatrixHandle mh2, const char *filename)
{

  DenseMatrix *dp_pos = dynamic_cast<DenseMatrix *>(mh1.get_rep());
  if (!dp_pos) {
    pr->error("Error -- input field wasn't a DenseMatrix");
    return false;
  }

  int nd=dp_pos->nrows();
  
  DenseMatrix *dp_vec = dynamic_cast<DenseMatrix *>(mh2.get_rep());
  if (!dp_vec) {
    pr->error("Error -- input field wasn't a DenseMatrix");
    return false;
  }

  int nr=dp_vec->nrows();
  int nc=dp_vec->ncols();

  pr->remark("Number of rows = " + to_string(nr) + ".");
  pr->remark("Number of columns = " + to_string(nc) + ".");	     


  FILE *f = fopen(filename, "wt");
  if (!f) {
    pr->error(std::string("Error -- couldn't open output file: ") + filename );
    return false;
  }

  // write header
  fprintf(f,"# 1 dipole(s) at %d time steps\n", nr);
  fprintf(f,"UnitPosition\tmm\n");
  fprintf(f,"UnitMoment\tnAm\n");
  fprintf(f,"UnitTime\tms\n");  
  fprintf(f,"NumberPositions=\t%d\n",nd);
  fprintf(f,"NumberTimeSteps=\t%d\n",nr);
  fprintf(f,"TimeSteps\t0(1)%d\n",(unsigned)nr-1);
  fprintf(f,"FirstTimeStep\t0\n");
  fprintf(f,"LastTimeStep\t%d\n",(unsigned)(nr-1));

  // write dipole moment
  fprintf(f,"MomentsFixed\n");
  double sum = 0.0;
  for (int c=0; c<nc; c++){
    sum += ((*dp_vec)[0][c])*((*dp_vec)[0][c]);
  }

  if(sum <= 0.0){
    std::cerr << "Dipole vectors are all zeros" << "\n";
    return false;
  }
 
  sum = sqrt(sum)/sqrt(2.0);;
  fprintf(f,"%f\t%f\t%f\n",((*dp_vec)[0][0])/sum,((*dp_vec)[0][1])/sum, ((*dp_vec)[0][2])/sum);

  // write dipole position
  fprintf(f,"PositionsFixed\n");
  for (int i=0; i<3; i++)
    fprintf(f,"%7.3f ",(*dp_pos)[0][i]);
  fprintf(f,"\n");

  //  write dipole magnitude
  fprintf(f,"Magnitudes\n");
  for (int r=0; r<nr; r++)
  {
    double sum = 0.0;
    for (int c=0; c<nc; c++)
    {
      sum += ((*dp_vec)[r][c])*((*dp_vec)[r][c]);
     }
    fprintf(f,"%e ",sqrt(sum)/sqrt(2.0));
  }
   
  fprintf(f,"\n");
  fprintf(f,"NoLabels\n");
  fclose(f);

  return true;
}


static bool
send_result_file(ProgressReporter *pr, Module *module,
		 std::string filename)
{
  std::string extfilename = filename + ".msr"; 
  std::ifstream matstream(extfilename.c_str(),std::ios::in);
  if (matstream.fail())
  {
    pr->error("Could not open results file " + filename + ".");
    return false;
  }

  std::string tmp;
  int nc, nr;
  matstream >> tmp >> nc;
  matstream >> tmp >> nr;
  pr->remark("Number of Electrodes = " + to_string(nc) + ".");
  pr->remark("Number of Time Steps = " + to_string(nr) + ".");

  // Skip text lines.
  for (int i=0; i<4; ++i)
  {
    getline(matstream, tmp);
   }

  MatrixHandle dm = new DenseMatrix(nr, nc);

  // Write number of row & column.
  //dm->put(0, 0, nr);
  //dm->put(0, 1, nc);

  // Write potentials on electrodes for each time step.
  int r, c;
  for (r=0; r < nr; r++)
  {
    matstream >> tmp;
    for (c=0; c < nc; c++)
    {
      long d;
      matstream >> d;
      dm->put(r, c, d);
    }
  }
  pr->remark("Done building output matrix.");
  
  module->send_output_handle("Forward Potential",dm,true);
  return true;
}


bool
InterfaceWithNeuroFEMForward::write_par_file(const std::string& filename)
{
  FILE *f = fopen(filename.c_str(), "w");
  if (f == NULL)
  {
    error("Unable to open file '" + filename + "' for writing.");
    return false;
  }

  fprintf(f, "#Parameter file: FEM for source simulation\n");
  fprintf(f, "\n");
  fprintf(f, "[ReferenceData]\n");
  fprintf(f, "#Number of channels\n");
  fprintf(f, "numchan= 1\n");
  fprintf(f, "numsamples= 1\n");
  fprintf(f, "startsample= 0\n");
  fprintf(f, "stopsample=  1\n");
  fprintf(f, "\n");
  fprintf(f, "[RegularizationTikhonov]\n");
  fprintf(f, "# estimated signal to noise ratio\n");
  fprintf(f, "estimatedSNR= 10000.0\n");
  fprintf(f, "\n");
  fprintf(f, "[STR]\n");
  fprintf(f, "# spatial regularization parameter\n");
  fprintf(f, "lambda= 0.1\n");
  fprintf(f, "# temporal regularization parameter\n");
  fprintf(f, "mue= 0.1\n");
  fprintf(f, "# degree of temporal smoothness\n");
  fprintf(f, "degree= 2\n");
  fprintf(f, "# fast modus 1 leads to approximated solution, much faster.\n");
  fprintf(f, "fast= 1\n");
  fprintf(f, "\n");
  fprintf(f, "[NeuroFEMSimulator]\n");
  fprintf(f, "# GENERAL PARAMETERS\n");
  fprintf(f, "# residuum of the forward solution \n");
  fprintf(f, "toleranceforwardsolution= .10000E-08\n");
  fprintf(f, "# degree of Gaussian integration (2 is recommended) \n");
  fprintf(f, "degreeofintegration= 2\n");
  fprintf(f, "# different analytical solutions for the EEG problem \n");
  fprintf(f, "analyticalsolution= 1\n");
  fprintf(f, "\n");
  fprintf(f, "# should METIS repartition for parallel run?\n");
  fprintf(f, "metisrepartitioning= 1\n");
  fprintf(f, "\n");
  fprintf(f, "# SOLVER (1:Jakobi-CG, 2:IC(0)-CG, 3:AMG-CG, 4:PILUTS(ILDLT)-CG) \n");
  fprintf(f, "solvermethod= 3\n");
  fprintf(f, "# use or use not lead field basis approach\n");
  fprintf(f, "associativity= %d\n", (unsigned int) associativityTCL_.get());
  fprintf(f, "# NeuroFEM Solver\n");
  fprintf(f, "# parameter file for Pebbles solver \n");
  fprintf(f, "pebbles= pebbles.inp\n");
  fprintf(f, "# ONLY for MultiRHS-Pebbles: number of right-hand-sides which are solved simultaneously\n");
  fprintf(f, "pebblesnumrhs= 1\n");
  fprintf(f, "\n");
  fprintf(f, "# SOURCE SIMULATION\n");
  fprintf(f, "# threshold (percentage of the greatest dipole strength) of all dipoles to appear in the result files \n");
  fprintf(f, "dipolethreshold= .10000E+01\n");
  fprintf(f, "# blurring single loads to adjacent nodes by means of the Gaussian (normal) distribution \n");
  fprintf(f, "sourcesimulationepsilondirac= .10000E-08\n");
  fprintf(f, "\n");
  fprintf(f, "# DIPOLE MODELING\n");
  fprintf(f, "# dipole modeling; weighting of the source distribution with the power of the declared value \n");
  fprintf(f, "dipolemodelingsmoothness= 2\n");
  fprintf(f, "# power of the dipole moments to be considered \n");
  fprintf(f, "dipolemodelingorder= 2\n");
  fprintf(f, "# necessary internal scaling factor; should be larger than twice the element edge length \n");
  fprintf(f, "dipolemodelingscale= 20.000\n");
  fprintf(f, "# Lagrangian multiplier for the (inverse) dipole modeling \n");
  fprintf(f, "dipolemodelinglambda= .10000E-05\n");
  fprintf(f, "# source-sink separation of the analytical dipole model \n");
  fprintf(f, "dipolemodelingdistance= 1.0000\n");
  fprintf(f, "# use rango dipole model \n");
  fprintf(f, "dipolemodelingrango= 0\n");
  fprintf(f, "\n");
  fprintf(f, "# Monopole/Dipole\n");
  fprintf(f, "# calculation of the potential distribution for spherical, homogeneous  structures by means of a analytical description \n");
  fprintf(f, "analyticaldipole= 0\n");
  fprintf(f, "# forward solution computation with monopoles or dipoles\n");
  fprintf(f, "dipolesources= 1\n");
  fprintf(f, "# spread a nodal load to adjacent nodes\n");
  fprintf(f, "sourcesimulation= 1\n");
  fprintf(f, "#to compare analytical and numerical solutions, an integral average value of the potential is subtracted\n");
  fprintf(f, "#from the analytical results because they are not related to a mean potential (in contrast to the numerical solution)\n");
  fprintf(f, "averagecorrection= 0\n");
  fprintf(f, "\n");
  fprintf(f, "[NeuroFEMGridGenerator] \n");
  fprintf(f, "# Number of materials \n");
  fprintf(f, "nummaterials= 7\n");
  fprintf(f, "# Conductivities of fem head model\n");
  fprintf(f, "\n");
  fprintf(f, "#define EXTRA           3       // tissue codes (see D1.2b)\n");
  fprintf(f, "#define SKULL           1\n");
  fprintf(f, "#define CSF             8\n");
  fprintf(f, "#define GREY_MATTER     7\n");
  fprintf(f, "#define WHITE_MATTER    6\n");
  fprintf(f, "\n");
  fprintf(f, "#skin_cond       0.33;\n");
  fprintf(f, "#skull_cond      0.0042;\n");
  fprintf(f, "#csf_cond        1.79;\n");
  fprintf(f, "#grey_cond       0.337;\n");
  fprintf(f, "#white_cond_iso  0.14;\n");
  fprintf(f, "\n");
  fprintf(f, "# The first value corresponds to the electrode which will be added to the mesh \n");
  fprintf(f, "conductivities= \n");
  fprintf(f, "1.0 0.33 0.0042 0.33 0.33 0.33 0.033\n");
  fprintf(f, "# Labels in the head model corresponding to the different tissues\n");
  fprintf(f, "# The first value corresponds to the electrode which will be added to the mesh \n");
  fprintf(f, "labels=\n");
  fprintf(f, "1000 3 1 8 7 6 9\n");
  fprintf(f, "\n");
  fprintf(f, "# Tissue labels in the head model for which the tensor valued coductivity \n");
  fprintf(f, "# should be used if available\n");
  fprintf(f, "tensorlabels=\n");
  fprintf(f, "#0 0 1 0 0 1 0\n");
  fprintf(f, "0 0 1 0 0 1 0\n");
  fprintf(f, "\n");
  fprintf(f, "# -1- for sorted conductivities\n");
  fprintf(f, "sortedconductivities= 0\n");
  fprintf(f, "\n");
  fprintf(f, "[InitialGuess]\n");
  fprintf(f, "numinitialguesspos= 1\n");
  fprintf(f, "posx=\n");
  fprintf(f, "#0.044175999 \n");
  fprintf(f, "0.097\n");
  fprintf(f, "posy=\n");
  fprintf(f, "#0.077307999\n");
  fprintf(f, "0.154\n");
  fprintf(f, "posz=\n");
  fprintf(f, "#0.135539993\n");
  fprintf(f, "0.128\n");
  fprintf(f, "numinitialguessdir= 1\n");
  fprintf(f, "dirx=\n");
  fprintf(f, "#0.0\n");
  fprintf(f, "1.0  \n");
  fprintf(f, "diry=\n");
  fprintf(f, "#1.0\n");
  fprintf(f, "0.0\n");
  fprintf(f, "dirz=\n");
  fprintf(f, "#1.0\n");
  fprintf(f, "0.0 \n");
  fprintf(f, "\n");
  fprintf(f, "[FileFormats]   \n");
  fprintf(f, "#ReferenceData file: 1= Vista, 2= ASA, 3= ASCII\n");
  fprintf(f, "ReferenceDataIn= 2\n");
  fprintf(f, "#LeadfieldMatrix input file: 1= Vista, 2= ASA \n");
  fprintf(f, "LeadfieldIn= 1\n");
  fprintf(f, "#SensorConfiguration file: 1= Vista, 2= ASA \n");
  fprintf(f, "SensorconfigurationIn= 2 \n");
  fprintf(f, "#Source Space Grid: 1= Vista, 2= ASA, 3= CAUCHY\n");
  fprintf(f, "SourceSpaceIn= 2\n");
  fprintf(f, "#FEM HeadGrid: 1= Vista, 2= ASA, 3= CAUCHY\n");
  fprintf(f, "HeadGridIn= 3\n");
  fprintf(f, "#Output File; Format: 1= Vista, 2 = ASA, 3 = Vista + ASA, 4 = Vista + Ascii, 5 = Cauchy, 6 = SCIRun \n");
  fprintf(f, "ResultOut= 4\n");
  
  fclose(f);
  return true;
}


bool
InterfaceWithNeuroFEMForward::write_pebbles_file(const std::string& filename)
{
  FILE *f = fopen(filename.c_str(), "wt");
  if (f == NULL)
  {
    std::cerr << "Unable to open file"<< filename <<" for writing."<<"\n";
    return false;
  }

  fprintf(f,"EPS_PCCG\n");
  fprintf(f,"1e-8\n");
  fprintf(f,"EPS_ITER\n");
  fprintf(f,"1e-8\n");
  fprintf(f,"EPS_INT\n");
  fprintf(f,"1e-6\n");
  fprintf(f,"EPS_MAT\n");
  fprintf(f,"%.1e\n", (double) eps_matTCL_.get());
  fprintf(f,"ALPHA\n");
  fprintf(f,"0.01\n");
  fprintf(f,"BETA\n");
  fprintf(f,"0.0\n");
  fprintf(f,"SOLVER:_AMG=1,_Jacobi=2,_ILU=3\n");
  fprintf(f,"1\n");
  fprintf(f,"SOLUTION_STRATEGY:_ITER=1,_PCCG=2,_BiCGStab=3\n");
  fprintf(f,"2\n");
  fprintf(f,"PRECOND_STEP\n");
  fprintf(f,"1\n");
  fprintf(f,"COARSE_SYSTEM\n");
  fprintf(f,"1000\n");
  fprintf(f,"COARSE_SOLVER:_LR=1,_LL=2,_DIAG=3,_NAG=4\n");
  fprintf(f,"2\n");
  fprintf(f,"MAX_NUMBER_ITER\n");
  fprintf(f,"200\n");
  fprintf(f,"MAX_NUMBER_SIZE\n");
  fprintf(f,"2000000\n");
  fprintf(f,"MAX_NEIGHBOUR\n");
  fprintf(f,"100\n");
  fprintf(f,"COARSENING:_strong=1,_easy=2,_vmb=3\n");
  fprintf(f,"1\n");
  fprintf(f,"INTERPOLATION:_sophisticated=1,_easy=2,_foolproof=3\n");
  fprintf(f,"3\n");
  fprintf(f,"NORM:_energy=1,_L2=2,_Max=3\n");
  fprintf(f,"1\n");
  fprintf(f,"CYCLE:_V=1,_W=2,_VV=3\n");
  fprintf(f,"1\n");
  fprintf(f,"SMOOTHING_TYPE\n");
  fprintf(f,"1\n");
  fprintf(f,"SMOOTHING_STEP\n");
  fprintf(f,"1\n");
  fprintf(f,"ELEMENT_PRECOND\n");
  fprintf(f,"0\n");
  fprintf(f,"ELEMENT_KAPPA\n");
  fprintf(f,"10\n");
  fprintf(f,"MATRIX_DIFF\n");
  fprintf(f,"1e-1\n");
  fprintf(f,"SHOW_DISPLAY:_nothing=0,_medium=1,_full=2\n");
  fprintf(f,"1\n");
  fprintf(f,"CONV_FACTOR\n");
  fprintf(f,"0\n");

  fclose(f);
  return true;
}

void
InterfaceWithNeuroFEMForward::execute()
{
  FieldHandle condmesh;
  get_input_handle("CondMesh", condmesh,true);

  FieldInformation fi(condmesh);
  
  if (!(fi.is_tetvolmesh() || fi.is_hexvolmesh()))
  {
    error("CondMesh must contain a TVMesh or HVMesh.");
    return;
  }
  if (!(condmesh->vfield()->is_tensor()))
  {
    error("CondMesh must be a tensor field.");
    return;
  }
  if (condmesh->vfield()->basis_order() != 0)
  {
    error("CondMesh must have conductivities at the cell centers.");
    return;
  }

  FieldHandle electrodes;
  get_input_handle("Electrodes", electrodes,true);

  MatrixHandle dipole1;
  get_input_handle("Dipole Positions", dipole1,true);

  MatrixHandle dipole2;
  get_input_handle("Dipole Moments", dipole2,true);

  // Make our tmp directory
  //const std::string tmpdir = "/tmp/InterfaceWithNeuroFEMForward" + to_string(getpid()) +"/";
  const std::string tmpdir = std::string(sci_getenv("SCIRUN_TMP_DIR")) + "/InterfaceWithNeuroFEMForward/";
  const std::string tmplog = tmpdir + "forward.log";
  const std::string resultfile = tmpdir + "fwd_result";
  const std::string condmeshfile = tmpdir + "ca_head.geo";
  const std::string condtensfile = tmpdir + "ca_perm.knw";
  const std::string electrodefile = tmpdir + "electrode.elc";
  const std::string dipolefile = tmpdir + "dipole.dip";
  const std::string parafile = tmpdir + "forward.par";

  try 
  {
    // Make our tmp directory.
    mode_t umsk = umask(00);
    if (MKDIR(tmpdir.c_str(), 0700) == -1)
    {
      //error("Unable to open a temporary working directory.");
      //throw false;
      warning("Unable to open a temporary working directory.");
    }
    umask(umsk);

    // Write out condmesh
    if (!write_geo_file(this, condmesh, condmeshfile.c_str()))
    {
      error("Unable to export CondMesh geo file.");
      throw false;
    }
    // Write out condmesh
    if (!write_knw_file(this, condmesh, condtensfile.c_str()))
    {
      error("Unable to export CondMesh knw file.");
      throw false;
    }

    // Write out Electrodes
    if (!write_elc_file(this, electrodes, electrodefile.c_str()))
    {
      error("Unable to export electrode file.");
      throw false;
    }

    // Write out Dipole
    if (!write_dip_file(this, dipole1, dipole2, dipolefile.c_str()))
    {
      error("Unable to export dipole file.");
      throw false;
    }
    
    // Write out our parameter file.
    if (!write_par_file(parafile)) { throw false; }

    // write out pebble.inp
    if(!write_pebbles_file(tmpdir+"pebbles.inp")) {throw false; }

    // Construct our command line.
    const std::string ipmfile = ipm_pathTCL_.get();
    const std::string command = "(cd " + tmpdir + ";" +
      ipmfile + " -i sourcesimulation" +
      " -h " + condmeshfile + " -p " + parafile + " -s " + electrodefile +
      " -dip " + dipolefile + " -o " + resultfile + " -fwd FEM -sens EEG)";

    
    // Execute the command.  Last arg just a tmp file for logging.
    if (!Exec_execute_command(this, command, tmplog))
    {
      error("The ipm program failed to run for some unknown reason.");
      throw false;
    }

    // Read in the results and send them along.
    if (!send_result_file(this, this, resultfile))
    {
      error("Unable to send output matrix.");
      throw false;
    }

    throw true; // cleanup.
  }
  catch (...)
  {
#if 0
    // Clean up our temporary files.
    unlink(condmeshfile.c_str());
    unlink(condtensfile.c_str());
    unlink(parafile.c_str());
    unlink(electrodefile.c_str());
    unlink(dipolefile.c_str());
    unlink(resultfile.c_str());
    unlink(tmplog.c_str());
    rmdir(tmpdir.c_str());
#endif
  }
}


} // End namespace SCIRun
