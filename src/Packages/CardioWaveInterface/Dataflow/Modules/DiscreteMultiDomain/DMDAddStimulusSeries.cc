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


#include <Dataflow/Network/Module.h>


#include <Core/Datatypes/Bundle.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/String.h>
#include <Core/Datatypes/Matrix.h>
#include <Dataflow/Network/Ports/BundlePort.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Ports/StringPort.h>
#include <Packages/CardioWaveInterface/Core/XML/SynapseXML.h>
#include <Core/Algorithms/Converter/ConverterAlgo.h>


#include <sstream>
#include <vector>
#include <string>
 

namespace CardioWaveInterface {

using namespace SCIRun;

class DMDAddStimulusSeries : public Module {
public:
  DMDAddStimulusSeries(GuiContext*);
  virtual void execute();

private:  
  // TCL tools

  GuiDouble  guidomain_;
  GuiDouble  guidomainmin_;
  GuiDouble  guidomainmax_;
  GuiInt     guidomainrange_;
  GuiDouble  guicurrent_;
  GuiDouble  guistart_;
  GuiDouble  guiend_;
  GuiDouble  guiduration_;
  GuiDouble  guicurrentdensity_;
  GuiInt     guiuseelements_;
};


DECLARE_MAKER(DMDAddStimulusSeries)

DMDAddStimulusSeries::DMDAddStimulusSeries(GuiContext* ctx)
  : Module("DMDAddStimulusSeries", ctx, Source, "DiscreteMultiDomain", "CardioWaveInterface"),
    guidomain_(get_ctx()->subVar("stim-domain"),0.0),
    guidomainmin_(get_ctx()->subVar("stim-domain-min"),0.0),
    guidomainmax_(get_ctx()->subVar("stim-domain-max"),0.0),
    guidomainrange_(get_ctx()->subVar("stim-domain-range"),0),
    guicurrent_(get_ctx()->subVar("stim-current"),400.0),
    guistart_(get_ctx()->subVar("stim-start"),0.1),
    guiend_(get_ctx()->subVar("stim-end"),0.4),
    guiduration_(get_ctx()->subVar("stim-duration"),500.0),
    guicurrentdensity_(get_ctx()->subVar("stim-is-current-density"),400.0),
    guiuseelements_(get_ctx()->subVar("stim-useelements"),0)
{
}

void DMDAddStimulusSeries::execute()
{
  BundleHandle StimulusBundle;
  FieldHandle Geometry;
  MatrixHandle DomainMin, DomainMax;
  MatrixHandle Current;
  MatrixHandle Start;
  MatrixHandle End;
  MatrixHandle Duration;
  MatrixHandle FieldDensity;
  
  // required ones
  if(!(get_input_handle("Geometry",Geometry,true))) return;
  // optional ones
  get_input_handle("StimulusBundle",StimulusBundle,false);
  get_input_handle("Domain",DomainMin,false);
  get_input_handle("Current",Current,false);
  get_input_handle("Start",Start,false);
  get_input_handle("End",End,false);
  get_input_handle("Repeat",Duration,false);

  // If nothing changed, do nothing:
  if (inputs_changed_ || guidomain_.changed() || guidomainmin_.changed() || 
      guidomainmax_.changed() || guidomainrange_.changed() || guicurrent_.changed() || 
      guiduration_.changed() || guistart_.changed() || guiend_.changed() || 
      guicurrentdensity_.changed() || guiuseelements_.changed() ||
      !oport_cached("StimulusBundle"))
  {
    // Get entry point to SCIRun core:
    SCIRunAlgo::ConverterAlgo mc(this);
    double val;

    // Get data from GUI:
    if (mc.MatrixToDouble(DomainMin,val)) guidomain_.set(val);
    if (mc.MatrixToDouble(DomainMin,val)) guidomainmin_.set(val);
    if (mc.MatrixToDouble(DomainMin,val)) guidomainmax_.set(val);
    if (mc.MatrixToDouble(Current,val)) guicurrent_.set(val);
    if (mc.MatrixToDouble(Start,val)) guistart_.set(val);      
    if (mc.MatrixToDouble(End,val)) guiend_.set(val);  
    if (mc.MatrixToDouble(Duration,val)) guiduration_.set(val);  
    
    if (StimulusBundle.get_rep() == 0)
    {
      StimulusBundle = new Bundle();
      if (StimulusBundle.get_rep() == 0)
      {
        error("Could not allocate new stimulus bundle");
        return;
      }
    }
    else
    {
      StimulusBundle.detach();
    }

    
    int stimulus_num = 0;
    std::string fieldname;
    
    {
      std::ostringstream oss;
      oss << "Stimulus_" << stimulus_num;
      fieldname = oss.str(); 
    }
    
    while (StimulusBundle->isBundle(fieldname))
    {
      stimulus_num++;
      {
        std::ostringstream oss;
        oss << "Stimulus_" << stimulus_num;
        fieldname = oss.str(); 
      }
    }

    if (stimulus_num > 0)
    {
      std::string oldsourcefile;    
      BundleHandle OldStimulus;
      StringHandle OldSourceFile;
      
      if (StimulusBundle->isBundle("Stimulus_0"))
      {
        OldStimulus = StimulusBundle->getBundle("Stimulus_0");
        OldSourceFile = OldStimulus->getString("SourceFile");
        if (OldSourceFile.get_rep())
        {
          oldsourcefile = OldSourceFile->get();
          if (oldsourcefile != "StimSeriesFile.c ")
          {
            error("CardioWave does not allow for different stimulus models to be mixed together");
            return;
          }
        }  
      }
    }



    BundleHandle Stimulus;
    Stimulus = new Bundle();
    if (Stimulus.get_rep() == 0)
    {
      error("Could not allocate new stimulus bundle");
      return;
    }
    
    {
      std::ostringstream oss;
      oss << "Stimulus_" << stimulus_num;
      fieldname = oss.str(); 
    }
    StimulusBundle->setBundle(fieldname,Stimulus);
      
    Stimulus->setField("Geometry",Geometry);

    if (guidomainrange_.get())
    {
      val = guidomainmin_.get(); mc.DoubleToMatrix(val,DomainMin); 
      val = guidomainmax_.get(); mc.DoubleToMatrix(val,DomainMax); 
    }
    else
    {
      val = guidomain_.get(); mc.DoubleToMatrix(val,DomainMin); 
      val = guidomain_.get(); mc.DoubleToMatrix(val,DomainMax); 
    } 
    val = guicurrent_.get(); mc.DoubleToMatrix(val,Current); 
    val = guistart_.get(); mc.DoubleToMatrix(val,Start); 
    val = guiend_.get(); mc.DoubleToMatrix(val,End); 
    val = guiduration_.get(); mc.DoubleToMatrix(val,Duration); 
   
    double fielddensity = guicurrentdensity_.get();
    mc.DoubleToMatrix(fielddensity,FieldDensity); 

    MatrixHandle UseElements;
    int useelements = guiuseelements_.get();
    mc.IntToMatrix(useelements,UseElements);

    Stimulus->setMatrix("DomainMin",DomainMin);
    Stimulus->setMatrix("DomainMax",DomainMax);
    Stimulus->setMatrix("Current",Current);
    Stimulus->setMatrix("Start",Start);
    Stimulus->setMatrix("End",End);
    Stimulus->setMatrix("Duration",Duration);
    Stimulus->setMatrix("FieldDensity",FieldDensity);
    Stimulus->setMatrix("UseElements",UseElements);

    StringHandle SourceFile = new String("StimSeriesFile.c ");
    Stimulus->setString("SourceFile",SourceFile);

    StringHandle Parameters = new String("");
    Stimulus->setString("Parameters",Parameters);

    // Send data downstream:
    send_output_handle("StimulusBundle",StimulusBundle,true);
  }
}

} // End namespace CardioWave

