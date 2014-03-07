/*
 *  ComputeActivation.cc:
 *
 *  Written by:
 *   jeroen
 *   TODAY'S DATE HERE
 *
 */

#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/MatrixPort.h>

#include <float.h>

namespace CardioWaveInterface {

using namespace SCIRun;

class ComputeActivation : public Module {
public:
  ComputeActivation(GuiContext*);

  virtual void execute();

private:
  vector<double> prev_;
  
  vector<double> vmax_;
  vector<double> dvdtmax_;
  vector<double> activation_;


  vector<double> threshold_60_;
  vector<double> threshold_40_;
  vector<double> threshold_20_;
  vector<double> threshold0_;
  vector<double> threshold20_;
  
  double            old_time_;
  Matrix::size_type data_size_;
};


DECLARE_MAKER(ComputeActivation)

ComputeActivation::ComputeActivation(GuiContext* ctx) :
  Module("ComputeActivation", ctx, Source, "Tools", "CardioWaveInterface")
{
  old_time_ = DBL_MAX;
}

void
ComputeActivation::execute()
{
  MatrixHandle Data;
  MatrixHandle Time;
  
  get_input_handle("Data",Data,true);
  get_input_handle("Time",Time,true);
  
  if (Data.get_rep() == 0  || Time.get_rep() == 0)
  {
    error("ComputeActivation: Invalid inputs");
    return;
  }

  double newtime;
  
  if (Time->ncols() != 1 || Time->nrows() != 1)
  {
    error("ComputeActivation: Time needs to be a scalar");
    return;
  }

  double* data;
  Matrix::size_type size;
  
  data = Time->get_data_pointer();
  size = Time->get_data_size();
  
  if (size == 0)
  {
    error("No time value specified");
    return;
  }

  newtime = data[0];
  
  if (newtime == old_time_ && (data_size_ == Data->get_data_size())) return;
  
  if (newtime < old_time_ || (data_size_ != Data->get_data_size()))
  {
    // reset module
  
    old_time_ = newtime;
    data_size_ = Data->get_data_size();
    data = Data->get_data_pointer();

    // initialize
    vmax_.resize(data_size_);
    for (Matrix::size_type j=0; j < data_size_; j++ ) vmax_[j] = -90.0;
    
    // initialize
    dvdtmax_.resize(data_size_);
    for (Matrix::size_type j=0; j < data_size_; j++ ) dvdtmax_[j] = 0.0;
    
    // initialize
    activation_.resize(data_size_);
    for (Matrix::size_type j=0; j < data_size_; j++ ) activation_[j] = 0.0;

    // initialize
    prev_.resize(data_size_);
    for (Matrix::size_type j=0; j < data_size_; j++ ) prev_[j] = data[j];    
    
    
    //thresholds
    threshold_60_.resize(data_size_);
    for (Matrix::size_type j=0; j < data_size_; j++ ) threshold_60_[j] = 0.0; 

    threshold_40_.resize(data_size_);
    for (Matrix::size_type j=0; j < data_size_; j++ ) threshold_40_[j] = 0.0; 

    threshold_20_.resize(data_size_);
    for (Matrix::size_type j=0; j < data_size_; j++ ) threshold_20_[j] = 0.0; 
    
    threshold0_.resize(data_size_);
    for (Matrix::size_type j=0; j < data_size_; j++ ) threshold0_[j] = 0.0; 

    threshold20_.resize(data_size_);
    for (Matrix::size_type j=0; j < data_size_; j++ ) threshold20_[j] = 0.0; 

  }
  else
  {
    data = Data->get_data_pointer();

    if (Data->get_data_size() != data_size_)
    {
      error("input vector of incorrect length");
      return;
    }
  }
  
  // Compute Vmax
  for (Matrix::size_type j=0; j < data_size_; j++ ) if (vmax_[j] < data[j]) vmax_[j] = data[j];
  
  // Compute dV/dtmax and activation
  for (Matrix::size_type j=0; j < data_size_; j++ )
  {
    double dvdt = (data[j]-prev_[j])/(newtime-old_time_);
    if (dvdtmax_[j] <  dvdt) 
    {
      dvdtmax_[j] = dvdt;
      activation_[j] = newtime;
    }
  }
  
  // Compute thresholds
  for (Matrix::size_type j=0; j < data_size_; j++ )
  {
    if (threshold_60_[j] == 0.0 && data[j] > -60.0) threshold_60_[j] = newtime;
  }

  for (Matrix::size_type j=0; j < data_size_; j++ )
  {
    if (threshold_40_[j] == 0.0 && data[j] > -40.0) threshold_40_[j] = newtime;
  }

  for (Matrix::size_type j=0; j < data_size_; j++ )
  {
    if (threshold_20_[j] == 0.0 && data[j] > -20.0) threshold_20_[j] = newtime;
  }
  
  for (Matrix::size_type j=0; j < data_size_; j++ )
  {
    if (threshold0_[j] == 0.0 && data[j] > 0.0) threshold0_[j] = newtime;
  }
  
  for (Matrix::size_type j=0; j < data_size_; j++ )
  {
    if (threshold20_[j] == 0.0 && data[j] > 20.0) threshold20_[j] = newtime;
  }
  
  
  for (Matrix::size_type j=0; j < data_size_; j++ )
  { 
    prev_[j] = data[j];
  }
  old_time_ = newtime;


  MatrixHandle Activation = new DenseMatrix(data_size_,1);
  if (Activation.get_rep() == 0)
  {
    error("Could not allocate Activation Matrix");
  }

  MatrixHandle Vmax = new DenseMatrix(data_size_,1);
  if (Vmax.get_rep() == 0)
  {
    error("Could not allocate Vmax Matrix");
  }

  MatrixHandle Dvdtmax = new DenseMatrix(data_size_,1);
  if (Dvdtmax.get_rep() == 0)
  {
    error("Could not allocate Vmax Matrix");
  }

  MatrixHandle Threshold = new DenseMatrix(data_size_,5);
  if (Dvdtmax.get_rep() == 0)
  {
    error("Could not allocate Vmax Matrix");
  }

  double *dst = Activation->get_data_pointer();
  for (Matrix::size_type j=0; j < data_size_; j++ )
  {
    dst[j] =activation_[j];
  }

  dst = Vmax->get_data_pointer();
  for (Matrix::size_type j=0; j < data_size_; j++ )
  {
    dst[j] =vmax_[j];
  }

  dst = Dvdtmax->get_data_pointer();
  for (Matrix::size_type j=0; j < data_size_; j++ )
  {
    dst[j] =dvdtmax_[j];
  }

  size_type k = 0;
  dst = Threshold->get_data_pointer();
  for (Matrix::size_type j=0; j < data_size_; j++ )
  {
    dst[k] = threshold_60_[j]; k++;
    dst[k] = threshold_40_[j]; k++;
    dst[k] = threshold_20_[j]; k++;
    dst[k] = threshold0_[j]; k++;
    dst[k] = threshold20_[j]; k++;
  }

  send_output_handle("Activation",Activation,true);
  send_output_handle("Vmax",Vmax,true);
  send_output_handle("dVdtmax",Dvdtmax,true);
  send_output_handle("Threshold",Threshold,true);

}


} // End namespace CardioWaveInterface


