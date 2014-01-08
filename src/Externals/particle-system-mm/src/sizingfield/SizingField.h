//----------------------------------------------------------------------//
// FILE        : SizingField.h                                           
// DESCRIPTION : This class creates the initial and final sizing fields
//               for a set of input multiple materials (or, just one
//               material). 
//----------------------------------------------------------------------//

#ifndef __SIZING_FIELD_H__
#define __SIZING_FIELD_H__

#include <mtxlib.h>
#include <multiDarrays.h>
#include <defines.h>
#include <IOScalarField.h>

class SizingField
{
public:
  SizingField(const char *file_basename, const char *volume_names,
              int main_indicator, bool useNrrdIO = false);
  SizingField(const char *file_basename, const char *volume_names,
              const char* redo_name, bool useNrrdIO = false);
  ~SizingField();

  void generateInitialSizingField();
  void generateSizingField();
	void generateConstantSizingField( float sf_value );
  void redoInitialSizingField();

  inline void delta(float d)
  { _delta = d; };
  inline float delta() const
  { return _delta; };

  inline void epsilon(float e) { _epsilon = e; };
  inline float epsilon() const { return _epsilon; };
  inline void useNrrdIO(bool useNrrd) { _useNrrdIO = useNrrd; }
  
private:
  char *_base_filename, *_volume_files, *_redo_file;
  int _main_indicator, _xdim, _ydim, _zdim;

  float _scaling[3];
  DataCenter _centers[3];

  array3D<bool> _crossing;
  array3D<float> _h0, _h;
  IOScalarField *_io_field;

  float _delta, _epsilon, _max_size;

  float _min_sf;

  void createCrossingVolume();
  void initializeh0ToMaxk();
  void compareMaxkToLFS();


  void outputNrrdFile(char* h_filename, array3D<float> & _h_pass);
  void outputVolFile(char* h_filename, array3D<float> & _h_pass);
  
   void outputNrrdFile(char* h_filename, array3D<bool> & _h_pass);
  void outputVolFile(char* h_filename, array3D<bool> & _h_pass);
  
  
  float gradLimUpdate(float grad_lim);
  bool _useNrrdIO;
  bool _restart_files;
  


};

#endif // __SIZING_FIELD_H__
