#include <iostream>

//VISPACK INCLUDES
#include <mat/matrix.h>
#include <vol/volume.h>
#include <vol/volumefile.h>
#include <volutils.h>

// TEEM INCLUDES
#include <teem/nrrd.h>

using namespace std;

int main(int argc, char **argv)
{
  if ( argc != 2 )
  {
    cout << "Usage: <file>"
         << endl;
    exit( 0 );
  }

  string filename(argv[1]);
  VISVolumeFile vf;
  VISVolume<float> in = VISVolume<float>(vf.read(filename.c_str()));  
  
  cout << "min: " << in.min() << endl;
  cout << "max: " << in.max() << endl;

  if (writeScalarVolumeFile("out.nrrd", in)) {
    cerr << "Succesfully converted. wrote out.nrrd" << endl;
  } else {
    cerr << "writing nrrd failed." << endl;
  }

  return 0;
}
