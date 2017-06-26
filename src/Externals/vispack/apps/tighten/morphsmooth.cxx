#include <volutils.h>
#include <volumefile.h>
#include <morphology.h>

int main(int argc, char** argv)
{

  if (argc < 3)
    {
      cout << "usage: morphsmooth radius infile outfile" << endl;
      exit(-1);
    }

  float radius;
	float vol_min, vol_max;
  int morph_size;
  VolumeScalar vol_in, vol_out;
  VISImageFile imfile;
  VISVolumeFile volfile;
  int args = 1;
  sscanf(argv[args++], "%f", &radius);

	cout << "radius is " << radius << endl;

  VolScale scale;

  vol_in = readScalarVolumeFile(argv[args++],scale);

	cout << "min input " << (vol_min = vol_in.min()) << " and max " << (vol_max = vol_in.max()) << endl;
  cout << "adjust " << (2.0f*(vol_min/(vol_max - vol_min)) + 1.0f) << endl;
  cout << "mult " << (2.0f/(vol_max - vol_min)) << endl;

	vol_in = vol_in*(2.0f/(vol_max - vol_min)) - (2.0f*(vol_min/(vol_max - vol_min)) + 1.0f);

  cout << "min input " << (vol_min = vol_in.min()) << " and max " << (vol_max = vol_in.max()) << endl;

  tighten(vol_in, radius);

  writeScalarVolumeFile(argv[args++], vol_in, scale);

  cout << "min out " << vol_in.min() << " and max " << vol_in.max() << endl;
  cout << "done" << endl;
}
