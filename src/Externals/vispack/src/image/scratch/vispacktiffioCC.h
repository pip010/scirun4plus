#include <tiff.h>
#include "vispacktiffioC.h"

extern "C"
{
int WriteTIFFRGB(const byte*, int, int, const char*, int);
int WriteTIFFfloat(const float*, int, int, int, const char*, int);
int WriteTIFFint(const int*, int, int, int, const char*, int);
int WriteTIFFbyte(const byte*, int, int, int, const char*, int);
int WriteTIFFshort(const short*, int, int, int, const char*, int);
int LoadTIFF(const char* fname, PICINFO* pinfo);
}
 
