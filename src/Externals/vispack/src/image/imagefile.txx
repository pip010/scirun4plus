#include <stdio.h>

template <class T>
int write_raw(const char* name, const VISImage<T>& im)
{
    FILE* out_file;
    char filename[80];
    sprintf(filename, "%s", name);
    if ((out_file = fopen(filename, "w")) == NULL)
	{
	    printf("ERROR: output file open failed\n");
	    return(-1);
	}
    fwrite((const void*)(im.rep()->buffer()), sizeof(T),
	   im.height()*im.width(), out_file);
    return(0);
}

template <class T>
int read_raw(const char* name, VISImage<T>& im)
{
    FILE* out_file;
    char filename[80];
    if ((out_file = fopen(name, "r")) == NULL)
	{
	    printf("ERROR: output file open failed\n");
	    return(-1);
	}
    fread((void*)(im.repRef()->bufferRef()), sizeof(T),
	   im.height()*im.width(), 
	   out_file);
    return(0);
}

template <class T>
VISImage<T> convertBL(const VISImage<T> &im)
{
  int w = im.width(), h = im.height();
  int j, k;
  VISImage<T> r(w, h);
  const char *buf_in = (const char*)(im.rep())->buffer();
  char *buf_out  = (char*)(r.repRef())->bufferRef();

  for (k = 0; k < w*h; k++)
    {
      for (j = 0; j < sizeof(T); j++)
	buf_out[(sizeof(T) - 1) - j] = buf_in[j];
      buf_out += sizeof(T);
      buf_in += sizeof(T);
    }

  return(r);
}
