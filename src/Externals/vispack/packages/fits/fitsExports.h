#undef fitsio_SHARE

#if defined(_WIN32) && defined(_USRDLL)
# ifdef BUILD_fitsio
#  define fitsio_SHARE __declspec(dllexport)
# else
#  define fitsio_SHARE __declspec(dllimport)
# endif
#else
#  define fitsio_SHARE
#endif
