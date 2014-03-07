#undef tiff_SHARE

#if defined(_WIN32) && defined(_USRDLL)
# ifdef BUILD_tiff
#  define tiff_SHARE __declspec(dllexport)
# else
#  define tiff_SHARE __declspec(dllimport)
# endif
#else
#  define tiff_SHARE
#endif
