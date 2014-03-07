#undef jpeg_SHARE

#if defined(_WIN32) && defined(_USRDLL)
# ifdef BUILD_jpeg
#  define jpeg_SHARE __declspec(dllexport)
# else
#  define jpeg_SHARE __declspec(dllimport)
# endif
#else
#  define jpeg_SHARE
#endif
