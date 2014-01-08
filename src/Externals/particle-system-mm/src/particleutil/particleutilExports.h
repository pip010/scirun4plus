#undef Particleutil_SHARE

#if defined(_WIN32) && defined(_USRDLL)
#  ifdef BUILD_Particleutil
#    define Particleutil_SHARE __declspec(dllexport)
#  else
#    define Particleutil_SHARE __declspec(dllimport)
#  endif
#else
#  define Particleutil_SHARE
#endif