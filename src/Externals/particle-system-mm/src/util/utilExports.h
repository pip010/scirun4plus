#undef MM_util_SHARE

#if defined(_WIN32) && defined(_USRDLL)
# ifdef BUILD_MM_util
#  define MM_util_SHARE __declspec(dllexport)
# else
#  define MM_util_SHARE __declspec(dllimport)
# endif
#else
# define MM_util_SHARE
#endif