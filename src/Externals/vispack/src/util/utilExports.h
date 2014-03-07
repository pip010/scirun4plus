#undef util_SHARE

#if defined(_WIN32) && defined(_USRDLL)
# ifdef BUILD_util
#  define util_SHARE __declspec(dllexport)
# else
#  define util_SHARE __declspec(dllimport)
# endif
#else
#  define util_SHARE
#endif
