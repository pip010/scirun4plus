#undef System_SHARE

#if defined(_WIN32) && defined(_USRDLL)
# ifdef BUILD_system
#  define System_SHARE __declspec(dllexport)
# else
#  define System_SHARE __declspec(dllimport)
# endif
#else
# define System_SHARE
#endif