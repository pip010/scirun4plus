#undef Features_SHARE

#if defined(_WIN32) && defined(_USRDLL)
# ifdef BUILD_features
#  define Features_SHARE __declspec(dllexport)
# else
#  define Features_SHARE __declspec(dllimport)
# endif
#else
# define Features_SHARE
#endif
