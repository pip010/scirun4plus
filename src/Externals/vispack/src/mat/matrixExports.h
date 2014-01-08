#undef matrix_SHARE

#if defined(_WIN32) && defined(_USRDLL)
# ifdef BUILD_matrix
#  define matrix_SHARE __declspec(dllexport)
# else
#  define matrix_SHARE __declspec(dllimport)
# endif
#else
#  define matrix_SHARE
#endif
