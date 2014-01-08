#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2010 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software. 
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.

message(STATUS "Configure libf2c library...")

set(libf2c_srcs_MISC
  #f77vers.c
  #i77vers.c
  #main.c
  #s_rnge.c
  #abort_.c
  exit_.c
  #getarg_.c
  #iargc_.c
  #getenv_.c
  #signal_.c
  #s_stop.c
  #s_paus.c
  #system_.c
  #cabs.c
  ctype.c
  #derf_.c
  #derfc_.c
  #erf_.c
  #erfc_.c
  sig_die.c
  #uninit.c
)

set(libf2c_srcs_POW
  pow_ci.c
  pow_dd.c
  pow_di.c
  pow_hh.c
  pow_ii.c
  pow_ri.c
  pow_zi.c
  pow_zz.c)

set(libf2c_srcs_CX
  c_abs.c
  c_cos.c
  c_div.c
  c_exp.c
  c_log.c
  c_sin.c
  c_sqrt.c)

set(libf2c_srcs_DCX
  z_abs.c
  z_cos.c
  z_div.c
  z_exp.c
  z_log.c
  z_sin.c
  z_sqrt.c)

set(libf2c_srcs_REAL
  r_abs.c
  r_acos.c
  r_asin.c
  r_atan.c
  r_atn2.c
  r_cnjg.c
  r_cos.c
  r_cosh.c
  r_dim.c
  r_exp.c
  r_imag.c
  r_int.c
  r_lg10.c
  r_log.c
  r_mod.c
  r_nint.c
  r_sign.c
  r_sin.c
  r_sinh.c
  r_sqrt.c
  r_tan.c
  r_tanh.c)

set(libf2c_srcs_DBL
  d_abs.c
  d_acos.c
  d_asin.c
  d_atan.c
  d_atn2.c
  d_cnjg.c
  d_cos.c
  d_cosh.c
  d_dim.c
  d_exp.c
  d_imag.c
  d_int.c
  d_lg10.c
  d_log.c
  d_mod.c
  d_nint.c
  d_prod.c
  d_sign.c
  d_sin.c
  d_sinh.c
  d_sqrt.c
  d_tan.c
  d_tanh.c)

set(libf2c_srcs_INT
  i_abs.c
  i_dim.c
  i_dnnt.c
  i_indx.c
  i_len.c
  i_len_trim.c
  i_mod.c
  i_nint.c
  i_sign.c
  lbitbits.c
  lbitshft.c
  i_ceiling.c)

set(libf2c_srcs_HALF
  h_abs.c
  h_dim.c
  h_dnnt.c
  h_indx.c
  h_len.c
  h_mod.c
  h_nint.c
  h_sign.c)

set(libf2c_srcs_CMP
  l_ge.c
  l_gt.c
  l_le.c
  l_lt.c
  hl_ge.c
  hl_gt.c
  hl_le.c
  hl_lt.c)

set(libf2c_srcs_EFL
  ef1asc_.c
  ef1cmc_.c)

set(libf2c_srcs_CHAR
  f77_aloc.c
  s_cat.c
  s_cmp.c
  s_copy.c)

set(libf2c_srcs_I77
  backspac.c
  close.c
  dfe.c
  dolio.c
  due.c
  endfile.c
  err.c
  fmt.c
  fmtlib.c
  ftell_.c
  iio.c
  ilnw.c
  inquire.c
  lread.c
  lwrite.c
  open.c
  rdfmt.c
  rewind.c
  rsfe.c
  rsli.c
  rsne.c
  sfe.c
  sue.c
  typesize.c
  uio.c
  util.c
  wref.c
  wrtfmt.c
  wsfe.c
  wsle.c
  wsne.c
  xwsne.c)

set(libf2c_srcs_QINT
  pow_qq.c
  qbitbits.c
  qbitshft.c
  ftell64_.c)

set(libf2c_srcs_TIME
  dtime_.c
  etime_.c)

# For INTEGER*8 support (which requires system-dependent adjustments to
# f2c.h), add ${libf2c_srcs_QINT} library
add_library(${SCI_LIBF2C_LIBRARY} STATIC
  ${libf2c_srcs_MISC}
  ${libf2c_srcs_POW}
  ${libf2c_srcs_CX}
  ${libf2c_srcs_DCX}
  ${libf2c_srcs_REAL}
  ${libf2c_srcs_DBL}
  ${libf2c_srcs_INT}
  ${libf2c_srcs_HALF}
  ${libf2c_srcs_CMP}
  ${libf2c_srcs_EFL}
  ${libf2c_srcs_CHAR}
  ${libf2c_srcs_I77}
  ${libf2c_srcs_TIME})

add_definitions("-DSkip_f2c_Undefs -DNO_ISATTY")

check_function_exists("onexit" HAVE_ONEXIT)
if(NOT HAVE_ONEXIT)
  add_definitions("-DNO_ONEXIT")
endif(NOT HAVE_ONEXIT)

if(UNIX AND NOT APPLE AND NOT CYGWIN)
  set_target_properties(${SCI_LIBF2C_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC")
endif(UNIX AND NOT APPLE AND NOT CYGWIN)
