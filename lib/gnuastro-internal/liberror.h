/*********************************************************************
error - similar error functions for all librarys.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2023-2023 Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef __GAL_LIBERROR_H__
#define __GAL_LIBERROR_H__

/* C++ Preparations */
#undef __BEGIN_C_DECLS
#undef __END_C_DECLS
#ifdef __cplusplus
# define __BEGIN_C_DECLS extern "C" {
# define __END_C_DECLS }
#else
# define __BEGIN_C_DECLS                /* empty */
# define __END_C_DECLS                  /* empty */
#endif
/* End of C++ preparations */

/* Actual header contants (the above were for the Pre-processor). */
__BEGIN_C_DECLS  /* From C++ preparations */


/* Library codes: To re-generate this list, run the following command:

      $ cd lib/gnuastro
      $ ls *.h \
           | sed 's/\.h//' \
           | awk '{printf "GAL_LIBERROR_CODE_%s,\n", toupper($1)}'

   You can then simply copy-paste the output below (in the specified
   region). */
enum gal_error_library_codes{
  GAL_LIBERROR_CODE_INVALID,     /* ==0: accoring to the C standard. */

  /*-------------------- Output of command above --------------------*/
  GAL_LIBERROR_CODE_ARITHMETIC,
  GAL_LIBERROR_CODE_ARRAY,
  GAL_LIBERROR_CODE_BINARY,
  GAL_LIBERROR_CODE_BLANK,
  GAL_LIBERROR_CODE_BOX,
  GAL_LIBERROR_CODE_COLOR,
  GAL_LIBERROR_CODE_CONVOLVE,
  GAL_LIBERROR_CODE_COSMOLOGY,
  GAL_LIBERROR_CODE_DATA,
  GAL_LIBERROR_CODE_DIMENSION,
  GAL_LIBERROR_CODE_DS9,
  GAL_LIBERROR_CODE_EPS,
  GAL_LIBERROR_CODE_ERROR,
  GAL_LIBERROR_CODE_ERRORINPROGRAM,
  GAL_LIBERROR_CODE_FIT,
  GAL_LIBERROR_CODE_FITS,
  GAL_LIBERROR_CODE_GIT,
  GAL_LIBERROR_CODE_INTERPOLATE,
  GAL_LIBERROR_CODE_JPEG,
  GAL_LIBERROR_CODE_KDTREE,
  GAL_LIBERROR_CODE_LABEL,
  GAL_LIBERROR_CODE_LIST,
  GAL_LIBERROR_CODE_MATCH,
  GAL_LIBERROR_CODE_PDF,
  GAL_LIBERROR_CODE_PERMUTATION,
  GAL_LIBERROR_CODE_POINTER,
  GAL_LIBERROR_CODE_POLYGON,
  GAL_LIBERROR_CODE_POOL,
  GAL_LIBERROR_CODE_PYTHON,
  GAL_LIBERROR_CODE_QSORT,
  GAL_LIBERROR_CODE_SPECLINES,
  GAL_LIBERROR_CODE_STATISTICS,
  GAL_LIBERROR_CODE_TABLE,
  GAL_LIBERROR_CODE_THREADS,
  GAL_LIBERROR_CODE_TIFF,
  GAL_LIBERROR_CODE_TILE,
  GAL_LIBERROR_CODE_TXT,
  GAL_LIBERROR_CODE_TYPE,
  GAL_LIBERROR_CODE_UNITS,
  GAL_LIBERROR_CODE_WARP,
  GAL_LIBERROR_CODE_WCS,
  /*-----------------------------------------------------------------*/

  GAL_LIBERROR_CODE_NUMLIBS /* Total number of libraies */
};





/* This macros contains the contents of the 'LIB_error()' function (where
   'LIB' is the name of the library). This is essentially a wrapper to
   avoid having to repeat the long library code macro within the analysis
   functions. It is defined as a macro to simplify the variable-arguments
   (va) interface (which needs to be within the function and is error-prone
   if we have to do it for all the library functions. */
#define GAL_LIBERROR_ADD_CONTENTS(LIBCODE) {                            \
    int status;                                                         \
    va_list args;                                                       \
    va_start(args, format);                                             \
    status=gal_error_add_va(err, LIBCODE, code,                         \
                            is_warning, format, args);                  \
    va_end(args);                                                       \
    return status;                                                      \
  }

__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_LIBERROR_H__ */
