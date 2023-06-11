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





/* This macros contains the contents of the 'LIB_error()' function (where
   'LIB' is the name of the library). This is essentially a wrapper to
   avoid having to repeat the long library code macro within the analysis
   functions. It is defined as a macro to simplify the variable-arguments
   (va) interface (which needs to be within the function and is error-prone
   if we have to do it for all the library functions. */
#define GAL_LIBERROR_ADD_CONTENTS(LIBCODE) {                            \
    int status;                                                         \
    va_list args;                                                       \
    va_start(args, template);                                           \
    status=gal_error_add_va(err, code, is_warning, LIBCODE, template,   \
                            args);                                      \
    va_end(args);                                                       \
    return status;                                                      \
  }





__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_LIBERROR_H__ */
