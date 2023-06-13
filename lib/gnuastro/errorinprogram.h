/*********************************************************************
Function to abort programs with a complete error message.

!!! ONLY INCLUDE IN COMPILED PROGRAMS, NO THE LIBRARY !!!

Original author:
     Jash Shah <jash28582@gmail.com>
Contributing author(s):
     Fathma Mehnoor <fathmamehnoor@gmail.com>
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Copyright (C) 2015-2023 Free Software Foundation, Inc.

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
#ifndef __ERRORINPROGRAM_H__
#define __ERRORINPROGRAM_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <gnuastro/error.h>

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





/* This function will call the 'error(EXIT_FAILURE, ...)' function, which
   should not be placed in the 'libgnuastro' library file. This is because
   the library should not contain any calls to 'exit(EXIT_FAILURE,
   ...)'. Only the executable programs should contain this. Therefore this
   function's full definition is in this header file and the user'sSo it is
   loaded into each source file that needs it separately. */
void
gal_errorinprogram(gal_error_t *error, int verbose)
{
  /* Find the last breaking error code. */
  int failcode=gal_error_breaking_last_code(error);

  /* If there is an error, this function will just print the reversed
     list, it will keep the input untouched. */
  if(gal_error_write_all_stderr_reverse(error, verbose))
    {
      gal_error_free(error);
      exit(failcode);
    }
}





__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_ERRORPROGRAM_H__ */
