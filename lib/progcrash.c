/*********************************************************************
Function to abort programs with a complete error message.

!!! ONLY INCLUDE IN COMPILED PROGRAMS, NO THE LIBRARY !!!

Authors:
     2022-2022 Jash Shah <jash28582@gmail.com>
     2022-2023 Mohammad Akhlaghi <mohammad@akhlaghi.org>
     2023-2023 Fathma Mehnoor <fathmamehnoor@gmail.com>
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
#include <config.h>

#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/error.h>





/* This function will call the 'error(EXIT_FAILURE, ...)' function, which
   should not be placed in the 'libgnuastro' library file. This is because
   the library should not contain any calls to 'exit(EXIT_FAILURE,
   ...)'. Only the executable programs should contain this. Therefore this
   function's full definition is in this header file and the user'sSo it is
   loaded into each source file that needs it separately. */
void
gal_progcrash_list(gal_error_t *error, int verbose)
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





/* This function gets a single error message and error code, and will print
   it and abort. */
void
gal_progcrash_one(int code, int is_warning, const char *func,
                  int verbose, char *template, ...)
{
  int status=0;
  va_list args;
  gal_error_t *err=NULL;

  /* Start reading the variable arguments ("va"). */
  va_start(args, template);

  /* Define the error. */
  status=gal_error_add_va(&err, code, is_warning, func, template, args);

  /* If the error structure couldn't be allocated, inform the user. Since
     we couldn't actually allocate an error structure, we can't use
     ('gal_error_write_all_stderr_reverse').*/
  if(status==GAL_ERROR_CODE_ERRNOTALLOC)
    error(EXIT_FAILURE, 0, "%s: couldn't allocate 'err' structure! "
          "This should not regularly happen, unless there is severe "
          "memory/kernel consumption on your system", __func__);

  /* Print the error and abort. */
  gal_progcrash_list(err, verbose);
}
