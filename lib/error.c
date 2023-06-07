/*********************************************************************
error - error handling throughout the Gnuastro library
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Jash Shah <jash28582@gmail.com>
Contributing author(s):
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
     Pedram Ashofteh-Ardakani <pedramardakani@pm.me>
Copyright (C) 2022 Free Software Foundation, Inc.

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

#include <gnuastro/error.h>
#include <gnuastro-internal/checkset.h>

#include <error.h>




/* Print the name of the library that the error belongs to. To re-generate
   the 'case' parts, run the following command:

   ls *.h \
      | sed 's/\.h//' \
      | awk '{printf "    case GAL_ERROR_LIB_%s: return \"%s.h\"; break;\n", \
                     toupper($1), $1}'
*/
char *
gal_error_print_lib_name(int lib_code)
{
  switch(lib_code)
    {
    /*-------------------- Output of command above --------------------*/
    case GAL_ERROR_LIB_ARITHMETIC: return "arithmetic.h"; break;
    case GAL_ERROR_LIB_ARRAY: return "array.h"; break;
    case GAL_ERROR_LIB_BINARY: return "binary.h"; break;
    case GAL_ERROR_LIB_BLANK: return "blank.h"; break;
    case GAL_ERROR_LIB_BOX: return "box.h"; break;
    case GAL_ERROR_LIB_COLOR: return "color.h"; break;
    case GAL_ERROR_LIB_CONVOLVE: return "convolve.h"; break;
    case GAL_ERROR_LIB_COSMOLOGY: return "cosmology.h"; break;
    case GAL_ERROR_LIB_DATA: return "data.h"; break;
    case GAL_ERROR_LIB_DIMENSION: return "dimension.h"; break;
    case GAL_ERROR_LIB_DS9: return "ds9.h"; break;
    case GAL_ERROR_LIB_EPS: return "eps.h"; break;
    case GAL_ERROR_LIB_ERROR: return "error.h"; break;
    case GAL_ERROR_LIB_ERRORINPROGRAM: return "errorinprogram.h"; break;
    case GAL_ERROR_LIB_FIT: return "fit.h"; break;
    case GAL_ERROR_LIB_FITS: return "fits.h"; break;
    case GAL_ERROR_LIB_GIT: return "git.h"; break;
    case GAL_ERROR_LIB_INTERPOLATE: return "interpolate.h"; break;
    case GAL_ERROR_LIB_JPEG: return "jpeg.h"; break;
    case GAL_ERROR_LIB_KDTREE: return "kdtree.h"; break;
    case GAL_ERROR_LIB_LABEL: return "label.h"; break;
    case GAL_ERROR_LIB_LIST: return "list.h"; break;
    case GAL_ERROR_LIB_MATCH: return "match.h"; break;
    case GAL_ERROR_LIB_PDF: return "pdf.h"; break;
    case GAL_ERROR_LIB_PERMUTATION: return "permutation.h"; break;
    case GAL_ERROR_LIB_POINTER: return "pointer.h"; break;
    case GAL_ERROR_LIB_POLYGON: return "polygon.h"; break;
    case GAL_ERROR_LIB_POOL: return "pool.h"; break;
    case GAL_ERROR_LIB_PYTHON: return "python.h"; break;
    case GAL_ERROR_LIB_QSORT: return "qsort.h"; break;
    case GAL_ERROR_LIB_SPECLINES: return "speclines.h"; break;
    case GAL_ERROR_LIB_STATISTICS: return "statistics.h"; break;
    case GAL_ERROR_LIB_TABLE: return "table.h"; break;
    case GAL_ERROR_LIB_THREADS: return "threads.h"; break;
    case GAL_ERROR_LIB_TIFF: return "tiff.h"; break;
    case GAL_ERROR_LIB_TILE: return "tile.h"; break;
    case GAL_ERROR_LIB_TXT: return "txt.h"; break;
    case GAL_ERROR_LIB_TYPE: return "type.h"; break;
    case GAL_ERROR_LIB_UNITS: return "units.h"; break;
    case GAL_ERROR_LIB_WARP: return "warp.h"; break;
    case GAL_ERROR_LIB_WCS: return "wcs.h"; break;
    /*-----------------------------------------------------------------*/

    default:
      return "NOT-DEFINED! A bug! Please contact us at "PACKAGE_BUGREPORT;
    }
}





char *
gal_error_to_string(gal_error_t *err, int verbose)
{
  char *out, *stat=NULL;

  /* If an error is found which is NOT a warning. */
  if(err->is_warning==0) stat="[BREAKING]";
  else                   stat="[WARNING]";

  /* Print the message. */
  if(err->front_msg)
    asprintf(&out, "%s: %s: %d: %s %s", err->front_msg,
             gal_error_print_lib_name(err->lib_code),
             err->code, err->back_msg, stat);
  else
    asprintf(&out, "%s: %d: %s %s",
             gal_error_print_lib_name(err->lib_code),
             err->code, err->back_msg, stat);

  /* Return the final string. */
  return out;
}





/* Prints all the error messages in the given structure in the following
   format for verbose mode (when 'verbose!=0')

     Frontend msg: code: Backend msg: [BREAKING]

   "Frontend msg" and "[BREAKING]" is printed only if they exist. Returns
   an int denoting if number of breaking errors is more than 1, thus giving
   the user the option to EXIT_FAILURE themselves. */
int
gal_error_to_stderr_all(gal_error_t *err, int verbose)
{
  char *errstr;
  int ncritical=0;
  gal_error_t *tmperr = NULL;

  /* If error structure is empty, everything is fine (there was no error to
     report), so simply return 0. */
  if(!err) return 0;

  /* Go over each component and print the message. */
  for(tmperr = err; tmperr!=NULL; tmperr = tmperr->next)
    {
      if(err->is_warning==0) ncritical++;
      errstr=gal_error_to_string(tmperr, verbose);
      error(EXIT_SUCCESS, 0, errstr);
      free(errstr);
    }

  /* Return the number of critical errors. */
  return ncritical;
}





/* Returns 0 or 1 depending on whether a breaking error occurred on the
   last function call i.e if the error on top of the stack is a breaking
   error. */
uint8_t
gal_error_occurred(gal_error_t *err)
{
  /* Return 0 when empty or a warning, return 1 otherwise. */
  return (!err || err->is_warning) ? 0 : 1;
}





/* Allocate an error data structure based on the given parameters.
   While allocating no frontend error message should be given. The
   frontend error should only be added using gal_error_add_front_msg. */
gal_error_t *
gal_error_allocate(uint8_t lib_code, uint8_t code, char *back_msg,
                   uint8_t is_warning)
{
  gal_error_t *outerr;

  /* Allocate the space for the structure.
     We use calloc here so that the error code and is_warning flags
     are set to 0 indicating generic error type and a breaking error
     by default. */
  outerr = calloc(1, sizeof *outerr);
  if(outerr == NULL)
    error(EXIT_FAILURE, 0, "%s: %zu bytes for gal_error_t,",
          __func__, sizeof *outerr);

  /* Initialize the allocated error data */
  outerr->lib_code = lib_code;
  outerr->code = code;
  outerr->is_warning = is_warning;
  gal_checkset_allocate_copy(back_msg, &outerr->back_msg);

  /* Return the final structure. */
  return outerr;
}





/* Takes in a 32-bit integer (value of an error macro) and extracts the
   error 'lib_code', 'code' and 'is_warning' flag values. */
void
gal_error_bits_extract(uint32_t bitstring, uint8_t *lib_code,
                       uint8_t *code, uint8_t *is_warning)
{
  /* The value of an error macro is a 32-bit integer. The first(starting
     from the LSB) 8 bits denote the 'is_warning' flag status. Since the
     status is either 0 or 1, if the macro value is odd then 'is_warning'
     flag is true.

      00000000 00000000 00000000 00000000
              |      | |      | |      |
              -------  -------  -------
                |         |        |
            {lib_code} {code} {is_warning}
*/
  if(bitstring % 2 != 0) *is_warning = 1;
  else *is_warning = 0;

  /* The next 8 bits (i.e. from bit 9-16) represent the 'code'. */
  *code = (bitstring >> 8) & 255;

  /* The next 8 bits (i.e. from bit 17-24) represent the 'code'. */
  *lib_code = (bitstring >> 16) & 255;
}





/* Adds a new error to the top of the given 'err' structure given the error
   'bitstring' as the third argument. From it extract the 'code' and the
   'is_warning' flags and save them in the structure. */
void
gal_error_add_back_msg(gal_error_t **err, char *back_msg,
                       uint32_t bitstring)
{
  uint8_t code = 0;
  gal_error_t *newerr;
  uint8_t lib_code = 0;
  uint8_t is_warning = 0;

  /* If no back_msg has been provided then return NULL. */
  if (back_msg == NULL) return;

  /* Extract the flags. */
  gal_error_bits_extract(bitstring, &lib_code, &code, &is_warning);

  /* Allocate a new error to be added at the top of the error stack. */
  newerr = gal_error_allocate(lib_code, code, back_msg, is_warning);

  /* Push the new error to the top of the stack. */
  newerr->next = *err;
  *err = newerr;
}





/* Adds a frontend error message to the top error in the given 'err'
   structure. If the 'replace' flag is 1 then the front_msg of the top
   error in the given 'err' structure is replaced. */
void
gal_error_add_front_msg(gal_error_t **err, char *front_msg,
                        uint8_t replace)
{
  /* Don't do anything if error structure is empty or no message provided. */
  if (!*err || front_msg==NULL)
       return;

  /* Only allocate if an error already exists. */
  if ((*err)->front_msg && !replace)
    error(EXIT_FAILURE, 0, "%s: A frontend error message already exists "
          "for the given error %d. If you wish to replace it then pass "
          "'1' to the replace flag while calling the function.", __func__,
          (*err)->code);
  else
    gal_checkset_allocate_copy(front_msg,&(*err)->front_msg);
}





/* Reverse the errors in the list. This is needed since we are treating
   the gal_error_t structure like a stack. */
void
gal_error_reverse(gal_error_t **err)
{
  /* Structure which will store the correct/reversed order. */
  gal_error_t *correctorder = NULL;

  /* 'bistring' has to be constructed from 'code' & 'is_warning'. */
  uint32_t bitstring = 0;

  /* If error structure is empty, no need to do anything. */
  if (!*err) return;

  /* Only do the reversal if there is more than one element. */
  if( *err && (*err)->next )
    {
      while(*err!=NULL)
        {
          /* The least significant 8 bits represents the 'is_warning' flag,
             next 8 bits represent the 'code' and next significant 8 bits
             represent the library code('lib_code'). */
          bitstring = (*err)->lib_code;
          bitstring = (bitstring << 8) | (*err)->code;
          bitstring = (bitstring << 8) | (*err)->is_warning;

          /* Pop top element and add to new list */
          gal_error_add_back_msg(&correctorder,
                                 (*err)->back_msg, bitstring);
          gal_error_add_front_msg(&correctorder, (*err)->front_msg, 0);
          (*err) = (*err)->next;
        }
      *err = correctorder;
    }
}





/* Given an 'err' structure and a 'bitstring', return 1 or 0 based on
   whether an error of the given type exists within the structure. */
uint8_t
gal_error_check(gal_error_t **err, uint32_t bitstring)
{
  uint8_t code = 0;
  uint8_t lib_code = 0;
  uint8_t is_warning = -1;
  gal_error_t *tmperr = NULL;

  gal_error_bits_extract(bitstring, &lib_code, &code, &is_warning);

  for(tmperr = *err; tmperr != NULL; tmperr = tmperr->next)
    if(tmperr->code == code) return 1;

  return 0;
}





void
gal_error(gal_error_t **err, int lib_code, int error_code,
          int is_warning, char *format, ...)
{
  va_list args;
  char *errstr=NULL;
  uint32_t code=GAL_ERROR_BITSET(lib_code, error_code, is_warning);

  /* Put '...' in 'args', and give it to  */
  va_start(args, format);

  /* Allocate the error string and put it in the pointer. */
  if(vasprintf(&errstr, format, args)<0)
    errstr=gal_checkset_malloc_cat((char *)__func__,
                                   ": can not use 'vasprintf'" );

  /* Put the error string in the error structure. */
  gal_error_add_back_msg(err, errstr, code);

  /* Close the variable argument list. */
  va_end(args);
}





/* Function to call at the start of library functions. This will check the
   '*err' pointer and if NULL, it will add a standard string saying that
   the function was not executed. */
int
gal_error_exists_leave_func(gal_error_t **err, int lib_code,
                            int error_code, int is_warning,
                            const char *func)
{
  if(*err)
    {
      gal_error(err, lib_code, error_code, is_warning,
                "%s: previous %s, will not continue", func,
                (*err)->next ? "errors exist" : "error exists");
      return 1;
    }
  else return 0;
}
