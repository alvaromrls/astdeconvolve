/*********************************************************************
error - error handling throughout the Gnuastro library
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     2022-2022 Jash Shah <jash28582@gmail.com>
     2022-2023 Mohammad Akhlaghi <mohammad@akhlaghi.org>
     2022-2022 Pedram Ashofteh-Ardakani <pedramardakani@pm.me>
Copyright (C) 2022-2023 Free Software Foundation, Inc.

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

#include <gnuastro-internal/liberror.h>


/* Print the name of the library that the error belongs to. To re-generate
   the 'case' parts, run the following command:

   ls *.h \
      | sed 's/\.h//' \
      | awk '{printf "    case GAL_LIBERROR_CODE_%s: return \"%s.h\"; break;\n", \
                     toupper($1), $1}'
*/
char *
gal_error_write_lib_name(int lib_code)
{
  switch(lib_code)
    {
    /*-------------------- Output of command above --------------------*/
    case GAL_LIBERROR_CODE_ARITHMETIC: return "arithmetic.h"; break;
    case GAL_LIBERROR_CODE_ARRAY: return "array.h"; break;
    case GAL_LIBERROR_CODE_BINARY: return "binary.h"; break;
    case GAL_LIBERROR_CODE_BLANK: return "blank.h"; break;
    case GAL_LIBERROR_CODE_BOX: return "box.h"; break;
    case GAL_LIBERROR_CODE_COLOR: return "color.h"; break;
    case GAL_LIBERROR_CODE_CONVOLVE: return "convolve.h"; break;
    case GAL_LIBERROR_CODE_COSMOLOGY: return "cosmology.h"; break;
    case GAL_LIBERROR_CODE_DATA: return "data.h"; break;
    case GAL_LIBERROR_CODE_DIMENSION: return "dimension.h"; break;
    case GAL_LIBERROR_CODE_DS9: return "ds9.h"; break;
    case GAL_LIBERROR_CODE_EPS: return "eps.h"; break;
    case GAL_LIBERROR_CODE_ERROR: return "error.h"; break;
    case GAL_LIBERROR_CODE_ERRORINPROGRAM: return "errorinprogram.h"; break;
    case GAL_LIBERROR_CODE_FIT: return "fit.h"; break;
    case GAL_LIBERROR_CODE_FITS: return "fits.h"; break;
    case GAL_LIBERROR_CODE_GIT: return "git.h"; break;
    case GAL_LIBERROR_CODE_INTERPOLATE: return "interpolate.h"; break;
    case GAL_LIBERROR_CODE_JPEG: return "jpeg.h"; break;
    case GAL_LIBERROR_CODE_KDTREE: return "kdtree.h"; break;
    case GAL_LIBERROR_CODE_LABEL: return "label.h"; break;
    case GAL_LIBERROR_CODE_LIST: return "list.h"; break;
    case GAL_LIBERROR_CODE_MATCH: return "match.h"; break;
    case GAL_LIBERROR_CODE_PDF: return "pdf.h"; break;
    case GAL_LIBERROR_CODE_PERMUTATION: return "permutation.h"; break;
    case GAL_LIBERROR_CODE_POINTER: return "pointer.h"; break;
    case GAL_LIBERROR_CODE_POLYGON: return "polygon.h"; break;
    case GAL_LIBERROR_CODE_POOL: return "pool.h"; break;
    case GAL_LIBERROR_CODE_PYTHON: return "python.h"; break;
    case GAL_LIBERROR_CODE_QSORT: return "qsort.h"; break;
    case GAL_LIBERROR_CODE_SPECLINES: return "speclines.h"; break;
    case GAL_LIBERROR_CODE_STATISTICS: return "statistics.h"; break;
    case GAL_LIBERROR_CODE_TABLE: return "table.h"; break;
    case GAL_LIBERROR_CODE_THREADS: return "threads.h"; break;
    case GAL_LIBERROR_CODE_TIFF: return "tiff.h"; break;
    case GAL_LIBERROR_CODE_TILE: return "tile.h"; break;
    case GAL_LIBERROR_CODE_TXT: return "txt.h"; break;
    case GAL_LIBERROR_CODE_TYPE: return "type.h"; break;
    case GAL_LIBERROR_CODE_UNITS: return "units.h"; break;
    case GAL_LIBERROR_CODE_WARP: return "warp.h"; break;
    case GAL_LIBERROR_CODE_WCS: return "wcs.h"; break;
    /*-----------------------------------------------------------------*/

    default:
      return "NOT-DEFINED! A bug! Please contact us at "PACKAGE_BUGREPORT;
    }
}





char *
gal_error_write_string(gal_error_t *err, int verbose)
{
  char *out, *stat=NULL;

  /* If an error is found which is NOT a warning. */
  if(err->is_warning==0) stat="[BREAKING]";
  else                   stat="[WARNING]";

  /* Print the message. */
  asprintf(&out, "%s: %d: %s %s",
           gal_error_write_lib_name(err->lib_code),
           err->code, err->back_msg, stat);

  /* Return the final string. */
  return out;
}





/* Prints all the error messages in the given structure to the standard
   error. It returns the number of breaking errors that were found, thus
   giving the caller the option to 'EXIT_FAILURE' if necessary. */
int
gal_error_write_all_stderr(gal_error_t *err, int verbose)
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
      errstr=gal_error_write_string(tmperr, verbose);
      error(EXIT_SUCCESS, 0, errstr);
      free(errstr);
    }

  /* Return the number of critical errors. */
  return ncritical;
}




















/****************************************************************
 ********************   New gal_error_t   ***********************
 ****************************************************************/
/* Allocate an error data structure based on the given parameters. */
static gal_error_t *
error_allocate(uint8_t code, uint8_t is_warning, uint8_t lib_code,
               char *message, int *alloc_failed)
{
  gal_error_t *out;

  /* Allocate the space for the structure. We use 'calloc' here so that the
     error code and is_warning flags are set to 0 indicating generic error
     type and a breaking error by default.  */
  out = calloc(1, sizeof *out);
  if(out) *alloc_failed=0;
  else {  *alloc_failed=1; return NULL; }

  /* Initialize the allocated error data */
  out->code = code;
  out->lib_code = lib_code;
  out->is_warning = is_warning;
  gal_checkset_allocate_copy(message, &out->back_msg);

  /* Return the final structure. */
  return out;
}





int
gal_error_add(gal_error_t **err, int code, int is_warning,
              int lib_code, char *format, ...)
{
  int status;
  va_list args;

  /* Start reading the variable arguments ("va"). */
  va_start(args, format);

  /* Add this error to the queue. */
  status=gal_error_add_va(err, lib_code, code,
                          is_warning, format, args);

  /* Close the variable arguments. */
  va_end(args);

  /* Return the error allocation status. */
  return status;
}





/* Add a new error node at the top of the error list. If this function
   returns with zero, then everything is fine. If the new error structure
   could not be allocated, this function will return with
   'GAL_ERROR_CODE_ERRNOTALLOC'. */
int
gal_error_add_va(gal_error_t **err, int code, int is_warning,
                 int lib_code, char *format, va_list args)
{
  int alloc_failed;
  char *message=NULL;
  gal_error_t *new=NULL;

  /* Allocate the error string and put it in the pointer. */
  if(vasprintf(&message, format, args)<0)
    message=gal_checkset_malloc_cat((char *)__func__,
                                    ": can not use 'vasprintf'" );

  /* Allocate the new error structure. */
  new=error_allocate(code, is_warning, lib_code, message,
                     &alloc_failed);

  /* Close the variable argument list and return the status. */
  if(alloc_failed) return GAL_ERROR_CODE_ERRNOTALLOC;
  else        { new->next=*err; *err=new; return 0; }
}





/* Function to call at the start of library functions. This will check the
   '*err' pointer and if it is NULL (empty), it will return a '0'. If the
   'err' structure was not empty a new error is added on the error stack
   with a fixed string to be clear and it will return the integer
   'GAL_ERROR_CODE_ERRLISTFULL'. */
int
gal_error_has_leave(gal_error_t **err, int lib_code, const char *func)
{
  if(*err)
    {
      gal_error_add(err, lib_code, GAL_ERROR_CODE_ERRLISTFULL, 0,
                    "%s: previous %s, will not continue", func,
                    (*err)->next ? "errors exist" : "error exists");
      return (*err)->code;
    }
  else return 0;
}
