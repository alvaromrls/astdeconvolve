/*********************************************************************
error - error handling throughout the Gnuastro library
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Jash Shah <jash28582@gmail.com>
Contributing author(s):
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
     Labeeb Asari <asari.r.labeeb7@@gmail.com>
     Pedram Ashofteh-Ardakani <pedramardakani@pm.me>
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
#ifndef __GAL_ERROR_H__
#define __GAL_ERROR_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <stdint.h>
#include <stdarg.h>
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





/* Actual header contants (the above were for the Pre-processor). */
__BEGIN_C_DECLS  /* From C++ preparations */





/* Given the `lib_code`,`code` and the `is_warning` flag of an error,
   returns the value whose least significant 8 bits represents the
   `is_warning` flag, next 8 bits represent the `code` and next
   significant 8 bits represent the library code(`lib_code`).


                        ┌──────────────────┐
                        │                  │
                        │32 Bit Macro Value│
                        │                  │
                        └─────────┬────────┘
                                  │
                                  │
             ┌────────────────────┼───────────────────┐
             │                    │                   │
        Bits 16-25           Bits 8-15           Bits 0-7
             │                    │                   │
     ┌───────▼────────┐  ┌────────▼────────┐ ┌────────▼───────┐
     │   lib_code     │  │      code       │ │   is_warning   │
     │                │  │                 │ │                │
     │   0000 0000    │  │    0000 0000    │ │    0000 0000   │
     └────────────────┘  └─────────────────┘ └────────────────┘
*/
#define GAL_ERROR_BITSET(lib_code, code, is_warning) ((lib_code << 16) | (code << 8) | is_warning)







/************************************************************
 **************        Error Structure        ***************
 ************************************************************/
/* Data type for storing errors */
typedef struct gal_error_t
{
  uint8_t code;             /* From 'GAL_ERROR_CODE_*' defined below.    */
  uint8_t lib_code;         /* Library which created the error.          */
  uint8_t is_warning;       /* Defines if the error is only a warning.   */
  char *back_msg;           /* Message of backend (library).             */
  char *front_msg;          /* Message of front end (caller of library). */
  struct gal_error_t *next; /* Next error message.                       */
} gal_error_t;





/* Library codes: To re-generate this list, run the following command:

      $ cd lib/gnuastro
      $ ls *.h \
           | sed 's/\.h//' \
           | awk '{printf "GAL_ERROR_LIB_%s,\n", toupper($1)}'

   You can then simply copy-paste the output below (in the specified
   region). */
enum gal_error_library_codes{
  GAL_ERROR_LIB_INVALID,     /* ==0: accoring to the C standard. */

  /*-------------------- Output of command above --------------------*/
  GAL_ERROR_LIB_ARITHMETIC,
  GAL_ERROR_LIB_ARRAY,
  GAL_ERROR_LIB_BINARY,
  GAL_ERROR_LIB_BLANK,
  GAL_ERROR_LIB_BOX,
  GAL_ERROR_LIB_COLOR,
  GAL_ERROR_LIB_CONVOLVE,
  GAL_ERROR_LIB_COSMOLOGY,
  GAL_ERROR_LIB_DATA,
  GAL_ERROR_LIB_DIMENSION,
  GAL_ERROR_LIB_DS9,
  GAL_ERROR_LIB_EPS,
  GAL_ERROR_LIB_ERROR,
  GAL_ERROR_LIB_ERRORINPROGRAM,
  GAL_ERROR_LIB_FIT,
  GAL_ERROR_LIB_FITS,
  GAL_ERROR_LIB_GIT,
  GAL_ERROR_LIB_INTERPOLATE,
  GAL_ERROR_LIB_JPEG,
  GAL_ERROR_LIB_KDTREE,
  GAL_ERROR_LIB_LABEL,
  GAL_ERROR_LIB_LIST,
  GAL_ERROR_LIB_MATCH,
  GAL_ERROR_LIB_PDF,
  GAL_ERROR_LIB_PERMUTATION,
  GAL_ERROR_LIB_POINTER,
  GAL_ERROR_LIB_POLYGON,
  GAL_ERROR_LIB_POOL,
  GAL_ERROR_LIB_PYTHON,
  GAL_ERROR_LIB_QSORT,
  GAL_ERROR_LIB_SPECLINES,
  GAL_ERROR_LIB_STATISTICS,
  GAL_ERROR_LIB_TABLE,
  GAL_ERROR_LIB_THREADS,
  GAL_ERROR_LIB_TIFF,
  GAL_ERROR_LIB_TILE,
  GAL_ERROR_LIB_TXT,
  GAL_ERROR_LIB_TYPE,
  GAL_ERROR_LIB_UNITS,
  GAL_ERROR_LIB_WARP,
  GAL_ERROR_LIB_WCS,
  /*-----------------------------------------------------------------*/

  GAL_ERROR_LIB_NUMLIBS /* Total number of libraies */
};





/* The error codes allow high-level classification of the various
   situations and should be specified when adding a new error on the
   stack. The low-level details of each error are given in their
   message/description.

   Gnuastry's error codes are primarily taken from the relevant GNU C
   Library's error codes that you can access from [1] or [2] below. The GNU
   C Library documentation contains a more complete description of each
   error type. Note that we are a little more loose in interpreting the
   errors compared to the GNU C library because the context is larger here
   (data analysis).

   [1] This command: $ info libc "error codes"
   [2] https://www.gnu.org/software/libc/manual/html_node/Error-Codes.html

   If a necessary error code is not present in the GNU C Library, we can
   add new codes here. The source that the code was inspired from should
   just be clear in the description with the abbreviations below.

     [CU]: Custom (not from any known source).
     [PY]: https://docs.python.org/3/library/exceptions.html
   */
enum gal_error_codes{

  /* Invalid (should be first!). */
  GAL_ERROR_CODE_INVALID,       /* ==0 in the C standard.                */


  /* File/directory or Input/output issues. */
  GAL_ERROR_CODE_EIO,           /* Generic I/O (only when not in below). */
  GAL_ERROR_CODE_EACCESS,       /* Cannot access the given location.     */
  GAL_ERROR_CODE_ENOENT,        /* No such file or directory.            */
  GAL_ERROR_CODE_ENXIO,         /* No such device or address.            */
  GAL_ERROR_CODE_EFTYPE,        /* Bad file format for operation.        */
  GAL_ERROR_CODE_EEXIST,        /* File exists (and we can't overwrite). */
  GAL_ERROR_CODE_ENOTDIR,       /* Not a directory (but expects dir).    */
  GAL_ERROR_CODE_EISDIR,        /* Is a directory (but expects file).    */
  GAL_ERROR_CODE_EFBIG,         /* File is too large.                    */
  GAL_ERROR_CODE_EOF,           /* [PY] Reached end-of-file, no content. */
  GAL_ERROR_CODE_ENOTEMPTY,     /* Directory not empty.                  */
  GAL_ERROR_CODE_ENODATA,       /* No data available in input file.      */

  /* Input values (usually checked at the start of a function). */
  GAL_ERROR_CODE_E2BIG,         /* Argument list is too long.            */
  GAL_ERROR_CODE_EINVAL,        /* Invalid argument (if not in below).   */
  GAL_ERROR_CODE_EDOM,          /* Numerical input value out of range.   */
  GAL_ERROR_CODE_INDEX,         /* [PY] An array index is out of range.  */
  GAL_ERROR_CODE_ENAMETOOLONG,  /* Given name is too long.               */
  GAL_ERROR_CODE_NAME,          /* [PY] Given name is not found.         */
  GAL_ERROR_CODE_TYPE,          /* [PY] Given type is not expected.      */

  /* Requested operation. */
  GAL_ERROR_CODE_EPERM,         /* Requested operation not permitted.    */
  GAL_ERROR_CODE_ENOTSUPP,      /* Requested operation not supported.    */
  GAL_ERROR_CODE_ENOSYS,        /* Requested operation not implemented.  */
  GAL_ERROR_CODE_ESRCH,         /* No such process/function.             */
  GAL_ERROR_CODE_EGREGIOUS,     /* The requested operation is not clear. */
  GAL_ERROR_CODE_ENOPKG,        /* Necessary lib (package) not installed.*/

  /* Output values. */
  GAL_ERROR_CODE_ZERODIVISION,  /* Division or modulo by zero, all types.*/
  GAL_ERROR_CODE_ERANGE,        /* Numerical output value out of range.  */
  GAL_ERROR_CODE_EOVERFLOW,     /* Output value has overflowed.          */

  /* Operational errors (in the middle of the function). */
  GAL_ERROR_CODE_ENOMEM,        /* Cannot allocate memory.               */
  GAL_ERROR_CODE_ETIMEDOUT,     /* Operation has taken too long.         */
  GAL_ERROR_CODE_RECURSION,     /* [PY] Maximum depth of recursion.      */
  GAL_ERROR_CODE_SYSTEMEXIT,    /* [PY] system()' function crashed.      */

  /* External interruptions. */
  GAL_ERROR_CODE_EINTR,         /* Interrupted system call or by signal. */
  GAL_ERROR_CODE_ENETDOWN,      /* Network is down.                      */
  GAL_ERROR_CODE_ENETUNREACH,   /* Network is not reachable.             */
  GAL_ERROR_CODE_KEYBOARD,      /* [PY] Keyboard interrupt, e.g., Ctrl+C)*/


  /* Total number of Gnuastro error types (should be last!) */
  GAL_ERROR_CODE_NTYPES
};










/****************************************************************
 ************************   Allocation   ************************
 ****************************************************************/
gal_error_t *
gal_error_allocate(uint8_t lib_code, uint8_t code, char *back_msg,
                   uint8_t is_warning);

void
gal_error_add_back_msg(gal_error_t **err, char *back_msg,
                       uint32_t macro_val);

void
gal_error_add_front_msg(gal_error_t **err, char *front_msg,
                        uint8_t replace);

void
gal_error(gal_error_t **err, int lib_code, int error_code,
          int is_warning, char *format, ...);

void
gal_error_reverse(gal_error_t **err);

/****************************************************************
 *************************   Checking   *************************
 ****************************************************************/
uint8_t
gal_error_check(gal_error_t **err, uint32_t macro_val);

int
gal_error_exists_leave_func(gal_error_t **err, int lib_code,
                            int error_code, int is_warning,
                            const char *func);

void
gal_error_parse_macro(uint32_t macro_val, uint8_t *lib_code, uint8_t *code,
                      uint8_t *is_warning);

uint8_t
gal_error_occurred(gal_error_t *err);


/****************************************************************
 *************************   Priting   **************************
 ****************************************************************/
char *
gal_error_to_string(gal_error_t *err, int verbose);

int
gal_error_to_stderr_all(gal_error_t *err, int verbose);

__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_ERROR_H__ */
