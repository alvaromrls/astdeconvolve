/*********************************************************************
error - error handling throughout the Gnuastro library
This is part of GNU Astronomy Utilities (Gnuastro) package.

Authors:
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

#include <error.h>

#include <gnuastro/error.h>

#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/liberror.h>





/****************************************************************
 ********************       Writing       ***********************
 ****************************************************************/
char *
gal_error_write_string(gal_error_t *err, int verbose)
{
  char *out, *stat=NULL;

  /* If an error is found which is NOT a warning. */
  if(err->is_warning==0) stat="BREAKING";
  else                   stat="WARNING";

  /* Print the message. */
  if(verbose)
    {
      if(err->func)             /* Within libraries. */
        asprintf(&out, "%s [code %d; %s]: %s", err->func,
                 err->code, stat, err->message);
      else                      /* Within programs. */
        asprintf(&out, "[code %d; %s]: %s",
                 err->code, stat, err->message);
    }
  else
    {
      if(err->func)             /* Within libraries. */
        asprintf(&out, "%s: %s", err->func, err->message);
      else                      /* Within programs. */
        asprintf(&out, "%s", err->message);
    }

  /* Return the final string. */
  return out;
}





/* Prints all the error messages in the given structure to the standard
   error. It returns the number of breaking errors that were found, thus
   giving the caller the option to 'EXIT_FAILURE' if necessary. */
int
gal_error_write_all_stderr_reverse(gal_error_t *err, int verbose)
{
  char *str;
  int ncritical=0;
  gal_error_t *rev=NULL, *tmp=NULL;

  /* If error structure is empty, everything is fine (there was no error to
     report), so simply return 0. */
  if(err==NULL) return 0;

  /* Reverse the errors */
  rev=gal_error_reverse_keep_in(err);

  /* Go over each component and print the message. */
  for(tmp=rev; tmp!=NULL; tmp=tmp->next)
    {
      if(tmp->is_warning==0) ncritical++;
      if(verbose || tmp->code!=GAL_ERROR_CODE_ERRLISTFULL)
        {
          str=gal_error_write_string(tmp, verbose);
          error(EXIT_SUCCESS, 0, str);
          free(str);
        }
    }

  /* Free the reversed list (the input is untouched) and return the number
     of critical errors. */
  gal_error_free(rev);
  return ncritical;
}




/* Put all the error information into a simple table. */
#define ERROR_INFO(NAME, CAT, INFO)               \
  gal_checkset_allocate_copy(CAT,  &catarr[i]);   \
  gal_checkset_allocate_copy(NAME, &namearr[i]);  \
  gal_checkset_allocate_copy(INFO, &infoarr[i]);
void *
gal_error_code_info(void *junk)
{
  uint8_t *carr;
  char **namearr, **infoarr, **catarr;
  size_t i, nerr=GAL_ERROR_CODE_NUMCODES;

  /* Allocate the two columns. */
  gal_data_t *codes=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &nerr, NULL,
                                   0, -1, 1, "CODE", "counter",
                                   "Code of error.");
  gal_data_t *names=gal_data_alloc(NULL, GAL_TYPE_STRING, 1, &nerr, NULL,
                                   0, -1, 1, "NAME", "counter",
                                   "Code of error.");
  gal_data_t *cat=gal_data_alloc(NULL, GAL_TYPE_STRING, 1, &nerr, NULL,
                                 0, -1, 1, "CATEGORY", "info",
                                 "Category of the error.");
  gal_data_t *info=gal_data_alloc(NULL, GAL_TYPE_STRING, 1, &nerr, NULL,
                                  0, -1, 1, "DESCRIPTION", "info",
                                  "Description of the error.");

  /* Go one by one and add all the error information. */
  carr=codes->array;
  catarr=cat->array;
  infoarr=info->array;
  namearr=names->array;
  for(i=0; i<nerr; ++i)
    {
      carr[i]=i;
      switch(i)
        {
        /* Success */
        case GAL_ERROR_CODE_INVALID:
          ERROR_INFO("INVALID", "invalid", "No error (Success!)."); break;

        /* File/directory. */
        case GAL_ERROR_CODE_EIO:
          ERROR_INFO("EIO", "file-dir",
                     "Generic I/O (only when not in below)."); break;

        case GAL_ERROR_CODE_EACCESS:
          ERROR_INFO("EACCESS", "file-dir",
                     "Cannot access the given location."); break;
        case GAL_ERROR_CODE_ENOENT:
          ERROR_INFO("ENOENT", "file-dir",
                     "No such file or directory."); break;
        case GAL_ERROR_CODE_ENXIO:
          ERROR_INFO("ENXIO", "file-dir",
                     "No such device or address."); break;
        case GAL_ERROR_CODE_EFTYPE:
          ERROR_INFO("EFTYPE", "file-dir",
                     "Bad file format for operation."); break;
        case GAL_ERROR_CODE_EEXIST:
          ERROR_INFO("EEXIST", "file-dir",
                     "File exists (and we can't overwrite)."); break;
        case GAL_ERROR_CODE_ENOTDIR:
          ERROR_INFO("ENOTDIR", "file-dir",
                     "Not a directory (but operation expects dir).");break;
        case GAL_ERROR_CODE_EISDIR:
          ERROR_INFO("EISDIR", "file-dir",
                     "Is a directory (but operation expects file)."); break;
        case GAL_ERROR_CODE_EFBIG:
          ERROR_INFO("EFBIG", "file-dir","File is too large."); break;
        case GAL_ERROR_CODE_EOF:
          ERROR_INFO("EOF", "file-dir",
                     "Reached end-of-file, no content."); break;
        case GAL_ERROR_CODE_ENOTEMPTY:
          ERROR_INFO("ENOTEMPTY", "file-dir",
                     "Directory not empty."); break;
        case GAL_ERROR_CODE_ENODATA:
          ERROR_INFO("ENODATA", "file-dir",
                     "No data available in input file."); break;

        /* Input values (usually checked at the start of a function). */
        case GAL_ERROR_CODE_E2BIG:
          ERROR_INFO("E2BIG", "inputs",
                     "Argument list is too long."); break;
        case GAL_ERROR_CODE_EINVAL:
          ERROR_INFO("EINVAL", "inputs",
                     "Invalid argument (if not in below)."); break;
        case GAL_ERROR_CODE_EDOM:
          ERROR_INFO("EDOM", "inputs",
                     "Numerical input value out of range."); break;
        case GAL_ERROR_CODE_INDEX:
          ERROR_INFO("INDEX", "inputs",
                     "An array index is out of range."); break;
        case GAL_ERROR_CODE_ENAMETOOLONG:
          ERROR_INFO("ENAMETOOLONG", "inputs",
                     "Given name is too long."); break;
        case GAL_ERROR_CODE_NAME:
          ERROR_INFO("NAME", "inputs",
                     "Given name is not found."); break;
        case GAL_ERROR_CODE_TYPE:
          ERROR_INFO("TYPE", "inputs",
                     "Given type is not expected."); break;

        /* Requested operation. */
        case GAL_ERROR_CODE_EPERM:
          ERROR_INFO("EPERM", "request-op",
                     "Requested operation not permitted"); break;
        case GAL_ERROR_CODE_ENOTSUPP:
          ERROR_INFO("ENOTSUPP", "request-op",
                     "Requested operation not supported."); break;
        case GAL_ERROR_CODE_ENOSYS:
          ERROR_INFO("ENOSYS", "request-op",
                     "Requested operation not implemented."); break;
        case GAL_ERROR_CODE_ESRCH:
          ERROR_INFO("ESRCH", "request-op",
                     "No such process/function."); break;
        case GAL_ERROR_CODE_EGREGIOUS:
          ERROR_INFO("EGREGIOUS", "request-op",
                     "The requested operation is not clear."); break;
        case GAL_ERROR_CODE_ENOPKG:
          ERROR_INFO("ENOPKG", "request-op",
                     "Necessary lib (package) not installed."); break;

        /* Output values. */
        case GAL_ERROR_CODE_ZERODIVISION:
          ERROR_INFO("ZERODIVISION", "output",
                     "Division or modulo by zero, all types."); break;
        case GAL_ERROR_CODE_ERANGE:
          ERROR_INFO("ERANGE", "output",
                     "Numerical output value out of range."); break;
        case GAL_ERROR_CODE_EOVERFLOW:
          ERROR_INFO("EOVERFLOW", "output",
                     "Output value has overflowed."); break;

        /* External interruptions. */
        case GAL_ERROR_CODE_EINTR:
          ERROR_INFO("EINTR", "external",
                     "Interrupted system call or by signal."); break;
        case GAL_ERROR_CODE_ENETDOWN:
          ERROR_INFO("ENETDOWN", "external", "Network is down."); break;
        case GAL_ERROR_CODE_ENETUNREACH:
          ERROR_INFO("ENETUNREACH", "external",
                     "Network is not reachable."); break;
        case GAL_ERROR_CODE_KEYBOARD:
          ERROR_INFO("KEYBOARD", "external",
                     "Keyboard interrupt, e.g., Ctrl+C)"); break;

        /* Operational errors (in the middle of the function). */
        case GAL_ERROR_CODE_ERRLISTFULL:
          ERROR_INFO("ERRLISTFULL", "operation",
                     "Error list not empty, not continuing."); break;
        case GAL_ERROR_CODE_ERRNOTALLOC:
          ERROR_INFO("ERRNOTALLOC", "operation",
                     "Couldn't allocate error struct."); break;
        case GAL_ERROR_CODE_BUG:
          ERROR_INFO("BUG", "operation",
                     "Unexpected situation, a bug!"); break;
        case GAL_ERROR_CODE_ENOMEM:
          ERROR_INFO("ENOMEM", "operation",
                     "Cannot allocate memory."); break;
        case GAL_ERROR_CODE_ETIMEDOUT:
          ERROR_INFO("ETIMEDOUT", "operation",
                     "Operation has taken too long."); break;
        case GAL_ERROR_CODE_RECURSION:
          ERROR_INFO("RECURSION", "operation",
                     "Maximum depth of recursion."); break;
        case GAL_ERROR_CODE_SYSTEMEXIT:
          ERROR_INFO("SYSTEMEXIT", "operation",
                     "system() function crashed."); break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
                "to fix the problem. The value '%zu' is not a recognized "
                "error code", __func__, PACKAGE_BUGREPORT, i);
        }
    }

  /* Put the information as a second column and return. */
  codes->next=names;
  names->next=cat;
  cat->next=info;
  return codes;
}


















/****************************************************************
 ********************   New gal_error_t   ***********************
 ****************************************************************/
/* Allocate an error data structure based on the given parameters. */
static gal_error_t *
error_allocate(uint8_t code, uint8_t is_warning, const char *func,
               char *message, int *alloc_failed)
{
  gal_error_t *out;

  /* Allocate the space for the structure. We use 'calloc' here so that the
     error code and is_warning flags are set to 0 indicating generic error
     type and a breaking error by default.  */
  out = calloc(1, sizeof *out);
  if(out) *alloc_failed=0;
  else {  *alloc_failed=1; return NULL; }

  /* Set the integer values. */
  out->code = code;
  out->is_warning = is_warning;

  /* The message was allocated by 'gal_error_add_va' (which calls this), so
     we don't need to re-allocate it here. */
  out->message = message;

  /* The function name is not allocated, so we need to allocate it as a
     string before continuing. */
  gal_checkset_allocate_copy(func, &out->func);

  /* Return the final structure. */
  return out;
}





int
gal_error_add(gal_error_t **err, int code, int is_warning,
              const char *func, char *template, ...)
{
  int status;
  va_list args;

  /* Start reading the variable arguments ("va"). */
  va_start(args, template);

  /* Add this error to the queue. */
  status=gal_error_add_va(err, code, is_warning, func, template, args);

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
                 const char *func, char *template, va_list args)
{
  int alloc_failed;
  char *message=NULL;
  gal_error_t *new=NULL;

  /* Allocate the error string and put it in the pointer. */
  if(vasprintf(&message, template, args)<0)
    message=gal_checkset_malloc_cat((char *)__func__,
                                    ": can not use 'vasprintf'" );

  /* Allocate the new error structure. */
  new=error_allocate(code, is_warning, func, message, &alloc_failed);

  /* Close the variable argument list and return the status. */
  if(alloc_failed) return GAL_ERROR_CODE_ERRNOTALLOC;
  else        { new->next=*err; *err=new; return 0; }
}





/* Reverse the input error list, without freeing the input. */
gal_error_t *
gal_error_reverse_keep_in(gal_error_t *err)
{
  gal_error_t *tmp, *out=NULL;
  if(err)
    for(tmp=err; tmp!=NULL; tmp=tmp->next)
      gal_error_add(&out, tmp->code, tmp->is_warning,
                    tmp->func, "%s", tmp->message);
  return out;
}





void
gal_error_free(gal_error_t *err)
{
  gal_error_t *tmp;
  while(err!=NULL)
    {
      tmp=err->next;
      free(err->func);
      free(err->message);
      free(err);
      err=tmp;
    }
}




















/****************************************************************
 ********************        Checks       ***********************
 ****************************************************************/

/* Function to call at the start of library functions. This will check the
   '*err' pointer and if it is NULL (empty), it will return a '0'. If the
   'err' structure was not empty a new error is added on the error stack
   with a fixed string to be clear and it will return the integer
   'GAL_ERROR_CODE_ERRLISTFULL'. */
int
gal_error_has_leave(gal_error_t **err, const char *func)
{
  if(*err)
    {
      gal_error_add(err, GAL_ERROR_CODE_ERRLISTFULL, 0, func,
                    "previous %s, will not continue",
                    (*err)->next ? "errors exist" : "error exists");
      return (*err)->code;
    }
  else return 0;
}





/* Return 1 if there is a breaking error. */
int
gal_error_breaking_present(gal_error_t *err)
{
  gal_error_t *tmp;

  /* Parse through the list of errors and return 1 if any are breaking. */
  for(tmp=err; tmp!=NULL; tmp=tmp->next)
    if(tmp->is_warning==0) return 1;

  /* If we got here, then there was no breaking errors. */
  return 0;
}





/* Return the last breaking error code; ignoring the ones that were due to
   an already existing error ('GAL_ERROR_CODE_ERRLISTFULL'). */
int
gal_error_breaking_last_code(gal_error_t *err)
{
  gal_error_t *tmp;

  /* Parse through the list of errors and return 1 if any are breaking. */
  for(tmp=err; tmp!=NULL; tmp=tmp->next)
    if(tmp->is_warning==0 && tmp->code!=GAL_ERROR_CODE_ERRLISTFULL)
      return tmp->code;

  /* If we got here, then there was no breaking errors. */
  return 0;
}





gal_error_t *
gal_error_breaking_last_err(gal_error_t *err)
{
  gal_error_t *tmp;

  /* Parse through the list of errors and return 1 if any are breaking. */
  for(tmp=err; tmp!=NULL; tmp=tmp->next)
    if(tmp->is_warning==0 && tmp->code!=GAL_ERROR_CODE_ERRLISTFULL)
      return tmp;

  /* If we got here, then there was no breaking errors. */
  return NULL;
}
