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
    asprintf(&out, "%s [code %d; %s]: %s", err->func,
             err->code, stat, err->message);
  else
    asprintf(&out, "%s: %s", err->func, err->message);

  /* Return the final string. */
  return out;
}





/* Prints all the error messages in the given structure to the standard
   error. It returns the number of breaking errors that were found, thus
   giving the caller the option to 'EXIT_FAILURE' if necessary. */
int
gal_error_write_all_stderr_reverse(gal_error_t **err, int verbose)
{
  char *str;
  int ncritical=0;
  gal_error_t *tmp=NULL;

  /* If error structure is empty, everything is fine (there was no error to
     report), so simply return 0. */
  if(*err==NULL) return 0;

  /* Reverse the errors */
  gal_error_reverse(err);

  /* Go over each component and print the message. */
  for(tmp=*err; tmp!=NULL; tmp=tmp->next)
    {
      if(tmp->is_warning==0) ncritical++;
      str=gal_error_write_string(tmp, verbose);
      error(EXIT_SUCCESS, 0, str);
      free(str);
    }

  /* Return the number of critical errors. */
  return ncritical;
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





void
gal_error_reverse(gal_error_t **err)
{
  gal_error_t *tmp, *out=NULL;

  if( *err && (*err)->next )
    {
      /* Parse the input and add them to the 'out' list. */
      for(tmp=*err; tmp!=NULL; tmp=tmp->next)
        gal_error_add(&out, tmp->code, tmp->is_warning,
                      tmp->func, "%s", tmp->message);

      /* Free the input.  */
      gal_error_free(*err);

      /* Put the output in the pointer of the input. */
      *err=out;
    }
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
gal_error_has_breaking(gal_error_t *err)
{
  gal_error_t *tmp;

  /* Parse through the list of errors and return 1 if any are breaking. */
  for(tmp=err; tmp!=NULL; tmp=tmp->next)
    if(tmp->is_warning==0) return 1;

  /* If we got here, then there was no breaking errors. */
  return 0;
}
