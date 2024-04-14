/*********************************************************************
Deconvolve - A minimal set of files and functions to define a program.
Deconvolve is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Alvaro Morales <alvaro96m@hotmail.com>
Contributing author(s):
Copyright (C) 2024-2025 Free Software Foundation, Inc.

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

#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <string.h>

#include <gnuastro/list.h>

#include <gnuastro/array.h>
#include <gnuastro/fits.h>
#include <gnuastro/threads.h>
#include <gnuastro/wcs.h>

#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/fixedstringmacros.h>
#include <gnuastro-internal/options.h>
#include <gnuastro-internal/timing.h>

#include "main.h"

#include "authors-cite.h"
#include "ui.h"

/**************************************************************/
/*********      Argp necessary global entities     ************/
/**************************************************************/
/* Definition parameters for the Argp: */
const char *argp_program_version = PROGRAM_STRING
    "\n" GAL_STRINGS_COPYRIGHT "\n\nWritten/developed by " PROGRAM_AUTHORS;

const char *argp_program_bug_address = PACKAGE_BUGREPORT;

static char args_doc[] = "ASTRdata";

const char doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME
    " Apply a deconvoluion given an image and a "
    "PSF.\n" GAL_STRINGS_MORE_HELP_INFO
    /* After the list of options: */
    "\v" PACKAGE_NAME " home page: " PACKAGE_URL;

/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options (struct deconvolve_params *p,
                       struct argp_option *program_options,
                       struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp = &p->cp;

  /* Set the necessary common parameters structure. */
  cp->program_struct = p;
  cp->poptions = program_options;
  cp->program_name = PROGRAM_NAME;
  cp->program_exec = PROGRAM_EXEC;
  cp->program_bibtex = PROGRAM_BIBTEX;
  cp->program_authors = PROGRAM_AUTHORS;
  cp->numthreads = gal_threads_number ();
  cp->coptions = gal_commonopts_options;

  /* Modify common options. */
  for (i = 0; !gal_options_is_last (&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch (cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_HDU:
        case GAL_OPTIONS_KEY_TYPE:
        case GAL_OPTIONS_KEY_MINMAPSIZE:
          cp->coptions[i].mandatory = GAL_OPTIONS_MANDATORY;
          break;
        }
    }
}

/* Parse a single option: */
error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  struct deconvolve_params *p = state->input;
  /* Pass 'gal_options_common_params' into the child parser. */
  state->child_inputs[0] = &p->cp;

  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then 'arg'
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if (arg && arg[0] == '=')
    argp_error (state,
                "incorrect use of the equal sign ('='). For short "
                "options, '=' should not be used and for long options, "
                "there should be no space between the option, equal sign "
                "and value");

  /* Set the key to this option. */
  switch (key)
    {
    /* Read the non-option tokens (arguments): */
    case ARGP_KEY_ARG:
      /* The user may give a shell variable that is empty! In that case
         'arg' will be an empty string! We don't want to account for such
         cases (and give a clear error that no input has been given). */
      if (p->filename)
        argp_error (state, "only one argument (input file) should be given");
      else if (arg[0] != '\0')
        p->filename = arg;
      break;

    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key (key, arg, p->cp.poptions, &p->cp);
    }

  return 0;
}

/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Check ONLY the options. When arguments are involved, do the check
   in 'ui_check_options_and_arguments'. */
static void
ui_check_only_options (struct deconvolve_params *p)
{
}

static void
ui_check_options_and_arguments (struct deconvolve_params *p)
{

  if (p->filename)
    {
      /* If input is FITS. */
      if ((p->isfits = gal_fits_file_recognized (p->filename)))
        {
          /* Check if a HDU is given. */
          if (p->cp.hdu == NULL)
            error (EXIT_FAILURE, 0,
                   "no HDU specified. When the input is a "
                   "FITS file, a HDU must also be specified, you can use "
                   "the '--hdu' ('-h') option and give it the HDU number "
                   "(starting from zero), extension name, or anything "
                   "acceptable by CFITSIO");

          /* If its an image, make sure column isn't given (in case the
             user confuses an image with a table). */
          p->hdu_type = gal_fits_hdu_format (p->filename, p->cp.hdu, "--hdu");
        }
    }

  if (p->kernelname)
    {
      /* If input is FITS. */
      if (gal_fits_file_recognized (p->kernelname))
        {
          /* Check if a HDU is given. */
          if (p->khdu == NULL)
            error (EXIT_FAILURE, 0,
                   "no HDU specified. When the kernel is a "
                   "FITS file, a HDU must also be specified, you can use "
                   "the '--khdu' ('-u') option and give it the HDU number "
                   "(starting from zero), extension name, or anything "
                   "acceptable by CFITSIO");

          /* If its an image, make sure column isn't given (in case the
             user confuses an image with a table). */
          gal_fits_hdu_format (p->kernelname, p->khdu, "--khdu");
        }
    }

  if (!strcmp (TIKHONOV_NAME, p->algorithmstr))
    {
      p->algorithm = DECONVOLUTION_ALGORITHM_TIKHONOV;
      if (p->lambda == 0)
        {
          error (EXIT_FAILURE, 0,
                 "Tikhonov algorithm requires a lambda value");
        }
    }

  else if (!strcmp (DIRECT_INVERSION_NAME, p->algorithmstr))
    {
      p->algorithm = DECONVOLUTION_ALGORITHM_DIRECT_INVERSION;
    }
  else if (!strcmp (RICHADSON_LUCY_NAME, p->algorithmstr))
    {
      p->algorithm = DECONVOLUTION_ALGORITHM_RL;
    }
  else if (!strcmp (AWMLE_NAME, p->algorithmstr))
    {
      p->algorithm = DECONVOLUTION_ALGORITHM_AWMLE;
    }
  else
    {
      error (EXIT_FAILURE, 0,
             "deconvolve algorithm not recognised (%s), "
             "please use a valid one: \n -%s \n -%s \n -%s\n",
             p->algorithmstr, TIKHONOV_NAME, DIRECT_INVERSION_NAME,
             RICHADSON_LUCY_NAME);
    }
}

/* Read the input dataset. */
static void
ui_read_input (struct deconvolve_params *p)
{
  /* If the input is a FITS image or any recognized array file format, then
     read it as an array, otherwise, as a table. */
  if (p->filename && gal_array_name_recognized (p->filename))
    if (p->isfits && p->hdu_type == IMAGE_HDU)
      {
        p->input = gal_array_read_one_ch_to_type (
            p->filename, p->cp.hdu, NULL, INPUT_USE_TYPE, p->cp.minmapsize,
            p->cp.quietmmap, "--hdu");
        p->input->wcs
            = gal_wcs_read (p->filename, p->cp.hdu, p->cp.wcslinearmatrix, 0,
                            0, &p->input->nwcs, "--hdu");
        p->input->ndim = gal_dimension_remove_extra (
            p->input->ndim, p->input->dsize, p->input->wcs);
      }
}

/* Read the kernel.  */
static void
ui_read_kernel (struct deconvolve_params *p)
{
  /* Read the image into file. */
  if (p->kernelname && gal_array_name_recognized (p->kernelname))
    {
      p->kernel = gal_array_read_one_ch_to_type (
          p->kernelname, p->khdu, NULL, INPUT_USE_TYPE, p->cp.minmapsize,
          p->cp.quietmmap, "--khdu");
      p->kernel->ndim = gal_dimension_remove_extra (
          p->kernel->ndim, p->kernel->dsize, p->kernel->wcs);
    }
}

/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/

static void
ui_preparations (struct deconvolve_params *p)
{
  struct gal_options_common_params *cp = &p->cp;
  char *outsuffix = "_deconvolved.fits";

  /* Read the input dataset. */
  ui_read_input (p);
  ui_read_kernel (p);

  /* Set the output name if the user hasn't set it. */
  if (cp->output == NULL)
    cp->output = gal_checkset_automatic_output (cp, p->filename, outsuffix);
  gal_checkset_writable_remove (cp->output, p->filename, 0, cp->dontdelete);
}

/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
ui_read_check_inputs_setup (int argc, char *argv[],
                            struct deconvolve_params *p)
{
  struct gal_options_common_params *cp = &p->cp;

  /* Include the parameters necessary for argp from this program ('args.h')
     and for the common options to all Gnuastro ('commonopts.h'). We want
     to directly put the pointers to the fields in 'p' and 'cp', so we are
     simply including the header here to not have to use long macros in
     those headers which make them hard to read and modify. This also helps
     in having a clean environment: everything in those headers is only
     available within the scope of this function. */
#include <gnuastro-internal/commonopts.h>

#include "args.h"

  /* Initialize the options and necessary information. */
  ui_initialize_options (p, program_options, gal_commonopts_options);

  /* Read the command-line options and arguments. */
  errno = 0;
  if (argp_parse (&thisargp, argc, argv, 0, 0, p))
    error (EXIT_FAILURE, errno, "parsing arguments");

  /* Read the configuration files and set the common values. */
  gal_options_read_config_set (&p->cp);

  /* Sanity check only on options. */
  ui_check_only_options (p);

  /* Print the option values if asked. Note that this needs to be done
     after the option checks so un-sane values are not printed in the
     output state. */
  gal_options_print_state (&p->cp);
  /* Prepare all the options as FITS keywords to write in output later. */
  gal_options_as_fits_keywords (&p->cp);

  /* Check that the options and arguments fit well with each other. Note
     that arguments don't go in a configuration file. So this test should
     be done after (possibly) printing the option values. */
  ui_check_options_and_arguments (p);

  /* Read/allocate all the necessary starting arrays. */
  ui_preparations (p);
}

/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report (struct deconvolve_params *p, struct timeval *t1)
{
  /* Free the allocated arrays. */
  free (p->cp.hdu);
  free (p->cp.output);

  /* Print the final message. */
  if (!p->cp.quiet)
    gal_timing_report (t1, PROGRAM_NAME " finished in: ", 0);
}
