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
#ifndef MAIN_H
#define MAIN_H

/* Include necessary headers */
#include <gnuastro/data.h>

#include <gnuastro-internal/options.h>

/* Progarm names.  */
#define PROGRAM_NAME "deconvolve"    /* Program full name.       */
#define PROGRAM_EXEC "astdeconvolve" /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME " (" PACKAGE_NAME ") " PACKAGE_VERSION

/* Macros */
#define CONVFLOATINGPOINTERR 1e-10
#define INPUT_USE_TYPE GAL_TYPE_FLOAT32

typedef enum
{
  DECONVOLUTION_ALGORITHM_TIKHONOV,
  DECONVOLUTION_ALGORITHM_RL,
  DECONVOLUTION_ALGORITHM_AWMLE,
  DECONVOLUTION_ALGORITHM_DIRECT_INVERSION
} DECONVOLUTION_ALGORITHM;

/* Main program parameters structure */
struct deconvolve_params
{
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.   */
  char *filename;     /* Input filename.                         */
  char *kernelname;   /* File name of kernel.                    */
  char *khdu;         /* HDU of kernel.                          */
  char *algorithmstr; /*String with the algorithm name*/
  double lambda;      /*Lambda value for Tikhonov algorithm      */
  size_t numberofitr; /*Number of iterations for RL (default 75)*/
  double alpha;       /* Rate for RL (default 1) */
  size_t waves;       /*AWMLE*/
  double tolerance;   /*AWMLE*/
  double sigma;       /*AWMLE*/
  /* Internal */
  int isfits;   /* Input is a FITS file.                   */
  int hdu_type; /* Type of HDU (image or table).           */
  DECONVOLUTION_ALGORITHM algorithm;
  gal_data_t *input;  /* Input image array.                      */
  gal_data_t *kernel; /* Input Kernel array.                     */
  /* Output: */
  time_t rawtime; /* Starting time of the program.             */
};

#endif
