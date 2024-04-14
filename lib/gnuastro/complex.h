/*********************************************************************
Functions to perform operations with complex arrays.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Alvaro Morales <alvaro96m@hotmail.com>
Contributing author(s):
Copyright (C) 2024-2024 Free Software Foundation, Inc.

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
#ifndef __GAL__COMPLEX__
#define __GAL__COMPLEX__

#include <config.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gnuastro/data.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>

// Enum with the different ways to convert a complex number into a real one.
typedef enum
{
  COMPLEX_TO_REAL_INVALID, // Error case
  COMPLEX_TO_REAL_SPEC,    // Take the module
  COMPLEX_TO_REAL_PHASE,   // Take the phase
  COMPLEX_TO_REAL_REAL,    // Take only the real part
} complex_to_real;

double *gal_complex_to_real (gsl_complex_packed_array complexarray,
                             size_t size, complex_to_real action);

gsl_complex_packed_array gal_complex_real_to_complex (double *realarray,
                                                      size_t size);

gsl_complex_packed_array gal_complex_multiply (gsl_complex_packed_array first,
                                               gsl_complex_packed_array second,
                                               size_t size);

gsl_complex_packed_array
gal_complex_substract (gsl_complex_packed_array first,
                       gsl_complex_packed_array second, size_t size);

void gal_complex_scale (gsl_complex_packed_array inout, double value,
                        size_t size);

gsl_complex_packed_array gal_complex_divide (gsl_complex_packed_array first,
                                             gsl_complex_packed_array second,
                                             size_t size, double minvalue);

gsl_complex_packed_array
gal_complex_conjugate (gsl_const_complex_packed_array input, size_t size);

gsl_complex_packed_array
gal_complex_add_scalar (gsl_const_complex_packed_array input, size_t size,
                        gsl_complex scalar);

void gal_complex_create_padding (const gal_data_t *image,
                                 const gal_data_t *kernel,
                                 gsl_complex_packed_array *outputimage,
                                 gsl_complex_packed_array *outputkernel,
                                 size_t *xdim, size_t *ydim);

void gal_complex_normalize (gsl_complex_packed_array inout, size_t size);

double gal_complex_cumulative_sum (gsl_complex_packed_array input,
                                   size_t size);

gsl_complex_packed_array gal_complex_power (gsl_complex_packed_array input,
                                            double exponent, size_t size);

#endif