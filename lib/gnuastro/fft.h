/*********************************************************************
Functions to perform Fast Fourier Transformations (fft).
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
#ifndef __GAL_FFT_H__
#define __GAL_FFT_H__
#include <config.h>

#include <gnuastro/complex.h>
#include <gsl/gsl_fft_complex.h>

gsl_complex_packed_array gal_fft_two_dimension_transformation (
    gsl_const_complex_packed_array input, size_t *dim, size_t numthreads,
    size_t minmapsize, gsl_fft_direction sign);

void gal_fft_shift_center (gsl_complex_packed_array kernel, size_t *dim);

#endif
