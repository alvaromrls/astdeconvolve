/*********************************************************************
Functions to perform Wavelet decompositions.
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
#ifndef __GAL_WAVELET_H__
#define __GAL_WAVELET_H__
#include <config.h>

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <gnuastro/complex.h>
#include <gnuastro/data.h>
#include <gsl/gsl_fft_complex.h>

gal_data_t *gal_wavelet_no_decimate (const gal_data_t *image,
                                     uint8_t numberplanes, size_t numthreads,
                                     size_t minmapsize);

void gal_wavelet_free (gal_data_t *wavelet);

double *gal_wavelet_add_padding (double *input, size_t *inputsize,
                                 size_t *outputsize);

#endif
