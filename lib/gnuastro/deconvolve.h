/*********************************************************************
Functions for deconvolution operations.
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
#ifndef __GAL_DECONVOLVE_H__
#define __GAL_DECONVOLVE_H__

#include <config.h>

#include <gnuastro/data.h>

gal_data_t *gal_deconvolve_tikhonov (const gal_data_t *image,
                                     const gal_data_t *PSF, double lambda,
                                     size_t numthreads, size_t minmapsize);

gal_data_t *gal_deconvolve_direct_inversion (const gal_data_t *image,
                                             const gal_data_t *PSF,
                                             size_t numthreads,
                                             size_t minmapsize);

gal_data_t *gal_deconvolve_richardson_lucy (const gal_data_t *image,
                                            const gal_data_t *PSF,
                                            size_t iterations, double alpha,
                                            size_t minmapsize,
                                            size_t numthreads);

gal_data_t *gal_deconvolve_AWMLE (const gal_data_t *image,
                                  const gal_data_t *PSF, size_t iterations,
                                  size_t waves, double tolerance, double sigma,
                                  double alpha, size_t minmapsize,
                                  size_t numthreads, size_t *early);

#endif