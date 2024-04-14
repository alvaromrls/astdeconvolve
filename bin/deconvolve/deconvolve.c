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

#include <error.h>
#include <stdio.h>
#include <stdlib.h>

#include <main.h>

#include <gnuastro/deconvolve.h>

/*******************************************************************/
/*************            Top-level function           *************/
/*******************************************************************/
void
deconvolve (struct deconvolve_params *p)
{
  gal_data_t *data = NULL;
  size_t awmle_iterations = 0;
  uint8_t quiet = p->cp.quiet;
  switch (p->algorithm)
    {
    case DECONVOLUTION_ALGORITHM_DIRECT_INVERSION:
      if (quiet == 0)
        {
          printf ("Executing direct inversion algorithm \n");
        }
      data = gal_deconvolve_direct_inversion (
          p->input, p->kernel, p->cp.numthreads, p->cp.minmapsize);
      break;
    case DECONVOLUTION_ALGORITHM_TIKHONOV:
      if (quiet == 0)
        {
          printf ("Executing tikhonov algorithm with lambda = %f \n",
                  p->lambda);
        }
      data = gal_deconvolve_tikhonov (p->input, p->kernel, p->lambda,
                                      p->cp.numthreads, p->cp.minmapsize);
      break;
    case DECONVOLUTION_ALGORITHM_RL:
      if (quiet == 0)
        {
          printf (
              "Executing richardson-lucy algorithm with alpha = %f and %zu "
              "iterations\n",
              p->alpha, p->numberofitr);
        }
      data = gal_deconvolve_richardson_lucy (
          p->input, p->kernel, p->numberofitr, p->alpha, p->cp.minmapsize,
          p->cp.numthreads);
      break;
    case DECONVOLUTION_ALGORITHM_AWMLE:
      if (quiet == 0)
        {
          printf ("Executing AWMLE algorithm with alpha = %f, waves = "
                  "%zu, tolerance = %f, sigma = %f and %zu "
                  "iterations\n",
                  p->alpha, p->waves, p->tolerance, p->sigma, p->numberofitr);
        }
      data = gal_deconvolve_AWMLE (p->input, p->kernel, p->numberofitr,
                                   p->waves, p->tolerance, p->sigma, p->alpha,
                                   p->cp.minmapsize, p->cp.numthreads,
                                   &awmle_iterations);
      if ((awmle_iterations > 0) && (quiet == 0))
        {
          printf ("-----Early stop at %zu -----\n", awmle_iterations);
        }
      break;
    default:
      error (EXIT_FAILURE, 0, "Not implemented algorithm");
      break;
    }

  gal_data_t *original_type = gal_data_copy_to_new_type (data, p->input->type);
  gal_fits_img_write (original_type, p->cp.output, NULL, 0); // output option
  free (data);
  free (original_type);
}
