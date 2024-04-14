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
#include <config.h>

#include <errno.h>
#include <error.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <gnuastro/convolve.h>
#include <gnuastro/data.h>
#include <gnuastro/fft.h>
#include <gnuastro/pointer.h>
#include <gnuastro/warp.h>
#include <gnuastro/wavelet.h>

/*****************************************************************/
/**************    Father function defines    ********************/
/*****************************************************************/

/****************** Final B3 SPLINE *******************
                  1    4   6   4   1
                  4   16  24  16   4
                  6   24  36  24   6
                  4   16  24  16   4
                  1    4   6   4   1
******************************************************/

#define B3_SPLINE_SIZE 5
#define B3_SPLINE_1ST_ROW                                                     \
  {                                                                           \
    1.0, 4.0, 6.0, 4.0, 1.0                                                   \
  }
#define B3_SPLINE_2ND_ROW                                                     \
  {                                                                           \
    4.0, 16.0, 24.0, 16.0, 4.0                                                \
  }
#define B3_SPLINE_3RD_ROW                                                     \
  {                                                                           \
    6.0, 24.0, 36.0, 24.0, 6.0                                                \
  }

double *wavelet_substract (double *first, double *second, size_t size);

gal_data_t *wavelet_expand_father_function (gal_data_t *father, size_t factor,
                                            size_t minmapsize);

gal_data_t *wavelet_init_father_function (size_t minmapsize);

/**
 * @brief Substract two double arrays. Needed for rest calculation.
 *
 * @param first First array in the substraction.
 * @param second Second array in the substraction.
 * @param size  The number of elements in the array.
 * @return double* New array with the difference.
 */
double *
wavelet_substract (double *first, double *second, size_t size)
{
  double *out = gal_pointer_allocate (GAL_TYPE_FLOAT64, size, 1, __func__,
                                      "substract");
  for (size_t i = 0; i < size; i++)
    {
      out[i] = first[i] - second[i];
    }
  return out;
}

/**
 * @brief Generates the first father function.
 *
 * @param minmapsize
 * @return gal_data_t* father function (5x5)
 */
gal_data_t *
wavelet_init_father_function (size_t minmapsize)
{
  gal_data_t *data = NULL;
  double *fatherp;
  size_t dsize[] = { B3_SPLINE_SIZE, B3_SPLINE_SIZE };
  fatherp = gal_pointer_allocate (GAL_TYPE_FLOAT64,
                                  B3_SPLINE_SIZE * B3_SPLINE_SIZE, 1, __func__,
                                  "initfather");
  double father[B3_SPLINE_SIZE][B3_SPLINE_SIZE]
      = { B3_SPLINE_1ST_ROW, B3_SPLINE_2ND_ROW, B3_SPLINE_3RD_ROW,
          B3_SPLINE_2ND_ROW, B3_SPLINE_1ST_ROW };
  for (size_t i = 0; i < B3_SPLINE_SIZE; i++)
    {
      for (size_t j = 0; j < B3_SPLINE_SIZE; j++)
        {
          fatherp[i * B3_SPLINE_SIZE + j] = father[i][j];
        }
    }
  data = gal_data_alloc (fatherp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1,
                         minmapsize, 1, NULL, NULL, NULL);
  return data;
}
/**
 * @brief Given a father function, this function scales it by a factor. Uses a
 * 2D interpolation.
 *
 * @param father The father function to be scaled [NxN].
 * @param factor The scale factor (S).
 * @param minmapsize
 * @return gal_data_t* A father function in a bigger scale [N*SxN*S].
 */
gal_data_t *
wavelet_expand_father_function (gal_data_t *father, size_t factor,
                                size_t minmapsize)
{
  size_t dsize[] = { father->dsize[0] * factor, father->dsize[1] * factor };
  size_t outsize = dsize[0] * dsize[1];
  double *fatherp = (double *)father->array;
  double *nextfather = gal_pointer_allocate (GAL_TYPE_FLOAT64, outsize, 1,
                                             __func__, "nextfather");
  double xgrid[father->dsize[0]];
  double ygrid[father->dsize[1]];
  for (size_t i = 0; i < father->dsize[0]; i++)
    {
      xgrid[i] = (double)(i * (dsize[0] - 1)) / (double)(father->dsize[0] - 1);
      ygrid[i] = (double)(i * (dsize[1] - 1)) / (double)(father->dsize[1] - 1);
    }
  // Init 2D interp
  gsl_interp2d *interp = gsl_interp2d_alloc (
      gsl_interp2d_bilinear, father->dsize[0], father->dsize[1]);
  gsl_interp2d_init (interp, xgrid, ygrid, fatherp, father->dsize[0],
                     father->dsize[1]);

  // Init accel
  gsl_interp_accel *x_acc = gsl_interp_accel_alloc ();
  gsl_interp_accel *y_acc = gsl_interp_accel_alloc ();

  for (size_t y = 0; y < dsize[1]; y++)
    {
      for (size_t x = 0; x < dsize[0]; x++)
        {
          // Target coordinates
          double xi = (double)x;
          double yi = (double)y;
          // Interpolate
          double zi = gsl_interp2d_eval (interp, xgrid, ygrid, fatherp, xi, yi,
                                         x_acc, y_acc);
          nextfather[x + y * dsize[1]] = zi;
        }
    }
  // Free resources
  gsl_interp2d_free (interp);
  gsl_interp_accel_free (x_acc);
  gsl_interp_accel_free (y_acc);
  return gal_data_alloc (nextfather, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1,
                         minmapsize, 1, NULL, NULL, NULL);
}

/**
 * @brief Apply a wavelet decomposition to a given imagen. Uses no decimate
 * algorithm.
 *
 * @param image The image to be decomposed into wavelet planes.
 * @param numberplanes The number of planes (N).
 * @param numthreads Number fo threads for FFT transformations.
 * @param minmapsize
 * @return gal_data_t* N+1 images. Use ->next to navigate.
 */
gal_data_t *
gal_wavelet_no_decimate (const gal_data_t *image, uint8_t numberplanes,
                         size_t numthreads, size_t minmapsize)
{
  gal_data_t *current = NULL;
  gal_data_t *out = NULL;
  gal_data_t *father = NULL;
  gal_data_t *nfather;
  double *fatherpadding;
  father = wavelet_init_father_function (minmapsize);
  gsl_complex_packed_array fcomplex = NULL;

  gal_data_t *input_float
      = gal_data_copy_to_new_type (image, GAL_TYPE_FLOAT64);

  gsl_complex_packed_array input
      = gal_complex_real_to_complex (input_float->array, image->size);
  double *input_real = input_float->array;

  gsl_complex_packed_array rest = NULL;
  double *rest_real;
  double *wavelet;

  for (size_t iteration = 0; iteration < numberplanes; iteration++)
    {
      fatherpadding = gal_wavelet_add_padding (father->array, father->dsize,
                                               image->dsize);
      fcomplex = gal_complex_real_to_complex (fatherpadding, image->size);
      gal_complex_normalize (fcomplex, image->size);
      free (fatherpadding);
      gal_fft_shift_center (fcomplex, image->dsize);

      /* CONVOLUTION */
      rest = gal_convolve_frequency (fcomplex, input, image->dsize[0],
                                     image->dsize[1], numthreads, minmapsize);
      rest_real
          = gal_complex_to_real (rest, image->size, COMPLEX_TO_REAL_REAL);

      /*Wavelet*/
      wavelet = wavelet_substract (input_real, rest_real, image->size);

      if (current == NULL)
        { // 1st iteration where current is still empty
          current = gal_data_alloc (wavelet, GAL_TYPE_FLOAT64, 2, image->dsize,
                                    NULL, 1, minmapsize, 1, NULL, NULL, NULL);
          out = current;
        }
      else
        {
          current->next
              = gal_data_alloc (wavelet, GAL_TYPE_FLOAT64, 2, image->dsize,
                                NULL, 1, minmapsize, 1, NULL, NULL, NULL);
          current = current->next;
        }

      /* Update params */
      if (iteration + 1 == numberplanes)
        { // clean data
          free (father->array);
          free (father);
          free (rest);
          free (input);
          free (input_real);
          current->next
              = gal_data_alloc (rest_real, GAL_TYPE_FLOAT64, 2, image->dsize,
                                NULL, 1, minmapsize, 1, NULL, NULL, NULL);
        }
      else
        { // expand father
          nfather = wavelet_expand_father_function (father, 2, minmapsize);
          free (father->array);
          free (father);
          father = nfather;
          free (input);
          free (input_real);
          input = rest;
          input_real = rest_real;
        }
    }

  return out;
}

/**
 * @brief Generate 0's around an image until it fixes a desired size.
 *
 * @param input The image to be filled.
 * @param inputsize Original size [X,Y].
 * @param outputsize Desired size [X,y].
 * @return double* New array with input in the middle and outputsize
 * dimensions.
 */
double *
gal_wavelet_add_padding (double *input, size_t *inputsize, size_t *outputsize)
{
  /* Calculate usefull sizes */
  size_t input_size_x = inputsize[0];
  size_t input_size_y = inputsize[1];
  size_t output_size_x = outputsize[0];
  size_t output_size_y = outputsize[1];
  size_t half_input_size_x = input_size_x / 2;
  size_t half_input_size_y = input_size_y / 2;
  size_t half_output_size_x = output_size_x / 2;
  size_t half_output_size_y = output_size_y / 2;

  /*Specific cases*/
  size_t fix_mismacht = 0; // This variable fixes cases with odd-pair sizes
  if ((output_size_y % 2 == 0) && (input_size_y % 2 == 1))
    {
      fix_mismacht = 1;
    }
  else if ((output_size_y % 2 == 1) && (input_size_y % 2 == 0))
    {
      fix_mismacht = -1;
    }

  /* Allocate resources */
  double *out = gal_pointer_allocate (
      GAL_TYPE_FLOAT64, output_size_x * output_size_y, 1, __func__, "out");

  /* Fill the array */
  for (int i = 0; i < input_size_x; i++)
    {
      for (int j = 0; j < input_size_y; j++)
        {
          size_t outindex
              = (half_output_size_x - half_input_size_x + i - fix_mismacht)
                    * output_size_y
                + half_output_size_y - half_input_size_y + j;
          size_t inputindex = i * input_size_y + j;
          out[outindex] = input[inputindex];
        }
    }
  return out;
}

/**
 * @brief Free wavelet resources.
 *
 * @param wavelet A wavelet array to be cleaned.
 */
void
gal_wavelet_free (gal_data_t *wavelet)
{
  if (wavelet->next != NULL)
    {
      gal_wavelet_free (wavelet->next);
    }
  free (wavelet->array);
  free (wavelet);
}