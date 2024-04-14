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
#include <config.h>

#include <errno.h>
#include <error.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <string.h>

#include <gnuastro/complex.h>
#include <gnuastro/deconvolve.h>
#include <gnuastro/fft.h>
#include <gnuastro/fits.h>
#include <gnuastro/list.h>
#include <gnuastro/pointer.h>
#include <gnuastro/wavelet.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>

// AWMLE defines
#define DECONVOLVE_AWMLE_MASK_FACTOR 1.5
#define DECONVOLVE_AWMLE_NOISE_LEVEL_ADJUST 1.0
#define DECONVOLVE_AWMLE_SIGMA_MIN 1.0e-6

void deconvolve_AWMLE_update_object (double *correctionterm, double energy,
                                     gal_data_t *modified_image_res,
                                     gal_data_t *projection,
                                     gsl_complex_packed_array psf_fft_conj,
                                     double *object, double alpha,
                                     size_t numthreads, size_t minmapsize);

gal_data_t *deconvolve_AWMLE_estimate_p_prime (gal_data_t *image,
                                               gal_data_t *projection,
                                               double sigma, double *likehood,
                                               size_t minmapsize);

gal_data_t *
deconvolve_calcule_AWMLE_noise_factor (size_t planes, size_t *dsize,
                                       gal_data_t *residue, double sigma,
                                       size_t minmapsize, size_t numthreads);

void deconvolve_AWMLE_increment_correction_term (
    double *correctionterm_A, gal_data_t *projection_wave, gal_data_t *mask,
    gal_data_t *modifiedimage_wave, gal_data_t *projection);

gal_data_t *deconvolve_calcule_AWMLE_mask (gal_data_t *wavelet,
                                           gal_data_t *projection,
                                           gal_data_t *noise, size_t ampl,
                                           size_t numthreads,
                                           size_t minmapsize);

void deconvolve_richardson_lucy_calculate_next_solution (
    gsl_complex_packed_array *solution, gsl_complex_packed_array h,
    gsl_complex_packed_array h_f, gsl_complex_packed_array image, double alpha,
    size_t *dsize, size_t minmapsize, size_t numthreads);

gsl_complex_packed_array
deconvolve_richardson_lucy_init_solution (size_t size);

/**
 * @brief Implement a deconvolution using the Wiener / Tikhonov method
 * described at https://ui.adsabs.harvard.edu/abs/2002PASP..114.1051S/abstract
 * (eq 7)
 *
 * @param image The image to be deconvolved.
 * @param PSF The distorsion kernel.
 * @param lambda The Tikhonov main parameter.
 * @param numthreads Number of threads in FFT.
 * @param minmapsize
 * @return gal_data_t* The image deconvolved.
 */
gal_data_t *
gal_deconvolve_tikhonov (const gal_data_t *image, const gal_data_t *PSF,
                         double lambda, size_t numthreads, size_t minmapsize)
{
  gsl_complex_packed_array imagepadding;  // original image after padding
  gsl_complex_packed_array psfpadding;    // kernel after padding
  gsl_complex_packed_array psffreq;       // PSF(u,v)
  gsl_complex_packed_array imagefreq;     // I(u,v)
  gsl_complex_packed_array psffconj;      // PSF*(u,v)
  gsl_complex_packed_array psffreqsquare; // |PSF(u,v)|^2
  gsl_complex_packed_array numerator;     // PSF*(u,v)I(u,v)
  gsl_complex_packed_array denominator;   // |PSF(u,v)|^2 + λ
  gsl_complex_packed_array
      deconvolutionfreq; // PSF*(u,v)I(u,v)/ (|PSF(u,v)|^2 + λ)
  gsl_complex_packed_array deconvolution; // Deconvolution in time domine

  gal_data_t *data = NULL;
  size_t dsize[2];
  size_t size;
  double *tmp;

  /* Check image type. */
  if (image->type != GAL_TYPE_FLOAT32)
    error (EXIT_FAILURE, 0, "%s: input data must be float 32", __func__);

  /* Process the image and kernel to have same size and be complex numbers. */
  gal_complex_create_padding (image, PSF, &imagepadding, &psfpadding, dsize,
                              &dsize[1]);
  size = dsize[0] * dsize[1]; // Total number of elements.

  /* Normalize and rearange the kernel. */
  gal_complex_normalize (psfpadding, size);
  gal_fft_shift_center (psfpadding, dsize);

  /* Convert to frequency domain. */
  psffreq = gal_fft_two_dimension_transformation (
      psfpadding, dsize, numthreads, minmapsize, gsl_fft_forward);
  imagefreq = gal_fft_two_dimension_transformation (
      imagepadding, dsize, numthreads, minmapsize, gsl_fft_forward);

  /* Calculate numerator PSF*(u,v)I(u,v) */
  psffconj = gal_complex_conjugate (psffreq, size);
  numerator = gal_complex_multiply (psffconj, imagefreq, size);

  /* Caculate denominator |PSF(u,v)|^2 + λ */
  psffreqsquare = gal_complex_multiply (psffreq, psffconj, size);
  denominator = gal_complex_add_scalar (psffreqsquare, size, lambda + I * 0);

  /* Calculate the deconvolve image (in frequency domain).*/
  deconvolutionfreq
      = gal_complex_divide (numerator, denominator, size, lambda);

  /* Go back to time domain. */
  deconvolution = gal_fft_two_dimension_transformation (
      deconvolutionfreq, dsize, numthreads, minmapsize, gsl_fft_backward);

  /* Convert to Real number and convert it to GAL TYPE.*/
  tmp = gal_complex_to_real (deconvolution, dsize[0] * dsize[1],
                             COMPLEX_TO_REAL_REAL);
  data = gal_data_alloc (tmp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1, minmapsize,
                         1, NULL, NULL, NULL); // data has to be 32

  /* Free resources. */
  free (imagepadding);
  free (psfpadding);
  free (psffreq);
  free (psffconj);
  free (imagefreq);
  free (psffreqsquare);
  free (numerator);
  free (denominator);
  free (deconvolutionfreq);
  free (deconvolution);
  return data;
}

/**
 * @brief Implements a deconvolution using the direct inversion method:
 * O(u,v) = I(u,v)/PSF(u,v)
 *
 * @param image The image to be deconvolved.
 * @param PSF The distorsion kernel.
 * @param numthreads Number of threads in FFT.
 * @param minmapsize
 * @return gal_data_t* The image deconvolved.
 */
gal_data_t *
gal_deconvolve_direct_inversion (const gal_data_t *image,
                                 const gal_data_t *PSF, size_t numthreads,
                                 size_t minmapsize)
{
  gsl_complex_packed_array imagepadding;      // original image after padding
  gsl_complex_packed_array psfpadding;        // kernel after padding
  gsl_complex_packed_array psffreq;           // PSF(u,v)
  gsl_complex_packed_array imagefreq;         // I(u,v)
  gsl_complex_packed_array deconvolutionfreq; // I(u,v) / PSF(u,v)
  gsl_complex_packed_array deconvolution;     // Deconvolution in time domine

  gal_data_t *data = NULL;
  size_t dsize[2];
  size_t size;
  double *tmp;

  /* Check image type. */
  if (image->type != GAL_TYPE_FLOAT32)
    error (EXIT_FAILURE, 0, "%s: input data must be float 32", __func__);

  /* Process the image and kernel to have same size and be complex numbers. */
  gal_complex_create_padding (image, PSF, &imagepadding, &psfpadding, dsize,
                              &dsize[1]);
  size = dsize[0] * dsize[1]; // Total number of elements.

  /* Normalize and rearange the kernel. */
  gal_complex_normalize (psfpadding, size);
  gal_fft_shift_center (psfpadding, dsize);

  /* Convert to frequency domain. */
  psffreq = gal_fft_two_dimension_transformation (
      psfpadding, dsize, numthreads, minmapsize, gsl_fft_forward);
  imagefreq = gal_fft_two_dimension_transformation (
      imagepadding, dsize, numthreads, minmapsize, gsl_fft_forward);

  /* Calculate the deconvolve image (in frequency domain).*/
  deconvolutionfreq = gal_complex_divide (imagefreq, psffreq, size, 1e-6);

  /* Go back to time domain. */
  deconvolution = gal_fft_two_dimension_transformation (
      deconvolutionfreq, dsize, numthreads, minmapsize, gsl_fft_backward);

  /* Convert to Real number and convert it to GAL TYPE.*/
  tmp = gal_complex_to_real (deconvolution, dsize[0] * dsize[1],
                             COMPLEX_TO_REAL_REAL);
  data = gal_data_alloc (tmp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1, minmapsize,
                         1, NULL, NULL, NULL); // data has to be 32

  /* Free resources. */
  free (imagepadding);
  free (psfpadding);
  free (psffreq);
  free (deconvolutionfreq);
  free (deconvolution);

  return data;
}

/**
 * @brief Inits an array with 1's as first solution for iterative methods.
 *
 * @param size Number of elements in the generated array.
 * @return gsl_complex_packed_array An array with 1's.
 */
gsl_complex_packed_array
deconvolve_richardson_lucy_init_solution (size_t size)
{
  gsl_complex_packed_array out;

  /* Allocate the space for the output array. */
  out = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size, 1, __func__,
                              "rlobject");

  for (size_t index = 0; index < size; index++)
    {
      out[index * 2] = 1;
    }
  return out;
}

/**
 * @brief Performs one iteration in Richardson-Lucy algorithm.
 *
 * @param solution Last iteration image and next image (in-out)
 * @param h The PSF in Fourier domain.
 * @param h_f The PSF conjugated in Fourier domain.
 * @param image The original image.
 * @param alpha The growth factor.
 * @param dsize 2D dimension size (X,Y).
 * @param minmapsize
 * @param numthreads Number of threads.
 */
void
deconvolve_richardson_lucy_calculate_next_solution (
    gsl_complex_packed_array *solution, gsl_complex_packed_array h,
    gsl_complex_packed_array h_f, gsl_complex_packed_array image, double alpha,
    size_t *dsize, size_t minmapsize, size_t numthreads)
{
  gsl_complex_packed_array solution_f;    // FFT(x)
  gsl_complex_packed_array yest_f;        // h×FFT(x(k))
  gsl_complex_packed_array yest;          // y_est := FFT−1(h×FFT(x))
  gsl_complex_packed_array division;      // y/y_est
  gsl_complex_packed_array divisionfreq;  // FFT(y/y_est)
  gsl_complex_packed_array bracket_f;     // h*xFFT(y/y_est)
  gsl_complex_packed_array bracket;       // FFT^-1(bracket_f)
  gsl_complex_packed_array bracket_alpha; // FFT^-1(bracket_f)^alpha
  gsl_complex_packed_array next_solution; // FFT^-1(bracket_f)^alpha·x
  size_t size = dsize[0] * dsize[1];

  /* Convert solution to frequency domain. */
  solution_f = gal_fft_two_dimension_transformation (
      *solution, dsize, numthreads, minmapsize, gsl_fft_forward);

  yest_f = gal_complex_multiply (solution_f, h, size);

  /*Calculate y_est*/
  yest = gal_fft_two_dimension_transformation (yest_f, dsize, numthreads,
                                               minmapsize, gsl_fft_backward);

  /*Calculate y/y_est and its FFT*/
  division = gal_complex_divide (image, yest, size,
                                 1e-6); // check min value issues

  divisionfreq = gal_fft_two_dimension_transformation (
      division, dsize, numthreads, minmapsize, gsl_fft_forward);

  /* Calculate bracket value in Freq and Space*/
  bracket_f = gal_complex_multiply (divisionfreq, h_f, size);

  bracket = gal_fft_two_dimension_transformation (
      bracket_f, dsize, numthreads, minmapsize, gsl_fft_backward);
  // alpha value
  bracket_alpha = gal_complex_power (bracket, alpha, size);
  /* Calculate next object */
  next_solution = gal_complex_multiply (bracket_alpha, *solution, size);

  /* Free all internal variables */
  free (*solution);
  free (solution_f);
  free (yest_f);
  free (yest);
  free (division);
  free (divisionfreq);
  free (bracket_f);
  free (bracket);
  free (bracket_alpha);

  *solution = next_solution;
}

/**
 * @brief Perform a deconvolution with Richardson-Lucy algorithm.
 *
 * @param image The image be deconvolved.
 * @param PSF Kernel used to take the image.
 * @param iterations The number of iterations.
 * @param alpha The growth factor.
 * @param minmapsize
 * @param numthreads Number of threads.
 * @return gal_data_t* The deconvolved image.
 */
gal_data_t *
gal_deconvolve_richardson_lucy (const gal_data_t *image, const gal_data_t *PSF,
                                size_t iterations, double alpha,
                                size_t minmapsize, size_t numthreads)
{
  gsl_complex_packed_array imagepadding; // original image after padding I(x,y)
  gsl_complex_packed_array psfpadding;   // kernel after padding PSF(x,y)
  gsl_complex_packed_array psffreq;      // PSF(u,v)
  gsl_complex_packed_array psffconj;     // PSF*(u,v)
  gsl_complex_packed_array object;       // O(x,y)

  gal_data_t *data = NULL;
  size_t dsize[2];
  size_t size;
  double *tmp;

  /* Check image type. */
  if (image->type != GAL_TYPE_FLOAT32)
    error (EXIT_FAILURE, 0, "%s: input data must be float 32", __func__);

  /* Process the image and kernel to have same size and be complex numbers. */
  gal_complex_create_padding (image, PSF, &imagepadding, &psfpadding, dsize,
                              &dsize[1]);
  size = dsize[0] * dsize[1]; // Total number of elements.

  /* Normalize and rearange the kernel. */
  gal_complex_normalize (psfpadding, size);
  gal_fft_shift_center (psfpadding, dsize);

  /* Init the solution */
  object = deconvolve_richardson_lucy_init_solution (size);

  /* Convert to frequency domain. */
  psffreq = gal_fft_two_dimension_transformation (
      psfpadding, dsize, numthreads, minmapsize, gsl_fft_forward);

  /* Calculate  PSF*(u,v) */
  psffconj = gal_complex_conjugate (psffreq, size);

  for (size_t iteration = 0; iteration < iterations; iteration++)
    {
      deconvolve_richardson_lucy_calculate_next_solution (
          &object, psffreq, psffconj, imagepadding, alpha, dsize, minmapsize,
          numthreads);
    }

  /* Convert to Real number and convert it to GAL TYPE.*/
  tmp = gal_complex_to_real (object, dsize[0] * dsize[1],
                             COMPLEX_TO_REAL_REAL);
  data = gal_data_alloc (tmp, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1, minmapsize,
                         1, NULL, NULL, NULL); // data has to be 32

  free (imagepadding);
  free (psfpadding);
  free (psffreq);
  free (psffconj);
  free (object);

  return data;
}

/**
 * @brief Subfunction to estimate p prime in AWMLE. Only for internal use.
 *
 * @param image The image to estimate p prime. Must be float 64.
 * @param projection The image projection in the psf (convolution).
 * @param sigma The gaussian noise estimation. If < DECONVOLVE_AWMLE_SIGMA_MIN,
 * only poisson noise will be considered.
 * @param likehood Output parameter, value used to early stop the algorithm.
 * @param minmapsize
 * @return gal_data_t* The p prime estimation
 */
gal_data_t *
deconvolve_AWMLE_estimate_p_prime (gal_data_t *image, gal_data_t *projection,
                                   double sigma, double *likehood,
                                   size_t minmapsize)
{

  /*Check image is F64*/
  if (image->type != GAL_TYPE_FLOAT64)
    error (EXIT_FAILURE, 0, "%s: input data must be float 64", __func__);

  size_t size = image->size;
  double *output
      = gal_pointer_allocate (GAL_TYPE_FLOAT64, size, 1, __func__, "pprime");
  double *imagearray = image->array;
  double *projectionarray = projection->array;
  *likehood = 0.0;

  /* Low level values */
  double v = sigma;
  double vv = sigma * sigma;
  double vv2 = 2. * vv;
  double exp2 = exp (-1. / vv);
  double expantor = exp (1. / vv2);
  double lc = log (sqrt (M_PI * 2 * v));
  double max = 2.5 * v;

  if (sigma < DECONVOLVE_AWMLE_SIGMA_MIN)
    // Sigma so small we can consider poisson noise
    {
      memcpy (output, imagearray, sizeof (double) * size);
      for (size_t i = 0; i < size; i++)
        {
          if (imagearray[i] != 0.0)
            {
              *likehood += (imagearray[i]
                                * (1 + log (projectionarray[i])
                                   - log (fabs (imagearray[i])))
                            - projectionarray[i]);
            }
          else
            {
              *likehood -= projectionarray[i];
            }
        }
    }
  else
    {
      for (size_t i = 0; i < size; i++)
        {
          double p = imagearray[i];
          double hh = projectionarray[i];
          int lower = (int)(p - 2.5 * v);
          int upper = (int)(p + 2.5 * v);
          double hp = (hh - p) * (hh - p) + 1.;
          double lhp = log (hp);
          double new = 0;
          double serieantden, serieantnum, termeantnum, termeantden, expant;
          double termeantnumsup, termeantdensup, termeantnuminf,
              termeantdeninf;

          if (p <= max)
            {
              lower = 0;
              if (upper < 1)
                {
                  upper = 1;
                }
              serieantden = exp (-(1. - p) * (1. - p) / vv2) * hh / hp;
              serieantnum = serieantden;
              termeantnum = serieantden;
              termeantden = serieantden;

              expant = exp (-(1. - 2. * p) / vv2);
              for (int k = 2; k <= upper; ++k)
                {
                  expant = expant * exp2;
                  termeantnum = termeantnum * hh / (k - 1) * expant;
                  termeantden = termeantden * hh / k * expant;
                  serieantnum += termeantnum;
                  serieantden += termeantden;
                }
              serieantden += exp (-p * p / vv2) / hp;
              new = -lc - hh + log (serieantden) + lhp;
            }
          else
            {
              double lp = log (p);
              double lh = log (hh);
              double l2 = log (sqrt (2 * M_PI * p));

              serieantden = 1. / hp;
              serieantnum = p / hp;
              termeantnumsup = serieantnum;
              termeantdensup = serieantden;
              termeantnuminf = serieantnum;
              termeantdeninf = serieantden;
              expant = expantor;

              int k = (int)(p) + 1;
              int l = (int)(p)-1;
              for (; k <= upper && l >= lower; k++, l--)
                {
                  expant = expant * exp2;
                  termeantnumsup = termeantnumsup * hh / (k - 1) * expant;
                  termeantdensup = termeantdensup * hh / k * expant;
                  termeantnuminf = termeantnuminf * l / hh * expant;
                  termeantdeninf = termeantdeninf * (l + 1) / hh * expant;
                  serieantnum += termeantnumsup + termeantnuminf;
                  serieantden += termeantdensup + termeantdeninf;
                }
              new = -lc - hh + log (serieantden) + p *lh - l2 - p *lp + p
                    + lhp;
            }
          double new_output = serieantnum / serieantden;
          if (!isnan (new_output) && !isinf (new_output) && !isnan (new)
              && !isinf (new))
            {
              output[i] = new_output;
              (*likehood) += new;
            }
        }
    }
  return gal_data_alloc (output, GAL_TYPE_FLOAT64, 2, image->dsize, NULL, 1,
                         minmapsize, 1, NULL, NULL, NULL);
}

/**
 * @brief Subfunction to calculate the stadistical mask (mj).
 *
 * @param wavelet The wavelet image.
 * @param projection The projection in that plane.
 * @param noise The noise factor in that wavelet plane.
 * @param ampl The size of the window.
 * @param numthreads The number of threads.
 * @param minmapsize
 * @return gal_data_t* The stadistical mask. Must be in range [0,1].
 */
gal_data_t *
deconvolve_calcule_AWMLE_mask (gal_data_t *wavelet, gal_data_t *projection,
                               gal_data_t *noise, size_t ampl,
                               size_t numthreads, size_t minmapsize)
{
  size_t *dsize = wavelet->dsize;
  size_t size = wavelet->size;
  size_t ampl2 = ampl * ampl;
  size_t dsize_ampl[] = { ampl, ampl };

  double *window0
      = gal_pointer_allocate (GAL_TYPE_FLOAT64, ampl2, 1, __func__, "window");
  for (size_t i = 0; i < ampl2; i++)
    {
      window0[i] = 1.0 / ((double)ampl2);
    }

  // Convert window to frequency domain
  double *window
      = gal_wavelet_add_padding (window0, dsize_ampl, wavelet->dsize);
  free (window0);
  gsl_complex_packed_array windowc
      = gal_complex_real_to_complex (window, size);
  free (window);
  gal_fft_shift_center (windowc, wavelet->dsize);
  gsl_complex_packed_array fft_window = gal_fft_two_dimension_transformation (
      windowc, dsize, numthreads, minmapsize, gsl_fft_forward);
  free (windowc);

  // Calculate and transform square error to wavelet and projection
  gsl_complex_packed_array waveletc
      = gal_complex_real_to_complex (wavelet->array, size);
  gsl_complex_packed_array projectionc
      = gal_complex_real_to_complex (projection->array, size);

  gsl_complex_packed_array error
      = gal_complex_substract (waveletc, projectionc, size);
  free (waveletc);
  free (projectionc);

  gsl_complex_packed_array error_square = gal_complex_power (error, 2, size);
  free (error);

  gsl_complex_packed_array fft_errorsq = gal_fft_two_dimension_transformation (
      error_square, dsize, numthreads, minmapsize, gsl_fft_forward);
  free (error_square);

  // Calculate mask
  gsl_complex_packed_array convolution_freq
      = gal_complex_multiply (fft_errorsq, fft_window, size);
  free (fft_errorsq);
  free (fft_window);

  gsl_complex_packed_array convolutionc
      = gal_fft_two_dimension_transformation (
          convolution_freq, dsize, numthreads, minmapsize, gsl_fft_backward);
  free (convolution_freq);

  double *convolution
      = gal_complex_to_real (convolutionc, size, COMPLEX_TO_REAL_REAL);
  free (convolutionc);

  double *aux
      = gal_pointer_allocate (GAL_TYPE_FLOAT64, size, 1, __func__, "aux");
  for (size_t i = 0; i < size; i++)
    {
      aux[i] = sqrt (convolution[i]);
    }
  free (convolution);

  double *mask
      = gal_pointer_allocate (GAL_TYPE_FLOAT64, size, 1, __func__, "mask");
  double *noisearray = noise->array;
  for (size_t i = 0; i < size; i++)
    {
      double a = DECONVOLVE_AWMLE_MASK_FACTOR * (aux[i] - noisearray[i]);
      mask[i]
          = 1.0
            - exp ((-1.0) * (a * a) / (2.0 * noisearray[i] * noisearray[i]));

      // Safety check!
      if (mask[i] > 1.0)
        {
          mask[i] = 1.0;
        }
      if (isnan (mask[i]) || (mask[i] <= DBL_EPSILON))
        {
          mask[i] = DBL_EPSILON;
        }
    }
  free (aux);
  return gal_data_alloc (mask, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1, minmapsize,
                         1, NULL, NULL, NULL);
}

/**
 * @brief Subfunction to calculate the noise factor. Only for interal use.
 *
 * @param planes Number of planes in the algorithm.
 * @param dsize Dimensions [x,y].
 * @param residue The residual from the image wavelet descomposition.
 * @param sigma Gaussian noise estimation.
 * @param minmapsize
 * @param numthreads Number of threads.
 */
gal_data_t *
deconvolve_calcule_AWMLE_noise_factor (size_t planes, size_t *dsize,
                                       gal_data_t *residue, double sigma,
                                       size_t minmapsize, size_t numthreads)
{
  // Allocate resources
  size_t size = dsize[0] * dsize[1];
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  // Generate a Gauss Noise Matrix
  double *gaussian
      = gal_pointer_allocate (GAL_TYPE_FLOAT64, size, 1, __func__, "gaussian");
  for (size_t i = 0; i < size; i++)
    {
      gaussian[i] = gsl_ran_gaussian (r, 1.0);
    }

  gsl_rng_free (r); // free random resources

  // Wavelet transformation to the Gaussian distr
  gal_data_t *gaussian_noise
      = gal_data_alloc (gaussian, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1,
                        minmapsize, 1, NULL, NULL, NULL);

  gal_data_t *gaussian_waves = gal_wavelet_no_decimate (
      gaussian_noise, planes, numthreads, minmapsize);
  gal_data_free (gaussian_noise);

  // Noise 0 estimation
  double *imagew = residue->array;
  double *noise0
      = gal_pointer_allocate (GAL_TYPE_FLOAT64, size, 1, __func__, "noise0");
  for (size_t i = 0; i < size; i++)
    {
      noise0[i] = sqrt (abs (imagew[i] + sigma * sigma));
    }

  // Noise Factor
  gal_data_t *p = gaussian_waves;
  gal_data_t *noise = NULL;
  for (size_t plane = 0; plane < planes; plane++)
    {
      double noiseFactor = gsl_stats_sd_m (p->array, 1, size, 0.0)
                           / DECONVOLVE_AWMLE_NOISE_LEVEL_ADJUST;
      p = p->next;
      double *noise_array = gal_pointer_allocate (GAL_TYPE_FLOAT64, size, 1,
                                                  __func__, "noise");
      for (size_t i = 0; i < size; i++)
        {
          noise_array[i] = noise0[i] * noiseFactor;
        }
      gal_data_t *new_noise
          = gal_data_alloc (noise_array, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1,
                            minmapsize, 1, NULL, NULL, NULL);
      gal_list_data_add (&noise, new_noise);
    }
  gal_list_data_reverse (&noise);
  gal_list_data_free (gaussian_waves);
  free (noise0);
  return noise;
}

/**
 * @brief Perform a deconvolution with AWMLE algorithm.
 *
 * @param image The image to be deconvolved.
 * @param PSF  The distorsion kernel.
 * @param iterations Number of iterations.
 * @param waves Number of planes in wavelet decomposition.
 * @param tolerance Threshold value to perform an early stop.
 * @param sigma Estimation of gaussian noise (standard desviation).
 * @param alpha The growth factor.
 * @param minmapsize
 * @param numthreads Number of threads.
 * @param early Out parameter: If > 0 an early stop took place and shows the
 * number of iterations. Early < Iterations in all cases.
 * @return gal_data_t* The deconvolved image.
 */
gal_data_t *
gal_deconvolve_AWMLE (const gal_data_t *image, const gal_data_t *PSF,
                      size_t iterations, size_t waves, double tolerance,
                      double sigma, double alpha, size_t minmapsize,
                      size_t numthreads, size_t *early)
{
  size_t size = image->size;
  size_t *dsize = image->dsize;
  *early = 0;
  /* Check and adapt data types */
  gal_data_t *image64;
  if (image->type == GAL_TYPE_FLOAT64)
    {
      image64 = gal_data_copy (image);
    }
  else
    {
      image64 = gal_data_copy_to_new_type (image, GAL_TYPE_FLOAT64);
    }

  if (PSF == NULL)
    {
      error (EXIT_FAILURE, 0, "%s: Error loading PSF fits file", __func__);
    }

  gal_data_t *psf64;
  if (PSF->type == GAL_TYPE_FLOAT64)
    {
      psf64 = gal_data_copy (PSF);
    }
  else
    {
      psf64 = gal_data_copy_to_new_type (PSF, GAL_TYPE_FLOAT64);
    }

  /* Prepare the PSF in freq domain*/
  double *psf_padding
      = gal_wavelet_add_padding (psf64->array, PSF->dsize, dsize);
  gsl_complex_packed_array psf_complex
      = gal_complex_real_to_complex (psf_padding, size);
  free (psf_padding);

  gal_fft_shift_center (psf_complex, dsize);
  gal_complex_normalize (psf_complex, size);

  gsl_complex_packed_array psf_fft = gal_fft_two_dimension_transformation (
      psf_complex, dsize, numthreads, minmapsize, gsl_fft_forward);
  free (psf_complex);
  gsl_complex_packed_array psf_fft_conj
      = gal_complex_conjugate (psf_fft, size);

  /* Decompose the image into wavelets*/
  gal_data_t *imagewaves
      = gal_wavelet_no_decimate (image, waves, numthreads, minmapsize);

  gal_data_t *imageres = gal_list_data_last (imagewaves);

  /* Get noise factor */
  gal_data_t *noise = deconvolve_calcule_AWMLE_noise_factor (
      waves, dsize, imageres, sigma, minmapsize, numthreads);

  double fx_new = DBL_MAX;
  // energy ?

  /* Initialize data*/
  gsl_complex_packed_array object_c
      = deconvolve_richardson_lucy_init_solution (size);
  double *object = gal_complex_to_real (object_c, size, COMPLEX_TO_REAL_REAL);

  double energy = 0;
  double *image_array = image64->array;
  for (size_t i = 0; i < size; i++)
    {
      energy += image_array[i];
    }

  double likehood;
  double last_likehood = DBL_EPSILON;
  for (size_t iteration = 0; iteration < iterations; iteration++)
    {
      double *correctionterm = gal_pointer_allocate (GAL_TYPE_FLOAT64, size, 1,
                                                     __func__, "correction");

      gsl_complex_packed_array object_fft
          = gal_fft_two_dimension_transformation (object_c, dsize, numthreads,
                                                  minmapsize, gsl_fft_forward);

      gsl_complex_packed_array projection_fft
          = gal_complex_multiply (object_fft, psf_fft, size);
      free (object_fft);

      gsl_complex_packed_array projection_c
          = gal_fft_two_dimension_transformation (
              projection_fft, dsize, numthreads, minmapsize, gsl_fft_backward);
      free (projection_fft);

      double *projection
          = gal_complex_to_real (projection_c, size, COMPLEX_TO_REAL_REAL);

      /* Secury check */
      for (size_t i = 0; i < size; i++)
        {
          if (projection[i] < DBL_EPSILON)
            {
              projection[i] = DBL_EPSILON;
            }
        }

      gal_data_t *projection_gal
          = gal_data_alloc (projection, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1,
                            minmapsize, 1, NULL, NULL, NULL);
      free (projection_c);

      gal_data_t *modifiedimage = deconvolve_AWMLE_estimate_p_prime (
          image64, projection_gal, sigma, &likehood, minmapsize);

      /* Decompose in wavelets */
      gal_data_t *modifiedimage_wavelet = gal_wavelet_no_decimate (
          modifiedimage, waves, numthreads, minmapsize);
      gal_data_t *projection_wavelet = gal_wavelet_no_decimate (
          projection_gal, waves, numthreads, minmapsize);

      /* For each wavelet */
      gal_data_t *imagewaves_p = imagewaves;
      gal_data_t *noise_p = noise;
      gal_data_t *projection_p = projection_wavelet;
      gal_data_t *modifiedimage_p = modifiedimage_wavelet;

      for (size_t wave = 0; wave < waves; wave++)
        {
          uint ampl = 2 * ((uint)pow (2, wave) + 1) + 1;
          gal_data_t *mask = deconvolve_calcule_AWMLE_mask (
              imagewaves_p, projection_p, noise_p, ampl, numthreads,
              minmapsize);

          deconvolve_AWMLE_increment_correction_term (
              correctionterm, projection_p, mask, modifiedimage_p,
              projection_gal);

          /* Increment wavelet pointers*/
          imagewaves_p = imagewaves_p->next;
          noise_p = noise_p->next;
          projection_p = projection_p->next;
          modifiedimage_p = modifiedimage_p->next;
          gal_data_free (mask);
        }

      deconvolve_AWMLE_update_object (correctionterm, energy, modifiedimage_p,
                                      projection_gal, psf_fft_conj, object,
                                      alpha, numthreads, minmapsize);

      free (object_c);
      object_c = gal_complex_real_to_complex (object, size);

      free (correctionterm);
      gal_list_data_free (projection_wavelet);
      gal_list_data_free (modifiedimage_wavelet);

      if (iteration == 0)
        {
          last_likehood = likehood;
        }
      else
        {
          double diff_likehood = likehood - last_likehood;
          double median_likehood = (last_likehood + likehood) / 2.0;
          if (fabs (diff_likehood / median_likehood) <= tolerance)
            {
              *early = iteration;
              iteration = iterations;
            }
          else
            {
              last_likehood = likehood;
            }
        }
    }

  free (object_c);
  free (psf_fft);
  free (psf_fft_conj);
  gal_data_free (image64);
  gal_data_free (psf64);
  gal_list_data_free (imagewaves);
  gal_list_data_free (noise);

  return gal_data_alloc (object, GAL_TYPE_FLOAT64, 2, dsize, NULL, 1,
                         minmapsize, 1, NULL, NULL, NULL);
}

/**
 * @brief Subfunction to update the object at the end of each iteration.
 * 1st it updates the correction term based on the residue from pprime.
 * Then it updates the object.
 * Last, it adjust the energy in the object.
 *
 * @param correctionterm The correction term.
 * @param energy The original energy value (sum of all elements).
 * @param modified_image_res Pprime residue.
 * @param projection The image projection.
 * @param psf_fft_conj psf conjugate.
 * @param object The current image reconstruction (in/out parameter).
 * @param alpha The growing factor.
 * @param numthreads The number of threads.
 * @param minmapsize
 */
void
deconvolve_AWMLE_update_object (double *correctionterm, double energy,
                                gal_data_t *modified_image_res,
                                gal_data_t *projection,
                                gsl_complex_packed_array psf_fft_conj,
                                double *object, double alpha,
                                size_t numthreads, size_t minmapsize)
{
  size_t size = projection->size;
  size_t *dsize = projection->dsize;
  double *projection_array = projection->array;
  double *res_array = modified_image_res->array;
  double newenergy = 0;
  /*Update correction term with the residue*/
  for (size_t i = 0; i < size; i++)
    {
      if (projection_array[i] > DBL_EPSILON)
        {
          correctionterm[i] += res_array[i] / projection_array[i];
          // printf ("Correction term at %zu is %f\n", i, correctionterm[i]);
          if (isnan (correctionterm[i]) || isinf (correctionterm[i]))
            {
              correctionterm[i] = DBL_EPSILON;
            }
        }
    }

  /* Calculate new object (in fft space)*/
  gsl_complex_packed_array correction_complex
      = gal_complex_real_to_complex (correctionterm, size);
  gsl_complex_packed_array correction_fft
      = gal_fft_two_dimension_transformation (
          correction_complex, dsize, numthreads, minmapsize, gsl_fft_forward);
  free (correction_complex);

  gsl_complex_packed_array new_object_fft
      = gal_complex_multiply (correction_fft, psf_fft_conj, size);

  free (correction_fft);

  gsl_complex_packed_array new_object_c
      = gal_fft_two_dimension_transformation (
          new_object_fft, dsize, numthreads, minmapsize, gsl_fft_backward);
  free (new_object_fft);

  double *new_object
      = gal_complex_to_real (new_object_c, size, COMPLEX_TO_REAL_REAL);
  free (new_object_c);

  /* Update object */
  for (size_t i = 0; i < size; i++)
    {
      object[i] *= (pow (new_object[i], alpha));
      if ((object[i] < DBL_EPSILON) || (isnan (object[i]))
          || (isinf (object[i])))
        {
          object[i] = DBL_EPSILON;
        }
      if (object[i] > DBL_MAX)
        {
          object[i] = DBL_MAX;
        }

      newenergy += object[i];
    }
  for (size_t i = 0; i < size; i++)
    {
      object[i] *= energy / newenergy;
    }
  free (new_object);
}

/**
 * @brief Subfunction to update the correction term. Only for
 * internal use.
 *
 * @param correctionterm In/out parameter.
 * @param projection_wave The projection in a specific wave.
 * @param mask The statistical mask in a specific wave.
 * @param modifiedimage_wave Pprime in a specific wave.
 * @param projection The object and psf convolution.
 */
void
deconvolve_AWMLE_increment_correction_term (double *correctionterm,
                                            gal_data_t *projection_wave,
                                            gal_data_t *mask,
                                            gal_data_t *modifiedimage_wave,
                                            gal_data_t *projection)
{
  size_t size = projection_wave->size;
  double *pw_array = projection_wave->array;
  double *mask_array = mask->array;
  double *mi_array = modifiedimage_wave->array;
  double *projection_array = projection->array;
  for (size_t i = 0; i < size; i++)
    {
      if (projection_array[i] > DBL_EPSILON)
        {
          double aux = (mi_array[i] - pw_array[i]) * mask_array[i];
          double aux2 = (pw_array[i] + aux) / projection_array[i];
          if (!isnan (aux2) || !isinf (aux2))
            {
              correctionterm[i] += aux2;
            }
        }
    }
}