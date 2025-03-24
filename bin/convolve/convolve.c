/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2025 Free Software Foundation, Inc.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>

#include <gnuastro/wcs.h>
#include <gnuastro/tile.h>
#include <gnuastro/fits.h>
#include <gnuastro/pointer.h>
#include <gnuastro/threads.h>
#include <gnuastro/convolve.h>
#include <gnuastro/complex.h>

#include <gnuastro-internal/timing.h>

#include "main.h"
#include "convolve.h"

/******************************************************************/
/*************      Padding and initializing      *****************/
/******************************************************************/
void
frequency_make_padded_complex(struct convolveparams *p)
{
  size_t i, ps0, ps1;
  double *o, *op, *pimg, *pker;
  size_t is0=p->input->dsize[0],  is1=p->input->dsize[1];
  size_t ks0=p->kernel->dsize[0], ks1=p->kernel->dsize[1];
  float *f, *ff, *input=p->input->array, *kernel=p->kernel->array;


  /* Find the sizes of the padded image, note that since the kernel
     sizes are always odd, the extra padding on the input image is
     always going to be an even number (clearly divisable). */
  ps0=p->ps0 = p->makekernel ? is0 : is0 + ks0 - 1;
  ps1=p->ps1 = p->makekernel ? is1 : is1 + ks1 - 1;


  /* The Discrete Fourier transforms operate faster on even-sized
     arrays. So if the padded sides are not even, make them so: */
  if(ps0%2) ps0=p->ps0=ps0+1;
  if(ps1%2) ps1=p->ps1=ps1+1;


  /* Allocate the space for the padded input image and fill it. */
  pimg=p->pimg=gal_pointer_allocate(GAL_TYPE_FLOAT64, 2*ps0*ps1, 0,
                                    __func__, "pimg");
  for(i=0;i<ps0;++i)
    {
      op=(o=pimg+i*2*ps1)+2*ps1; /* pimg is complex.            */
      if(i<is0)
        {
          ff=(f=input+i*is1)+is1;
          do {*o++=*f; *o++=0.0f;} while(++f<ff);
        }
      do *o++=0.0f; while(o<op);
    }


  /* Allocate the space for the padded Kernel and fill it. */
  pker=p->pker=gal_pointer_allocate(GAL_TYPE_FLOAT64, 2*ps0*ps1, 0,
                                    __func__, "pker");
  for(i=0;i<ps0;++i)
    {
      op=(o=pker+i*2*ps1)+2*ps1; /* pker is complex.            */
      if(i<ks0)
        {
          ff=(f=kernel+i*ks1)+ks1;
          do {*o++=*f; *o++=0.0f;} while(++f<ff);
        }
      do *o++=0.0f; while(o<op);
    }
}





/*  Remove the padding from the final convolved image and also correct for
    roundoff errors.

    NOTE: The padding to the input image (on the first axis for example)
          was 'p->kernel->dsize[0]-1'. Since 'p->kernel->dsize[0]' is
          always odd, the padding will always be even.  */
void
removepaddingcorrectroundoff(struct convolveparams *p)
{
  size_t ps1=p->ps1;
  size_t *isize=p->input->dsize;
  float *o, *input=p->input->array;
  double *d, *df, *start, *rpad=p->rpad;
  size_t i, hi0, hi1, mkwidth=2*p->makekernel-1;

  /* Set all the necessary parameters to crop the desired region. hi0 and
     hi1 are the coordinates of the first pixel in the output image. In the
     case of deconvolution, if the maximum radius is larger than the input
     image, we will also only be using region that contains non-zero rows
     and columns. */
  if(p->makekernel)
    {
      hi0      = mkwidth < isize[0] ? p->ps0/2-p->makekernel : 0;
      hi1      = mkwidth < isize[1] ? p->ps1/2-p->makekernel : 0;
      isize[0] = mkwidth < isize[0] ? 2*p->makekernel-1 : isize[0];
      isize[1] = mkwidth < isize[1] ? 2*p->makekernel-1 : isize[1];
    }
  else
    {
      hi0 = ( p->kernel->dsize[0] - 1 )/2;
      hi1 = ( p->kernel->dsize[1] - 1 )/2;
    }

  /* To start with, 'start' points to the first pixel in the final
     image: */
  start=&rpad[hi0*ps1+hi1];
  for(i=0;i<isize[0];++i)
    {
      o = &input[ i * isize[1] ];

      df = ( d = start + i * ps1 ) + isize[1];
      do
        *o++ = ( *d<-CONVFLOATINGPOINTERR || *d>CONVFLOATINGPOINTERR )
          ? *d
          : 0.0f;
      while (++d<df);
    }
}

/* Unfortunately I don't understand why the division operation in
   deconvolution (makekernel) does not produce a centered image, the
   image is translated by half the input size in both dimensions. So I
   am correcting this in the spatial domain here. */
void
correctdeconvolve(struct convolveparams *p, double **spatial)
{
  double r, *s, *n, *d, *df, sum=0.0f;
  size_t i, j, ps0=p->ps0, ps1=p->ps1;
  int ii, jj, ci=p->ps0/2-1, cj=p->ps1/2-1;

  /* Check if the image has even sides. */
  if(ps0%2 || ps1%2)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s. The padded "
          "image sides are not an even number", __func__, PACKAGE_BUGREPORT);

  /* First convert the complex image to a real image: */
  s = gal_complex_to_real (p->pimg, ps0 * ps1, COMPLEX_TO_REAL_REAL);

  /* Allocate the array to keep the new values. */
  errno=0;
  n=malloc(ps0*ps1*sizeof *n);
  if(n==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'n'",
          __func__, ps0*ps1*sizeof *n);


  /* Put the elements in their proper place: For example in one
     dimension where the values are actually the true distances:

        s[0]=0, s[1]=1, s[2]=2, s[3]=3, s[4]=4, s[5]=5

     We want the value 0 to be in the 'center'. Note that 's' is
     periodic, for example the next 6 elements have distances:

        s[6]=0, s[7]=1, s[8]=2, s[9]=3, s[10]=4, s[11]=5

     So a 'center'ed array would be like:

        s[0]=4, s[1]=5, s[2]=0, s[3]=1, s[4]=2, s[5]=3

     The relations between the old (i and j) and new (ii and jj) come
     from something like the above line.
   */
  for(i=0;i<ps0;++i)
    {
      ii= i>ps0/2 ? i-(ps0/2+1) : i+ps0/2-1;
      for(j=0;j<ps1;++j)
        {
          jj = j>ps1/2 ? j-(ps1/2+1) : j+ps1/2-1;

          r=sqrt( (ii-ci)*(ii-ci) + (jj-cj)*(jj-cj) );
          sum += n[ii*ps1+jj] = r < p->makekernel ? s[i*ps1+j] : 0;

          /*printf("(%zu, %zu) --> (%zu, %zu)\n", i, j, ii, jj); */
        }
    }


  /* Divide all elements by the sum so the kernel is normalized: */
  df=(d=n)+ps0*ps1; do *d++/=sum; while(d<df);


  /* Clean up. */
  free(s);
  *spatial=n;
}

void
convolve_frequency(struct convolveparams *p)
{
  double *tmp;
  size_t dsize[2];
  size_t total_size;
  struct timeval t1;
  gal_data_t *data = NULL;

  // Pointer to the image in freq domain
  gsl_complex_packed_array image_frequency;
  // Pointer to the kernel in freq domain
  gsl_complex_packed_array kernel_frequency;
  // Pointer to the result (multiplication or division) in freq domain
  gsl_complex_packed_array result_frequency;
  // Pointer to the result (multiplication or division) in space domain
  gsl_complex_packed_array result;

  /* Make the padded arrays. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);

  frequency_make_padded_complex(p);
  dsize[0] = p->ps0;
  dsize[1] = p->ps1;
  total_size = dsize[0] * dsize[1];

  if(!p->cp.quiet)
    gal_timing_report(&t1, "Input and Kernel images padded.", 1);
  if(p->checkfreqsteps)
    {
    /* Prepare the data structure for viewing the steps, note that we
       don't need the array that is initially made. */
    data = gal_data_alloc (NULL, GAL_TYPE_FLOAT64, 2, dsize, NULL, 0,
                           p->cp.minmapsize, p->cp.quietmmap, NULL, NULL, NULL);
    free (data->array);

    /* Save the padded input image. */
    tmp = gal_complex_to_real (p->pimg, total_size, COMPLEX_TO_REAL_REAL);
    data->array = tmp;
    data->name = "input padded";
    gal_fits_img_write (data, p->freqstepsname, NULL, 0);
    free (tmp);
    data->name = NULL;

    /* Save the padded kernel image. */
    tmp = gal_complex_to_real (p->pker, total_size, COMPLEX_TO_REAL_REAL);
    data->array = tmp;
    data->name = "kernel padded";
    gal_fits_img_write (data, p->freqstepsname, NULL, 0);
    free (tmp);
    data->name = NULL;
    }

    /* Forward 2D FFT on each image. */
    if (!p->cp.quiet)
      gettimeofday (&t1, NULL);

    image_frequency = gal_fft_two_dimension_transformation (
        p->pimg, dsize, p->cp.numthreads, p->cp.minmapsize, gsl_fft_forward);
    free (p->pimg);
    p->pimg = image_frequency;

    kernel_frequency = gal_fft_two_dimension_transformation (
        p->pker, dsize, p->cp.numthreads, p->cp.minmapsize, gsl_fft_forward);
    free (p->pker);
    p->pker = kernel_frequency;

    if (!p->cp.quiet)
      gal_timing_report (&t1, "Images converted to frequency domain.", 1);

    if (p->checkfreqsteps) {
      tmp = gal_complex_to_real (p->pimg, total_size, COMPLEX_TO_REAL_REAL);
      data->array = tmp;
      data->name = "input transformed";
      gal_fits_img_write (data, p->freqstepsname, NULL, 0);
      free (tmp);
      data->name = NULL;

      tmp = gal_complex_to_real (p->pker, total_size, COMPLEX_TO_REAL_REAL);
      data->array = tmp;
      data->name = "kernel transformed";
      gal_fits_img_write (data, p->freqstepsname, NULL, 0);
      free (tmp);
      data->name = NULL;
    }

  /* Multiply or divide the two arrays and save them in the output.*/
  if(!p->cp.quiet) gettimeofday(&t1, NULL);

  if(p->makekernel)
    {
    // Deconvolution
    result_frequency
        = gal_complex_divide (p->pimg, p->pker, total_size, p->minsharpspec);
    free (p->pimg);
    p->pimg = result_frequency;

    if (!p->cp.quiet)
      gal_timing_report (&t1, "Divided in the frequency domain.", 1);
    }
  else
    {
      // Convolution
      result_frequency = gal_complex_multiply (p->pimg, p->pker, total_size);
      free (p->pimg);
      p->pimg = result_frequency;

      if(!p->cp.quiet)
        gal_timing_report(&t1, "Multiplied in the frequency domain.", 1);
    }
  if(p->checkfreqsteps)
    {
    tmp = gal_complex_to_real (p->pimg, total_size, COMPLEX_TO_REAL_REAL);
    data->array = tmp;
    data->name = p->makekernel ? "Divided" : "Multiplied";
    gal_fits_img_write (data, p->freqstepsname, NULL, 0);
    free (tmp);
    data->name = NULL;
    }

  /* Forward (in practice inverse) 2D FFT on each image. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);

  result = gal_fft_two_dimension_transformation (
      p->pimg, dsize, p->cp.numthreads, p->cp.minmapsize, gsl_fft_backward);
  free (p->pimg);
  p->pimg = result;

  if (p->makekernel) {
    correctdeconvolve (p, &p->rpad);
  } else {
    p->rpad = gal_complex_to_real (p->pimg, total_size, COMPLEX_TO_REAL_REAL);
  }

  if(!p->cp.quiet)
    gal_timing_report(&t1, "Converted back to the spatial domain.", 1);
  if(p->checkfreqsteps)
    {
      data->array=p->rpad; data->name="padded output";
      gal_fits_img_write(data, p->freqstepsname, NULL, 0);
      data->name=NULL; data->array=NULL;
    }

  /* Free the padded arrays (they are no longer needed) and put the
     converted array (that is real, not complex) in p->pimg. */
  gal_data_free(data);
  free(p->pimg);
  free(p->pker);

  /* Crop out the center, numbers smaller than 10^{-17} are errors,
     remove them. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  removepaddingcorrectroundoff(p);
  if (!p->cp.quiet)
    gal_timing_report (&t1, "Padded parts removed.", 1);
}





void
convolve_spatial(struct convolveparams *p)
{
  gal_data_t *out, *check;
  int multidim=p->input->ndim>1;
  struct gal_options_common_params *cp=&p->cp;


  /* Prepare the mesh structure. */
  if(multidim) gal_tile_full_two_layers(p->input, &cp->tl);

  /* Save the tile IDs if they are requested. */
  if(multidim && cp->tl.tilecheckname)
    {
      check=gal_tile_block_check_tiles(cp->tl.tiles);
      gal_fits_img_write(check, cp->tl.tilecheckname, NULL, 0);
      gal_data_free(check);
    }

  /* Do the spatial convolution. One of the main reason someone would
     want to do spatial domain convolution with this Convolve program
     is edge correction. So by default we assume it and will only
     ignore it if the user asks. */
  out=gal_convolve_spatial(multidim ? cp->tl.tiles : p->input,
                           p->kernel,
                           cp->numthreads,
                           multidim ? !p->noedgecorrection : 1,
                           multidim ? cp->tl.workoverch : 1,
                           p->conv_on_blank);

  /* Clean up: free the actual input and replace it's pointer with the
     convolved dataset to save as output. */
  gal_tile_full_free_contents(&cp->tl);
  gal_data_free(p->input);
  p->input=out;
}




















/******************************************************************/
/*************         Top-level function         *****************/
/******************************************************************/
void
convolve(struct convolveparams *p)
{
  struct gal_options_common_params *cp=&p->cp;

  /* Do the convolution. */
  if(p->domain==CONVOLVE_DOMAIN_SPATIAL) convolve_spatial(p);
  else                                   convolve_frequency(p);

  /* Write Convolve's parameters as keywords into the first extension of
     the output. */
  if( gal_fits_name_is_fits(p->cp.output) )
    {
      gal_fits_key_write_filename("input", p->filename, &cp->ckeys, 1,
                                  cp->quiet);
      gal_fits_key_write(cp->ckeys, cp->output, "0", "NONE", 1, 1);
    }

  /* Save the output (which is in p->input) array. */
  if(p->input->ndim==1)
    gal_table_write(p->input, NULL, NULL, p->cp.tableformat, p->cp.output,
                    "CONVOLVED", 0, 0);
  else
    gal_fits_img_write_to_type(p->input, cp->output, NULL, cp->type, 0);

  /* Inform the user that the job is done. */
  if(!p->cp.quiet)
    printf("  - Output: %s\n", p->cp.output);
}
