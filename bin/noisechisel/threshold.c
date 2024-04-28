/*********************************************************************
NoiseChisel - Detect signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2024 Free Software Foundation, Inc.

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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
#include <gnuastro/statistics.h>
#include <gnuastro/permutation.h>
#include <gnuastro/interpolate.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/tile-internal.h>

#include "main.h"

#include "ui.h"
#include "threshold.h"










/**********************************************************************/
/***************        Apply a given threshold.      *****************/
/**********************************************************************/
struct threshold_apply_p
{
  float               *value1;
  float               *value2;
  int                    kind;
  struct noisechiselparams *p;
};




/* Apply the threshold on the tiles given to this thread. */
static void *
threshold_apply_on_thread(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct threshold_apply_p *taprm=(struct threshold_apply_p *)(tprm->params);
  struct noisechiselparams *p=taprm->p;

  size_t i, tid;
  void *tarray=NULL;
  gal_data_t *tile, *tblock=NULL;
  float *value1=taprm->value1, *value2=taprm->value2;

  /* Go over all the tiles assigned to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading. */
      tid=tprm->indexs[i];
      tile=&p->cp.tl.tiles[tid];

      /* Based on the kind of threshold. */
      switch(taprm->kind)
        {

        /* This is a quantile threshold. */
        case THRESHOLD_QUANTILES:
          /* Correct the tile's pointers to apply the threshold on the
             convolved image. */
          if(p->conv)
            {
              tarray=tile->array; tblock=tile->block;
              tile->array=gal_tile_block_relative_to_other(tile, p->conv);
              tile->block=p->conv;
            }

          /* Apply the threshold: When the '>' comparison fails, it can be
             either because the pixel was actually smaller than the
             threshold, or that it was a NaN value. In the first case,
             return 0, in the second, return a blank value.

             We already know if a tile contains a blank value (which is a
             constant over the whole loop). So before checking if the value
             is blank, see if the tile actually has a blank value. This
             will help in efficiency, because the compiler can move this
             check out of the loop and only check for NaN values when we
             know the tile has blank pixels. */
          GAL_TILE_PO_OISET(float, uint8_t, tile, p->binary, 1, 0, {
              *o = ( *i > value1[tid]
                     ? ( *i > value2[tid] ? THRESHOLD_NO_ERODE_VALUE : 1 )
                     : ( (tile->flag & GAL_DATA_FLAG_HASBLANK) && !(*i==*i)
                         ? GAL_BLANK_UINT8 : 0 ) );
            });

          /* Revert the tile's pointers back to what they were. */
          if(p->conv) { tile->array=tarray; tile->block=tblock; }
          break;


        /* This is a Sky and Sky STD threshold. */
        case THRESHOLD_SKY_STD:

          /* See the explanation above the same step in the quantile
             threshold for an explanation. */
          GAL_TILE_PO_OISET(float, uint8_t, tile, p->binary, 1, 0, {
              *o = ( ( *i - value1[tid] > p->dthresh * value2[tid] )
                     ? 1
                     : ( (tile->flag & GAL_DATA_FLAG_HASBLANK) && !(*i==*i)
                         ? GAL_BLANK_UINT8 : 0 ) );
            });
          break;


        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so we "
                "can address the problem. A value of %d had for "
                "'taprm->kind' is not valid", __func__, PACKAGE_BUGREPORT,
                taprm->kind);
        }
    }

  /* Wait until all the other threads finish. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* Apply a given threshold threshold on the tiles. */
void
threshold_apply(struct noisechiselparams *p, float *value1,
                float *value2, int kind)
{
  struct threshold_apply_p taprm={value1, value2, kind, p};
  gal_threads_spin_off(threshold_apply_on_thread, &taprm,
                       p->cp.tl.tottiles, p->cp.numthreads,
                       p->cp.minmapsize, p->cp.quietmmap);
}



















/**********************************************************************/
/***************            Write S/N values          *****************/
/**********************************************************************/
void
threshold_write_sn_table(struct noisechiselparams *p, gal_data_t *insn,
                         gal_data_t *inind, char *filename,
                         gal_list_str_t *comments, char *extname)
{
  gal_data_t *sn, *ind, *cols;

  /* Remove all blank elements. The index and sn values must have the same
     set of blank elements, but checking on the integer array is faster. */
  if( gal_blank_present(inind, 1) )
    {
      /* Remove blank elements. */
      ind=gal_data_copy(inind);
      sn=gal_data_copy(insn);
      gal_blank_remove(ind);
      gal_blank_remove(sn);
    }
  else
    {
      sn  = insn;
      ind = inind;
    }

  /* Set the columns. */
  cols       = ind;
  cols->next = sn;


  /* Prepare the comments. */
  gal_table_comments_add_intro(&comments, PROGRAM_STRING, &p->rawtime);


  /* Write the table. Note that we'll set the 'dontdelete' argument to 0
     because when the output is a FITS table, we want all the tables in one
     FITS file. We have already deleted any existing file with the same
     name in 'ui_set_output_names'.*/
  gal_table_write(cols, NULL, comments, p->cp.tableformat, filename,
                  extname, 0, 0);


  /* Clean up (if necessary). */
  if(sn!=insn) gal_data_free(sn);
  if(ind==inind) ind->next=NULL; else gal_data_free(ind);
}



















/**********************************************************************/
/***************     Interpolation and smoothing      *****************/
/**********************************************************************/
/* Interpolate and smooth the values for each tile over the whole image. */
void
threshold_interp_smooth(struct noisechiselparams *p, gal_data_t **first,
                        gal_data_t **second, gal_data_t **third,
                        char *filename)
{
  gal_data_t *tmp;
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;

  /* A small sanity check. */
  if( (*first)->next )
    error(EXIT_FAILURE, 0, "%s: 'first' must not have any 'next' pointer.",
          __func__);
  if( (*second)->next )
    error(EXIT_FAILURE, 0, "%s: 'second' must not have any 'next' pointer.",
          __func__);
  if( third && (*third)->next )
    error(EXIT_FAILURE, 0, "%s: 'third' must not have any 'next' pointer.",
          __func__);

  /* Do the interpolation of both arrays. */
  (*first)->next = *second;
  if(third) (*second)->next = *third;
  tmp=gal_interpolate_neighbors(*first, tl, cp->interpmetric,
                                cp->interpnumngb, cp->numthreads,
                                cp->interponlyblank, 1,
                                GAL_INTERPOLATE_NEIGHBORS_FUNC_MEDIAN);
  gal_data_free(*first);
  gal_data_free(*second);
  if(third) gal_data_free(*third);
  *first=tmp;
  *second=(*first)->next;
  if(third)
    {
      *third=(*second)->next;
      (*third)->next=NULL;
    }
  (*first)->next=(*second)->next=NULL;
  if(filename)
    {
      (*first)->name="THRESH1_INTERP";
      (*second)->name="THRESH2_INTERP";
      if(third) (*third)->name="THRESH3_INTERP";
      gal_tile_full_values_write(*first, tl, !p->ignoreblankintiles,
                                 filename, NULL, 0);
      gal_tile_full_values_write(*second, tl, !p->ignoreblankintiles,
                                 filename, NULL, 0);
      if(third)
        gal_tile_full_values_write(*third, tl, !p->ignoreblankintiles,
                                   filename, NULL, 0);
      (*first)->name = (*second)->name = NULL;
      if(third) (*third)->name=NULL;
    }

  /* Smooth the threshold if requested. */
  if(p->smoothwidth>1)
    {
      /* Smooth the first. */
      tmp=gal_tile_full_values_smooth(*first, tl, p->smoothwidth,
                                      p->cp.numthreads);
      gal_data_free(*first);
      *first=tmp;

      /* Smooth the second */
      tmp=gal_tile_full_values_smooth(*second, tl, p->smoothwidth,
                                      p->cp.numthreads);
      gal_data_free(*second);
      *second=tmp;

      /* Smooth the third */
      if(third)
        {
          tmp=gal_tile_full_values_smooth(*third, tl, p->smoothwidth,
                                          p->cp.numthreads);
          gal_data_free(*third);
          *third=tmp;
        }

      /* Add them to the check image. */
      if(filename)
        {
          (*first)->name="THRESH1_SMOOTH";
          (*second)->name="THRESH2_SMOOTH";
          if(third) (*third)->name="THRESH3_SMOOTH";
          gal_tile_full_values_write(*first, tl, !p->ignoreblankintiles,
                                     filename, NULL, 0);
          gal_tile_full_values_write(*second, tl, !p->ignoreblankintiles,
                                     filename, NULL, 0);
          if(third)
            gal_tile_full_values_write(*third, tl, !p->ignoreblankintiles,
                                       filename, NULL, 0);
          (*first)->name = (*second)->name = NULL;
          if(third) (*third)->name=NULL;
        }
    }
}




















/****************************************************************
 ************           Quantile threshold           ************
 ****************************************************************/
/* Prepare the 'usage' array from the tile */
static size_t
qthresh_on_tile_usage_prepare(struct noisechiselparams *p,
                              gal_data_t *usage, gal_data_t *tile,
                              gal_data_t *meanconv)
{
  void *tarray=NULL;
  gal_data_t *tblock=NULL;
  size_t i, ndim=p->input->ndim;

  /* Re-initialize the usage array's space (will be changed in
     'gal_data_copy_to_allocated' for each tile). */
  usage->ndim=ndim;
  usage->size=p->maxtcontig;
  memcpy(usage->dsize, p->maxtsize, ndim*sizeof *p->maxtsize);

  /* Temporarily change the tile's pointers so we can do the work on
     the convolved image, then copy the desired contents into the
     already allocated 'usage' array and set the pointers back to what
     they were. */
  tarray=tile->array;
  tblock=tile->block;
  tile->array=gal_tile_block_relative_to_other(tile, meanconv);
  tile->block=meanconv;
  gal_data_copy_to_allocated(tile, usage);
  tile->array=tarray;
  tile->block=tblock;

  /* Adjust the 'usage' array pointers that should not be inherited
     from the tile. We are setting it to 1D because after the clipping,
     the dimensionality is going to be lost anyway. */
  usage->size=1;
  usage->next=NULL;
  for(i=0;i<usage->ndim;++i) usage->size*=usage->dsize[i];
  usage->dsize[0]=usage->size; /* Must be after the loop above. */
  usage->ndim=1;               /* Must be after the loop above. */

  /* Return the size of the final 'usage' array (based on this particular
     tile). */
  return usage->size;
}





/* Calculate the MAD-clipped mean quantile. */
static double
qthresh_on_tile_mean_quant(gal_data_t *usage)
{
  size_t one=1;
  double meanquant;
  gal_data_t *clip, *mean, *meanq;
  uint8_t extrastats=GAL_STATISTICS_CLIP_OUTCOL_OPTIONAL_MEAN;

  /* Do the MAD-clipping in-place. */
  clip=gal_statistics_clip_mad(usage, 4.5, 0.01, extrastats, 1, 1);
  mean=gal_data_alloc(NULL, clip->type, 1, &one, NULL, 0, -1, 1,
                      NULL, NULL, NULL);
  memcpy(mean->array,
         gal_pointer_increment(clip->array,
                               GAL_STATISTICS_CLIP_OUTCOL_MEAN,
                               clip->type),
         gal_type_sizeof(clip->type));
  meanq = ( usage->size
            ? gal_statistics_quantile_function(usage, mean, 1)
            : NULL );
  meanquant=meanq ? *(double *)(meanq->array) : NAN;

  /* Clean up and return. */
  gal_data_free(meanq);
  gal_data_free(clip);
  gal_data_free(mean);
  return meanquant;
}





/* Bring 'usage' in the range from 0 to 1; then find the difference between
   the values at quantile 0.25 to 0.75. In a Uniform distribution, this is
   0.5 (0.75-0.25), but in a Gaussian distribution this is about 0.16
   (because most values are concentrated in the center of the
   distribution). */
static int
qthresh_on_tile_concentrated(gal_data_t *usage, double q_value_diff,
                             size_t tind)
{
  int out=0;
  size_t i, one=1;
  double qf1, qf2;
  float *uarr=usage->array;
  gal_data_t *qf1_d, *qf2_d, *qfinput;
  float umin=uarr[0], umax=uarr[usage->size-1];

  /* Small sanity check. */
  if(usage->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to fix "
          "the problem. The type of the 'usage' array should be float32, "
          "but it is '%s'", __func__, PACKAGE_BUGREPORT,
          gal_type_name(usage->type, 1));

  /* Allocate the dataset to hold the quantile value. */
  qfinput=gal_data_alloc(NULL, usage->type, 1, &one, NULL, 0, -1, 1,
                         NULL, NULL, NULL);

  /* Bring all values between 0 to 1 and estimate the quantiles. */
  for(i=0;i<usage->size;++i) uarr[i]=(uarr[i]-umin)/(umax-umin);
  ((float *)(qfinput->array))[0]=0.25;
  qf1_d=gal_statistics_quantile_function(usage, qfinput, 1);
  ((float *)(qfinput->array))[0]=0.75;
  qf2_d=gal_statistics_quantile_function(usage, qfinput, 1);
  qf1=((double *)(qf1_d->array))[0];
  qf2=((double *)(qf2_d->array))[0];

  /**********************************************/
  /* You are considering dividing by the given width to produce a value of
     1 and larger, like below. When we decrease the central difference (for
     example '--quantile=0.4,0.6', the tile-3 of the 'isdss' image
     separates better from the normal distribution. */

  /* for f in junk-sigma.fits junk-uniform.fits tile-3-i-sdss.txt; do printf "$f: "; aststatistics $f --quantfunc=0.25,0.75 | asttable -c'arith $2 $1 - 0.5 /' -Y; done */
  /**********************************************/


  /* For a check. */
  if(tind!=-1)
    {
      printf("Tile %zu values:\n", tind);
      for(i=0;i<usage->size;++i) printf("\t%f\n", uarr[i]);
    }

  /* Final returned value (placed here so we can check it if necessary):
     also because we need to correct the values. */
  out = qf2-qf1>q_value_diff;
  if(out)
    for(i=0;i<usage->size;++i) uarr[i]=uarr[i]*(umax-umin)+umin;

  /* For a check (to disable the print, set 'tind' to -1). */
  if(tind!=-1)
    {
      printf("%s:%zu: umin: %f, umax: %f\n", __func__, tind, umin, umax);
      printf("%s:%zu: qf1: %f,  qf2: %f, qdiff=%f\n",  __func__, tind,
             qf1, qf2, qf2-qf1);
      printf("%s:%zu: Thresh: %f. Result: %s\n", __func__, tind,
             q_value_diff, out?"PASSED":"FAILED");
      gal_table_write(usage, NULL, NULL, GAL_TABLE_FORMAT_BFITS,
                      "table-usage.fits", NULL, 0, 0);
      exit(0);
    }

  /* Clean up and return. */
  gal_data_free(qfinput);
  gal_data_free(qf1_d);
  gal_data_free(qf2_d);
  return out;
}





struct qthreshparams
{
  gal_data_t        *erode_th;
  gal_data_t      *noerode_th;
  gal_data_t       *expand_th;
  void                 *usage;
  struct noisechiselparams *p;
};





static void
qthresh_on_tile_write(struct qthreshparams *qprm, gal_data_t *usage,
                      size_t tind)
{
  gal_data_t *qvalue;
  int type=qprm->erode_th->type;
  size_t twidth=gal_type_sizeof(type);
  struct noisechiselparams *p=qprm->p;

  /* Get the erosion quantile for this tile and save it. Note that
     the type of 'qvalue' is the same as the input dataset. */
  qvalue=gal_statistics_quantile(usage, p->qthresh, 1);
  memcpy(gal_pointer_increment(qprm->erode_th->array, tind, type),
         qvalue->array, twidth);
  gal_data_free(qvalue);

  /* Same for the no-erode quantile. */
  qvalue=gal_statistics_quantile(usage, p->noerodequant, 1);
  memcpy(gal_pointer_increment(qprm->noerode_th->array, tind, type),
         qvalue->array, twidth);
  gal_data_free(qvalue);

  /* Same for the expansion quantile. */
  if(qprm->expand_th)
    {
      qvalue=gal_statistics_quantile(usage, p->detgrowquant, 1);
      memcpy(gal_pointer_increment(qprm->expand_th->array, tind,
                                   type),
             qvalue->array, twidth);
      gal_data_free(qvalue);
    }
}





static void *
qthresh_on_tile(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct qthreshparams *qprm=(struct qthreshparams *)tprm->params;
  struct noisechiselparams *p=qprm->p;

  size_t initsize;
  double meanquant;
  void *tarray=NULL;
  int type=qprm->erode_th->type;
  size_t i, tind, ndim=p->input->ndim;
  gal_data_t *tile, *usage, *tblock=NULL;
  gal_data_t *meanconv = p->wconv ? p->wconv : p->conv;

  /* Put the temporary usage space for this thread into a data set for easy
     processing. */
  usage=gal_data_alloc(gal_pointer_increment(qprm->usage,
                                             tprm->id*p->maxtcontig, type),
                       type, ndim, p->maxtsize, NULL, 0, p->cp.minmapsize,
                       p->cp.quietmmap, NULL, NULL, NULL);


  /* Go over all the tiles given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* Copy this tile's data into the already allocated 'usage' array
         which we can comfortably (without editing the original data)
         change, reorder and etc. */
      tind = tprm->indexs[i];
      tile=&p->cp.tl.tiles[tind];
      initsize=qthresh_on_tile_usage_prepare(p, usage, tile, meanconv);

      /* For a check...
      char *plab=NULL;
      switch(tind)
        {
        case -1: plab="BAD"; break;
      //case 11031: case 14951: case 27470: case 27673: plab="GOOD"; break;
        }
      if(plab)
        {
          char tname[100];
          sprintf(tname, "table-%zu-before.fits", tind);
          gal_table_write(usage, NULL, NULL, GAL_TABLE_FORMAT_BFITS, tname,
                          NULL, 0, 0);
        }
      */

      /* Find the mean's quantile after clipping inplace. */
      meanquant=qthresh_on_tile_mean_quant(usage);

      /* For a check: to have each tile that you want to study in a
         separate FITS file, you can use the commands below:

         $ astnoisechisel img.fits --checkqthresh --checktiles \
                          --continueaftercheck
         $ ids="20904 21324 21527 11031 14951 27470 27673"
         $ for i in $ids; do
             astarithmetic img_tiles.fits   -h2 set-l \
                           img_qthresh.fits -h1 set-i \
                           i l $i ne nan where trim -o crop-$i.fits;
           done

      if(plab)
        {
          char tname[100];
          printf("%s:%s:%zu: %f (before: %zu, after: %zu). result: %d\n",
                 __func__, plab, tind, meanq, beforeclip, usage->size,
                 (meanq<0.5f && meanq>0.5-p->meanmedqdiff) );
          sprintf(tname, "table-%zu-after.fits", tind);
          gal_table_write(usage, NULL, NULL, GAL_TABLE_FORMAT_BFITS, tname,
                          NULL, 0, 0);
        }
      */

      /* Only continue if the mean's quantile is close enough to the
         median.  */
      if( meanquant
          && meanquant<0.5f
          && meanquant>0.5f-p->meanmedqdiff
          && (float)usage->size/(float)initsize > 0.75
          && qthresh_on_tile_concentrated(usage, 0.8, tind) )
        {

          /* The mean was found on the wider convolved image, but the
             qthresh values have to be found on the sharper convolved
             images. This is because the distribution becomes more skewed
             with a wider kernel, helping us find tiles with no data more
             easily. But for the quantile threshold, we want to use the
             sharper convolved image to loose less of the spatial
             information. */
          if(meanconv!=p->conv)
            {
              /* Corrections with MAD clipping have not been implemented
                 yet. */
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at "
                    "'%s' to fix the problem. The newly added "
                    "corrections to the quantile threshold have not yet "
                    "been implemented with the '--widekernel' option",
                    __func__, PACKAGE_BUGREPORT);
              tarray=tile->array; tblock=tile->block;
              tile->array=gal_tile_block_relative_to_other(tile, p->conv);
              tile->block=p->conv;
              usage->ndim=ndim;           /* Since usage was modified in */
              usage->size=p->maxtcontig;  /* place, it needs to be       */
              gal_data_copy_to_allocated(tile, usage);/* re-initialized. */
              tile->array=tarray; tile->block=tblock;
            }

          /* Calculate the desired quantile and write them in the
             output. */
          qthresh_on_tile_write(qprm, usage, tind);
        }
      else
        {
          gal_blank_write(gal_pointer_increment(qprm->erode_th->array,
                                                 tind, type), type);
          gal_blank_write(gal_pointer_increment(qprm->noerode_th->array,
                                                 tind, type), type);
          if(qprm->expand_th)
            gal_blank_write(gal_pointer_increment(qprm->expand_th->array,
                                                   tind, type), type);
        }
    }

  /* Clean up and wait for the other threads to finish, then return. */
  usage->array=NULL;  /* Not allocated here. */
  gal_data_free(usage);
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





static void
threshold_good_error(size_t number, int before0_after1,
                     size_t interpnumngb)
{
  before0_after1=1;

  /* Set the differing strings. */
  char *in1 = ( before0_after1
                ? "after removing outliers"
                : "for defining a quantile threshold" );
  char *in2 = ( before0_after1
                ? ""
                : "NOTE that this is happening *BEFORE* outlier rejection "
                  "(where the number may decrease even further).");
  char *in3 = ( before0_after1
                ? "\n"
                  "  - (slightly) Decrease '--outliernumngb' to use less "
                  "tiles to find outliers.\n"
                  "  - (slightly) Increase '--outliersclip' to reject less "
                  "as outliers.\n"
                  "  - (slightly) Increase '--outliersigma' to reject less "
                  "as outliers.\n"
                : "\n");

  /* Print the error message and abort. */
  error(EXIT_FAILURE, 0, "%zu tiles usable %s!\n\n"

        "This is smaller than the requested minimum value of %zu (value to "
        "the '--interpnumngb' option). %s\n\n"

        "There are several ways to address the problem. The best and most "
        "highly recommended is to use a larger input if possible (when the "
        "input is a crop from a larger dataset). If this is not the case, "
        "or it doesn't solve the problem, you need to loosen the "
        "parameters mentioned below in the respective order (and therefore "
        "cause scatter/inaccuracy in the final result). Hence its best to "
        "not loosen them too much (recall that you can see all the option "
        "values to Gnuastro's programs by appending '-P' to the end of your "
        "command).\n"
        "  - (slightly) Decrease '--tilesize' so your tile-grid has more "
        "tiles.\n"
        "  - (slightly) Increase '--meanmedqdiff' to accept more tiles.%s"
        "  - (slightly) Decrease '--interpnumngb' to be less than %zu.\n\n"

        "---- Tip ----\n"
        "Append your command with '--checkqthresh' to see the "
        "successful tiles in relation with this dataset's contents "
        "before this crash. A visual inspection will greatly help in "
        "finding the cause/solution for this particular dataset (note "
        "that the output of '--checkqthresh' is a multi-extension FITS "
        "file).\n\n"
        "To better understand this important step, please run the "
        "following command (press 'SPACE'/arrow-keys to navigate and "
        "'Q' to return back to the command-line):\n\n"
        "    $ info gnuastro \"Quantifying signal in a tile\"\n", number,
        in1, interpnumngb, in2, in3, number);
}





void
threshold_quantile_find_apply(struct noisechiselparams *p)
{
  char *msg;
  size_t nval;
  gal_data_t *num;
  struct timeval t1;
  struct qthreshparams qprm;
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;


  /* Get the starting time if necessary. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);


  /* Add image to check image if requested. If the user has asked for
     'oneelempertile', then the size of values is not going to be the same
     as the input, making it hard to inspect visually. So we'll only put
     the full input when 'oneelempertile' isn't requested. */
  if(p->qthreshname && !tl->oneelempertile)
    {
      gal_fits_img_write(p->conv ? p->conv : p->input, p->qthreshname,
                         NULL, 0);
      if(p->wconv)
        gal_fits_img_write(p->wconv ? p->wconv : p->input, p->qthreshname,
                           NULL, 0);
    }


  /* Allocate space for the quantile threshold values. */
  qprm.erode_th=gal_data_alloc(NULL, p->input->type, p->input->ndim,
                               tl->numtiles, NULL, 0, cp->minmapsize,
                               p->cp.quietmmap, NULL, p->input->unit,
                               NULL);
  qprm.noerode_th=gal_data_alloc(NULL, p->input->type, p->input->ndim,
                                 tl->numtiles, NULL, 0, cp->minmapsize,
                                 p->cp.quietmmap, NULL, p->input->unit,
                                 NULL);
  qprm.expand_th = ( p->detgrowquant!=1.0f
                     ? gal_data_alloc(NULL, p->input->type, p->input->ndim,
                                      tl->numtiles, NULL, 0, cp->minmapsize,
                                      p->cp.quietmmap, NULL, p->input->unit,
                                      NULL)
                     : NULL );


  /* Allocate temporary space for processing in each tile. */
  qprm.usage=gal_pointer_allocate(p->input->type,
                                  cp->numthreads * p->maxtcontig, 0,
                                  __func__, "qprm.usage");

  /**********************************
  cp->numthreads=1;
  **********************************/

  /* Find the threshold on each tile, free the temporary processing space
     and set the blank flag on both. Since they have the same blank
     elements, it is only necessary to check one (with the 'updateflag'
     value set to 1), then update the next. */
  qprm.p=p;
  gal_threads_spin_off(qthresh_on_tile, &qprm, tl->tottiles,
                       cp->numthreads, cp->minmapsize,
                       cp->quietmmap);
  free(qprm.usage);
  if( gal_blank_present(qprm.erode_th, 1) )
    {
      qprm.noerode_th->flag |= GAL_DATA_FLAG_HASBLANK;
      if(qprm.expand_th) qprm.expand_th->flag  |= GAL_DATA_FLAG_HASBLANK;
    }
  qprm.noerode_th->flag |= GAL_DATA_FLAG_BLANK_CH;
  if(qprm.expand_th) qprm.expand_th->flag  |= GAL_DATA_FLAG_BLANK_CH;
  if(p->qthreshname)
    {
      qprm.erode_th->name="QTHRESH_ERODE";
      qprm.noerode_th->name="QTHRESH_NOERODE";
      gal_tile_full_values_write(qprm.erode_th, tl,
                                 !p->ignoreblankintiles,
                                 p->qthreshname, NULL, 0);
      gal_tile_full_values_write(qprm.noerode_th, tl,
                                 !p->ignoreblankintiles,
                                 p->qthreshname, NULL, 0);
      qprm.erode_th->name=qprm.noerode_th->name=NULL;

      if(qprm.expand_th)
        {
          qprm.expand_th->name="QTHRESH_EXPAND";
          gal_tile_full_values_write(qprm.expand_th, tl,
                                     !p->ignoreblankintiles,
                                     p->qthreshname, NULL, 0);
          qprm.expand_th->name=NULL;
        }
    }


  /* Remove the outliers. */
  if(p->outliernumngb)
    gal_tileinternal_no_outlier_local(qprm.erode_th, qprm.noerode_th,
                                      qprm.expand_th, &cp->tl,
                                      cp->interpmetric, p->outliernumngb,
                                      cp->numthreads, p->outliersclip,
                                      p->outliersigma, p->qthreshname,
                                      "--outliernumngb");


  /* Use the no-outlier grid as a basis for later estimating the sky. To
     see this array on the image, use 'gal_tile_full_values_write'. */
  p->noskytiles=gal_blank_flag(qprm.erode_th);
  /* For a check:
  gal_tile_full_values_write(p->noskytiles, &cp->tl, 1,
                             "noskytiles.fits", NULL, NULL);
  */


  /* Check if the number of acceptable tiles is more than the minimum
     interpolated number. Since this is a common problem for users, it is
     much more useful to do the check here rather than printing multiple
     errors in parallel. */
  num=gal_statistics_number(qprm.erode_th);
  nval=((size_t *)(num->array))[0];
  if( nval < cp->interpnumngb )
    threshold_good_error(nval, 1, cp->interpnumngb);
  gal_data_free(num);


  /* Interpolate and smooth the derived values. */
  threshold_interp_smooth(p, &qprm.erode_th, &qprm.noerode_th,
                          qprm.expand_th ? &qprm.expand_th : NULL,
                          p->qthreshname);


  /* We now have a threshold for all tiles, apply it. */
  threshold_apply(p, qprm.erode_th->array, qprm.noerode_th->array,
                  THRESHOLD_QUANTILES);


  /* Write the binary image if check is requested. */
  if(p->qthreshname && !tl->oneelempertile)
    {
      p->binary->name="QTHRESH-APPLIED";
      gal_fits_img_write(p->binary, p->qthreshname, NULL, 0);
      p->binary->name=NULL;
    }


  /* Set the expansion quantile if necessary. */
  p->expand_thresh = qprm.expand_th ? qprm.expand_th : NULL;


  /* Clean up and report duration if necessary. */
  gal_data_free(qprm.erode_th);
  gal_data_free(qprm.noerode_th);
  if(!p->cp.quiet)
    {
      if( asprintf(&msg, "%.2f & %0.2f quantile thresholds applied.",
                   p->qthresh, p->noerodequant)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }


  /* If the user wanted to check the threshold and hasn't called
     'continueaftercheck', then stop NoiseChisel. */
  if(p->qthreshname && !p->continueaftercheck)
    ui_abort_after_check(p, p->qthreshname, NULL, "quantile threshold check");
}
