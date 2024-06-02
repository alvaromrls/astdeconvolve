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
      gal_blank_remove(ind, 0);
      gal_blank_remove(sn, 0);
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
  tmp=gal_interpolate_neighbors(*first, tl, p->interpmetric,
                                p->interpnumngb, cp->numthreads,
                                p->interponlyblank, 1,
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
      gal_tile_full_values_write(*first, tl, !p->ignoreblankintiles,
                                 filename, NULL, 0);
      (*first)->name=NULL;
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
          gal_tile_full_values_write(*first, tl, !p->ignoreblankintiles,
                                     filename, NULL, 0);
          (*first)->name=NULL;
        }
    }
}




















/****************************************************************
 ************           Quantile threshold           ************
 ****************************************************************/
struct qthreshparams
{
  gal_data_t        *erode_th;
  gal_data_t      *noerode_th;
  gal_data_t       *expand_th;
  void                 *usage;
  struct noisechiselparams *p;
};





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





/* See if the tile's distribution is concentrated or not. */
static int
qthresh_on_tile_concentrated(gal_data_t *usage, double width,
                             double thresh, size_t tind)
{
  int out=0;
  double *m;
  gal_data_t *measured;

  /* Small sanity check. */
  if(usage->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
          "fix the problem. The type of the 'usage' array should be "
          "float32, but it is '%s'", __func__, PACKAGE_BUGREPORT,
          gal_type_name(usage->type, 1));

  /* Measure the concentration. */
  measured=gal_statistics_concentration(usage, width, 1);
  m=measured->array;

  /* See if it is above the threshold or not. */
  out = m[0] > thresh;

  /* Clean up and return. */
  gal_data_free(measured);
  return out;
}





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
  void *tarray=NULL;
  int type=qprm->erode_th->type;
  size_t i, tind, ndim=p->input->ndim;
  gal_data_t *tile, *usage, *tblock=NULL;
  double meanquant, *concent=p->concentration->array;
  gal_data_t *meanconv = p->wconv ? p->wconv : p->conv;

  /* Put the temporary usage space for this thread into a dataset for easy
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

      /* Find the mean's quantile after clipping inplace. */
      meanquant=qthresh_on_tile_mean_quant(usage);

      /* Only continue when: 1) the mean's quantile is below the median,
         but not too much (close enough to the median). 2) The faction of
         usable pixels is not too small. 3) the flux distribution is
         concentrated. */
      if(    meanquant<0.5f+p->meanmedqdiff
          && meanquant>0.5f-p->meanmedqdiff
          && (float)usage->size/(float)initsize > p->minskyfrac
          && qthresh_on_tile_concentrated(usage, concent[0], concent[1],
                                          tind) )
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
threshold_quantilf_prepare(struct noisechiselparams *p,
                           struct qthreshparams *qprm)
{
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;

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
  qprm->p=p;
  qprm->erode_th=gal_data_alloc(NULL, p->input->type, p->input->ndim,
                                tl->numtiles, NULL, 0, cp->minmapsize,
                                p->cp.quietmmap, NULL, p->input->unit,
                                NULL);
  qprm->noerode_th=gal_data_alloc(NULL, p->input->type, p->input->ndim,
                                  tl->numtiles, NULL, 0, cp->minmapsize,
                                  p->cp.quietmmap, NULL, p->input->unit,
                                  NULL);
  qprm->expand_th = ( p->detgrowquant!=1.0f
                      ? gal_data_alloc(NULL, p->input->type, p->input->ndim,
                                       tl->numtiles, NULL, 0, cp->minmapsize,
                                       p->cp.quietmmap, NULL, p->input->unit,
                                       NULL)
                      : NULL );


  /* Allocate temporary space for processing in each tile. */
  qprm->usage=gal_pointer_allocate(p->input->type,
                                  cp->numthreads * p->maxtcontig, 0,
                                  __func__, "qprm.usage");

}





static int
threshold_quantilf_find_tiles_good(struct noisechiselparams *p,
                                   struct qthreshparams *qprm)
{
  size_t tc, cc, numgood;
  float *f=qprm->erode_th->array;
  struct gal_tile_two_layer_params *tl=&p->cp.tl;
  int out=1; /* Assume the output is good and change it if necessary. */

  /* The input should be float32, if not, it is a bug! */
  if(p->input->type!=GAL_TYPE_FLOAT32)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to fix "
          "the problem. The input dataset to this function should have "
          "a float32 type, but its type is %s", __func__,
          PACKAGE_BUGREPORT, gal_type_name(p->input->type, 1));

  /* It is mandatory that all channels have at least one non-blank
     tile. */
  for(cc=0; cc<tl->totchannels; ++cc)
    {
      /* Initialize the number of good tiles to zero. */
      numgood=0;

      /* Go over all the tiles in this chanel. */
      for(tc=0; tc<tl->tottilesinch; ++tc)
        {
          /* Increment the number of good tiles if the tile's value is not
             NaN. */
          numgood += isnan(f[tc+cc*tl->tottilesinch])==0;

          /* For a check
          printf("%s: (%zu,%zu) %f; numgood: %zu\n", __func__, cc, tc,
                 f[tc+cc*tl->tottilesinch], numgood);
          */
        }

      /* If there weren't any good tiles, then set the output and abort
         (even if a single channel doesn't have enough tiles, the over-all
         result is bad). */
      if(numgood==0) {out=0; break;}
    }

  /* Return the output value. */
  return out;
}





/* Find the threshold on each tile, free the temporary processing space and
   set the blank flag on both. Since they have the same blank elements, it
   is only necessary to check one (with the 'updateflag' value set to 1),
   then update the next. */
static void
threshold_quantilf_find_tiles(struct noisechiselparams *p,
                              struct qthreshparams *qprm)
{
  char *msg;
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;

  /* Find the good tiles. */
  while(1)
    {
      /* Find the good tiles. */
      gal_threads_spin_off(qthresh_on_tile, qprm, tl->tottiles,
                           cp->numthreads, cp->minmapsize,
                           cp->quietmmap);

      /* If we have a non-zero number of good tiles, we can break out of
         this loop. Otherwise, we need to decrease the quantile threshold
         and re-calculate the good tiles. */
      if( threshold_quantilf_find_tiles_good(p, qprm) ) break;
      else
        {
          if(p->meanmedqdiff>0.4)
            error(EXIT_FAILURE, 0, "no tiles could be found to estimate "
                  "the quantile threshold! Decrease '--tilesize' "
                  "to check for smaller (more numerous) regions within "
                  "the image that are not affected significantly by "
                  "signal and re-run NoiseChisel");
          else
            {
              /* Increment the quantile threshold. */
              p->meanmedqdiff *= 1.5;

              /* Let the user know. */
              if( asprintf(&msg, "WARNING: changed "
                           "'--meanmedqdiff=%.3g'.", p->meanmedqdiff)<0 )
                error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
              gal_timing_report(NULL, msg, 2);
            }
        }
    }
  free(qprm->usage);

  /* Set the flags accordingly.  */
  if( gal_blank_present(qprm->erode_th, 1) )
    {
      qprm->noerode_th->flag |= GAL_DATA_FLAG_HASBLANK;
      if(qprm->expand_th) qprm->expand_th->flag  |= GAL_DATA_FLAG_HASBLANK;
    }
  qprm->noerode_th->flag |= GAL_DATA_FLAG_BLANK_CH;
  if(qprm->expand_th) qprm->expand_th->flag  |= GAL_DATA_FLAG_BLANK_CH;

  /* Generate the check image if requested. */
  if(p->qthreshname)
    {
      qprm->erode_th->name="QTHRESH_ERODE";
      gal_tile_full_values_write(qprm->erode_th, tl,
                                 !p->ignoreblankintiles,
                                 p->qthreshname, NULL, 0);
      qprm->erode_th->name=NULL;
    }
}





static void
threshold_quantilf_find_apply_outlier(struct noisechiselparams *p,
                                      struct qthreshparams *qprm)
{
  char *msg, *reason;
  size_t numchg=0, ostat;
  struct gal_options_common_params *cp=&p->cp;
  size_t stat=1; /* Must initialize to non-zero value. */

  /* Break the loop only when 'stat==0' or if we have made two changes
     already (ideally, no more than two should be necessary). */
  while(stat && numchg<2)
    {
      /* Attempt the outlier rejection and save the status. */
      stat=gal_tileinternal_no_outlier_local(qprm->erode_th,
                                             qprm->noerode_th,
                                             qprm->expand_th,
                                             &cp->tl,
                                             p->interpmetric,
                                             p->outliernumngb,
                                             cp->numthreads,
                                             p->outliermclip,
                                             p->outliermad,
                                             p->qthreshname,
                                             "--outliernumngb");

      /* Decide on what to do. */
      switch(stat)
        {
        /* Minimum requested number of tiles were used successfully. */
        case 0: break;

        /* Number of tiles too small for a meaningful outlier rejection; so
           there is no need to continue and we'll just print a message and
           set stat to zero (so the loop doesn't continue). */
        case 1: case 2: case 3: case 4:
          if( asprintf(&msg, "WARNING: not enough tiles for outlier "
                       "rejection.")<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_timing_report(NULL, msg, 2);
          free(msg);
          stat=0;
          break;

        /* We need to change '--outliernumngb'. */
        default:

          /* For the printed message. */
          ostat=stat;

          /* When'stat' (the number of non-blank elements in the input), is
             larger than the minimum number of outliers, there is a problem
             (otherwise, 'stat' should have been zero). The problem happens
             when the maximum number is more than the total number of
             available points. As a result: the "measure" of all tiles
             becomes the same and there is no outliers by definition! In
             this case, we should decrease the maximum number to about half
             of the total numer so we get a good distribution. We will also
             decrease the minimum, so the minimum is smaller than the
             maximum (has no effect otherwise!). */
          if(stat>=p->outliernumngb[0])
            {
              reason="max was too large";
              p->outliernumngb[1]=stat/2;
              p->outliernumngb[0]=stat/2-1;

              /* Increment the number of changes we are making (we don't
                 want too many of these changes). */
              ++numchg;
            }

          /* When 'stat' is less than the minimum number of outliers, we
             just have to set the minimum number to half the total number
             of available tiles (similar to above: to have sufficient
             diversity). */
          else
            {
              reason="min was too large";
              p->outliernumngb[0]-=1;
            }

          /* In case a check image was made, delete the previously added
             HDU. */
          if(p->qthreshname)
            gal_fits_hdu_delete(p->qthreshname,
                                GAL_TILEINTERNAL_OUTLIER_LOCAL_HDUNAME,
                                NULL);

          /* Print a warning message. */
          if( asprintf(&msg, "WARNING: changed "
                       "'--outliernumngb=%zu,%zu' (all: %zu, %s)",
                       p->outliernumngb[0], p->outliernumngb[1],
                       ostat, reason)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_timing_report(NULL, msg, 2);
          free(msg);
        }
    }
}





static void
threshold_quantilf_interpolate(struct noisechiselparams *p,
                               struct qthreshparams *qprm)
{
  /* Use the no-outlier grid as a basis for later estimating the sky. To
     see this array on the image, use 'gal_tile_full_values_write'. */
  p->noskytiles=gal_blank_flag(qprm->erode_th);
  /* For a check:
  gal_tile_full_values_write(p->noskytiles, &cp->tl, 1,
                             "noskytiles.fits", NULL, NULL);
  */

  /* Interpolate and smooth the derived values. */
  threshold_interp_smooth(p, &qprm->erode_th, &qprm->noerode_th,
                          qprm->expand_th ? &qprm->expand_th : NULL,
                          p->qthreshname);
}




void
threshold_quantile_find_apply(struct noisechiselparams *p)
{
  char *msg;
  struct timeval t1;
  struct qthreshparams qprm;
  struct gal_options_common_params *cp=&p->cp;
  struct gal_tile_two_layer_params *tl=&cp->tl;

  /* Get the starting time if necessary. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);

  /* Allocate necessary datasets. */
  threshold_quantilf_prepare(p, &qprm);

  /* Find the good tiles. */
  threshold_quantilf_find_tiles(p, &qprm);

  /* Remove outliers. */
  if(p->outliernumngb)
    threshold_quantilf_find_apply_outlier(p, &qprm);

  /* Interpolate to fill the grid. */
  threshold_quantilf_interpolate(p, &qprm);

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
    ui_abort_after_check(p, p->qthreshname, NULL,
                         "quantile threshold check");
}
