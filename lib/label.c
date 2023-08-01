/*********************************************************************
label -- Work on labeled (integer valued) datasets.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2024 Free Software Foundation, Inc.

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

#include <gnuastro/list.h>
#include <gnuastro/qsort.h>
#include <gnuastro/label.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/checkset.h>







/****************************************************************
 *****************         Internal          ********************
 ****************************************************************/
static void
label_check_type(gal_data_t *in, uint8_t needed_type, char *variable,
                 const char *func)
{
  if(in->type!=needed_type)
    error(EXIT_FAILURE, 0, "%s: the '%s' dataset has '%s' type, but it "
          "must have a '%s' type.\n\n"
          "You can use 'gal_data_copy_to_new_type' or "
          "'gal_data_copy_to_new_type_free' to convert your input dataset "
          "to this type before calling this function", func, variable,
          gal_type_name(in->type, 1), gal_type_name(needed_type, 1));
}




















/****************************************************************
 *****************          Indexs           ********************
 ****************************************************************/
/* Put the indexs of each labeled region into an array of 'gal_data_t's
   (where each element is a dataset containing the respective label's
   indexs). */
gal_data_t *
gal_label_indexs(gal_data_t *labels, size_t numlabs, size_t minmapsize,
                 int quietmmap)
{
  size_t i, *areas;
  int32_t *a, *l, *lf;
  gal_data_t *max, *labindexs;

  /* Sanity check. */
  label_check_type(labels, GAL_TYPE_INT32, "labels", __func__);

  /* If the user hasn't given the number of labels, find it (maximum
     label). */
  if(numlabs==0)
    {
      max=gal_statistics_maximum(labels);
      numlabs=*((int32_t *)(max->array));
      gal_data_free(max);
    }
  labindexs=gal_data_array_calloc(numlabs+1);

  /* Find the area in each detected object (to see how much space we need
     to allocate). If blank values are present, an extra check is
     necessary, so to get faster results when there aren't any blank
     values, we'll also do a check. */
  areas=gal_pointer_allocate(GAL_TYPE_SIZE_T, numlabs+1, 1, __func__,
                             "areas");
  lf=(l=labels->array)+labels->size;
  do
    if(*l>0)  /* Only labeled regions: *l==0 (undetected), *l<0 (blank). */
      ++areas[*l];
  while(++l<lf);

  /* For a check.
  for(i=0;i<numlabs+1;++i)
    printf("detection %zu: %zu\n", i, areas[i]);
  exit(0);
  */

  /* Allocate/Initialize the dataset containing the indexs of each
     object. We don't want the labels of the non-detected regions
     (areas[0]). So we'll set that to zero. */
  for(i=1;i<numlabs+1;++i)
    gal_data_initialize(&labindexs[i], NULL, GAL_TYPE_SIZE_T, 1,
                        &areas[i], NULL, 0, minmapsize, quietmmap,
                        NULL, NULL, NULL);

  /* Put the indexs into each dataset. We will use the areas array again,
     but this time, use it as a counter. */
  memset(areas, 0, (numlabs+1)*sizeof *areas);
  lf=(a=l=labels->array)+labels->size;
  do
    if(*l>0)  /* No undetected regions (*l==0), or blank (<0). */
      ((size_t *)(labindexs[*l].array))[ areas[*l]++ ] = l-a;
  while(++l<lf);

  /* Clean up and return. */
  free(areas);
  return labindexs;
}




















/****************************************************************
 *****************   Over segmentation       ********************
 ****************************************************************/
static void
label_watershed_azimuth(double azimuth_rad, uint8_t *bytes)
{
  uint8_t thisbyte=0;
  int reset=1, bitbin=(int)( (M_PI+azimuth_rad) * 8 / (2*M_PI) );

  //printf("%s: %f, %f\n", __func__, azimuth_rad, azimuth_rad*180/M_PI);

  /* Fill in the respective byte. */
  switch(bitbin)
    {
    case 0: thisbyte |= 0x01; break;
    case 1: thisbyte |= 0x02; break;
    case 2: thisbyte |= 0x04; break;
    case 3: thisbyte |= 0x08; break;
    case 4: thisbyte |= 0x10; break;
    case 5: thisbyte |= 0x20; break;
    case 6: thisbyte |= 0x40; break;
    case 7:
    case 8:/* when azimuth_rad==M_PI (impossible to be larger) */
      thisbyte |= 0x80; break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us "
            "at '%s' to fix the problem. The azimuthal bin "
            "value of %d is not recognized", __func__,
            PACKAGE_BUGREPORT, bitbin);
    }

  /* If this bit is touching the previous bit (equal, or just one above or
     below), then we should increment the counter (not reset the
     counter). The goal of the checks above is therefore to see if the
     counter should be reset or not. We will first check the equality with
     the previous bit (which is easier to check), then the before/after
     cases (which involve a 'switch' due to the first/last bit). Recall
     that 'reset' has been initialized above to '1' at the start of this
     function. */
  if(thisbyte==bytes[0]) reset=0;
  else
    switch(thisbyte)
      {
      case 0x01:              /* 00000001 */
        if( (bytes[0] & 0x02) || (bytes[0] & 0x80) ) {reset=0;} break;
      case 0x80:              /* 10000000 */
        if( (bytes[0] & 0x01) || (bytes[0] & 0x40) ) {reset=0;} break;
      default:                /* Other (a middle bit). */
        if(thisbyte & (bytes[0]>>1 | bytes[0]<<1) )  {reset=0;} break;
      }

  /* If we need to reset the counter ('thisbyte' is more than one bit away
     from the previous byte), set the counter to zero. */
  bytes[1] = reset ? 0 : bytes[1]+1;

  /* Set 'thisbyte' to be the "previous" of next pixel in this clump. */
  bytes[0]=thisbyte;
}





/* Over-segment the region specified by its indexs into peaks and their
   respective regions (clumps). This is very similar to the immersion
   method of Vincent & Soille(1991), but here, we will not separate the
   image into layers, instead, we will work based on the ordered flux
   values. If a certain pixel (at a certain level) has no neighbors, it is
   a local maximum and will be assigned a new label. If it has a labeled
   neighbor, it will take that label and if there is more than one
   neighboring labeled region that pixel will be a 'river' pixel.

   DON'T FORGET: SET THE FLAGS FOR CONV EQUAL TO INPUT IN SEGMENT.

*/
size_t
gal_label_watershed(gal_data_t *values, gal_data_t *indexs,
                    gal_data_t *labels, size_t *topinds,
                    int min0_max1, uint8_t num_similar_azimuth)
{
  int hasblank;
  double azimuth_rad;
  gal_data_t *azcheck=NULL;
  uint8_t *azcheckarr=NULL;
  float *azarr=NULL, *arr=values->array;
  size_t azchecksize=0, ndim=values->ndim;
  gal_list_sizet_t *Q=NULL, *cleanup=NULL;
  size_t *a, *af, ind, *dsize=values->dsize;
  int32_t n1, nlab, rlab, curlab=1, *labs=labels->array;
  size_t *dinc=gal_dimension_increment(values->ndim, dsize);
  float xu, xd, yu, yd; /* xu: X-axis-up, xd: X-axis-down. */

  /* Sanity checks. */
  label_check_type(values, GAL_TYPE_FLOAT32, "values", __func__);
  label_check_type(indexs, GAL_TYPE_SIZE_T,  "indexs", __func__);
  label_check_type(labels, GAL_TYPE_INT32,   "labels", __func__);
  if( gal_dimension_is_different(values, labels) )
    error(EXIT_FAILURE, 0, "%s: the 'values' and 'labels' arguments must "
          "have the same size", __func__);
  if(indexs->ndim!=1)
    error(EXIT_FAILURE, 0, "%s: 'indexs' has to be a 1D array, but it is "
          "%zuD", __func__, indexs->ndim);


  /* See if there are blank values in the input dataset. */
  hasblank=gal_blank_present(values, 0);


  /*********************************************
   For checks and debugging:*
  int32_t labtocheck=117;
  size_t checkdsize[2]={20,20};    // Width of crop to check.
  size_t checkstart[2]={326,167};  // Starting coordinate (in C).

  gal_data_t *crop;
  size_t extcount=1;
  int32_t *cr, *crf;
  char *azfilename="azimuth.fits";
  char *filename="clumpbuild.fits";
  size_t checkstartind=gal_dimension_coord_to_index(2, dsize, checkstart);
  gal_data_t *tile=gal_data_alloc(gal_pointer_increment(arr, checkstartind,
                                                        values->type),
                                  GAL_TYPE_INVALID, 2, checkdsize, NULL,
                                  0, values->minmapsize,
                                  values->quietmmap, NULL, NULL, NULL);
  gal_data_t *azim=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, values->ndim,
                                 values->dsize, values->wcs, 1,
                                 values->minmapsize, values->quietmmap,
                                 NULL, NULL, NULL);
  azarr=azim->array;
  tile->block=values;
  gal_checkset_writable_remove(filename, NULL, 0, 0);
  gal_checkset_writable_remove(azfilename, NULL, 0, 0);
  crop=gal_data_copy(tile);
  gal_fits_img_write(crop, filename, NULL, NULL);
  gal_data_free(crop);
  printf("blank: %d\nriver: %d\ntmpcheck: %d\ninit: %d\n",
         (int32_t)GAL_BLANK_INT32, (int32_t)GAL_LABEL_RIVER,
         (int32_t)GAL_LABEL_TMPCHECK, (int32_t)GAL_LABEL_INIT);
  tile->array=gal_tile_block_relative_to_other(tile, labels);
  tile->block=labels;
  **********************************************/


  /* If the size of the indexs is zero, then this function is pointless. */
  if(indexs->size==0) return 0;


  /* If the indexs aren't already sorted (by the value they correspond to),
     sort them given indexs based on their flux ('gal_qsort_index_arr' is
     defined as static in 'gnuastro/qsort.h'). */
  if( !( (indexs->flag & GAL_DATA_FLAG_SORT_CH)
        && ( indexs->flag
             & (GAL_DATA_FLAG_SORTED_I
                | GAL_DATA_FLAG_SORTED_D) ) ) )
    {
      gal_qsort_index_single=values->array;
      qsort(indexs->array, indexs->size, sizeof(size_t),
            ( min0_max1
              ? gal_qsort_index_single_float32_d
              : gal_qsort_index_single_float32_i) );
    }


  /* Allocate an array to keep the azimuthal coverage. This array will
     contain 2 bytes for every clump.

       - The first byte (8 bits) for each clump will flag the azimuthal
         slice of the last pixel of this clump. The azimuthal slice is
         defined in the units of 1/8th of the azimuthal angle circle (45
         degree slices), so every bit in a byte can act as a flag for one
         azimuthal slice.

       - The second byte (8 bits) for each clump will store the number of
         previous pixels in similar azimuthal slice.

     These will allow us to stop the growth of a clump that is growing due
     to the wing of a larger clump: when the added pixels to a clump are
     fixed to a certain azimuthal range, it means that they belong to
     another clump.

     One technicality: we do not know the number of clumps a-priori! In the
     worst case that every pixel of the image is a local maximum surrounded
     by 8 minima, the number of clumps will be 1/9th the total number of
     pixels. So we'll use that. We also need to add one to this because we
     are counting clumps from 1.*/
  if(num_similar_azimuth)
    {
      azchecksize=(indexs->size/9+1)*2;
      azcheck=gal_data_alloc(NULL, GAL_TYPE_UINT8, 1, &azchecksize, NULL,
                             1, values->minmapsize, values->quietmmap,
                             NULL, NULL, NULL);
      azcheckarr=azcheck->array;
    }


  /* Initialize the region we want to over-segment. */
  af=(a=indexs->array)+indexs->size;
  do labs[*a]=GAL_LABEL_INIT; while(++a<af);


  /* Go over all the given indexs and pull out the clumps. */
  af=(a=indexs->array)+indexs->size;
  do

    /* When regions of a constant flux or masked regions exist, some later
       indexs (although they have same flux) will be filled before hand. If
       they are done, there is no need to do them again. */
    if(labs[*a]==GAL_LABEL_INIT)
      {
        /* Make sure that the total number of clumps does not exceed the
           available space for their azimuthal coverage. */
        if(azcheckarr && curlab>azchecksize/2)
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
                "explore and fix the problem. The number of clumps has "
                "exceed the assumed maximum", __func__, PACKAGE_BUGREPORT);


        /* It might happen where one or multiple regions of the pixels
           under study have the same flux. So two equal valued pixels of
           two separate (but equal flux) regions will fall immediately
           after each other in the sorted list of indexs and we have to
           account for this.

           Therefore, if we see that the next pixel in the index list has
           the same flux as this one, it does not guarantee that it should
           be given the same label. Similar to the breadth first search
           algorithm for finding connected components, we will search all
           the neighbours and the neighbours of those neighbours that have
           the same flux of this pixel to see if they touch any label or
           not and to finally give them all the same label. */
        if( (a+1)<af && arr[*a]==arr[*(a+1)] )
          {
            /* Label of first neighbor found. */
            n1=0;

            /* A small sanity check. */
            if(Q!=NULL || cleanup!=NULL)
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s "
                    "so we can fix this problem. 'Q' and 'cleanup' should "
                    "be NULL but while checking the equal flux regions "
                    "they aren't", __func__, PACKAGE_BUGREPORT);

            /* Add this pixel to a queue. */
            gal_list_sizet_add(&Q, *a);
            gal_list_sizet_add(&cleanup, *a);
            labs[*a] = GAL_LABEL_TMPCHECK;

            /* Find all the pixels that have the same flux and are
               connected. */
            while(Q!=NULL)
              {
                /* Pop an element from the queue. */
                ind=gal_list_sizet_pop(&Q);

                /* Look at the neighbors and see if we already have a
                   label. */
                GAL_DIMENSION_NEIGHBOR_OP(ind, ndim, dsize, ndim, dinc,
                   {
                     /* If it is already decided to be a river, then stop
                        looking at the neighbors. */
                     if(n1!=GAL_LABEL_RIVER)
                       {
                         /* For easy reading. */
                         nlab=labs[ nind ];

                         /* This neighbor's label isn't zero. */
                         if(nlab)
                           {
                             /* If this neighbor has not been labeled yet
                                and has an equal flux, add it to the queue
                                to expand the studied region.*/
                             if( nlab==GAL_LABEL_INIT&&arr[nind]==arr[*a] )
                               {
                                 labs[nind]=GAL_LABEL_TMPCHECK;
                                 gal_list_sizet_add(&Q, nind);
                                 gal_list_sizet_add(&cleanup, nind);
                               }
                             else
                               n1=( nlab>0

                                    /* If this neighbor has a positive
                                       nlab, it belongs to another object,
                                       so if 'n1' has not been set for the
                                       whole region (n1==0), put 'nlab'
                                       into 'n1'. If 'n1' has been set and
                                       is different from 'nlab' then this
                                       whole equal flux region should be a
                                       wide river because it is connecting
                                       two connected regions. */
                                    ? ( n1
                                        ? (n1==nlab ? n1 : GAL_LABEL_RIVER)
                                        : nlab )

                                    /* If the data has blank pixels, see if
                                       the neighbor is blank. If so, set
                                       the label to a river. Checking for
                                       the presence of blank values in the
                                       dataset can be done outside this
                                       loop (or even outside this function
                                       if flags are set). So to help the
                                       compiler optimize the program, we'll
                                       first use the pre-checked value. */
                                    : ( ( hasblank && isnan(arr[nind]) )
                                        ? GAL_LABEL_RIVER
                                        : n1 ) );
                           }

                         /* If this neigbour has a label of zero, then we
                            are on the edge of the indexed region (the
                            neighbor is not in the initial list of pixels
                            to segment). When over-segmenting the noise and
                            the detections, 'label' is zero for the parts
                            of the image that we are not interested in
                            here. */
                         else labs[*a]=GAL_LABEL_RIVER;
                       }
                   } );
              }

            /* Set the label that is to be given to this equal flux
               region. If 'n1' was set to any value, then that label should
               be used for the whole region. Otherwise, this is a new
               label, see the case for a non-flat region. */
            if(n1) rlab = n1;
            else
              {
                rlab = curlab++;
                if( topinds )           /* This is a local maximum of   */
                  topinds[rlab]=*a;     /* this region, save its index. */
              }

            /* Give the same label to the whole connected equal flux
               region, except those that might have been on the side of
               the image and were a river pixel. */
            while(cleanup!=NULL)
              {
                ind=gal_list_sizet_pop(&cleanup);
                /* If it was on the sides of the image, it has been
                   changed to a river pixel. */
                if( labs[ ind ]==GAL_LABEL_TMPCHECK ) labs[ ind ]=rlab;
              }
          }

        /* The flux of this pixel is not the same as the next sorted
           flux, so simply find the label for this object. */
        else
          {
            /* 'n1' is the label of the first labeled neighbor found, so
               we'll initialize it to zero. */
            n1=0;

            /* Go over all the fully connected neighbors of this pixel and
               see if all the neighbors (with maximum connectivity: the
               number of dimensions) that have a non-macro value belong to
               one label or not. If the pixel is neighboured by more than
               one label, set it as a river pixel. Also if it is touching a
               zero valued pixel (which does not belong to this object),
               set it as a river pixel. */
            xu=xd=yu=yd=NAN;    /* Initialize all the gradient values. */
            GAL_DIMENSION_NEIGHBOR_OP(*a, ndim, dsize, ndim, dinc,
               {
                 /* For the calculation of azimuthal angle. */
                 if(azcheckarr)
                   {
                     if     (nind==*a+1       ) xu=arr[nind];
                     else if(nind==*a-1       ) xd=arr[nind];
                     else if(nind==*a+dsize[1]) yu=arr[nind];
                     else if(nind==*a-dsize[1]) yd=arr[nind];
                   }

                 /* When 'n1' has already been set as a river, there is no
                    point in looking at the other neighbors. */
                 if(n1!=GAL_LABEL_RIVER)
                   {
                     /* For easy reading. */
                     nlab=labs[ nind ];

                     /* If this neighbor is out of the bounds of this
                        detection (sky or detections) it will be zero and
                        we have a river. Otherwise, set the first neighbor
                        accordingly. Note that we also want the zero valued
                        neighbors (detections if working on sky, and sky if
                        working on detection): we want rivers between the
                        two domains. */
                     n1 = ( nlab

                            /* nlab is usable (0: outside of the
                               detection. negative: non-usable). */
                            ? ( nlab>0

                                /* Neighbor is a previously labeled
                                   clump. */
                                ? ( n1

                                    /* The previously checked neighbor
                                       ('n1') had found a usable (positive)
                                       label. */
                                    ? ( nlab==n1

                                        /* The previously checked neighbor
                                           had the same label as this
                                           neighbor. */
                                        ? n1

                                        /* The previously checked neighbor
                                           had a different label to this
                                           neighbor */
                                        : GAL_LABEL_RIVER )

                                    /* No previously checked neighbor had a
                                       label ('n1==0'). */
                                    : ( azcheckarr
                                        && azcheckarr[nlab*2]==UINT8_MAX

                                        /* It has been previously
                                           discovered that new pixels are
                                           not adding azimuthal
                                           coverage. We need to stop the
                                           growth of this clump by defining
                                           this pixel as a river. */
                                        ? GAL_LABEL_RIVER

                                        /* New pixels are adding azimuthal
                                           coverage, so this pixel can be
                                           the same label as the found
                                           label. */
                                        : nlab) )

                                /* If the data has blank pixels, see if the
                                   neighbor is blank. If so, set the label
                                   to a river. Checking for the presence of
                                   blank values in the dataset can be done
                                   outside this loop (or even outside this
                                   function if flags are set). So to help
                                   the compiler optimize the program, we'll
                                   first use the pre-checked value. */
                                : ( ( hasblank && isnan(arr[nind]) )
                                    ? GAL_LABEL_RIVER
                                    : n1 ) )

                            /* 'nlab==0' (the neighbor lies in the other
                               domain (sky or detections). To avoid the
                               different domains touching, this pixel
                               should be a river. */
                            : GAL_LABEL_RIVER );
                   }
               });

            /* Either assign a new label to this pixel, or give it the one
               of its neighbors. If n1 equals zero, then this is a new
               peak, and a new label should be created.  But if n1!=0, it
               is either a river pixel (has more than one labeled neighbor
               and has been set to 'GAL_LABEL_RIVER' before) or all its
               neighbors have the same label. In both such cases, rlab
               should be set to n1. */
            if(n1) rlab = n1;
            else
              {
                rlab = curlab++;
                if( topinds )
                  topinds[ rlab ]=*a;
              }

            /* Calculate the azimuthal angle if we are on a relevant pixel
               (with a positive rlab, not a river for example). We are
               calculating the azimuth here to avoid copying all the
               up-down values into 'label_watershed_azimuth' and keep it
               simple. */
            if(rlab>0 && azcheckarr)
              {
                /* Find the azimuthal angle and use bit-wise flags to count
                   how similar it is to previous azimuthal angles. The
                   'azimuth_rad' is calculated here (and not within
                   'label_watershed_azimuth') to avoid having to copy all
                   the 'yu', 'yd', 'xu' and 'xd' variables into there.*/
                azimuth_rad=atan2( (double)yu-(double)yd,
                                   (double)xu-(double)xd );
                if( !isnan(azimuth_rad) )
                  {
                    /* Check the azimuthal region bits. */
                    label_watershed_azimuth(azimuth_rad,
                                            &azcheckarr[rlab*2]);

                    /* If the counter for similar azimuths exceeds the
                       threshold, then set the first byte to 'UINT8_MAX'
                       (all bits active: a sign that it should no longer be
                       grown). */
                    if(azcheckarr[rlab*2+1]>num_similar_azimuth)
                      azcheckarr[rlab*2]=UINT8_MAX;

                    /* For a check (make sure it is defined above) */
                    if(azarr) azarr[*a]=azimuth_rad*180.0f/M_PI;
                  }
              }

            /* Put the label in the pixel. */
            labs[ *a ] = rlab;
          }

        /*********************************************
         For checks and debugging: *
        if(    *a / dsize[1] >= checkstart[0]
            && *a / dsize[1] <  checkstart[0] + checkdsize[0]
            && *a % dsize[1] >= checkstart[1]
            && *a % dsize[1] <  checkstart[1] + checkdsize[1]
            && labs[*a]==labtocheck // Only see values for one clump.
          )
          {
            size_t L=labs[*a];
            char *bits=gal_type_bit_string(azcheckarr + L*2, 1);
            printf("%-4zu %zu %-4zu %-4zu %-3d %-10.2f %s (%u)\n",
                   ++extcount, *a, (*a%dsize[1])-checkstart[1],
                   (*a/dsize[1])-checkstart[0], labs[*a],
                   azimuth_rad*180.0f/M_PI, bits, azcheckarr[L*2+1]);
            crop=gal_data_copy(tile);
            crf=(cr=crop->array)+crop->size;
            do if(*cr==GAL_LABEL_RIVER) *cr=0; while(++cr<crf);
            gal_fits_img_write(crop, filename, NULL, NULL);
            gal_data_free(crop);
            free(bits);
          }
         **********************************************/
      }
  while(++a<af);

  /*********************************************
   For checks and debugging: *
  gal_fits_img_write(azim, azfilename, NULL, NULL);
  tile->array=NULL;
  gal_data_free(tile);
  gal_data_free(azim);
  printf("Total number of clumps: %u\n", curlab-1);
  exit(0);
  **********************************************/

  /* Clean up and return the total number of clumps*/
  gal_data_free(azcheck);
  free(dinc);
  return curlab-1;
}




















/**********************************************************************/
/*************             Clump significance             *************/
/**********************************************************************/
static int
label_clump_significance_sanity(gal_data_t *values,
                                gal_data_t *index_values, gal_data_t *std,
                                gal_data_t *label, gal_data_t *indexs,
                                struct gal_tile_two_layer_params *tl,
                                gal_data_t *sig, const char *func)
{
  size_t *a, *af;
  float *iv, first=NAN, second=NAN;

  /* Type of values. */
  if( values->type!=GAL_TYPE_FLOAT32 )
    error(EXIT_FAILURE, 0, "%s: the 'values' dataset must have a 'float' "
          "type, but it has a '%s' type", func,
          gal_type_name(values->type, 1));
  if( index_values->type!=GAL_TYPE_FLOAT32 )
    error(EXIT_FAILURE, 0, "%s: the 'index_values' dataset must have a "
          "'float' type, but it has a '%s' type", func,
          gal_type_name(values->type, 1));

  /* Type of standard deviation. */
  if( std->type!=GAL_TYPE_FLOAT32 )
    error(EXIT_FAILURE, 0, "%s: the standard deviation dataset must have a "
          "'float' ('float32') type, but it has a '%s' type", func,
          gal_type_name(std->type, 1));

  /* Type of labels image. */
  if( label->type!=GAL_TYPE_INT32 )
    error(EXIT_FAILURE, 0, "%s: the labels dataset must have an 'int32' "
          "type, but it has a '%s' type", func,
          gal_type_name(label->type, 1));

  /* Dimentionality of the values dataset. */
  if( values->ndim>3 )
    error(EXIT_FAILURE, 0, "%s: currently only supports 1, 2 or 3 "
          "dimensional datasets, but a %zu-dimensional dataset is given",
          func, values->ndim);

  /* Type of indexs image. */
  if( indexs->type!=GAL_TYPE_SIZE_T )
    error(EXIT_FAILURE, 0, "%s: the indexs dataset must have a 'size_t' "
          "type, but it has a '%s' type", func,
          gal_type_name(label->type, 1));

  /* Dimensionality of indexs (must be 1D). */
  if( indexs->ndim!=1 )
    error(EXIT_FAILURE, 0, "%s: the indexs dataset must be a 1D dataset, "
          "but it has %zu dimensions", func, indexs->ndim);

  /* Similar sizes between values, index_values, and label images.. */
  if( gal_dimension_is_different(values, index_values) )
    error(EXIT_FAILURE, 0, "%s: the values and label arrays don't have the "
          "same size.", func);
  if( gal_dimension_is_different(values, label) )
    error(EXIT_FAILURE, 0, "%s: the values and label arrays don't have the "
          "same size.", func);

  /* Size of the standard deviation. */
  if( !( std->size==1
         || std->size==values->size
         || (tl && std->size==tl->tottiles) ) )
    error(EXIT_FAILURE, 0, "%s: the standard deviation dataset has %zu "
          "elements. But it can only have one of these sizes: 1) a "
          "single value (used for the whole dataset), 2) The size of "
          "the values dataset (%zu elements, one value for each "
          "element), 3) The size of the number of tiles in the input "
          "tessellation (when a tessellation is given)",
          func, std->size, values->size);

  /* If the 'array' and 'dsize' elements of 'sig' have already been set. */
  if(sig->array)
    error(EXIT_FAILURE, 0, "%s: the dataset that will contain the "
          "significance values must have NULL pointers for its 'array' "
          "and 'dsize' pointers (they will be allocated here)", func);

  /* See if the clumps are to be built starting from local maxima or local
     minima. */
  iv=index_values->array;
  af=(a=indexs->array)+indexs->size;
  do
    /* A label may have NAN values. */
    if( !isnan(iv[*a]) )
      {
        if( isnan(first) ) /* 'first' is not set yet. */
          first=iv[*a];
        else
          {
            if( isnan(second) ) /* 'second' is not set yet. */
              {
                /* Note that the elements may have equal values, so for
                   'second', we want the first non-blank AND different
                   value. */
                if( iv[*a]!=first )
                  second=iv[*a];
              }
            else
              break;
          }
      }
  while(++a<af);

  /* Note that if all the values are blank or there is only one value
     covered by all the indexs, then both (or one) of 'first' or 'second'
     will be NAN. In either case, the significance measure is not going to
     be meaningful if we assume the clumps start from the maxima or
     minima. So we won't check if they are NaN or not. */
  return first>second ? 1 : 0;
}





/* This function is like a minimal MakeCatalog parsing function for clumps
   ('parse_clumps' in MakeCatalogs's 'parse.c'). To help in debugging and
   reading, the column names are defined with the same macros. There is no
   chance of collision because MakeCatalog is program, not a library. */
enum infocols
  {
    CCOL_NUM,          /* Total area within clump.                      */
    CCOL_SUM,          /* Sum of values within clump.                   */
    CCOL_STD,          /* Sum of values within clump.                   */
    CCOL_SUM_VAR,      /* Total area within clump.                      */
    CCOL_RIV_NUM,      /* Total area within rivers around clump.        */
    CCOL_RIV_SUM,      /* Sum of values within rivers around each clump.*/
    CCOL_RIV_SUM_VAR,  /* Sum of values within rivers around each clump.*/

    CCOL_NCOLS,        /* Total number of columns in the 'info' table.  */
  };
static void
label_clump_significance_raw(gal_data_t *values_d, gal_data_t *std_d,
                             gal_data_t *label_d, gal_data_t *indexs,
                             struct gal_tile_two_layer_params *tl,
                             double *info, int variance, int sky0_det1)
{
  size_t ndim=values_d->ndim, *dsize=values_d->dsize;

  double skystd, skyvar, *row;
  size_t i, *a, *af, ii, coord[3];
  size_t nngb=gal_dimension_num_neighbors(ndim);
  int32_t nlab, *ngblabs, *label=label_d->array;
  float *values=values_d->array, *std=std_d->array;
  size_t *dinc=gal_dimension_increment(ndim, dsize);

  /* Allocate the array to keep the neighbor labels of river pixels. */
  ngblabs=gal_pointer_allocate(GAL_TYPE_INT32, nngb, 0, __func__,
                               "ngblabs");

  /* Go over all the pixels in this region. */
  af=(a=indexs->array)+indexs->size;
  do
    if( !isnan(values[ *a ]) )
      {
        /* Extract the sky variance (which we need for both clump or river
           pixels). First, we'll get the sky standard deviation, then take
           it to the power of two to have variance (unless 'variance' is 1,
           showing that the input std dataset was actually variance). */
        if( std_d->size==1 || std_d->size==values_d->size)
          skystd = std_d->size==1 ? std[0] : std[*a];
        else
          { /* The standard deviation is in a tessellation. */
            gal_dimension_index_to_coord(*a, ndim, dsize, coord);
            skystd = std[ gal_tile_full_id_from_coord(tl, coord)];
          }
        skyvar = variance ? skystd : skystd*skystd;


        /* This pixel belongs to a clump. */
        if( label[ *a ]>0 )
          {
            /* For easy reading. */
            row = &info [ label[*a] * CCOL_NCOLS ];

            /* Add the raw measurements of this pixel into the clump
               statistics. */
            ++row[ CCOL_NUM ];
            row[ CCOL_STD ] = skystd;
            row[ CCOL_SUM ] += values[*a];
            row[ CCOL_SUM_VAR ] += skyvar + fabs( values[*a] );
          }

        /* This pixel belongs to a river (has a value of zero and isn't
           blank). */
        else
          {
            /* We are on a river pixel. So the value of this pixel has to
               be added to the rivers of any of the clumps it touches. But
               since it might touch a labeled region more than once, we use
               'ngblabs' to keep track of which label we have already added
               its value to. 'ii' is the number of different labels this
               river pixel has already been considered for. 'ngblabs' will
               keep the list labels. */
            ii=0;
            memset(ngblabs, 0, nngb*sizeof *ngblabs);

            /* Look into the 8-connected neighbors (recall that a
               connectivity of 'ndim' means all pixels touching it (even on
               one vertice). */
            GAL_DIMENSION_NEIGHBOR_OP(*a, ndim, dsize, ndim, dinc, {

                /* This neighbor's label. */
                nlab=label[ nind ];

                /* We only want those neighbors that are not rivers (>0) or
                   any of the flag values. */
                if(nlab>0)
                  {
                    /* Go over all already checked labels and make sure
                       this clump hasn't already been considered. */
                    for(i=0;i<ii;++i) if(ngblabs[i]==nlab) break;

                    /* This neighbor clump hasn't been considered yet: */
                    if(i==ii)
                      {
                        /* For easy reading. */
                        ngblabs[ii++] = nlab;
                        row = &info[ nlab * CCOL_NCOLS ];

                        /* Write the values. */
                        ++row[CCOL_RIV_NUM];
                        row[CCOL_RIV_SUM]+=values[*a];
                        row[CCOL_RIV_SUM_VAR ]+=skyvar+fabs( values[*a] );
                      }
                  }
              } );
          }
      }
  while(++a<af);

  /* Clean up. */
  free(dinc);
  free(ngblabs);
}




/* Make an S/N table for the clumps in a given region. */
void
gal_label_clump_significance(gal_data_t *values, gal_data_t *index_values,
                             gal_data_t *std, gal_data_t *label,
                             gal_data_t *indexs,
                             struct gal_tile_two_layer_params *tl,
                             size_t numclumps, size_t minarea,
                             int variance, int sky0_det1,
                             gal_data_t *sig, gal_data_t *sigind,
                             float cpscorr)
{
  double *info;
  int max1_min0;
  float *sigarr;
  int32_t *indarr=NULL;
  size_t i, ind, counter=0;
  size_t tablen=numclumps+1;
  double I, O, V, OV, N, NR, D, *row;

  /* If there were no initial clumps, then ignore this function. */
  if(numclumps==0) { sig->size=0; return; }

  /* Basic sanity checks. */
  max1_min0=label_clump_significance_sanity(values, index_values, std,
                                            label, indexs, tl, sig,
                                            __func__);

  /* Allocate the arrays to keep the final significance measure (and
     possibly the indexs). */
  sig->ndim  = 1;                        /* Depends on 'cltprm->sn' */
  sig->type  = GAL_TYPE_FLOAT32;
  if(sig->dsize==NULL)
    sig->dsize = gal_pointer_allocate(GAL_TYPE_SIZE_T, 1, 0, __func__,
                                      "sig->dsize");
  sig->array = gal_pointer_allocate(sig->type, tablen, 0, __func__,
                                    "sig->array");
  sig->size  = sig->dsize[0] = tablen;  /* MUST BE AFTER dsize. */
  info=gal_pointer_allocate(GAL_TYPE_FLOAT64, tablen*CCOL_NCOLS, 1,
                            __func__, "info");
  if( sigind )
    {
      sigind->ndim  = 1;
      sigind->type  = GAL_TYPE_INT32;
      sigind->dsize = gal_pointer_allocate(GAL_TYPE_SIZE_T, 1, 0,
                                           __func__, "sigind->dsize");
      sigind->size  = sigind->dsize[0] = tablen;/* After dsize */
      sigind->array = gal_pointer_allocate(sigind->type, tablen, 0,
                                           __func__, "sigind->array");
    }


  /* First, get the raw information necessary for making the S/N table. */
  label_clump_significance_raw(values, std, label, indexs, tl, info,
                               variance, sky0_det1);

  /* Calculate the signficance value for successful clumps. */
  sigarr=sig->array;
  if(sky0_det1) sigarr[0]=NAN;
  if(sigind) indarr=sigind->array;
  for(i=1;i<tablen;++i)
    {
      /* For readability. */
      row = &info[ i * CCOL_NCOLS ];

      /* If we have a sufficient area and any rivers were actually found
         for this clump, then do the measurement. */
      if( row[ CCOL_NUM ]>minarea && row[ CCOL_RIV_NUM ])
        {
          /* Set the index to write the values. If 'sky0_det1' is not
             called, we don't care about the IDs of the clumps anymore, so
             store the signal-to-noise ratios contiguously. Note that
             counter will always be smaller and equal to i. */
          ind = sky0_det1 ? i : counter++;

          /* For easy reading (as in 'columns_sn' of 'columns.c' in
             MakeCatalog). */
          N = row[ CCOL_NUM ];
          I = row[ CCOL_SUM ];
          V = row[ CCOL_SUM_VAR ];
          NR = row[ CCOL_RIV_NUM ];
          OV = row[ CCOL_RIV_SUM_VAR ];
          O = N*row[ CCOL_RIV_SUM ]/NR;

          /* Write the significance measure into the array. */
          if(sigind) indarr[ind]=i;
          D = max1_min0 ? I-O : O-I;
          sigarr[ind] = N>0.0f ? D/sqrt( (V+OV)*cpscorr) : NAN;

          /* For a check
          if(ind==30019)
            {
              printf("max1_min0: %d\n", max1_min0);
              printf("N: %.0f, I: %f, V: %f, NR=%.0f, OV: %f, O: %f\n",
                     N, I, V, NR, OV, O);
              printf("check: %f\n", sigarr[ind]);
              exit(0);
            }
          */
        }

      /* The clump is smaller than the minimum area or doesn't have any
         rivers. In this case, over detections, we should put a NaN when
         the S/N isn't calculated (for clumps, it is irrelevant). */
      else
        {
          if(sky0_det1)
            {
              sigarr[i]=NAN;
              if(sigind) indarr[i]=i;
            }
        }
    }

  /* If we don't want to keep the small clumps, the size of the S/N table
     has to be corrected. */
  if(sky0_det1==0)
    {
      sig->dsize[0] = sig->size = counter;
      if(sigind) sigind->dsize[0] = sigind->size = counter;
    }


  /* Clean up. */
  free(info);
}




















/**********************************************************************/
/*************               Growing labels               *************/
/**********************************************************************/
/* Grow the given labels without creating new ones. */
void
gal_label_grow_indexs(gal_data_t *labels, gal_data_t *indexs, int withrivers,
                      int connectivity)
{
  int searchngb;
  size_t *iarray=indexs->array;
  int32_t n1, nlab, *olabel=labels->array;
  size_t *s, *sf, thisround, ninds=indexs->size;
  size_t *dinc=gal_dimension_increment(labels->ndim, labels->dsize);

  /* Some basic sanity checks: */
  label_check_type(indexs, GAL_TYPE_SIZE_T, "indexs", __func__);
  label_check_type(labels, GAL_TYPE_INT32,  "labels", __func__);
  if(indexs->ndim!=1)
    error(EXIT_FAILURE, 0, "%s: 'indexs' has to be a 1D array, but it is "
          "%zuD", __func__, indexs->ndim);

  /* The basic idea is this: after growing, not all the blank pixels are
     necessarily filled, for example the pixels might belong to two regions
     above the growth threshold. So the pixels in between them (which are
     below the threshold will not ever be able to get a label, even if they
     are in the indexs list). Therefore, the safest way we can terminate
     the loop of growing the objects is to stop it when the number of
     pixels left to fill in this round (thisround) equals the number of
     blanks.

     To start the loop, we set 'thisround' to one more than the number of
     indexed pixels. Note that it will be corrected immediately after the
     loop has started, it is just important to pass the 'while'. */
  thisround=ninds+1;
  while( thisround > ninds )
    {
      /* 'thisround' will keep the number of pixels to be inspected in this
         round. 'ninds' will count the number of pixels left without a
         label by the end of this round. Since 'ninds' comes from the
         previous loop (or outside, for the first round) it has to be saved
         in 'thisround' to begin counting a fresh. */
      thisround=ninds;
      ninds=0;

      /* Go over all the available indexs. NOTE: while the 'indexs->array'
         pointer remains unchanged, 'indexs->size' can/will change (get
         smaller) in every loop. */
      sf = (s=indexs->array) + indexs->size;
      do
        {
          /* We'll begin by assuming the nearest neighbor of this pixel
             has no label (has a value of 0). */
          n1=0;

          /* Check the neighbors of this pixel. Note that since this
             macro has multiple loops within it, we can't use
             break. We'll use the 'searchngb' variable instead. */
          searchngb=1;
          GAL_DIMENSION_NEIGHBOR_OP(*s, labels->ndim, labels->dsize,
            connectivity, dinc,
            {
              if(searchngb)
                {
                  /* For easy reading. */
                  nlab = olabel[nind];

                  /* This neighbor's label is meaningful. */
                  if(nlab>0)                   /* This is a real label. */
                    {
                      if(n1)       /* A prev. ngb label has been found. */
                        {
                          if( n1 != nlab )    /* Different label from */
                            {    /* prevously found ngb for this pixel. */
                              n1=GAL_LABEL_RIVER;
                              searchngb=0;
                            }
                        }
                      else
                        {   /* This is the first labeld neighbor found. */
                          n1=nlab;

                          /* If we want to completely fill in the region
                             ('withrivers==0'), then there is no point in
                             looking in other neighbors, the first
                             neighbor we find, is the one we'll use. */
                          if(!withrivers) searchngb=0;
                        }
                    }
                }
            } );

          /* The loop over neighbors (above) finishes with three
             possibilities:

             n1==0                    --> No labeled neighbor was found.
             n1==GAL_LABEL_RIVER      --> Connecting two labeled regions.
             n1>0                     --> Only has one neighbouring label.

             The first one means that no neighbors were found and this
             pixel should be kept for the next loop (we'll be growing the
             objects pixel-layer by pixel-layer). In the other two cases,
             we just need to write in the value of 'n1'. */
          if(n1)
            {
              /* Set the label. */
              olabel[*s]=n1;

              /* If this pixel is a river (can only happen when
                 'withrivers' is zero), keep it in the loop, because we
                 want the 'indexs' dataset to contain all non-positive
                 (non-labeled) pixels, including rivers. */
              if(n1==GAL_LABEL_RIVER)
                iarray[ ninds++ ] = *s;
            }
          else
            iarray[ ninds++ ] = *s;

          /* Correct the size of the 'indexs' dataset. */
          indexs->size = indexs->dsize[0] = ninds;
        }
      while(++s<sf);
    }

  /* Clean up. */
  free(dinc);
}
