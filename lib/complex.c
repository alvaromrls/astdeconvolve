/*********************************************************************
Functions to perform operations with complex arrays.
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

#include <error.h>

#include <gnuastro/complex.h>
#include <gnuastro/data.h>
#include <gnuastro/pointer.h>

/**
 * @brief Function to convert a complex number array into a real number array.
 *
 * @param complexarray The complex array to be converted to real.
 * @param size The number of elements in the array.
 * @param action How to convert the complex numbers.
 * @return double* A pointer to the new data.
 */
double *
gal_complex_to_real (gsl_complex_packed_array complexarray, size_t size,
                     complex_to_real action)
{
  double *out;

  /* Allocate the space for the real array. */
  out = gal_pointer_allocate (GAL_TYPE_FLOAT64, size, 1, __func__, "out");

  /* Fill the real array with the derived value from the complex array. */
  switch (action)
    {
    case COMPLEX_TO_REAL_SPEC:
      for (size_t index = 0; index < size; index++)
        {
          size_t complexIndex = 2 * index;
          double real = complexarray[complexIndex];
          double im = complexarray[complexIndex + 1];
          out[index] = sqrt ((real * real) + (im * im));
        }
      break;
    case COMPLEX_TO_REAL_PHASE:
      for (size_t index = 0; index < size; index++)
        {
          size_t complexIndex = 2 * index;
          double real = complexarray[complexIndex];
          double im = complexarray[complexIndex + 1];
          out[index] = atan2 (real, im);
        }
      break;
    case COMPLEX_TO_REAL_REAL:
      for (size_t index = 0; index < size; index++)
        {
          size_t complexIndex = 2 * index;
          out[index] = complexarray[complexIndex];
        }
      break;
    default:
      error (EXIT_FAILURE, 0,
             "%s: a bug! Please contact us at %s so we can "
             "correct it. The 'action' code %d is not recognized",
             __func__, PACKAGE_BUGREPORT, action);
    }
  return out;
}

/**
 * @brief Converts from double to complex (only copying real part).
 *
 * @param realarray Real array to be converted to complex.
 * @param size The number of elements in the array.
 * @return gsl_complex_packed_array A pointer to the new data.
 */
gsl_complex_packed_array
gal_complex_real_to_complex (double *realarray, size_t size)
{
  gsl_complex_packed_array out;
  /* Allocate the space for the real array. */
  out = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size, 1, __func__, "out");

  for (size_t index = 0; index < size; index++)
    {
      size_t complexIndex = 2 * index;
      out[complexIndex] = realarray[index];
    }
  return out;
}

/**
 * @brief Perform an element-wise multiplication of two complex arrays.
 * @note Both arrays are assumed to have the same dimensions.
 *
 * @param first First element in the multiplication.
 * @param second Second element in the multiplication.
 * @param size The number of elements in the array.
 * @return gsl_complex_packed_array A pointer to the new data.
 */
gsl_complex_packed_array
gal_complex_multiply (gsl_complex_packed_array first,
                      gsl_complex_packed_array second, size_t size)
{
  gsl_complex_packed_array out;

  /* Allocate the space for the output array. */
  out = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size, 1, __func__,
                              "product");

  /* Iterate over the arrays: multiply the elements */
  for (size_t index = 0; index < size; index++)
    {
      gsl_complex x = first[index * 2] + I * first[index * 2 + 1];
      gsl_complex y = second[index * 2] + I * second[index * 2 + 1];
      gsl_complex result = gsl_complex_mul (x, y);
      out[index * 2] = GSL_REAL (result);
      out[index * 2 + 1] = GSL_IMAG (result);
    }
  return out;
}

/**
 * @brief Perform an element-wise substraction of two complex arrays.
 * @note Both arrays are assumed to have the same dimensions.
 *
 * @param first First element in the substraction
 * @param second second element in the substraction
 * @param size The number of elements in the array.
 * @return gsl_complex_packed_array A pointer to the new data
 */
gsl_complex_packed_array
gal_complex_substract (gsl_complex_packed_array first,
                       gsl_complex_packed_array second, size_t size)
{
  gsl_complex_packed_array out;

  /* Allocate the space for the output array. */
  out = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size, 1, __func__,
                              "product");

  /* Iterate over the arrays: substract the elements */
  for (size_t index = 0; index < size; index++)
    {
      out[index * 2] = first[index * 2] - second[index * 2];
      out[index * 2 + 1] = first[index * 2 + 1] - second[index * 2 + 1];
    }
  return out;
}

/**
 * @brief Function to multiply each element by a given double value (inplace).
 *
 * @param inout Array to be scaled.
 * @param value The scale factor.
 * @param size The number of elements in the array.
 */
void
gal_complex_scale (gsl_complex_packed_array inout, double value, size_t size)
{
  for (size_t index = 0; index < size; index++)
    {
      inout[index * 2] *= value;
      inout[index * 2 + 1] *= value;
    }
}

/**
 * @brief Perform an element-wise division of two complex arrays.
 * @note Both arrays are assumed to have the same dimensions.
 *
 * @param first Dividend array.
 * @param second Divisor array.
 * @param size The number of elements in the array.
 * @param minValue If divisor is less than this value, the result will be 0.
 * @return gsl_complex_packed_array A pointer to the new data
 */
gsl_complex_packed_array
gal_complex_divide (gsl_complex_packed_array first,
                    gsl_complex_packed_array second, size_t size,
                    double minvalue)
{
  gsl_complex_packed_array out;

  /* Allocate the space for the output array. */
  out = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size, 1, __func__,
                              "product");

  /* Iterate over the arrays: calculate division if element > minValue or set
   * to 0 otherwise*/
  for (size_t index = 0; index < size; index++)
    {
      gsl_complex x = first[index * 2] + I * first[index * 2 + 1];
      gsl_complex y = second[index * 2] + I * second[index * 2 + 1];
      if (gsl_complex_abs (y) > minvalue)
        {
          gsl_complex result = gsl_complex_div (x, y);
          out[index * 2] = GSL_REAL (result);
          out[index * 2 + 1] = GSL_IMAG (result);
        }
      else
        {
          out[index * 2] = 0;
          out[index * 2 + 1] = 0;
        }
    }
  return out;
}

/**
 * @brief Calculate the conjugate of a given complex array.
 * Conjugate(a+bi) = a-bi.
 *
 * @param input The array to conjugate.
 * @param size The number of elements in the array.
 * @return gsl_complex_packed_array A pointer to the new data
 */
gsl_complex_packed_array
gal_complex_conjugate (gsl_const_complex_packed_array input, size_t size)
{
  gsl_complex_packed_array out;
  /* Allocate the space for the output array. */
  out = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size, 1, __func__,
                              "product");

  /* Iterate over the array: copy real part and negate the imaginary */
  for (size_t index = 0; index < size; index++)
    {
      out[index * 2] = input[index * 2];
      out[index * 2 + 1] = -input[index * 2 + 1];
    }
  return out;
}

/**
 * @brief Add a scalar number to each of the elements of a given complex array.
 *
 * @param input The input array.
 * @param size The number of elements in the array.
 * @param scalar A complex number to add to each element.
 * @return gsl_complex_packed_array A pointer to the new data
 */
gsl_complex_packed_array
gal_complex_add_scalar (gsl_const_complex_packed_array input, size_t size,
                        gsl_complex scalar)
{
  gsl_complex_packed_array out;
  /* Allocate the space for the output array. */
  out = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size, 1, __func__,
                              "product");

  /* Iterate over the array: add the scalar to each element */
  for (size_t index = 0; index < size; index++)
    {
      out[index * 2] = input[index * 2] + GSL_REAL (scalar);
      out[index * 2 + 1] = input[index * 2 + 1] + GSL_IMAG (scalar);
    }
  return out;
}

/**
 * @brief Function to generate a padding in the kernel (adapt the image if
 * needed) and converts them from real images to complex ones
 *
 * @param image The base image.
 * @param kernel The kernel to be scale to the base image dimensions.
 * @param outputimage The base image in complex data.
 * @param outputkernel The kernel in complex data with base image dimensions.
 * @param xdim Base image X size if odd, Base image  X size - 1 if even.
 * @param ydim Base image Y size if odd, Base image  Y size - 1 if even.
 */
void
gal_complex_create_padding (const gal_data_t *image, const gal_data_t *kernel,
                            gsl_complex_packed_array *outputimage,
                            gsl_complex_packed_array *outputkernel,
                            size_t *xdim, size_t *ydim)
{
  size_t padsizex, padsizey, totalsize;
  size_t imgsizeX = image->dsize[0];
  size_t imgsizeY = image->dsize[1];
  size_t kernelsizeX = kernel->dsize[0];
  size_t kernelsizeY = kernel->dsize[1];
  float *imagepointer = image->array;
  float *kernelpointer = kernel->array;
  gsl_complex_packed_array pimg;    // pointer to the new image
  gsl_complex_packed_array pkernel; // pointer to the new kernel

  // It's better to have a odd dimension x odd dimension space (so kernel will
  // be centered)
  if (imgsizeX % 2 == 0)
    {
      padsizex = imgsizeX - 1;
    }
  else
    {
      padsizex = imgsizeX;
    }

  if (imgsizeY % 2 == 0)
    {
      padsizey = imgsizeY - 1;
    }
  else
    {
      padsizey = imgsizeY;
    }
  totalsize = padsizex * padsizey;

  // Allocate the image and fill it
  pimg = gal_pointer_allocate (GAL_TYPE_COMPLEX64, totalsize, 1, __func__,
                               "pimg");
  for (size_t x = 0; x < padsizex; x++)
    {
      for (size_t y = 0; y < padsizey; y++)
        {
          // base image index
          size_t index = imgsizeY * x + y;
          // output image index
          size_t indexcomplex = (padsizey * x + y) * 2;
          pimg[indexcomplex] = (double)imagepointer[index];
          pimg[indexcomplex + 1] = 0.0;
        }
    }

  // Allocate the kernel in the center and fill it
  pkernel = gal_pointer_allocate (GAL_TYPE_COMPLEX64, totalsize, 1, __func__,
                                  "pkernel");
  for (size_t x = 0; x < kernelsizeX; x++)
    {
      for (size_t y = 0; y < kernelsizeY; y++)
        {
          // index for the kernel
          size_t kindex = y + x * kernelsizeX;
          // index for the kernel position in the image
          size_t index = (x + padsizex / 2 - kernelsizeX / 2) * padsizex
                         + (y + padsizey / 2 - kernelsizeY / 2);
          pkernel[index * 2] = kernelpointer[kindex];
        }
    }
  // Assign values
  *xdim = padsizex;
  *ydim = padsizey;
  *outputimage = pimg;
  *outputkernel = pkernel;
}

/**
 * @brief Normalize an array to have a total sum of 1, this function works
 * inplace.
 *
 * @param inout The array to be normalized.
 * @param size The number of elements in the array.
 */
void
gal_complex_normalize (gsl_complex_packed_array inout, size_t size)
{
  double cumulativesume = gal_complex_cumulative_sum (inout, size);
  if (cumulativesume == 0.0)
    {
      error (EXIT_FAILURE, 0, "%s: error: the module can't be 0", __func__);
    }
  gal_complex_scale (inout, 1.0 / cumulativesume, size);
}

/**
 * @brief Function to calculate cumulative sum of an array.
 *
 * @param input The target array.
 * @param size The number of elements in the array.
 * @return double Cumulative sum.
 */
double
gal_complex_cumulative_sum (gsl_complex_packed_array input, size_t size)
{
  double cumulativesume = 0.0;
  for (size_t index = 0; index < size; index++)
    {
      cumulativesume
          += gsl_complex_abs (input[index * 2] + I * input[index * 2 + 1]);
    }
  return cumulativesume;
}

/**
 * @brief Function to raise the complex array power x
 *
 * @param input Input array to be raised.
 * @param exponent The exponent (real number).
 * @param size The number of elements in the array.
 * @return gsl_complex_packed_array A pointer to the new data
 */
gsl_complex_packed_array
gal_complex_power (gsl_complex_packed_array input, double exponent,
                   size_t size)
{
  gsl_complex_packed_array out;

  /* Allocate the space for the output array. */
  out = gal_pointer_allocate (GAL_TYPE_COMPLEX64, size, 1, __func__, "power");

  /* Iterate over the arrays: multiply the elements */
  for (size_t index = 0; index < size; index++)
    {
      gsl_complex z = input[index * 2] + I * input[index * 2 + 1];
      gsl_complex result = gsl_complex_pow (z, exponent);
      out[index * 2] = GSL_REAL (result);
      out[index * 2 + 1] = GSL_IMAG (result);
    }
  return out;
}
