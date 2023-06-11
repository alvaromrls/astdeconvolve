/*********************************************************************
Cosmological calculations.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
     Jash Shah <jash28582@gmail.com>
     Pedram Ashofteh-Ardakani <pedramardakani@pm.me>
Copyright (C) 2015-2023 Free Software Foundation, Inc.

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

#include <time.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>

#include <gnuastro/error.h>
#include <gnuastro/cosmology.h>

#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/liberror.h>

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_integration.h>




/**************************************************************/
/************             Definitions             *************/
/**************************************************************/
/* These are basic definitions that commonly go into the header files. But
   because this is a library and the user imports the header file, it is
   easier to just have them here in the main C file to avoid filling up the
   user's name-space with junk. */
struct cosmology_integrand_t
{
  double o_lambda_0;
  double o_curv_0;
  double o_matter_0;
  double o_radiation_0;
};





/* For the GSL integrations */
#define GSLIEPSABS 0
#define GSLILIMIT  1000
#define GSLIEPSREL 1e-7





/**************************************************************/
/************      Constraint Check Function      *************/
/**************************************************************/
/* Check if input parameters are witihin the required constraints.  i.e All
   density fractions should be between 0 and 1 AND their sum should not
   exceed 1. */
static int
cosmology_sanity_check(double o_lambda_0, double o_matter_0,
                       double o_radiation_0, const char *func,
                       gal_error_t **err)
{
  int stat=0;
  double sum;

  /* If there is a prior error, signal a failure (by returning 1). */
  if( gal_error_has_leave(err, func)) return 1;

  /* Check if the density fractions are between 0 and 1. */
  if(o_lambda_0 > 1 || o_lambda_0 < 0)
    stat=gal_error_add(err, GAL_ERROR_CODE_EDOM, 0, func,
                       "value to option 'olambda' must be between zero "
                       "and one (inclusive), but the given value is "
                       "'%g'. Recall that 'olambda' is the current "
                       "cosmological constant density per critical "
                       "density", o_lambda_0);

  if(o_matter_0 > 1 || o_matter_0 < 0)
    stat=gal_error_add(err, GAL_ERROR_CODE_EDOM, 0, func, "value to "
                       "option 'omatter' must be between zero and "
                       "one (inclusive), but the given value is '%g'. "
                       "Recall that 'omatter' is 'Current matter "
                       "density per critical density'", o_matter_0);

  if(o_radiation_0 > 1 || o_radiation_0 < 0)
    stat=gal_error_add(err, GAL_ERROR_CODE_EDOM, 0, func, "value to "
                       "option 'oradiation' must be between zero and "
                       "one (inclusive), but the given value is '%g'. "
                       "Recall that 'oradiation' is 'Current radiation "
                       "density per critical density", o_radiation_0);

  /* Check if the density fractions add up to 1 (within floating point
     error). */
  sum = o_lambda_0 + o_matter_0 + o_radiation_0;
  if( sum > (1+1e-8) || sum < (1-1e-8) )
    stat=gal_error_add(err, GAL_ERROR_CODE_EDOM, 0, func, "sum of "
                       "fractional densities is not 1, but %g. The "
                       "cosmological constant ('olambda'), matter "
                       "('omatter') and radiation ('oradiation') "
                       "densities are given as %g, %g, %g.", sum,
                       o_lambda_0, o_matter_0, o_radiation_0);

  /* If a new error was found, return non-zero. */
  return stat ? GAL_ERROR_CODE_ERRNOTALLOC : (*err==NULL ? 0 : 1);
}




















/**************************************************************/
/************         Integrand functions         *************/
/**************************************************************/
/* These are integrands, they won't be giving the final value. */
static double
cosmology_integrand_Ez(double z, void *params)
{
  struct cosmology_integrand_t *p=(struct cosmology_integrand_t *)params;
  return sqrt( p->o_lambda_0
               + p->o_curv_0      * (1+z) * (1+z)
               + p->o_matter_0    * (1+z) * (1+z) * (1+z)
               + p->o_radiation_0 * (1+z) * (1+z) * (1+z) * (1+z));
}





static double
cosmology_integrand_age(double z, void *params)
{
  return 1 / ( (1.0 + z) * cosmology_integrand_Ez(z,params) );
}





static double
cosmology_integrand_proper_dist(double z, void *params)
{
  return 1 / ( cosmology_integrand_Ez(z,params) );
}





static double
cosmology_integrand_comoving_volume(double z, void *params)
{
  size_t neval;
  gsl_function F;
  double result, error;

  /* Set the GSL function parameters */
  F.params=params;
  F.function=&cosmology_integrand_proper_dist;

  gsl_integration_qng(&F, 0.0, z, GSLIEPSABS, GSLIEPSREL,
                      &result, &error, &neval);

  return result * result / ( cosmology_integrand_Ez(z,params) );
}




















/**************************************************************/
/************      Basic cosmology functions      *************/
/**************************************************************/
/* Age of the universe (in Gyrs). H0 is in units of (km/sec/Mpc) and the
   fractional densities must add up to 1. */
double
gal_cosmology_age(double z, double H0, double o_lambda_0,
                  double o_matter_0, double o_radiation_0,
                  gal_error_t **err)
{
  gsl_function F;
  double result, error;
  gsl_integration_workspace *w;
  double o_curv_0 = 1.0 - ( o_lambda_0 + o_matter_0 + o_radiation_0 );
  double H0s=H0/1000/GSL_CONST_MKSA_PARSEC;  /* H0 in units of seconds. */
  struct cosmology_integrand_t p={o_lambda_0, o_curv_0, o_matter_0,
                                  o_radiation_0};

  /* Sanity checks (no errors and good inputs). */
  if( cosmology_sanity_check(o_lambda_0, o_matter_0, o_radiation_0,
                             __func__, err) ) return NAN;

  /* Set the GSL function parameters. */
  F.params=&p;
  w=gsl_integration_workspace_alloc(GSLILIMIT);
  F.function=&cosmology_integrand_age;
  gsl_integration_qagiu(&F, z, GSLIEPSABS, GSLIEPSREL, GSLILIMIT, w,
                        &result, &error);

  return result / H0s / (365*GSL_CONST_MKSA_DAY) / 1e9;
}





/* Proper distance to z (Mpc). */
double
gal_cosmology_proper_distance(double z, double H0, double o_lambda_0,
                              double o_matter_0, double o_radiation_0,
                              gal_error_t **err)
{
  size_t neval;
  gsl_function F;
  double result, error, c=GSL_CONST_MKSA_SPEED_OF_LIGHT;
  double o_curv_0 = 1.0 - ( o_lambda_0 + o_matter_0 + o_radiation_0 );
  double H0s=H0/1000/GSL_CONST_MKSA_PARSEC;  /* H0 in units of seconds. */
  struct cosmology_integrand_t p={o_lambda_0, o_curv_0, o_matter_0,
                                  o_radiation_0};

  /* Sanity checks (no errors and good inputs). */
  if( cosmology_sanity_check(o_lambda_0, o_matter_0, o_radiation_0,
                             __func__, err) ) return NAN;

  /* Set the GSL function parameters */
  F.params=&p;
  F.function=&cosmology_integrand_proper_dist;

  /* Do the integration. */
  gsl_integration_qng(&F, 0.0f, z, GSLIEPSABS, GSLIEPSREL, &result,
                      &error, &neval);

  /* Return the result. */
  return result * c / H0s / (1e6 * GSL_CONST_MKSA_PARSEC);
}





/* Comoving volume over 4pi stradian to z (Mpc^3). */
double
gal_cosmology_comoving_volume(double z, double H0, double o_lambda_0,
                              double o_matter_0, double o_radiation_0,
                              gal_error_t **err)
{
  size_t neval;
  gsl_function F;
  double result, error;
  double c=GSL_CONST_MKSA_SPEED_OF_LIGHT;
  double H0s=H0/1000/GSL_CONST_MKSA_PARSEC; /* H0 in units of seconds. */
  double cH = c / H0s / (1e6 * GSL_CONST_MKSA_PARSEC);
  double o_curv_0 = 1.0 - ( o_lambda_0 + o_matter_0 + o_radiation_0 );
  struct cosmology_integrand_t p={o_lambda_0, o_curv_0, o_matter_0,
                                  o_radiation_0};

  /* Sanity checks (no errors and good inputs). */
  if( cosmology_sanity_check(o_lambda_0, o_matter_0, o_radiation_0,
                             __func__, err) ) return NAN;

  /* Set the GSL function parameters */
  F.params=&p;
  F.function=&cosmology_integrand_comoving_volume;

  /* Do the integration. */
  gsl_integration_qng(&F, 0.0f, z, GSLIEPSABS, GSLIEPSREL,
                      &result, &error, &neval);

  /* Return the result. */
  return result * 4 * M_PI * cH*cH*cH;
}





/* Critical density at redshift z in units of g/cm^3. */
double
gal_cosmology_critical_density(double z, double H0, double o_lambda_0,
                               double o_matter_0, double o_radiation_0,
                               gal_error_t **err)
{
  double H;
  double H0s=H0/1000/GSL_CONST_MKSA_PARSEC; /* H0 in units of seconds. */
  double o_curv_0 = 1.0 - ( o_lambda_0 + o_matter_0 + o_radiation_0 );
  struct cosmology_integrand_t p={o_lambda_0, o_curv_0, o_matter_0,
                                  o_radiation_0};

  /* Sanity checks (no errors and good inputs). */
  if( cosmology_sanity_check(o_lambda_0, o_matter_0, o_radiation_0,
                             __func__, err) ) return NAN;

  /* Set the place holder, then return the result. */
  H = H0s * cosmology_integrand_Ez(z, &p);
  return 3*H*H/(8*M_PI*GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT)/1000;
}





/* Angular diameter distance to z (Mpc). */
double
gal_cosmology_angular_distance(double z, double H0, double o_lambda_0,
                               double o_matter_0, double o_radiation_0,
                               gal_error_t **err)
{
  /* Sanity checks (no errors and good inputs). */
  if( cosmology_sanity_check(o_lambda_0, o_matter_0, o_radiation_0,
                             __func__, err) ) return NAN;

  /* Do the calculation. */
  return gal_cosmology_proper_distance(z, H0, o_lambda_0, o_matter_0,
                                       o_radiation_0, err) / (1+z);
}





/* Luminosity distance to z (Mpc). */
double
gal_cosmology_luminosity_distance(double z, double H0, double o_lambda_0,
                                  double o_matter_0, double o_radiation_0,
                                  gal_error_t **err)
{
  /* Sanity checks (no errors and good inputs). */
  if( cosmology_sanity_check(o_lambda_0, o_matter_0, o_radiation_0,
                             __func__, err) ) return NAN;

  /* Do the calculation. */
  return gal_cosmology_proper_distance(z, H0, o_lambda_0, o_matter_0,
                                       o_radiation_0, err) * (1+z);
}





/* Distance modulus at z (no units). */
double
gal_cosmology_distance_modulus(double z, double H0, double o_lambda_0,
                               double o_matter_0, double o_radiation_0,
                               gal_error_t **err)
{
  /* Sanity checks (no errors and good inputs). */
  if( cosmology_sanity_check(o_lambda_0, o_matter_0, o_radiation_0,
                             __func__, err) ) return NAN;

  /* Do the calculation. */
  double ld=gal_cosmology_luminosity_distance(z, H0, o_lambda_0, o_matter_0,
                                              o_radiation_0, err);
  return 5*(log10(ld*1000000)-1);
}





/* Convert apparent to absolute magnitude. */
double
gal_cosmology_to_absolute_mag(double z, double H0, double o_lambda_0,
                              double o_matter_0, double o_radiation_0,
                              gal_error_t **err)
{
  /* Sanity checks (no errors and good inputs). */
  if( cosmology_sanity_check(o_lambda_0, o_matter_0, o_radiation_0,
                             __func__, err) ) return NAN;

  /* Do the calculation. */
  double dm=gal_cosmology_distance_modulus(z, H0, o_lambda_0, o_matter_0,
                                           o_radiation_0, err);
  return dm-2.5*log10(1.0+z);
}





/* Velocity at given redshift in units of km/s. */
double
gal_cosmology_velocity_from_z(double z)
{
  double c=GSL_CONST_MKSA_SPEED_OF_LIGHT;
  return c * ( (1+z)*(1+z) - 1 ) / ( (1+z)*(1+z) + 1 ) / 1000;
}





/* Redshift at given velocity (in units of km/s). */
double
gal_cosmology_z_from_velocity(double v)
{
  double c=GSL_CONST_MKSA_SPEED_OF_LIGHT/1000;
  return sqrt( (c+v)/(c-v) ) - 1;
}
