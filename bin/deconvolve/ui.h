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
#ifndef UI_H
#define UI_H

/* For common options groups. */
#include <gnuastro-internal/options.h>

#define TIKHONOV_NAME "tikhonov"
#define DIRECT_INVERSION_NAME "direct"
#define RICHADSON_LUCY_NAME "richardson-lucy"
#define AWMLE_NAME "awmle"

/* Available letters for short options:

   b c d e f g j l n p r s v x y z
   B C E G H J Q R W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_MULTIVALUE = 'm',
  UI_KEY_ONOFF = 'O',
  UI_KEY_KERNEL = 'k',
  UI_KEY_KHDU = 'u',
  UI_KEY_LAMBDA = 'L',
  UI_KEY_ALGORITHM = 'A',
  UI_KEY_ITERATIONS = 'i',
  UI_KEY_APLHA = 'a',
  UI_KEY_WAVES = 'w',
  UI_KEY_TOLERANCE = 't',
  UI_KEY_SIGMA = 's'
  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};

void ui_read_check_inputs_setup (int argc, char *argv[],
                                 struct deconvolve_params *p);

void ui_free_report (struct deconvolve_params *p, struct timeval *t1);

#endif
