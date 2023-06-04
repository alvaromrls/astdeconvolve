/*********************************************************************
Arithmetic operations on data structures.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
     Siyang He <siyang.he.uw@gmail.com>
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
#ifndef __ARITHMETIC_BINARY_H__
#define __ARITHMETIC_BINARY_H__



#include <string.h>
#include <inttypes.h>

#include <gnuastro/type.h>
#include <gnuastro/pointer.h>
#include <gnuastro/arithmetic.h>





/************************************************************************/
/*************               Overflow checks            *****************/
/************************************************************************/
#define OVERFLOW_CHECK_PLUS(OP, OT_UPPER) \
  if(lval>=0 && rval>=0)                  \
      { if(rval>max-lval) overflows=1; }  \
  else if(lval<0 && rval<0)               \
      { if(lval<min-rval) overflows=1; }





#define OVERFLOW_CHECK_MINUS(OP, OT_UPPER)                       \
  switch(o->type)                                                \
    {                                                            \
    case GAL_TYPE_UINT8:                                         \
    case GAL_TYPE_UINT16:                                        \
    case GAL_TYPE_UINT32:                                        \
    case GAL_TYPE_UINT64:                                        \
      /* Unsigned subtraction underflows when LHS < RHS */       \
      if(lval<rval) overflows=1;                                 \
      break;                                                     \
    case GAL_TYPE_INT8:                                          \
    case GAL_TYPE_INT16:                                         \
    case GAL_TYPE_INT32:                                         \
    case GAL_TYPE_INT64:                                         \
    case GAL_TYPE_FLOAT32:                                       \
    case GAL_TYPE_FLOAT64:                                       \
      /* Signed subtraction can underflow or overflow */         \
      if(lval>0 && rval<0)                                       \
        { if(rval<lval-max) overflows=1; }                       \
      else if(lval<0 && rval>0)                                  \
        { if(lval<min+rval) overflows=1; }                       \
      break;                                                     \
    default:                                                     \
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",  \
            __func__, o->type);                                  \
    }





#define OVERFLOW_CHECK_MULTIPLY(OP, OT_UPPER)                       \
  switch(o->type)                                                   \
    {                                                               \
    case GAL_TYPE_INT8:                                             \
    case GAL_TYPE_INT16:                                            \
    case GAL_TYPE_INT32:                                            \
    case GAL_TYPE_INT64:                                            \
    case GAL_TYPE_FLOAT32:                                          \
    case GAL_TYPE_FLOAT64:                                          \
      /* Only perform the "<-1" checks for signed types.     */     \
      /* Reason: due to wrapping, almost all unsigned values */     \
      /* are less than -1, so there would be many false      */     \
      /* positives. */                                              \
      if(lval<-1)                                                   \
        {                                                           \
          /* When lval and rval are integers, max/lval and */       \
          /* min/lval both truncate towards zero (rather */         \
          /* than truncating down), which is exactly what */        \
          /* we need for an upper/lower limit on rval. */           \
          if(rval<max/lval || rval>min/lval) overflows=1;           \
        }                                                           \
      if(rval<-1)                                                   \
        { if(lval<max/rval || lval>min/rval) overflows=1; }         \
      /* Don't break here, because the ">1" checks below */         \
      /* also apply to signed types. */                             \
    case GAL_TYPE_UINT8:                                            \
    case GAL_TYPE_UINT16:                                           \
    case GAL_TYPE_UINT32:                                           \
    case GAL_TYPE_UINT64:                                           \
      if(lval>1)                                                    \
        {                                                           \
          if(rval>max/lval || rval<min/lval) overflows=1;           \
        }                                                           \
      if(rval>1)                                                    \
        {                                                           \
          if(lval>max/rval || lval<min/rval) overflows=1;           \
        }                                                           \
      break;                                                        \
    default:                                                        \
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",     \
            __func__, o->type);                                     \
    }





#define OVERFLOW_CHECK_DIVIDE(OP, OT_UPPER)                         \
  switch(o->type)                                                   \
    {                                                               \
    /* Integer divisions cannot overflow unless the divisor is */   \
    /* zero, and division by zero is already handled elsewhere.*/   \
    case GAL_TYPE_UINT8:                                            \
    case GAL_TYPE_UINT16:                                           \
    case GAL_TYPE_UINT32:                                           \
    case GAL_TYPE_UINT64:                                           \
    case GAL_TYPE_INT8:                                             \
    case GAL_TYPE_INT16:                                            \
    case GAL_TYPE_INT32:                                            \
    case GAL_TYPE_INT64: break;                                     \
    case GAL_TYPE_FLOAT32:                                          \
    case GAL_TYPE_FLOAT64:                                          \
      if(rval>0 && rval<1)                                          \
        { if(lval>max*rval || lval<min*rval) overflows=1; }         \
      if(rval<0 && rval>-1)                                         \
        { if(lval<max*rval || lval>min*rval) overflows=1; }         \
      break;                                                        \
    default:                                                        \
      error(EXIT_FAILURE, 0, "%s: type code %d not recognized",     \
            __func__, o->type);                                     \
    }




















/************************************************************************/
/*************             Low-level operators          *****************/
/************************************************************************/
#define BINARY_OP_OT_OVERFLOW_BLANK(OP, OT, OT_UPPER, OP_NAME) {  \
  gal_blank_write(&lb, l->type);                                  \
  gal_blank_write(&rb, r->type);                                  \
  gal_blank_write(&ob, o->type);                                  \
  do                                                              \
    {                                                             \
      if(lb==lb && rb==rb)/* Both are integers.                */ \
        {                                                         \
          if(*la!=lb  && *ra!=rb)                                 \
            {                                                     \
              lval=(OT)*la, rval=(OT)*ra;                         \
              if(checkoverflow && !overflows)                     \
                { OVERFLOW_CHECK_ ## OP_NAME( OP, OT_UPPER ); }   \
              *oa=*la OP *ra;                                     \
            }                                                     \
          else *oa=ob;                                            \
        }                                                         \
      else if(lb==lb)     /* Only left operand is an integer.  */ \
        {                                                         \
          if(*la!=lb  && *ra==*ra)                                \
            {                                                     \
              lval=(OT)*la, rval=(OT)*ra;                         \
              if(checkoverflow && !overflows)                     \
                { OVERFLOW_CHECK_ ## OP_NAME( OP, OT_UPPER ); }   \
              *oa=*la OP *ra;                                     \
            }                                                     \
          else *oa=ob;                                            \
        }                                                         \
      else                /* Only right operand is an integer. */ \
        {                                                         \
          if(*la==*la && *ra!=rb)                                 \
            {                                                     \
              lval=(OT)*la, rval=(OT)*ra;                         \
              if(checkoverflow && !overflows)                     \
                { OVERFLOW_CHECK_ ## OP_NAME( OP, OT_UPPER ); }   \
              *oa=*la OP *ra;                                     \
            }                                                     \
          else *oa=ob;                                            \
        }                                                         \
      if(l->size>1) ++la;                                         \
      if(r->size>1) ++ra;                                         \
    }                                                             \
  while(++oa<of);                                                 \
}





#define BINARY_OP_OT_OVERFLOW_NO_BLANK(OP, OT, OT_UPPER, OP_NAME) { \
  if(l->size==r->size)                                              \
    do                                                              \
      {                                                             \
        lval=(OT)*la, rval=(OT)*ra;                                 \
        if(checkoverflow && !overflows)                             \
          { OVERFLOW_CHECK_ ## OP_NAME( OP, OT_UPPER ); }           \
        *oa = *la++ OP *ra++;                                       \
      }                                                             \
    while(++oa<of);                                                 \
  else if(l->size==1)                                               \
    {                                                               \
      lval=(OT)*la;                                                 \
      do                                                            \
        {                                                           \
          rval=(OT)*ra;                                             \
          if(checkoverflow && !overflows)                           \
            { OVERFLOW_CHECK_ ## OP_NAME( OP, OT_UPPER ); }         \
          *oa = *la OP *ra++;                                       \
        }                                                           \
      while(++oa<of);                                               \
    }                                                               \
  else                                                              \
    {                                                               \
      rval=(OT)*ra;                                                 \
      do                                                            \
        {                                                           \
          lval=(OT)*la;                                             \
          if(checkoverflow && !overflows)                           \
            { OVERFLOW_CHECK_ ## OP_NAME( OP, OT_UPPER ); }         \
          *oa = *la++ OP *ra;                                       \
        }                                                           \
      while(++oa<of);                                               \
    }                                                               \
}





#define BINARY_OP_OT_RT_LT_SET_OVERFLOW_IMPOSSIBLE(OP, OT, OT_UPPER,    \
                                                   LT, RT, OP_NAME) {   \
    LT lb, *la=l->array;                                                \
    RT rb, *ra=r->array;                                                \
    OT ob, *oa=o->array, *of=oa + o->size;                              \
    if(checkblank)                                                      \
      {                                                                 \
        gal_blank_write(&lb, l->type);                                  \
        gal_blank_write(&rb, r->type);                                  \
        gal_blank_write(&ob, o->type);                                  \
        do                                                              \
          {                                                             \
            if(lb==lb && rb==rb)/* Both are integers.                */ \
              *oa = (*la!=lb  && *ra!=rb)  ? *la OP *ra : ob ;          \
            else if(lb==lb)     /* Only left operand is an integer.  */ \
              *oa = (*la!=lb  && *ra==*ra) ? *la OP *ra : ob;           \
            else                /* Only right operand is an integer. */ \
              *oa = (*la==*la && *ra!=rb)  ? *la OP *ra : ob;           \
            if(l->size>1) ++la;                                         \
            if(r->size>1) ++ra;                                         \
          }                                                             \
        while(++oa<of);                                                 \
      }                                                                 \
    else                                                                \
      {                                                                 \
        if(l->size==r->size) do *oa = *la++ OP *ra++; while(++oa<of);   \
        else if(l->size==1)  do *oa = *la   OP *ra++; while(++oa<of);   \
        else                 do *oa = *la++ OP *ra;   while(++oa<of);   \
      }                                                                 \
  }





/* Final step to be used by all operators and all types. */
#define BINARY_OP_OT_RT_LT_SET_OVERFLOW_POSSIBLE(OP, OT, OT_UPPER, LT,  \
                                                 RT, OP_NAME) {         \
    LT lb, *la=l->array;                                                \
    RT rb, *ra=r->array;                                                \
    void *type_max, *type_min;                                          \
    OT max, min, lval, rval, ob, *oa=o->array, *of=oa + o->size;        \
                                                                        \
    type_max=gal_pointer_allocate(o->type, 1, 0, __func__, "type_max"); \
    gal_type_max(o->type, type_max);                                    \
    max=*(OT*)type_max;                                                 \
    free(type_max);                                                     \
                                                                        \
    type_min=gal_pointer_allocate(o->type, 1, 0, __func__, "type_min"); \
    gal_type_min(o->type, type_min);                                    \
    min=*(OT*)type_min;                                                 \
    free(type_min);                                                     \
                                                                        \
    if(checkblank)                                                      \
      { BINARY_OP_OT_OVERFLOW_BLANK(    OP, OT, OT_UPPER, OP_NAME ); }  \
    else                                                                \
      { BINARY_OP_OT_OVERFLOW_NO_BLANK( OP, OT, OT_UPPER, OP_NAME ); }  \
}





/* Blank values aren't defined for integer operators. */
#define BINARY_INT_OP_OT_RT_LT_SET(OP, OT, LT, RT) {               \
    LT *la=l->array;                                               \
    RT *ra=r->array;                                               \
    OT *oa=o->array, *of=oa + o->size;                             \
    if(l->size==r->size) do *oa = *la++ OP *ra++; while(++oa<of);  \
    else if(l->size==1)  do *oa = *la   OP *ra++; while(++oa<of);  \
    else                 do *oa = *la++ OP *ra;   while(++oa<of);  \
  }





/* This is for operators like '&&' and '||', where the right operator is
   not necessarily read (and thus incremented). */
#define BINARY_OP_INCR_OT_RT_LT_SET(OP, OT, LT, RT) {                   \
    LT *la=l->array;                                                    \
    RT *ra=r->array;                                                    \
    OT *oa=o->array, *of=oa + o->size;                                  \
    if(l->size==r->size) do {*oa = *la++ OP *ra; ++ra;} while(++oa<of); \
    else if(l->size==1)  do {*oa = *la   OP *ra; ++ra;} while(++oa<of); \
    else                 do  *oa = *la++ OP *ra;        while(++oa<of); \
  }




















/************************************************************************/
/*************              Type specifiers             *****************/
/************************************************************************/

/* Flags for BINARY_OP_RT_LT_SET to identify the output type. */
enum arithmetic_binary_outtype_flags
{
  ARITHMETIC_BINARY_INVALID,                /* ==0 by C standard. */

  ARITHMETIC_BINARY_OUT_TYPE_LEFT,
  ARITHMETIC_BINARY_OUT_TYPE_RIGHT,
  ARITHMETIC_BINARY_OUT_TYPE_UINT8,
  ARITHMETIC_BINARY_OUT_TYPE_INCR_SEP,
};





/* For operators whose type may be any of the given inputs only for
   integers. For integer operators we have less options. */
#define BINARY_SET_OUT_INT(F, OP, LT, RT)                               \
  switch(F)                                                             \
    {                                                                   \
    case ARITHMETIC_BINARY_OUT_TYPE_LEFT:                               \
      BINARY_INT_OP_OT_RT_LT_SET(OP, LT, LT, RT);                       \
      break;                                                            \
    case ARITHMETIC_BINARY_OUT_TYPE_RIGHT:                              \
      BINARY_INT_OP_OT_RT_LT_SET(OP, RT, LT, RT);                       \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to "   \
            "address the problem. %d not recognized for 'F'",           \
            "BINARY_SET_OUT_INT", PACKAGE_BUGREPORT, F);                \
    }





/* For operators whose type may be any of the given inputs. */
#define BINARY_SET_OUT(F, OP, LT, LT_UPPER, RT, RT_UPPER,               \
                       IS_OVERFLOW_POSSIBLE, OP_NAME)                   \
  switch(F)                                                             \
    {                                                                   \
    case ARITHMETIC_BINARY_OUT_TYPE_LEFT:                               \
      BINARY_OP_OT_RT_LT_SET_ ## IS_OVERFLOW_POSSIBLE(OP, LT, LT_UPPER, \
                                                      LT, RT, OP_NAME); \
      break;                                                            \
    case ARITHMETIC_BINARY_OUT_TYPE_RIGHT:                              \
      BINARY_OP_OT_RT_LT_SET_ ## IS_OVERFLOW_POSSIBLE(OP, RT, RT_UPPER, \
                                                      LT, RT, OP_NAME); \
      break;                                                            \
    case ARITHMETIC_BINARY_OUT_TYPE_UINT8:                              \
      BINARY_OP_OT_RT_LT_SET_ ## IS_OVERFLOW_POSSIBLE(OP, uint8_t,      \
                                                      UINT8_T, LT,      \
                                                      RT, OP_NAME);     \
      break;                                                            \
    case ARITHMETIC_BINARY_OUT_TYPE_INCR_SEP:                           \
      BINARY_OP_INCR_OT_RT_LT_SET(OP, uint8_t, LT, RT);                 \
      break;                                                            \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to "   \
            "address the problem. %d not recognized for 'F'",           \
            "BINARY_SET_OUT", PACKAGE_BUGREPORT, F);                    \
    }





/* Set the right operator type only for integers (integer operators can't
   take floating point types). So floating point types must not be in the
   list of possibilities. */
#define BINARY_SET_RT_INT(F, OP, LT)                                        \
  switch(r->type)                                                           \
    {                                                                       \
    case GAL_TYPE_UINT8:   BINARY_SET_OUT_INT( F, OP, LT, uint8_t  ) break; \
    case GAL_TYPE_INT8:    BINARY_SET_OUT_INT( F, OP, LT, int8_t   ) break; \
    case GAL_TYPE_UINT16:  BINARY_SET_OUT_INT( F, OP, LT, uint16_t ) break; \
    case GAL_TYPE_INT16:   BINARY_SET_OUT_INT( F, OP, LT, int16_t  ) break; \
    case GAL_TYPE_UINT32:  BINARY_SET_OUT_INT( F, OP, LT, uint32_t ) break; \
    case GAL_TYPE_INT32:   BINARY_SET_OUT_INT( F, OP, LT, int32_t  ) break; \
    case GAL_TYPE_UINT64:  BINARY_SET_OUT_INT( F, OP, LT, uint64_t ) break; \
    case GAL_TYPE_INT64:   BINARY_SET_OUT_INT( F, OP, LT, int64_t  ) break; \
    default:                                                                \
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "       \
            "address the problem. %d is not a usable type code",            \
            "BINARY_SET_RT", PACKAGE_BUGREPORT, r->type);                   \
    }





/* Set the right operator type. */
#define BINARY_SET_RT(F, OP, LT, LT_UPPER, IS_OVERFLOW_POSSIBLE,      \
                      OP_NAME)                                        \
  switch(r->type)                                                     \
    {                                                                 \
    case GAL_TYPE_UINT8:                                              \
      BINARY_SET_OUT( F, OP, LT, LT_UPPER,  uint8_t,  UINT8_T,        \
                      IS_OVERFLOW_POSSIBLE, OP_NAME )                 \
      break;                                                          \
    case GAL_TYPE_INT8:                                               \
      BINARY_SET_OUT( F, OP, LT, LT_UPPER,  int8_t,   INT8_T,         \
                      IS_OVERFLOW_POSSIBLE, OP_NAME )                 \
      break;                                                          \
    case GAL_TYPE_UINT16:                                             \
      BINARY_SET_OUT( F, OP, LT, LT_UPPER,  uint16_t, UINT16_T,       \
                      IS_OVERFLOW_POSSIBLE, OP_NAME )                 \
      break;                                                          \
    case GAL_TYPE_INT16:                                              \
      BINARY_SET_OUT( F, OP, LT, LT_UPPER,  int16_t,  INT16_T,        \
                      IS_OVERFLOW_POSSIBLE, OP_NAME )                 \
      break;                                                          \
    case GAL_TYPE_UINT32:                                             \
      BINARY_SET_OUT( F, OP, LT, LT_UPPER,  uint32_t, UINT32_T,       \
                      IS_OVERFLOW_POSSIBLE, OP_NAME )                 \
      break;                                                          \
    case GAL_TYPE_INT32:                                              \
      BINARY_SET_OUT( F, OP, LT, LT_UPPER,  int32_t,  INT32_T,        \
                      IS_OVERFLOW_POSSIBLE, OP_NAME )                 \
      break;                                                          \
    case GAL_TYPE_UINT64:                                             \
      BINARY_SET_OUT( F, OP, LT, LT_UPPER,  uint64_t, UINT64_T,       \
                      IS_OVERFLOW_POSSIBLE, OP_NAME )                 \
      break;                                                          \
    case GAL_TYPE_INT64:                                              \
      BINARY_SET_OUT( F, OP, LT, LT_UPPER,  int64_t,  INT64_T,        \
                      IS_OVERFLOW_POSSIBLE, OP_NAME )                 \
      break;                                                          \
    case GAL_TYPE_FLOAT32:                                            \
      BINARY_SET_OUT( F, OP, LT, LT_UPPER,  float,    FLOAT,          \
                      IS_OVERFLOW_POSSIBLE, OP_NAME )                 \
      break;                                                          \
    case GAL_TYPE_FLOAT64:                                            \
      BINARY_SET_OUT( F, OP, LT, LT_UPPER,  double,   DOUBLE,         \
                      IS_OVERFLOW_POSSIBLE, OP_NAME )                 \
      break;                                                          \
    default:                                                          \
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to " \
            "address the problem. %d is not a usable type code",      \
            "BINARY_SET_RT", PACKAGE_BUGREPORT, r->type);             \
    }





/* Set the left operator type only for integers (integer operators can't
   take floating point types). So floating point types must not be in the
   list of possibilities. */
#define BINARY_SET_LT_INT(F, OP)                                        \
  switch(l->type)                                                       \
    {                                                                   \
    case GAL_TYPE_UINT8:   BINARY_SET_RT_INT( F, OP, uint8_t  ) break;  \
    case GAL_TYPE_INT8:    BINARY_SET_RT_INT( F, OP, int8_t   ) break;  \
    case GAL_TYPE_UINT16:  BINARY_SET_RT_INT( F, OP, uint16_t ) break;  \
    case GAL_TYPE_INT16:   BINARY_SET_RT_INT( F, OP, int16_t  ) break;  \
    case GAL_TYPE_UINT32:  BINARY_SET_RT_INT( F, OP, uint32_t ) break;  \
    case GAL_TYPE_INT32:   BINARY_SET_RT_INT( F, OP, int32_t  ) break;  \
    case GAL_TYPE_UINT64:  BINARY_SET_RT_INT( F, OP, uint64_t ) break;  \
    case GAL_TYPE_INT64:   BINARY_SET_RT_INT( F, OP, int64_t  ) break;  \
    default:                                                            \
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "   \
            "address the problem. %d is not a usable type code",        \
            "BINARY_SET_LT_INT", PACKAGE_BUGREPORT, l->type);           \
    }





/* Set the left operator type. */
#define BINARY_SET_LT(F, OP, IS_OVERFLOW_POSSIBLE, OP_NAME)           \
  switch(l->type)                                                     \
    {                                                                 \
    case GAL_TYPE_UINT8:                                              \
      BINARY_SET_RT( F, OP, uint8_t,  UINT8_T,  IS_OVERFLOW_POSSIBLE, \
                     OP_NAME )                                        \
      break;                                                          \
    case GAL_TYPE_INT8:                                               \
      BINARY_SET_RT( F, OP, int8_t,   INT8_T,   IS_OVERFLOW_POSSIBLE, \
                     OP_NAME )                                        \
      break;                                                          \
    case GAL_TYPE_UINT16:                                             \
      BINARY_SET_RT( F, OP, uint16_t, UINT16_T, IS_OVERFLOW_POSSIBLE, \
                     OP_NAME )                                        \
      break;                                                          \
    case GAL_TYPE_INT16:                                              \
      BINARY_SET_RT( F, OP, int16_t,  INT16_T,  IS_OVERFLOW_POSSIBLE, \
                     OP_NAME )                                        \
      break;                                                          \
    case GAL_TYPE_UINT32:                                             \
      BINARY_SET_RT( F, OP, uint32_t, UINT32_T, IS_OVERFLOW_POSSIBLE, \
                     OP_NAME )                                        \
      break;                                                          \
    case GAL_TYPE_INT32:                                              \
      BINARY_SET_RT( F, OP, int32_t,  INT32_T,  IS_OVERFLOW_POSSIBLE, \
                     OP_NAME )                                        \
      break;                                                          \
    case GAL_TYPE_UINT64:                                             \
      BINARY_SET_RT( F, OP, uint64_t, UINT64_T, IS_OVERFLOW_POSSIBLE, \
                     OP_NAME )                                        \
      break;                                                          \
    case GAL_TYPE_INT64:                                              \
      BINARY_SET_RT( F, OP, int64_t,  INT64_T,  IS_OVERFLOW_POSSIBLE, \
                     OP_NAME )                                        \
      break;                                                          \
    case GAL_TYPE_FLOAT32:                                            \
      BINARY_SET_RT( F, OP, float,    FLOAT,    IS_OVERFLOW_POSSIBLE, \
                     OP_NAME )                                        \
      break;                                                          \
    case GAL_TYPE_FLOAT64:                                            \
      BINARY_SET_RT( F, OP, double,   DOUBLE,   IS_OVERFLOW_POSSIBLE, \
                     OP_NAME )                                        \
      break;                                                          \
    default:                                                          \
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to " \
            "address the problem. %d is not a usable type code",      \
            "BINARY_SET_LT", PACKAGE_BUGREPORT, l->type);             \
    }



#endif
