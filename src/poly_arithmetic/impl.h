/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2015
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of either:
    - If the archive contains a file named toom-gpl.c (not a trivial
    placeholder), the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
    - If the archive contains a file named toom-gpl.c which is a trivial
    placeholder, the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
   
   You should have received a copy of the GNU General Public License as
   well as the GNU Lesser General Public License along with this program;
   see the files COPYING and COPYING.LIB.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301, USA.
*/

#ifndef IMPL_H_
#define IMPL_H_

/* first include the stuff that even users of the library have access to */
#include "export.h"

struct mul_pool {
	unsigned long * stk;
	size_t stk_size;
};
typedef struct mul_pool mul_pool_t[1];

/* then proceed to the really internal stuff */

/* These flags are for internal use. When a new routine is added, don't
   change the flags, otherwise the tuning in the already_tuned directory
   will become invalid. */
#define	GF2X_SELECT_KARA	0	/* do not change ! */
#define	GF2X_SELECT_TC3		1
#define	GF2X_SELECT_TC3W	2
#define	GF2X_SELECT_TC4		3
#define	GF2X_SELECT_KARAX	4
#define	GF2X_SELECT_TC3X	5

#define	GF2X_SELECT_UNB_DFLT	0
#define	GF2X_SELECT_UNB_TC3U	1	/* do not change ! */

#include <assert.h>
#ifndef ASSERT
#define ASSERT(x)	assert(x)
#endif

/* We use it for internal checking of the proper propagation of errors.
 */
/* https://gcc.gnu.org/onlinedocs/gcc-3.3.6/gcc/Function-Attributes.html#Function-Attributes
 * https://gcc.gnu.org/onlinedocs/gcc-3.4.6/gcc/Function-Attributes.html#Function-Attributes
 */
#define GF2X_LEXGE2(X,Y,A,B) (X>A || (X == A && Y >= B))
#define GF2X_LEXGE3(X,Y,Z,A,B,C) (X>A || (X == A && GF2X_LEXGE2(Y,Z,B,C)))
#define GF2X_LEXLE2(X,Y,A,B) GF2X_LEXGE2(A,B,X,Y)
#define GF2X_LEXLE3(X,Y,Z,A,B,C) GF2X_LEXGE3(A,B,C,X,Y,Z)

#define MAX(a,b)        ((a)<(b) ? (b) : (a))

extern void mul_base_case_inner(unsigned long * c, const unsigned long * a,
			 long na, const unsigned long * b, long nb);

#define gf_mul_basecase(c,a,na,b,nb) do {			        \
    mul_base_case_inner(c, a, na, b, nb);				\
} while (0)


#endif	/* IMPL_H_ */
