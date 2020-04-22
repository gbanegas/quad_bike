/* This file is part of the gf library.

 Copyright 2007, 2008, 2009, 2010, 2012, 2013, 2015
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

/* Main file for Karatsuba and Toom-Cook multiplication over GF(2)[x]. */

#ifndef SMALL_H_
#define SMALL_H_

#include "export.h"

void
gf_mul1(unsigned long *c, unsigned long a, unsigned long b);
unsigned long
gf_mul_1_n(unsigned long *cp, const unsigned long *bp, long sb, unsigned long a);
unsigned long
gf_addmul_1_n(unsigned long *dp, const unsigned long *cp,
		const unsigned long *bp, long sb, unsigned long a);
void
gf_mul2(unsigned long *c, const unsigned long *a, const unsigned long *b);
void
gf_mul3(unsigned long *c, const unsigned long *a, const unsigned long *b);
void
gf_mul4(unsigned long *c, const unsigned long *a, const unsigned long *b);
void
gf_mul5(unsigned long *c, const unsigned long *a, const unsigned long *b);
void
gf_mul6(unsigned long *c, const unsigned long *a, const unsigned long *b);
void
gf_mul7(unsigned long *c, const unsigned long *a, const unsigned long *b);
void
gf_mul8(unsigned long *c, const unsigned long *a, const unsigned long *b);
void
gf_mul9(unsigned long *c, const unsigned long *a, const unsigned long *b);


/* This file provides all the small-sized gf_mul1..gf_mul9 routines. It is
 * meant to be possibly included directly by applications. */



#endif  /* SMALL_H_ */
