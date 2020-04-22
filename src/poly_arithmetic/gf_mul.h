/*
 * gf_mul.h
 *
 *  Created on: Apr 15, 2020
 *      Author: yoda
 */

#ifndef GF_MUL_H_
#define GF_MUL_H_

#include <stdio.h>
#include "impl.h"
#include "export.h"
#include "small.h"

/*
#include "mul1.h"
#include "mul2.h"
#include "mul3.h"
#include "mul4.h"
#include "mul5.h"
#include "mul6.h"
#include "mul7.h"
#include "mul8.h"
#include "mul9.h"
*/

extern void gf_mul_kara(unsigned long * c, const unsigned long * a, const unsigned long * b,
	      long n, unsigned long * stk);

#endif /* GF_MUL_H_ */
