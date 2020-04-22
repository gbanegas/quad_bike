/*
 * poly_r_8096.c
 *
 *  Created on: Apr 14, 2020
 *      Author: yoda
 */

#include "poly_r_8192.h"

poly* create_polynomial() {
	poly *polynomial = (poly*) calloc(1, sizeof(poly));
	return polynomial;
}

void set_pos(int pos, poly *polynomial) {
	element_p one = 1;
	int idx_block = pos / 64;
	int idx_pos = 64 - ((idx_block * 64) - pos);
	polynomial->coeffs[idx_block] = polynomial->coeffs[idx_block]
			| (one << idx_pos);

}

void flip_pos(int pos, poly *polynomial) {
	element_p one = 1;
	int idx_block = pos / 64;
	int idx_pos = 64 - ((idx_block * 64) - pos);

	polynomial->coeffs[idx_block] = (polynomial->coeffs[idx_block]
			^ (one << idx_pos));

}

void add_poly(poly *result, const poly *p1, const poly *p2) {

	for (int i = 0; i < 128; i++) {
		result->coeffs[i] = p1->coeffs[i] ^ p2->coeffs[i];
	}
}

long toomspace(unsigned long n) {
	return 5 * n + 30;
}

long toomuspace(unsigned long sa) {
	return 2 * sa + 32 + toomspace(sa / 4 + 4);
}

int gf_mul_r(unsigned long *c, const unsigned long *a, unsigned long sa,
		const unsigned long *b, unsigned long sb) {
	int rc = 0;
	unsigned long sc = sa + sb;
	/* As a starting guess, assume that c may alias a or b */
	unsigned long *dst = c;
	mul_pool_t xxpool;
	struct mul_pool *xpool = NULL;

	unsigned long sp, sp2;
	sp = toomspace(sa); // Space for balanced TC routines

	if (sa != sb) {
		sp2 = toomuspace(2 * sa); // Space for unbalanced TC routines
		if (sp < sp2)
			sp = sp2; /* worst-case required */
		sp2 = 2 * sa + toomspace(sa); // Space for unbalanced TC routines w/ lazy cut
		if (sp < sp2)
			sp = sp2; /* worst-case required */
	}

	memset(xxpool, 0, sizeof(mul_pool_t));
	xpool = xxpool;

	if (xpool->stk_size < sp) {
		void *p = realloc(xpool->stk, sp * sizeof(element_p));
		if (p == NULL) {
			rc = -1;
			goto end_of_gf2x_mul_r;
		}
		xpool->stk = p;
		xpool->stk_size = sp;
		free(p);
	}

	// Avoid copy in common case
	gf_mul_kara(dst, a, b, sa, xpool->stk);

	end_of_gf2x_mul_r:
	//gf2x_mul_pool_clear(xxpool);

	if (dst && dst != c) {
		/* Then we have allocated a temp buffer */
		memcpy(c, dst, sc * sizeof(element_p));
		free(dst);
	}

	return rc;
}

void mul_poly(poly *result, const poly *p1, const poly *p2) {
	unsigned long r_tmp[256] = { 0 };
	gf_mul_r(r_tmp, p1->coeffs, 128, p2->coeffs, 128);
	for (int i = 0; i < 128; i++) {
		result->coeffs[i] = r_tmp[i] ^ r_tmp[i + 128];
	}

}

unsigned long get_deg(const poly *p) {
	unsigned long degree = 8191;
	while (degree > 0) {
		int idx_block = degree / 64;
		element_p one = 1;
		int idx_pos = 64 - ((idx_block * 64) - degree);
		if (p->coeffs[idx_block] & (one << idx_pos))
			return degree;
		degree--;
	}
	return 0;
}

void div_poly(poly *quo, poly *re, const poly *dividend, const poly *divisor) {

	unsigned long degree_dividend = get_deg(dividend);
	unsigned long degree_divisor = get_deg(divisor);
	poly temp = *divisor;
	poly d_copy = *dividend;
	print_polynomial(&d_copy);

	//quotient = create_polynomial(degree_quotient);
	while (degree_dividend >= degree_divisor) {
		righ_bit_shift_by_any(temp.coeffs, (degree_dividend - degree_divisor));

	//	print_polynomial(&temp);

		set_pos(degree_dividend - degree_divisor, quo);
		//print_polynomial(quo);

		add_poly(&d_copy, &temp, &d_copy);
	//	print_polynomial(&d_copy);
		degree_dividend = get_deg(&d_copy);

		temp = *divisor;

	}
	*re = d_copy;
	//polynomial_free(dividend);

}

__attribute__((optimize("unroll-loops")))
int pop_count(const poly *p) {
	int result = 0;
	for(int i = 0; i < 128; i ++){
		result +=__builtin_popcountl(p->coeffs[i]);
 	}
	return result;

}

void print_polynomial(poly *polynomial) {

	int idx = 0;
	element_p one = 1;
	for (int i = 0; i < 128; i++) {
		for (element_p j = 0; j < 64; j++) {
			element_p result = (one << j);
			if (polynomial->coeffs[i] & result) {
				if (idx != 0)
					printf("x^%d +", idx);
				else
					printf("1 +");
			}
			idx++;
		}
	}
	printf("\n");

}

