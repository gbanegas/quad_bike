/*
 * gf_bit_op.c
 *
 *  Created on: Apr 15, 2020
 *      Author: yoda
 */

#include "gf_bit_op.h"

/*----------------------------------------------------------------------------*/

void right_bit_shift(const int length, element_p in[]) {
	int j;
	__m256i a, b;
	for (j = length - 1; j > 4; j = j - 4) {

		a = _mm256_lddqu_si256((__m256i*) &in[j - 3]); //load from in[j-3] to in[j]
		b = _mm256_lddqu_si256((__m256i*) &in[j - 4]); //load from in[j-4] to in[j-1]

		a = _mm256_slli_epi64(a, 1);
		b = _mm256_srli_epi64(b, (64 - 1));

		_mm256_storeu_si256(((__m256i*) &in[j - 3]), _mm256_or_si256(a, b));

	}
	for (; j > 0; j--) {
		in[j] <<= 1;
		in[j] |= (in[j - 1] & (element_p) 0x01) << (64 - 1);
	}
	in[j] <<= 1;

} // end right_bit_shift

void left_bit_shift(const int length, element_p in[]) {
	int j;
	__m256i a, b;
	for (j = length - 1; j > 4; j = j - 4) {

		a = _mm256_lddqu_si256((__m256i*) &in[j - 3]); //load from in[j-3] to in[j]
		b = _mm256_lddqu_si256((__m256i*) &in[j - 4]); //load from in[j-4] to in[j-1]

		a = _mm256_srli_epi64(a, 1);
		b = _mm256_slli_epi64(b, (64 - 1));

		_mm256_storeu_si256(((__m256i*) &in[j - 3]), _mm256_or_si256(a, b));

	}
	for (; j > 0; j--) {
		in[j] >>= 1;
		in[j] |= (in[j - 1] & (element_p) 0x01) << (64 - 1);
	}
	in[j] >>= 1;

} // end right_bit_shift

void rotate_bit_right(element_p in[]) /*  equivalent to x * in(x) mod x^P+1 */
{

	element_p mask, rotated_bit;

	/* NUM_DIGITS_GF2X_MODULUS == 1 + NUM_DIGITS_GF2X_ELEMENT and
	 * MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS == 0
	 */
	mask = ((element_p) 0x1) << (64 - 1);
	rotated_bit = !!(in[127] & mask);
	in[127] &= ~mask; /* clear shifted bit */
	right_bit_shift(128, in);

	in[0] |= rotated_bit;
} // end rotate_bit_left

void rotate_bit_left(element_p in[]) /*  equivalent to x * in(x) mod x^P+1 */
{

	element_p mask, rotated_bit;

	/* NUM_DIGITS_GF2X_MODULUS == 1 + NUM_DIGITS_GF2X_ELEMENT and
	 * MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS == 0
	 */
	mask = ((element_p) 0x1);
	rotated_bit = !!(in[0] & mask);
	in[0] &= ~mask; /* clear shifted bit */
	left_bit_shift(128, in);

	in[127] |= (rotated_bit << (64 - 1));
} // end rotate_bit_left

 __m256i right_bit_shift_n(__m256i *data, int count) {
	__m256i innerCarry, carryOut, rotate;
	innerCarry = _mm256_srli_epi64(*data, 64 - count); // carry outs in bit 0 of each qword
	rotate = _mm256_permute4x64_epi64(innerCarry, 0x93); // rotate ymm left 64 bits
	innerCarry = _mm256_blend_epi32(_mm256_setzero_si256(), rotate, 0xFC); // clear lower qword
	*data = _mm256_slli_epi64(*data, count);        // shift all qwords left
	*data = _mm256_or_si256(*data, innerCarry); // propagate carrys from low qwords
	carryOut = _mm256_xor_si256(innerCarry, rotate); // clear all except lower qword
	return carryOut;
}
__attribute__((optimize("unroll-loops")))
inline void shift_n_bits_to_right(int count, element_p *in) {
	for (int i = 0; i < 124; i = i + 4) {
		__m256i tmp;
		__m256i a = _mm256_lddqu_si256((__m256i*) &in[i]);
		__m256i carryOut = right_bit_shift_n(&a, count);

		_mm256_storeu_si256(((__m256i*) &in[i]), a);

		__m256i a4 = _mm256_lddqu_si256((__m256i*) &in[i + 4]);
		tmp = _mm256_or_si256(a4, carryOut);
		_mm256_storeu_si256(((__m256i*) &in[i + 4]), tmp);
		//printf("i: %d -  i+4: %d\n", i, (i + 4));
	}
	__m256i a = _mm256_lddqu_si256((__m256i*) &in[124]);
	right_bit_shift_n(&a, count);

	_mm256_storeu_si256(((__m256i*) &in[124]), a);

}


__attribute__((optimize("unroll-loops")))
void tmp(unsigned long *in, int count) {

	__m256i tmp;
	__m256i as[32];
	__m256i carry[32] = { 0 };
	int j = 0;
	for (int i = 0; i < 32; i++) {
		as[i] = _mm256_lddqu_si256((__m256i*) &in[j]);
		j = j + 4;
	}
	//__m256i a = _mm256_lddqu_si256((__m256i*) &in[i]);
	carry[0] = right_bit_shift_n(&as[0], count);

	_mm256_storeu_si256(((__m256i*) &in[0]), as[0]);
	for (int i = 1; i < 32; i++) {
		carry[i] = right_bit_shift_n(&as[i], count);
	}

	j = 4;
	for (int i = 1; i < 32; i++) {
		tmp = _mm256_or_si256(as[i], carry[i - 1]);
		_mm256_storeu_si256(((__m256i*) &in[j]), tmp);
		j = j + 4;

	}

}


__attribute__((optimize("unroll-loops")))
void righ_bit_shift_by_any(element_p *in, int count) {
	int temp = count;
	int n = temp % 63;
	tmp(in, n);
	temp = temp - n;
	temp = temp / 63;
	for (int i = 0; i < temp; i++) {
		tmp(in, 63);
		//print_polynomial(in);
	}

}
