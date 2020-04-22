/*
 * gf_mul.c
 *
 *  Created on: Apr 15, 2020
 *      Author: yoda
 */

#include "gf_mul.h"

void mul_base_case_inner(unsigned long *c, const unsigned long *a, long na,
		const unsigned long *b, long nb) {
	assert(c != a);
	assert(c != b);

	if (na == nb) {
		switch (na) {
		case 0:		// This can occur in call from KarMul return;
		case 1:
			gf_mul1(c, a[0], b[0]);
			return;
		case 2:
			gf_mul2(c, a, b);
			return;
		case 3:
			gf_mul3(c, a, b);
			return;
		case 4:
			gf_mul4(c, a, b);
			return;
		case 5:
			gf_mul5(c, a, b);
			return;
		case 6:
			gf_mul6(c, a, b);
			return;
		case 7:
			gf_mul7(c, a, b);
			return;
		case 8:
			gf_mul8(c, a, b);
			return;
		case 9:
			gf_mul9(c, a, b);
			return;
		default:
			return;
		}
	} else if (na < nb) {
		/* FIXME -- this does not seem efficient */
		long i;
		/* No need to clear c[], it's filled up progressively */
		c[nb] = gf_mul_1_n(c, b, nb, a[0]);
		for (i = 1; i < na; i++) {
			c[nb + i] = gf_addmul_1_n(c + i, c + i, b, nb, a[i]);
		}
	} else {
		mul_base_case_inner(c, b, nb, a, na);
	}
}

void print_polynomial2(const unsigned long *temp, int size) {

	int idx = 0;
	unsigned long long one = 1;
	for (int i = 0; i < size; i++) {
		for (unsigned long long j = 0; j < 64; j++) {
			unsigned long long result = (one << j);
			if (temp[i] & result) {
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

void gf_mul_kara(unsigned long *c, const unsigned long *a,
		const unsigned long *b, long n, unsigned long *stk) {
	unsigned long t;
	unsigned long *aa, *bb, *cc;
	long j, d, n2;

	assert(c != a);
	assert(c != b);

#if 0
    if (n <= 0)
    {				/* if turned on this test shows that calls with n == 0 */
	/* do occur (e.g from tunefft, FFT(19683)), but don't  */
	/* seem to be harmful if mul_basecase_n just returns.  */
	printf("\nWarning: n %ld in call to KarMul\n", n);
	fflush(stdout);
    }
#endif

	if (n < 10) {
		gf_mul_basecase(c, a, n, b, n);
		return;
	}

	n2 = (n + 1) / 2; /* ceil(n/2) */
	d = n & 1; /* 2*n2 - n = 1 if n odd, 0 if n even */
	aa = stk; /* Size n2   */
	bb = aa + n2; /* Size n2   */
	cc = bb + n2; /* Size n2   */

	stk = cc + n2; /* sp(n) = 3*ceil(n/2)) + sp(ceil(n/2)) */

	const unsigned long *a1 = a + n2; /* a[n2] */
	const unsigned long *b1 = b + n2; /* b[n2] */
	unsigned long *c1 = c + n2; /* c[n2]   */
	unsigned long *c2 = c1 + n2; /* c[2*n2] */
	unsigned long *c3 = c2 + n2; /* c[3*n2] */

	gf_mul_kara(c, a, b, n2, stk); /* Low */
	gf_mul_kara(c2, a1, b1, n2 - d, stk); /* High */

//	print_polynomial2(a, n);
//	print_polynomial2(a1, n2);

	//print_polynomial2(aa, stk);

	for (j = 0; j < n2 - d; j++) {
		aa[j] = a[j] ^ a1[j];
		bb[j] = b[j] ^ b1[j];
		cc[j] = c1[j] ^ c2[j];
	}
	for (; j < n2; j++) { /* Only when n odd */
		aa[j] = a[j];
		bb[j] = b[j];
		cc[j] = c1[j] ^ c2[j];
	}

	gf_mul_kara(c1, aa, bb, n2, stk); /* Middle */

	for (j = 0; j < n2 - 2 * d; j++) {
		t = cc[j];
		c1[j] ^= t ^ c[j];
		c2[j] ^= t ^ c3[j];
	}
	for (; j < n2; j++) { /* Only when n odd */
		c1[j] ^= cc[j] ^ c[j];
		c2[j] ^= cc[j];
	}
}

void gf_mul1(unsigned long *c, unsigned long a, unsigned long b) {
	__m128i aa = _gf2x_mm_setr_epi64(a, 0);
	__m128i bb = _gf2x_mm_setr_epi64(b, 0);
	_mm_storeu_si128((__m128i*) c, _mm_clmulepi64_si128(aa, bb, 0));
}

unsigned long gf_mul_1_n(unsigned long *cp, const unsigned long *bp, long sb,
		unsigned long a) {
	long i;
	typedef union {
		__m128i s;
		unsigned long x[2];
	} __m128i_proxy;

	__m128i y = _gf2x_mm_set1_epi64(a);
	__m128i x;
	__m128i_proxy cc;

	// do two at a time
	for (i = 0; i + 2 < sb; i += 2) {
		x = _gf2x_mm_setr_epi64(bp[i], bp[i + 1]);
		cc.s = _mm_clmulepi64_si128(x, y, 0);
		if (i == 0)
			cp[i] = cc.x[0];
		else
			cp[i] ^= cc.x[0];
		cp[i + 1] = cc.x[1];
		cc.s = _mm_clmulepi64_si128(x, y, 1);
		cp[i + 1] ^= cc.x[0];
		cp[i + 2] = cc.x[1];
	}
	// last is different, to handle carry out
	unsigned long cy;
	if (i == sb - 2) {  // case bp is even
		x = _gf2x_mm_setr_epi64(bp[i], bp[i + 1]);
		cc.s = _mm_clmulepi64_si128(x, y, 0);
		if (i == 0)
			cp[i] = cc.x[0];
		else
			cp[i] ^= cc.x[0];
		cp[i + 1] = cc.x[1];
		cc.s = _mm_clmulepi64_si128(x, y, 1);
		cp[i + 1] ^= cc.x[0];
		cy = cc.x[1];
	} else { //case bp is odd
		x = _gf2x_mm_setr_epi64(bp[i], 0);
		cc.s = _mm_clmulepi64_si128(x, y, 0);
		if (i == 0)
			cp[i] = cc.x[0];
		else
			cp[i] ^= cc.x[0];
		cy = cc.x[1];
	}
	return cy;
}

unsigned long gf_addmul_1_n(unsigned long *dp, const unsigned long *cp,
		const unsigned long *bp, long sb, unsigned long a) {
	long i;
	typedef union {
		__m128i s;
		unsigned long x[2];
	} __m128i_proxy;

	__m128i y = _gf2x_mm_set1_epi64(a);
	__m128i x;
	__m128i_proxy dd;

	// do two at a time
	for (i = 0; i + 2 < sb; i += 2) {
		x = _gf2x_mm_setr_epi64(bp[i], bp[i + 1]);
		dd.s = _mm_clmulepi64_si128(x, y, 0);
		if (i == 0)
			dp[i] = cp[i] ^ dd.x[0];
		else
			dp[i] ^= dd.x[0];
		dp[i + 1] = cp[i + 1] ^ dd.x[1];
		dd.s = _mm_clmulepi64_si128(x, y, 1);
		dp[i + 1] ^= dd.x[0];
		dp[i + 2] = cp[i + 2] ^ dd.x[1];
	}
	unsigned long cy;
	if (i == sb - 2) {  // case bp is even
		x = _gf2x_mm_setr_epi64(bp[i], bp[i + 1]);
		dd.s = _mm_clmulepi64_si128(x, y, 0);
		if (i == 0)
			dp[i] = cp[i] ^ dd.x[0];
		else
			dp[i] ^= dd.x[0];
		dp[i + 1] = cp[i + 1] ^ dd.x[1];
		dd.s = _mm_clmulepi64_si128(x, y, 1);
		dp[i + 1] ^= dd.x[0];
		cy = dd.x[1];
	} else {
		x = _gf2x_mm_setr_epi64(bp[i], 0);
		dd.s = _mm_clmulepi64_si128(x, y, 0);
		if (i == 0)
			dp[i] = cp[i] ^ dd.x[0];
		else
			dp[i] ^= dd.x[0];
		cy = dd.x[1];
	}
	return cy;
}

/* Karatsuba with 3 multiplications */

void gf_mul2(unsigned long *t, unsigned long const *s1, unsigned long const *s2) {
#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
	__m128i ss1 = _mm_loadu_si128((__m128i*) s1);
	__m128i ss2 = _mm_loadu_si128((__m128i*) s2);

	__m128i t00 = _mm_clmulepi64_si128(ss1, ss2, 0);
	__m128i t11 = _mm_clmulepi64_si128(ss1, ss2, 0x11);

	ss1 = PXOR(ss1, _mm_shuffle_epi32(ss1, _MM_SHUFFLE(1,0,3,2)));
	ss2 = PXOR(ss2, _mm_shuffle_epi32(ss2, _MM_SHUFFLE(1,0,3,2)));

	__m128i tk = PXOR(t00, PXOR(t11, _mm_clmulepi64_si128(ss1, ss2, 0)));

	/* mul2cl.c is essentially identical, just replaces srli and srli by
	 * unpacklo and unpackhi */
	_mm_storeu_si128((__m128i*) t, PXOR(t00, _mm_slli_si128(tk, 8)));
	_mm_storeu_si128((__m128i*) (t + 2), PXOR(t11, _mm_srli_si128(tk, 8)));
#undef PXOR
}

static inline __m128i mul3cl_mul1(unsigned long a, unsigned long b) {
	__m128i aa = _gf2x_mm_setr_epi64(a, 0);
	__m128i bb = _gf2x_mm_setr_epi64(b, 0);
	return _mm_clmulepi64_si128(aa, bb, 0);
}

/* uses the variant of Karatsuba with 6 multiplications */

void gf_mul3(unsigned long *c, const unsigned long *a, const unsigned long *b) {
	unsigned long aa[3], bb[3];
	aa[0] = a[1] ^ a[2];
	aa[1] = a[0] ^ a[2];
	aa[2] = a[0] ^ a[1];
	bb[0] = b[1] ^ b[2];
	bb[1] = b[0] ^ b[2];
	bb[2] = b[0] ^ b[1];
	__m128i p0 = mul3cl_mul1(a[0], b[0]);
	__m128i p1 = mul3cl_mul1(a[1], b[1]);
	__m128i p2 = mul3cl_mul1(a[2], b[2]);
	__m128i pp0 = mul3cl_mul1(aa[0], bb[0]);
	__m128i pp1 = mul3cl_mul1(aa[1], bb[1]);
	__m128i pp2 = mul3cl_mul1(aa[2], bb[2]);

#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))

	__m128i ce0 = p0;
	__m128i ce2 = PXOR(p0, PXOR(p1, PXOR(p2, pp1)));
	__m128i ce4 = p2;

	__m128i co1 = PXOR(p0, PXOR(p1, pp2));
	__m128i co3 = PXOR(pp0, PXOR(p1, p2));

	_mm_storeu_si128((__m128i*) (c), PXOR(ce0, _mm_slli_si128(co1, 8)));
	_mm_storeu_si128((__m128i*) (c + 2),
			PXOR(ce2, PXOR(_mm_srli_si128(co1, 8), _mm_slli_si128(co3, 8))));
	_mm_storeu_si128((__m128i*) (c + 4), PXOR(ce4, _mm_srli_si128(co3, 8)));
#undef PXOR
}

#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PZERO    _mm_setzero_si128()
static inline void

mul4clk_mul2(__m128i *t, __m128i ss1, __m128i ss2) {
	__m128i t00 = _mm_clmulepi64_si128(ss1, ss2, 0);
	__m128i t11 = _mm_clmulepi64_si128(ss1, ss2, 0x11);
	ss1 = PXOR(ss1, _mm_shuffle_epi32(ss1, _MM_SHUFFLE(1,0,3,2)));
	ss2 = PXOR(ss2, _mm_shuffle_epi32(ss2, _MM_SHUFFLE(1,0,3,2)));
	__m128i tk = PXOR(PXOR(t00, t11), _mm_clmulepi64_si128(ss1, ss2, 0));
	t[0] = PXOR(t00, _mm_unpacklo_epi64(PZERO, tk));
	t[1] = PXOR(t11, _mm_unpackhi_epi64(tk, PZERO));
}
/* specialized Karatsuba with 3 calls to mul2, i.e., 9 multiplications */

void gf_mul4(unsigned long *c, const unsigned long *a, const unsigned long *b) {
	__m128i ab[2];
	__m128i lo[2], hi[2];
	__m128i a0 = _mm_loadu_si128((__m128i*) a);
	__m128i a2 = _mm_loadu_si128((__m128i*) (a + 2));
	__m128i b0 = _mm_loadu_si128((__m128i*) b);
	__m128i b2 = _mm_loadu_si128((__m128i*) (b + 2));
	mul4clk_mul2(lo, a0, b0);
	mul4clk_mul2(hi, a2, b2);
	__m128i middle = PXOR(lo[1], hi[0]);
	mul4clk_mul2(ab, PXOR(a0, a2), PXOR(b0, b2));
	_mm_storeu_si128((__m128i*) (c + 0), lo[0]);
	_mm_storeu_si128((__m128i*) (c + 2), PXOR(ab[0], PXOR(lo[0], middle)));
	_mm_storeu_si128((__m128i*) (c + 4), PXOR(ab[1], PXOR(hi[1], middle)));
	_mm_storeu_si128((__m128i*) (c + 6), hi[1]);
}
#undef PXOR
#undef PZERO

static inline __m128i mul5clk_c_mul1(unsigned long a, unsigned long b) {
	return _mm_clmulepi64_si128(_gf2x_mm_setr_epi64(a,0),
			_gf2x_mm_setr_epi64(b,0), 0);
}

void gf_mul5(unsigned long *c, const unsigned long *a, const unsigned long *b) {
#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
	/* Montgomery formulae with 13 multiplications, see
	 Five, Six, and Seven-Term {K}aratsuba-Like Formulae,
	 IEEE Transactions on Computers, volume 54, number 3, p. 362-369, 2005 */
	unsigned long ta[3], tb[3];

	ta[0] = a[0] ^ a[4];
	tb[0] = b[0] ^ b[4];
	ta[1] = a[1] ^ a[2];
	tb[1] = b[1] ^ b[2];
	ta[2] = a[3] ^ ta[0];
	tb[2] = b[3] ^ tb[0];

	unsigned long pa[8], pb[8];
	pa[0] = ta[1] ^ ta[2];
	pb[0] = tb[1] ^ tb[2];
	pa[1] = a[2] ^ ta[2];
	pb[1] = b[2] ^ tb[2];
	pa[2] = ta[0] ^ ta[1];
	pb[2] = tb[0] ^ tb[1];
	pa[3] = a[1] ^ ta[2];
	pb[3] = b[1] ^ tb[2];
	pa[4] = a[0] ^ a[2] ^ a[3];
	pb[4] = b[0] ^ b[2] ^ b[3];
	pa[5] = a[4] ^ ta[1];
	pb[5] = b[4] ^ tb[1];
	pa[6] = a[3] ^ a[4];
	pb[6] = b[3] ^ b[4];
	pa[7] = a[0] ^ a[1];
	pb[7] = b[0] ^ b[1];

	__m128i p0 = mul5clk_c_mul1(pa[0], pb[0]);
	__m128i p2 = mul5clk_c_mul1(pa[1], pb[1]);
	__m128i p4 = mul5clk_c_mul1(pa[2], pb[2]);
	__m128i p6 = mul5clk_c_mul1(pa[3], pb[3]);
	__m128i p8 = mul5clk_c_mul1(pa[4], pb[4]);
	__m128i p10 = mul5clk_c_mul1(pa[5], pb[5]);
	__m128i p12 = mul5clk_c_mul1(pa[6], pb[6]);
	__m128i p14 = mul5clk_c_mul1(pa[7], pb[7]);
	__m128i p16 = mul5clk_c_mul1(ta[0], tb[0]);
	__m128i p18 = mul5clk_c_mul1(a[4], b[4]);
	__m128i p20 = mul5clk_c_mul1(a[3], b[3]);
	__m128i p22 = mul5clk_c_mul1(a[1], b[1]);
	__m128i p24 = mul5clk_c_mul1(a[0], b[0]);
	__m128i t0 = PXOR(p14, p24);
	__m128i t2 = PXOR(p12, p18);
	__m128i t4 = PXOR(p2, p16);
	__m128i t6 = PXOR(p0, p6);
	__m128i t8 = PXOR(p4, p16);
	__m128i t10 = PXOR(p10, t0);
	__m128i t12 = PXOR(p8, t2);

	__m128i ce0 = p24;
	__m128i ce2 = PXOR(p18, PXOR(t8, t10));
	__m128i ce4 = PXOR(p0, PXOR(p20, PXOR(p22, PXOR(t10, t12))));
	__m128i ce6 = PXOR(p24, PXOR(t4, t12));
	__m128i ce8 = p18;

	__m128i co1 = PXOR(p22, t0);
	__m128i co3 = PXOR(t2, PXOR(t4, t6));
	__m128i co5 = PXOR(t0, PXOR(t6, t8));
	__m128i co7 = PXOR(p20, t2);

	_mm_storeu_si128((__m128i*) (c), PXOR(ce0, _mm_slli_si128(co1, 8)));
	_mm_storeu_si128((__m128i*) (c + 2),
			PXOR(ce2, PXOR(_mm_srli_si128(co1, 8), _mm_slli_si128(co3, 8))));
	_mm_storeu_si128((__m128i*) (c + 4),
			PXOR(ce4, PXOR(_mm_srli_si128(co3, 8), _mm_slli_si128(co5, 8))));
	_mm_storeu_si128((__m128i*) (c + 6),
			PXOR(ce6, PXOR(_mm_srli_si128(co5, 8), _mm_slli_si128(co7, 8))));
	_mm_storeu_si128((__m128i*) (c + 8), PXOR(ce8, _mm_srli_si128(co7, 8)));
#undef PXOR
}

#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PXOR3(op1, op2, op3) PXOR(op1, PXOR(op2, op3))
#define PXOR4(op1, op2, op3, op4) PXOR(op1, PXOR3(op2, op3, op4))
#define PZERO    _mm_setzero_si128()

static inline void mul6clk2_mul2(__m128i *t, __m128i ss1, __m128i ss2) {
	__m128i t00 = _mm_clmulepi64_si128(ss1, ss2, 0);
	__m128i t11 = _mm_clmulepi64_si128(ss1, ss2, 0x11);
	ss1 = PXOR(ss1, _mm_shuffle_epi32(ss1, _MM_SHUFFLE(1,0,3,2)));
	ss2 = PXOR(ss2, _mm_shuffle_epi32(ss2, _MM_SHUFFLE(1,0,3,2)));
	__m128i tk = PXOR(PXOR(t00, t11), _mm_clmulepi64_si128(ss1, ss2, 0));
	t[0] = PXOR(t00, _mm_unpacklo_epi64(PZERO, tk));
	t[1] = PXOR(t11, _mm_unpackhi_epi64(tk, PZERO));
}

/* variant with 6 calls to mul2, i.e., 18 multiplications */

void gf_mul6(unsigned long *c, const unsigned long *a, const unsigned long *b) {
	__m128i aa[3], bb[3];
	__m128i p0[2], p1[2], p2[2];
	__m128i pp0[2], pp1[2], pp2[2];
	__m128i a0 = _mm_loadu_si128((__m128i*) (a));
	__m128i a1 = _mm_loadu_si128((__m128i*) (a + 2));
	__m128i a2 = _mm_loadu_si128((__m128i*) (a + 4));
	__m128i b0 = _mm_loadu_si128((__m128i*) (b));
	__m128i b1 = _mm_loadu_si128((__m128i*) (b + 2));
	__m128i b2 = _mm_loadu_si128((__m128i*) (b + 4));
	aa[0] = PXOR(a1, a2);
	aa[1] = PXOR(a0, a2);
	aa[2] = PXOR(a0, a1);
	bb[0] = PXOR(b1, b2);
	bb[1] = PXOR(b0, b2);
	bb[2] = PXOR(b0, b1);
	mul6clk2_mul2(p0, a0, b0);
	mul6clk2_mul2(p1, a1, b1);
	mul6clk2_mul2(p2, a2, b2);
	mul6clk2_mul2(pp0, aa[0], bb[0]);
	mul6clk2_mul2(pp1, aa[1], bb[1]);
	mul6clk2_mul2(pp2, aa[2], bb[2]);
	_mm_storeu_si128((__m128i*) (c + 0), p0[0]);
	_mm_storeu_si128((__m128i*) (c + 2),
			PXOR(PXOR3(p0[0], p1[0], pp2[0]), p0[1]));
	_mm_storeu_si128((__m128i*) (c + 4),
			PXOR(PXOR4(p0[0], p1[0], p2[0], pp1[0]),
					PXOR3(p0[1], p1[1], pp2[1])));
	_mm_storeu_si128((__m128i*) (c + 6),
			PXOR(PXOR3(pp0[0], p1[0], p2[0]),
					PXOR4(p0[1], p1[1], p2[1], pp1[1])));
	_mm_storeu_si128((__m128i*) (c + 8),
			PXOR(p2[0], PXOR3(pp0[1], p1[1], p2[1])));
	_mm_storeu_si128((__m128i*) (c + 10), p2[1]);
}

#undef PXOR
#undef PXOR3
#undef PXOR4
#undef PZERO

static inline __m128i mul7cl_mul1(unsigned long a, unsigned long b) {
	return _mm_clmulepi64_si128(_gf2x_mm_setr_epi64(a,0),
			_gf2x_mm_setr_epi64(b,0), 0);
}

/* variant with 22 multiplications */

void gf_mul7(unsigned long *c, const unsigned long *a, const unsigned long *b) {
#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PXOR3(op1, op2, op3) PXOR(op1, PXOR(op2, op3))
#define PXOR4(op1, op2, op3, op4) PXOR(op1, PXOR3(op2, op3, op4))
#define PXOR5(op1, op2, op3, op4, op5) PXOR(op1, PXOR4(op2, op3, op4, op5))
#define PXOR6(op1, op2, op3, op4, op5, op6) PXOR(op1, PXOR5(op2, op3, op4, op5, op6))
#define PZERO    _mm_setzero_si128()
	/* Montgomery formulae with 22 multiplications, see
	 Five, Six, and Seven-Term {K}aratsuba-Like Formulae,
	 IEEE Transactions on Computers, volume 54, number 3, p. 362-369, 2005 */
	unsigned long ta[5], tb[5], pa[22], pb[22];
	ta[0] = a[0] ^ a[4];
	ta[1] = a[3] ^ a[5];
	ta[2] = a[2] ^ a[6];
	ta[3] = a[1] ^ ta[0];
	ta[4] = ta[1] ^ ta[2];
	pa[0] = a[6];
	pa[1] = a[5];
	pa[2] = a[5] ^ a[6];
	pa[3] = a[4];
	pa[4] = a[4] ^ a[6];
	pa[5] = a[3];
	pa[6] = ta[1];
	pa[7] = a[2];
	pa[8] = ta[2];
	pa[9] = a[1];
	pa[10] = a[1] ^ a[3];
	pa[11] = a[1] ^ a[2] ^ a[4] ^ a[5];
	pa[12] = a[1] ^ ta[4];
	pa[13] = a[0];
	pa[14] = ta[0];
	pa[15] = a[0] ^ a[2];
	pa[16] = a[0] ^ ta[4];
	pa[17] = a[3] ^ ta[0] ^ ta[2];
	pa[18] = a[0] ^ a[1];
	pa[19] = a[3] ^ a[6] ^ ta[3];
	pa[20] = ta[1] ^ ta[3];
	pa[21] = ta[3] ^ ta[4];

	tb[0] = b[0] ^ b[4];
	tb[1] = b[3] ^ b[5];
	tb[2] = b[2] ^ b[6];
	tb[3] = b[1] ^ tb[0];
	tb[4] = tb[1] ^ tb[2];
	pb[0] = b[6];
	pb[1] = b[5];
	pb[2] = b[5] ^ b[6];
	pb[3] = b[4];
	pb[4] = b[4] ^ b[6];
	pb[5] = b[3];
	pb[6] = tb[1];
	pb[7] = b[2];
	pb[8] = tb[2];
	pb[9] = b[1];
	pb[10] = b[1] ^ b[3];
	pb[11] = b[1] ^ b[2] ^ b[4] ^ b[5];
	pb[12] = b[1] ^ tb[4];
	pb[13] = b[0];
	pb[14] = tb[0];
	pb[15] = b[0] ^ b[2];
	pb[16] = b[0] ^ tb[4];
	pb[17] = b[3] ^ tb[0] ^ tb[2];
	pb[18] = b[0] ^ b[1];
	pb[19] = b[3] ^ b[6] ^ tb[3];
	pb[20] = tb[1] ^ tb[3];
	pb[21] = tb[3] ^ tb[4];

	__m128i p[22];

	p[0] = mul7cl_mul1(pa[0], pb[0]);
	p[1] = mul7cl_mul1(pa[1], pb[1]);
	p[2] = mul7cl_mul1(pa[2], pb[2]);
	p[3] = mul7cl_mul1(pa[3], pb[3]);
	p[4] = mul7cl_mul1(pa[4], pb[4]);
	p[5] = mul7cl_mul1(pa[5], pb[5]);
	p[6] = mul7cl_mul1(pa[6], pb[6]);
	p[7] = mul7cl_mul1(pa[7], pb[7]);
	p[8] = mul7cl_mul1(pa[8], pb[8]);
	p[9] = mul7cl_mul1(pa[9], pb[9]);
	p[10] = mul7cl_mul1(pa[10], pb[10]);
	p[11] = mul7cl_mul1(pa[11], pb[11]);
	p[12] = mul7cl_mul1(pa[12], pb[12]);
	p[13] = mul7cl_mul1(pa[13], pb[13]);
	p[14] = mul7cl_mul1(pa[14], pb[14]);
	p[15] = mul7cl_mul1(pa[15], pb[15]);
	p[16] = mul7cl_mul1(pa[16], pb[16]);
	p[17] = mul7cl_mul1(pa[17], pb[17]);
	p[18] = mul7cl_mul1(pa[18], pb[18]);
	p[19] = mul7cl_mul1(pa[19], pb[19]);
	p[20] = mul7cl_mul1(pa[20], pb[20]);
	p[21] = mul7cl_mul1(pa[21], pb[21]);

	__m128i t[13];

	t[0] = PXOR(p[0], p[1]);
	t[1] = PXOR(p[9], p[13]);
	t[2] = PXOR(p[3], p[6]);
	t[3] = PXOR(p[7], p[10]);
	t[4] = PXOR(p[11], p[18]);
	t[5] = PXOR(p[4], t[3]);
	t[6] = PXOR(p[15], t[2]);
	t[7] = PXOR(p[20], t[5]);
	t[8] = PXOR(p[5], p[14]);
	t[9] = PXOR(p[2], p[17]);
	t[10] = PXOR(p[5], p[8]);
	t[11] = PXOR(p[21], t[6]);
	t[12] = PXOR(p[16], t[4]);

	__m128i cc[13];

	cc[0] = p[13];
	cc[2] = PXOR3(p[7], p[15], t[1]);
	cc[4] = PXOR4(p[3], t[1], t[3], t[8]);
	cc[6] = PXOR5(p[0], p[12], p[13], t[7], t[11]);
	cc[8] = PXOR4(p[7], t[0], t[2], t[10]);
	cc[10] = PXOR3(p[3], p[4], t[0]);
	cc[12] = p[0];

	cc[1] = PXOR(p[18], t[1]);
	cc[3] = PXOR3(t[7], t[9], t[12]);
	cc[5] = PXOR6(p[2], p[11], p[19], t[1], t[10], t[11]);
	cc[7] = PXOR5(p[21], t[0], t[5], t[8], t[12]);
	cc[9] = PXOR5(p[12], p[19], t[4], t[6], t[9]);
	cc[11] = PXOR(p[2], t[0]);

	_mm_storeu_si128((__m128i*) (c), PXOR(cc[0], _mm_slli_si128(cc[1], 8)));
	_mm_storeu_si128((__m128i*) (c + 2),
			PXOR(cc[2],
					PXOR(_mm_srli_si128(cc[1], 8), _mm_slli_si128(cc[3], 8))));
	_mm_storeu_si128((__m128i*) (c + 4),
			PXOR(cc[4],
					PXOR(_mm_srli_si128(cc[3], 8), _mm_slli_si128(cc[5], 8))));
	_mm_storeu_si128((__m128i*) (c + 6),
			PXOR(cc[6],
					PXOR(_mm_srli_si128(cc[5], 8), _mm_slli_si128(cc[7], 8))));
	_mm_storeu_si128((__m128i*) (c + 8),
			PXOR(cc[8],
					PXOR(_mm_srli_si128(cc[7], 8), _mm_slli_si128(cc[9], 8))));
	_mm_storeu_si128((__m128i*) (c + 10),
			PXOR(cc[10],
					PXOR(_mm_srli_si128(cc[9], 8), _mm_slli_si128(cc[11], 8))));
	_mm_storeu_si128((__m128i*) (c + 12),
			PXOR(cc[12], _mm_srli_si128(cc[11], 8)));
#undef PXOR
#undef PXOR3
#undef PXOR4
#undef PXOR5
#undef PXOR6
#undef PZERO
}

void gf_mul8(unsigned long *c, const unsigned long *a, const unsigned long *b) {
	/* specialized Karatsuba, RPB 20070518 */
	/* slightly faster on bogong than version with loops */
	/* this version uses minimal temporary storage (12 = 3*n/2 words) */
	unsigned long aa[4], bb[4], cc[4];
	gf_mul4(c + 8, a + 4, b + 4);
	gf_mul4(c, a, b);
	cc[0] = c[4] ^ c[8];
	cc[1] = c[5] ^ c[9];
	cc[2] = c[6] ^ c[10];
	cc[3] = c[7] ^ c[11];
	aa[0] = a[0] ^ a[4];
	aa[1] = a[1] ^ a[5];
	aa[2] = a[2] ^ a[6];
	aa[3] = a[3] ^ a[7];
	bb[0] = b[0] ^ b[4];
	bb[1] = b[1] ^ b[5];
	bb[2] = b[2] ^ b[6];
	bb[3] = b[3] ^ b[7];
	gf_mul4(c + 4, aa, bb);
	c[4] ^= c[0] ^ cc[0];
	c[5] ^= c[1] ^ cc[1];
	c[6] ^= c[2] ^ cc[2];
	c[7] ^= c[3] ^ cc[3];
	c[8] ^= c[12] ^ cc[0];
	c[9] ^= c[13] ^ cc[1];
	c[10] ^= c[14] ^ cc[2];
	c[11] ^= c[15] ^ cc[3];
}

#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PXOR3(op1, op2, op3) PXOR(op1, PXOR(op2, op3))
#define PXOR4(op1, op2, op3, op4) PXOR(op1, PXOR3(op2, op3, op4))
#define PXOR5(op1, op2, op3, op4, op5) PXOR(op1, PXOR4(op2, op3, op4, op5))
#define PXOR6(op1, op2, op3, op4, op5, op6) PXOR(op1, PXOR5(op2, op3, op4, op5, op6))
#define PXOR7(op1, op2, op3, op4, op5, op6, op7) PXOR(op1, PXOR6(op2, op3, op4, op5, op6, op7))
#define PZERO    _mm_setzero_si128()

/* variant with 30 multiplications */
void gf_mul9(unsigned long *c, const unsigned long *a, const unsigned long *b) {
	/* Taken from Cenk & Ozbudak 2009 */
	/* We reserve one more to follow notations of the paper */
	__m128i ab[9] = { _gf2x_mm_setr_epi64(a[0], b[0]), _gf2x_mm_setr_epi64(a[1],
			b[1]), _gf2x_mm_setr_epi64(a[2], b[2]), _gf2x_mm_setr_epi64(a[3],
			b[3]), _gf2x_mm_setr_epi64(a[4], b[4]), _gf2x_mm_setr_epi64(a[5],
			b[5]), _gf2x_mm_setr_epi64(a[6], b[6]), _gf2x_mm_setr_epi64(a[7],
			b[7]), _gf2x_mm_setr_epi64(a[8], b[8]), };
	__m128i pab[30];

#if 0
    pab[ 0] = ab[0]^ab[1]^ab[2]^ab[3]^ab[4]^ab[5]^ab[6]^ab[7]^ab[8];
    pab[ 1] = ab[0]^      ab[2]^      ab[4]^      ab[6]^      ab[8];
    pab[ 2] =       ab[1]^ab[2]^ab[3]^      ab[5]^            ab[8];
    pab[ 3] = ab[0]^      ab[2]^ab[3]^      ab[5]^ab[6]^      ab[8];
    pab[ 4] = ab[0]^ab[1]^      ab[3]^ab[4]^      ab[6]^ab[7];
    pab[ 5] = ab[0]^            ab[3]^ab[4]^ab[5]^      ab[7];
    pab[ 6] =       ab[1]^ab[2]^ab[3]^            ab[6]^      ab[8];
    pab[ 7] =             ab[2]^      ab[4]^ab[5]^ab[6];
    pab[ 8] =             ab[2]^ab[3]^ab[4]^      ab[6];
    pab[ 9] =       ab[1]^      ab[3]^      ab[5]^      ab[7];
    pab[10] = ab[0]^ab[1]^            ab[4]^      ab[6]^ab[7]^ab[8];
    pab[11] = ab[0]^            ab[3]^      ab[5]^ab[6]^ab[7];
    pab[12] = ab[0]^ab[1]^            ab[4]^ab[5]^            ab[8];
    pab[13] =       ab[1]^ab[2]^      ab[4]^ab[5]^      ab[7]^ab[8];
    pab[14] = ab[0]^ab[1]^      ab[3]^            ab[6]^ab[7]^ab[8];
    pab[15] =       ab[1]^      ab[3]^ab[4]^ab[5]^            ab[8];
    pab[16] = ab[0]^      ab[2]^ab[3]^ab[4]^            ab[7];
    pab[17] =       ab[1]^            ab[4]^ab[5]^ab[6]^      ab[8];
    pab[18] = ab[0]^      ab[2]^            ab[5]^ab[6]^ab[7];
    pab[19] =             ab[2]^ab[3]^            ab[6]^ab[7];
    pab[20] =                                     ab[6]^      ab[8];
    pab[21] = ab[0]^      ab[2];
    pab[22] = ab[0]^ab[1];
    pab[23] = ab[0];
    pab[24] =       ab[1];
    pab[25] =                                           ab[7];
    pab[26] =                                           ab[7]^ab[8];
    pab[27] =                                     ab[6];
    pab[28] =                                                 ab[8];
    pab[29] =             ab[2];
#else
	/* same as above, but optimized with Maple's codegen[optimize] function
	 with 'tryhard' option: 89 XORs -> 46 XORs */
	__m128i t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63,
			t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75;
	t51 = ab[8];
	t55 = ab[4];
	t75 = PXOR(t51, t55);
	t54 = ab[5];
	t57 = ab[2];
	t74 = PXOR(t54, t57);
	t56 = ab[3];
	t58 = ab[1];
	t73 = PXOR(t56, t58);
	t59 = ab[0];
	t72 = PXOR(t59, t57);
	t71 = PXOR(t58, t75);
	t52 = ab[7];
	t70 = PXOR3(t52, t56, t59);
	t53 = ab[6];
	t69 = PXOR3(t53, t56, t57);
	t68 = PXOR(t53, t72);
	t67 = PXOR(t54, t71);
	t66 = PXOR(t59, t71);
	t65 = PXOR3(t51, t57, t73);
	t64 = PXOR(t53, t70);
	t63 = PXOR(t55, t70);
	t62 = PXOR(t54, t68);
	t61 = PXOR(t58, t64);
	t60 = PXOR(t55, t61);
	pab[0] = PXOR3(t51, t60, t74);
	pab[1] = PXOR(t68, t75);
	pab[2] = PXOR(t54, t65);
	pab[3] = PXOR3(t56, t51, t62);
	pab[4] = t60;
	pab[5] = PXOR(t54, t63);
	pab[6] = PXOR(t53, t65);
	pab[7] = PXOR3(t55, t53, t74);
	pab[8] = PXOR(t55, t69);
	pab[9] = PXOR3(t54, t52, t73);
	pab[10] = PXOR3(t53, t52, t66);
	pab[11] = PXOR(t54, t64);
	pab[12] = PXOR(t54, t66);
	pab[13] = PXOR3(t57, t52, t67);
	pab[14] = PXOR(t51, t61);
	pab[15] = PXOR(t56, t67);
	pab[16] = PXOR(t57, t63);
	pab[17] = PXOR(t53, t67);
	pab[18] = PXOR(t52, t62);
	pab[19] = PXOR(t52, t69);
	pab[20] = PXOR(t53, t51);
	pab[21] = t72;
	pab[22] = PXOR(t59, t58);
	pab[23] = t59;
	pab[24] = t58;
	pab[25] = t52;
	pab[26] = PXOR(t52, t51);
	pab[27] = t53;
	pab[28] = t51;
	pab[29] = t57;
#endif

	int i;
	for (i = 0; i < 30; ++i)
		pab[i] = _mm_clmulepi64_si128(pab[i], pab[i], 0x10);

	__m128i cc[17];

#if 0
    cc[0 ] = pab[23];
    cc[1 ] = pab[22]^pab[23]^pab[24];
    cc[2 ] = pab[21]^pab[23]^pab[24]^pab[29];
    cc[3 ] = pab[28]^pab[17]^pab[2]^pab[15]^pab[7]^pab[6]^pab[5]^pab[29]^pab[21]^pab[22]^pab[12]^pab[19]^pab[9]^pab[13]^pab[11]^pab[3]^pab[26]^pab[20]^pab[27];
    cc[4 ] = pab[4]^pab[3]^pab[10]^pab[11]^pab[6]^pab[2]^pab[8]^pab[14]^pab[9]^pab[22]^pab[23]^pab[24]^pab[1]^pab[20]^pab[27]^pab[28]^pab[25];
    cc[5 ] = pab[26]^pab[25]^pab[28]^pab[0]^pab[9]^pab[21]^pab[23]^pab[29]^pab[24]^pab[1]^pab[3]^pab[13]^pab[14]^pab[5]^pab[18]^pab[16]^pab[11]^pab[15];
    cc[6 ] = pab[26]^pab[12]^pab[19]^pab[21]^pab[23]^pab[29]^pab[4]^pab[3]^pab[14]^pab[5]^pab[18]^pab[22]^pab[1]^pab[20]^pab[27];
    cc[7 ] = pab[20]^pab[27]^pab[28]^pab[25]^pab[23]^pab[0]^pab[15]^pab[7]^pab[11]^pab[6]^pab[14]^pab[5]^pab[18];
    cc[8 ] = pab[0]^pab[23]^pab[24]^pab[10]^pab[15]^pab[7]^pab[2]^pab[18]^pab[14]^pab[17]^pab[22]^pab[26]^pab[25]^pab[28];
    cc[9 ] = pab[21]^pab[23]^pab[29]^pab[24]^pab[0]^pab[16]^pab[11]^pab[7]^pab[10]^pab[2]^pab[8]^pab[18]^pab[5]^pab[28];
    cc[10] = pab[12]^pab[0]^pab[19]^pab[9]^pab[21]^pab[29]^pab[22]^pab[3]^pab[13]^pab[16]^pab[11]^pab[7]^pab[10]^pab[20]^pab[27]^pab[28]^pab[26];
    cc[11] = pab[16]^pab[11]^pab[7]^pab[10]^pab[17]^pab[5]^pab[2]^pab[0]^pab[9]^pab[4]^pab[3]^pab[22]^pab[23]^pab[24]^pab[1]^pab[20]^pab[27]^pab[28]^pab[25];
    cc[12] = pab[26]^pab[25]^pab[28]^pab[8]^pab[14]^pab[5]^pab[17]^pab[10]^pab[6]^pab[16]^pab[15]^pab[3]^pab[13]^pab[1]^pab[9]^pab[21]^pab[23]^pab[29]^pab[24];
    cc[13] = pab[8]^pab[18]^pab[2]^pab[15]^pab[16]^pab[5]^pab[29]^pab[21]^pab[23]^pab[22]^pab[0]^pab[12]^pab[19]^pab[1]^pab[11]^pab[4]^pab[3]^pab[26]^pab[20]^pab[27];
    cc[14] = pab[20]^pab[27]^pab[28]^pab[25];
    cc[15] = pab[25]^pab[26]^pab[28];
    cc[16] = pab[28];
#else
	/* same as above, optimized with codegen[optimize] with 'tryhard' */
	__m128i t100, t101, t102, t103, t104, t105, t106, t107, t108, t109, t110,
			t111, t112, t113, t114, t115, t116, t117, t118, t119, t120, t121,
			t122, t123, t124, t125, t126, t127, t128, t129, t130, t77, t79, t80,
			t82, t83, t87, t88, t89, t90, t91, t92, t94, t95, t96, t97, t98,
			t99;
	t82 = pab[23];
	t87 = pab[18];
	t130 = PXOR(t82, t87);
	t77 = pab[28];
	t98 = pab[7];
	t129 = PXOR(t77, t98);
	t79 = pab[26];
	t83 = pab[22];
	t128 = PXOR(t79, t83);
	t90 = pab[15];
	t91 = pab[14];
	t127 = PXOR(t90, t91);
	t97 = pab[8];
	t99 = pab[6];
	t126 = PXOR(t97, t99);
	t100 = pab[5];
	t125 = PXOR(t100, t90);
	t117 = PXOR(pab[27], pab[20]);
	t80 = pab[25];
	t118 = PXOR(t77, t80);
	t112 = PXOR(t117, t118);
	t94 = pab[11];
	t124 = PXOR(t112, t94);
	t103 = pab[2];
	t105 = pab[0];
	t123 = PXOR(t103, t105);
	t89 = pab[16];
	t122 = PXOR3(t89, t94, t97);
	t121 = PXOR3(t100, t105, t98);
	t102 = pab[3];
	t104 = pab[1];
	t96 = pab[9];
	t120 = PXOR3(t102, t104, t96);
	t119 = PXOR(pab[29], pab[21]);
	t116 = PXOR(pab[24], t82);
	t115 = PXOR(t79, t118);
	t114 = PXOR(t83, t116);
	t113 = PXOR(t116, t119);
	t95 = pab[10];
	t111 = PXOR5(t87, t95, t116, t123, t129);
	t110 = PXOR6(t102, pab[19], pab[12], t117, t119, t128);
	t92 = pab[13];
	t109 = PXOR5(t92, t94, t96, t110, t129);
	t101 = pab[4];
	t108 = PXOR5(t100, t101, t104, t110, t130);
	t107 = PXOR6(t101, t103, t95, t114, t120, t124);
	t106 = PXOR7(t89, t91, t92, t113, t115, t120, t125);
	t88 = pab[17];
	cc[0] = t82;
	cc[1] = t114;
	cc[2] = t113;
	cc[3] = PXOR5(t88, t99, t103, t109, t125);
	cc[4] = PXOR3(t91, t107, t126);
	cc[5] = PXOR4(t87, t94, t105, t106);
	cc[6] = PXOR(t91, t108);
	cc[7] = PXOR5(t99, t121, t124, t127, t130);
	cc[8] = PXOR5(t88, t80, t111, t127, t128);
	cc[9] = PXOR4(t100, t111, t119, t122);
	cc[10] = PXOR4(t89, t95, t105, t109);
	cc[11] = PXOR4(t88, t89, t107, t121);
	cc[12] = PXOR4(t88, t95, t106, t126);
	cc[13] = PXOR4(t90, t108, t122, t123);
	cc[14] = t112;
	cc[15] = t115;
	cc[16] = t77;
#endif

	_mm_storeu_si128((__m128i*) (c), PXOR(cc[0], _mm_slli_si128(cc[1], 8)));
	_mm_storeu_si128((__m128i*) (c + 2),
			PXOR3(cc[2], _mm_srli_si128(cc[1], 8), _mm_slli_si128(cc[3], 8)));
	_mm_storeu_si128((__m128i*) (c + 4),
			PXOR3(cc[4], _mm_srli_si128(cc[3], 8), _mm_slli_si128(cc[5], 8)));
	_mm_storeu_si128((__m128i*) (c + 6),
			PXOR3(cc[6], _mm_srli_si128(cc[5], 8), _mm_slli_si128(cc[7], 8)));
	_mm_storeu_si128((__m128i*) (c + 8),
			PXOR3(cc[8], _mm_srli_si128(cc[7], 8), _mm_slli_si128(cc[9], 8)));
	_mm_storeu_si128((__m128i*) (c + 10),
			PXOR3(cc[10], _mm_srli_si128(cc[9], 8), _mm_slli_si128(cc[11], 8)));
	_mm_storeu_si128((__m128i*) (c + 12),
			PXOR3(cc[12], _mm_srli_si128(cc[11], 8),
					_mm_slli_si128(cc[13], 8)));
	_mm_storeu_si128((__m128i*) (c + 14),
			PXOR3(cc[14], _mm_srli_si128(cc[13], 8),
					_mm_slli_si128(cc[15], 8)));
	_mm_storeu_si128((__m128i*) (c + 16),
			PXOR(cc[16], _mm_srli_si128(cc[15], 8)));
#undef PXOR
#undef PXOR3
#undef PXOR4
#undef PXOR5
#undef PXOR6
#undef PXOR7
#undef PZERO
}

