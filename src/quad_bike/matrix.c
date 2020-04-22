/*
 * matrix.c
 *
 *  Created on: Apr 22, 2020
 *      Author: yoda
 */

#include "matrix.h"

matrix make_matrix(uint32_t n_rows, uint32_t n_cols) {
	matrix m;
	// set dimensions
	m.rows = n_rows;
	m.cols = n_cols;

	memset(m.data, 0, N_0*R_0*sizeof(poly));

	return m;

}

/*
 * free_matrix
 * 	Frees the memory that was allocated by make_matrix
 */
void free_matrix(matrix *mtx) {
	if (NULL != mtx) {
		free(mtx->data);
		free(mtx);
	}
}

matrix submatrix(const matrix *m, const int i, const int j, const int nr_row,
		const int nr_col) {

	matrix m_new = make_matrix(nr_row, nr_col);
	int j_temp = j;
	int i_temp = i;
	int new_row = 0;
	for (int row = 0; row < nr_row; row++) {
		int new_col = 0;
		j_temp = j;
		for (int col = 0; col < nr_col; col++) {
			m_new.data[new_row * nr_col + new_col] = m->data[i_temp * m->cols
					+ j_temp];
			new_col++;
			j_temp++;
		}
		i_temp++;
		new_row++;
	}
	return m_new;

}

matrix augment(const matrix *restrict a, const matrix *restrict b) {
	const int n_rows = a->rows;
	const int n_cols = a->cols + b->cols;
	matrix result = make_matrix(n_rows, n_cols);

	for (int i = 0; i < n_rows; i++) {
		memcpy(&result.data[result.cols * i], &a->data[a->cols * i],
				a->cols * sizeof(poly));
		memcpy(&result.data[result.cols * i + a->cols], &b->data[b->cols * i],
				(n_cols - a->cols) * sizeof(poly));
	}
	return result;
}

void echelon_form(matrix *a) {
	int nrows = a->rows;
	int ncols = a->cols;
	int c;
	int lead = 0;

	poly d, m;
	while (lead < nrows) {

		for (int r = 0; r < nrows; r++) { // for each row ...
			/* calculate divisor and multiplier */
			d = a->data[lead * ncols + lead];
			//m = gf_div(a->data[r * ncols + lead], a->data[lead * ncols + lead]);

			// for (int c = 0; c < ncols; c++) { // for each column ...
			// 	if (r == lead)
			// 		a->data[r * ncols + c] = gf_div(a->data[r * ncols + c], d); // make pivot = 1
			// 	else
			// 		a->data[r * ncols + c] ^= gf_mult(a->data[lead * ncols + c], m); // make other = 0
			// }

			if (r == lead) {
				for (c = 0; c < ncols; c++) {
					div_poly(&a->data[r * ncols + c], 0, &a->data[r * ncols + c], &d); // make pivot = 1
				}
			} else {
				for (c = 0; c < ncols; c++) {
				//	a->data[r * ncols + c] ^= gf_mult(a->data[lead * ncols + c],
				//			m); // make other = 0
					//TODO: Verify
				}
			}
		}
		lead++;
	}
}

matrix transpose_matrix(const matrix *restrict a) {

	int n_a_cols = a->cols;
	int n_a_rows = a->rows;
	matrix transposed = make_matrix(n_a_cols, n_a_rows);

	for (int c = 0; c < n_a_rows; c++)
		for (int d = 0; d < n_a_cols; d++)
			transposed.data[d * n_a_rows + c] = a->data[c * n_a_cols + d];

	return transposed;
}

void print_matrix(matrix *m) {
	int i, j;

	for (i = 0; i < m->rows; i++) {
		printf("| ");
		for (j = 0; j < m->cols; j++)
			print_polynomial(&m->data[i * m->cols + j]);
		printf(" |\n");
	}
	printf("\n");
}

void multiply_vector_matrix(poly *restrict u, matrix *G, poly *c) {

	int i, k;
	for (i = 0; i < G->cols; i++) {
		for (k = 0; k < G->rows; k++) {
			//c[i] ^= gf_mult(G->data[k * G->cols + i], u[k]);
		}
	}
}

/*
 void quasi_dyadic_bloc_matrix(matrix *M, poly *sig, const int ind_col, const int ind_rown) {
 int i, j;
 for (i = ind_rown; i < signature_block_size + ind_rown; i++) {
 for (j = ind_col; j < signature_block_size + ind_col; j++) {
 M->data[i * M->cols + j] = sig[(i ^ j) & (signature_block_size-1)];
 }
 }
 }
 */

