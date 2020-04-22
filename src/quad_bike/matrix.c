/*
 * matrix.c
 *
 *  Created on: Apr 22, 2020
 *      Author: yoda
 */

#include "matrix.h"

matrix* generate_matrix(uint32_t n_rows, uint32_t n_cols) {
	matrix *m = (matrix*) malloc(sizeof(matrix));
	// set dimensions
	m->rows = n_rows;
	m->cols = n_cols;

	// allocate a double array of length rows * cols
	m->data = (poly*) calloc(n_rows * n_cols, sizeof(poly));

	return m;

}

/*
 * free_matrix
 * 	Frees the memory that was allocated by make_matrix
 */
void free_matrix(matrix* mtx) {
	if (NULL != mtx){
		free(mtx->data);
		free(mtx);
	}
}
