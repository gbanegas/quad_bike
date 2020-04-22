/*
 * matrix.h
 *
 *  Created on: Apr 22, 2020
 *      Author: yoda
 */

#ifndef QUAD_BIKE_MATRIX_H_
#define QUAD_BIKE_MATRIX_H_

#include <stdint.h>
#include <stddef.h>

#include "../poly_arithmetic/export.h"
#include "../poly_arithmetic/poly_r_8192.h"

struct matrix {
    int rows; // number of rows
    int cols; // number of columns
    poly* data;
};
typedef struct matrix matrix;
matrix * generate_matrix(uint32_t n_rows, uint32_t n_cols);

void free_matrix(matrix *mat);

#endif /* QUAD_BIKE_MATRIX_H_ */
