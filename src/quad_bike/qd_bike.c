/*
 * qd_bike.c
 *
 *  Created on: Apr 22, 2020
 *      Author: yoda
 */

#include "qd_bike.h"

int key_gen(matrix *pk, matrix *sk) {
	for(int i = 0; i < N_0; i++){
		for(int j = 0; j  < R_0;j++){
			sk->data[i * sk->cols + j] = create_random_polynomial_with_weight(WEIGHT);
		}
	}

	copy_matrix(pk, sk);

	echelon_form(pk);




	return 0;
}
