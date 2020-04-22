/*
 ============================================================================
 Name        : qd_bike.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <sodium.h>

#include "poly_arithmetic/poly_r_8192.h"
#include "poly_arithmetic/export.h"

#include "quad_bike/matrix.h"

#include "quad_bike/qd_bike.h"

#include "params.h"

int main(void) {
	puts("!!!Hello World!!!"); /* prints !!!Hello World!!! */

	if (sodium_init() == -1) {
		return 1;
	}

	matrix pk;
	matrix sk;
	sk = make_matrix(N_0, R_0);

	key_gen(&pk, &sk);

	printf("SK: \n");
	print_matrix(&sk);

	puts("!!! Finish the fish !!!");
	return EXIT_SUCCESS;
}
