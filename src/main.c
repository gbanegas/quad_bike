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

int main(void) {
	puts("!!!Hello World!!!"); /* prints !!!Hello World!!! */

	if (sodium_init() == -1) {
		return 1;
	}

	poly *random = create_random_polynomial_with_weight(30);

	print_polynomial(random);

	int result = pop_count(random);
	printf("weight = %d \n", result);
	puts("!!! Finish the fish !!!");
	return EXIT_SUCCESS;
}
