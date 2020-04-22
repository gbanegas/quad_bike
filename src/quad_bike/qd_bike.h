/*
 * qd_bike.h
 *
 *  Created on: Apr 22, 2020
 *      Author: yoda
 */

#ifndef QUAD_BIKE_QD_BIKE_H_
#define QUAD_BIKE_QD_BIKE_H_

#include <stdlib.h>
#include <sodium.h>
#include <keccak/KangarooTwelve.h>

#include "../poly_arithmetic/export.h"
#include "../poly_arithmetic/poly_r_8192.h"

#include "../params.h"

#include "matrix.h"

int key_gen(matrix *pk, matrix *sk);



#endif /* QUAD_BIKE_QD_BIKE_H_ */
