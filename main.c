/*
 * main.c
 *
 *  	Date: June 10, 2017
 *      Author: Mert Yassı
 */
// QUADRATIC SIEVE

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gmp.h>
#include"quad_sieve.h"

//mpz_set_str(n, "794490252688143120079", 10);
#define M 10000000
#define FBBOUND 800
#define THRESHOLD 5

/*
//mpz_set_str(n, "51104184699520549779052870841", 10); // test it with local magma
#define M 15000000
#define FBBOUND 3800
#define THRESHOLD 4
*/

/*
//mpz_set_str(n, "6241019306901997355512619014111", 10); // test it with local magma
#define M 20000000
#define FBBOUND 6000
#define THRESHOLD 4
*/

int main() {
	mpz_t n;
	mpz_init(n);
	mpz_set_str(n, "794490252688143120079", 10);
	quadratic_sieve(n, FBBOUND, THRESHOLD, M);
	return 0;
}











