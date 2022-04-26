/*
 * quad_sieve.h
 *
 * 		Date: June 8, 2017
 *      Author: Mert Yassı
 */
// QUADRATIC SIEVE

#ifndef QUAD_SIEVE_H_
#define QUAD_SIEVE_H_
typedef unsigned long long int ull_int;

typedef struct ARRAYS *ARRAYS_PT;
typedef struct ARRAYS{
	ull_int numOfPrimes;
	ull_int *factorBase, *exponents;
}ARRAYS_T[1];

static ull_int ctrSize = 0;
static mpz_t *lArr, *kArr2;

void sieve_of_eratosthenes(ARRAYS_PT ARRS, int n);
void quadratic_sieve(mpz_t n, int FBBOUND, double THRESHOLD, ull_int M);
ull_int** trial_division(ARRAYS_PT ARRS, mpz_t *kValues, ull_int size, ull_int ct, ull_int **arr, mpz_t n);
void calc_nullspace(ARRAYS_PT ARRS, ull_int **arr);
void evaluate_nullspace(mpz_t *kValues, mpz_t n, ull_int ctk);

#endif /* QUAD_SIEVE_H_ */
