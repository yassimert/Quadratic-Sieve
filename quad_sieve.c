/*
 * quad_sieve.c
 *
 *  	Date: June 8, 2017
 *      Author: Mert Yassı
 */
// QUADRATIC SIEVE

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gmp.h>	
#include"quad_sieve.h"

void sieve_of_eratosthenes(ARRAYS_PT ARRS, int n) {
	ull_int *primes, *factorArr, i, j, count = 0, k = 0;
	primes = (ull_int*)malloc(sizeof(ull_int) * n);

	for(i = 2; i < n; i++){
		primes[i] = 1; // make all elements true
	}

	for(i = 2; i < n; i++){
		if(primes[i] == 1){
			for(j = i; (i * j) < n; j++){
				primes[i * j] = 0; // mark composites
			}
		}
	}

	for(i = 2; i < n; i++){ // loop for determining numOfPrimes
		if(primes[i] == 1){
			count++;
		}
	}

	factorArr = (ull_int*)malloc(sizeof(ull_int) * count);
	for(i = 2; i < n; i++){
		if(primes[i] == 1){
			factorArr[k] = i;
			k++;
		}
	}
	ARRS->factorBase = factorArr;
	ARRS->numOfPrimes = count;
	free(primes);
}

void quadratic_sieve(mpz_t n, int FBBOUND, double THRESHOLD, ull_int M) {
	ARRAYS_PT ARRS = (ARRAYS_PT)malloc(sizeof(ARRAYS_T));
	sieve_of_eratosthenes(ARRS, FBBOUND);
	double *logValues;
	ull_int size = ARRS->numOfPrimes + 3, ct = 0, i, j, index = 0, **arr;
	mpz_t k, temp_k, tn, temp, tk, *kArr, *kValues;

	mpz_init(k);
	mpz_init(tk);
	mpz_init(tn);
	mpz_init(temp);
	mpz_init(temp_k);

	mpz_sub_ui(tn, n, 1);
	mpz_sqrt(k, tn);
	mpz_add_ui(k, k, 1);
	mpz_set(temp_k, k);

	arr = (ull_int**)malloc(sizeof(ull_int*) * ARRS->numOfPrimes);
	for(i = 0; i < ARRS->numOfPrimes; i++){
		arr[i] = (ull_int*)malloc(sizeof(ull_int) * size);
	}
	ARRS->exponents = (ull_int*)malloc(sizeof(ull_int) * ARRS->numOfPrimes);
	logValues = (double*)malloc(sizeof(double) * M);
	kArr = (mpz_t*)malloc(sizeof(mpz_t) * M);
	for(i = 0; i < M; i++){
		mpz_init(kArr[i]);
	}

	printf("Factorbase primes: ");
	for(i = 0; i < ARRS->numOfPrimes; i++){
		printf("%llu ",ARRS->factorBase[i]);
	}
	printf("\n");

	for(i = 0; i < M; i++){
		logValues[i] = 0;
		mpz_set(kArr[i], temp_k);
		mpz_add_ui(temp_k, temp_k, 1);
	}

	for(i = 0; i < ARRS->numOfPrimes; ){
		mpz_mod_ui(temp, kArr[index], ARRS->factorBase[i]);
		if(mpz_sgn(temp) == 0){
			double mylog = log(ARRS->factorBase[i]);
			for(j = index; j < M; j = j + ARRS->factorBase[i]){
				logValues[j] += mylog;
			}
			i++;
			index = 0;
		}else{
			index++;
		}
	}

	for(i = 0; i < M; i++){
		if(logValues[i] > THRESHOLD){
			ct++;
		}
	}

	j = 0;
	kValues = (mpz_t*)malloc(sizeof(mpz_t) * ct);
	for(i = 0; i < ct; i++){
		mpz_init(kValues[i]);
	}

	for(i = 0; i < M; i++){
		if(logValues[i] > THRESHOLD){
			mpz_add_ui(tk, k, i);
			mpz_set(kValues[j], tk);
			j++;
		}
	}

	arr = trial_division(ARRS, kValues, size, ct, arr, n);
	calc_nullspace(ARRS, arr);
	evaluate_nullspace(kValues, n, ct);

	mpz_clear(k);
	mpz_clear(tk);
	mpz_clear(tn);
	mpz_clear(temp);
	mpz_clear(temp_k);

	for(i = 0; i < M; i++){
		mpz_clear(kArr[i]);
	}

	for(i = 0; i < ct; i++){
		mpz_clear(kValues[i]);
	}

	for(i = 0; i < size; i++){
		mpz_clear(lArr[i]);
		mpz_clear(kArr2[i]);
	}

	free(ARRS->factorBase);
	free(ARRS->exponents);
	for(i = 0; i < ARRS->numOfPrimes; i++){
		free(arr[i]);
	}
	free(arr);
	free(kArr);
	free(kValues);
	free(logValues);
	free(ARRS);
}

ull_int** trial_division(ARRAYS_PT ARRS, mpz_t *kValues, ull_int size, ull_int ct, ull_int **arr, mpz_t n) {
	ull_int r, i, ctr = 0, j = 0;
	mpz_t l, mulVal, modTemp, lTemp, kTemp;
	mpz_init(l);
	mpz_init(mulVal);
	mpz_init(modTemp);
	mpz_init(lTemp);
	mpz_init(kTemp);

	lArr = (mpz_t*)malloc(sizeof(mpz_t) * size);
	kArr2 = (mpz_t*)malloc(sizeof(mpz_t) * size);

	for(i = 0; i < size; i++){
		mpz_init(lArr[i]);
		mpz_init(kArr2[i]);
	}

	for(r = 0; r < size; ){
		for(i = 0; i < ARRS->numOfPrimes; i++){
			ARRS->exponents[i] = 0;
		}
		mpz_mul(mulVal, kValues[ctr], kValues[ctr]);
		mpz_mod(l, mulVal, n);
		mpz_set(lTemp, l);
		mpz_set(kTemp, kValues[ctr]);
		int control = 0; i = 0;
		while(control == 0){
			mpz_mod_ui(modTemp, l, ARRS->factorBase[i]);
			if((mpz_sgn(modTemp) == 0) && (mpz_sgn(l) != 0)){
				while((mpz_sgn(modTemp) == 0) && (mpz_sgn(l) != 0)){
					mpz_divexact_ui(l, l, ARRS->factorBase[i]);
					ARRS->exponents[i]++;
					mpz_mod_ui(modTemp, l, ARRS->factorBase[i]);
				}
			}else{
				if(i < ARRS->numOfPrimes - 1){
					i++;
				}else{
					control = 1;
				}
			}
		}
		if(mpz_cmp_ui(l, 1) == 0){
			mpz_set(lArr[j], lTemp);
			mpz_set(kArr2[j], kTemp);
			j++;
			printf("  [ ");
			for(i = 0; i < ARRS->numOfPrimes; i++){
				arr[i][r] = ARRS->exponents[i] % 2;
				printf("%llu ",arr[i][r]);
			}
			printf("]\n");
			r++;
			ctrSize++;
		}
		if(ctr < ct){
			ctr++;
		}else{
			r = size;
		}
	}

	mpz_clear(l);
	mpz_clear(mulVal);
	mpz_clear(modTemp);
	return arr;
}

void calc_nullspace(ARRAYS_PT ARRS, ull_int **arr) {
	ull_int i, j, count = 0, count2 = 0;

	FILE *f = fopen("matrix.txt", "w");
	if (f == NULL){
		printf("Error opening file!\n");
		exit(1);
	}

	fprintf(f,"M:=Matrix(GF(2),%llu,%llu,[", ctrSize, ARRS->numOfPrimes);
	for(j = 0; j < ctrSize; j++){
		count++;
		for(i = 0; i < ARRS->numOfPrimes; i++){
			count2++;
			fprintf(f,"%llu",arr[i][j]);
			if(count2 < ARRS->numOfPrimes){
				fprintf(f,",");
			}
		}
		if(count < ctrSize){
			fprintf(f,",");
		}
		count2 = 0;
	}
	fprintf(f,"]);");
	fprintf(f,"\nN:=NullSpace(M);");
	fprintf(f,"\nB:=Basis(N);");
	fprintf(f,"\nB;");
	fclose(f);
	//int ret = system("~/magma < matrix.txt > result.txt"); // you can factor larger numbers using local Magma
	int ret = system("python3 encoder.py matrix.txt > result.txt"); // will not work if the size of matrix.txt > 50 kB (Online Magma input limit)
}

void evaluate_nullspace(mpz_t *kValues, mpz_t n, ull_int ctk) {
	mpz_t old_x, old_y, x, y2, y, xpy, xmy, p, q;
	ull_int i, j, k;
	mpz_init(old_x);
	mpz_init(old_y);
	mpz_init(y);
	mpz_init(xpy);
	mpz_init(xmy);
	mpz_init(p);
	mpz_init(q);
	mpz_init_set_ui(x, 1);
	mpz_init_set_ui(y2, 1);


	FILE *f = fopen("result.txt", "rb");
	fseek(f, 0, SEEK_END);
	long fsize = ftell(f);
	fseek(f, 0, SEEK_SET);
	char *string = malloc(fsize + 1);
	int ret = fread(string, fsize, 1, f);
	fclose(f);
	string[fsize] = 0;

	int sct = 1; // sct = 0 when using local Magma
	i = 0;
	while(string[i] != '\0'){
		if(string[i] == 44){
			sct++;
		}
		i++;
	}

	char **arr2;
	arr2 = (char**)malloc(sizeof(char*) * sct);
	for(i = 0; i < sct; i++){
		arr2[i] = (char*)malloc(sizeof(char) * ctrSize);
	}

	i = 0, j = 0, k = 0;
	while(string[i] != '\0'){
		if(string[i] == '('){
			while(string[i] != ')'){
				if(((string[i] == 48) || (string[i] == 49))){
					arr2[k][j] = string[i];
					j++;
				}
				i++;
			}
			j = 0;
			k++;
		}else{
			i++;
		}
	}

	if(sct == 1){
		printf("Magma could not find NullSpace of M!\n");
	}else{
		printf("NullSpace:\n");
		for(i = 0; i < sct; i++){
			printf("  [ ");
			for(j = 0; j < ctrSize; j++){
				printf("%c ",arr2[i][j]);
			}
			printf("]\n");
		}
	}
	int ctrl = 0;
	for(i = 0; i < sct; i++){
		for(j = 0; j < ctrSize; j++){
			if(arr2[i][j] == 49){
				mpz_set(old_x, kArr2[j]);
				mpz_set(old_y, lArr[j]);
				mpz_mul(x, old_x, x);
				mpz_mul(y2, old_y, y2);
			}
		}
		mpz_sqrt(y, y2);
		mpz_add(xpy, x, y);
		mpz_sub(xmy, x, y);
		mpz_gcd(p, xpy, n);
		mpz_gcd(q, xmy, n);

		if(((mpz_cmp_ui(p, 1) != 0) && (mpz_cmp(p, n) != 0)) && ((mpz_cmp_ui(q, 1) != 0) && (mpz_cmp(q, n) != 0))){
			printf("1ST FACTOR P: ");
			mpz_out_str(stdout, 10, p);
			printf("\n");
			printf("2ND FACTOR Q: ");
			mpz_out_str(stdout, 10, q);
			printf("\n");
			ctrl = 1;
			i = sct;
		}else{
			mpz_set_ui(x, 1);
			mpz_set_ui(y2, 1);
		}
	}
	if(ctrl != 1){
		printf("The algorithm could not find proper factors!\n");
	}

}




