#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "gmp.h"
#include "utils.h"
// #include "refine.h"
#include "pminus1.h"

#define DEBUG 0
#define BUF_LENGTH 32


static void ReadDifFile(mpz_t p, FILE *file){
    mpz_t q;
    int dp;

    /* we have to read file while read_p <= p */
    mpz_init_set_ui(q, 1);
    while(fscanf(file, "%d", &dp) != EOF){
		mpz_add_ui(q, q, dp << 1);
		if(mpz_cmp(q, p) > 0)
			break;
    }
    mpz_set(p, q);
    mpz_clear(q);
}

/* Starting from nextprime(p) >= p+1. 
   On unsuccessful exit, p is the smallest prime > bound1.
*/
int PollardPminus1Step1(mpz_t f, const mpz_t N,	long bound1, FILE* file,
						mpz_t b, mpz_t p){
	char* line = (char*) malloc(BUF_LENGTH * sizeof(char));
	long k0 = 0;
	mpz_t curr_prime,temp,q;
	int ret = FACTOR_NOT_FOUND;
	
	mpz_set_ui(b, 2);	
	mpz_init_set_ui(q, 1);
	mpz_init(temp);
	mpz_init_set_ui(curr_prime, 2);
	
	while(mpz_cmp_ui(curr_prime,bound1) <= 0){

		mpz_set_ui(q, 1);
		mpz_set(temp, curr_prime);

		while(mpz_cmp_ui(temp, bound1) <= 0){
			mpz_set(q, temp);
			mpz_mul(temp, temp, curr_prime);
		}

		mpz_powm(b, b, q, N);
		mpz_sub_ui(temp, b, 1);
		mpz_gcd(f, temp, N);

		if(mpz_cmp_ui(f, 1) > 0){
			ret = FACTOR_FOUND;
			break;
		}

		fgets(line, BUF_LENGTH, file);
		if(line == NULL) {
			ret = FACTOR_FOUND;
			break;
		}
		
		mpz_set_str(temp, line, 10);
		if(k0 == 0)
			mpz_mul_ui(temp, temp, 2);
		mpz_add(curr_prime, curr_prime, temp);
		k0++;
	}

	mpz_clear(curr_prime);
	mpz_clear(temp);
	mpz_clear(q);
	return ret;
}

int PollardPminus1Step2(mpz_t f, const mpz_t N, long bound2, FILE* file,
						mpz_t b, mpz_t p){
    mpz_t bm1;
    unsigned long d;
    int dp, status = FACTOR_NOT_FOUND;
    int B = (int)log((double)bound2);
    B = B * B;

    mpz_init(bm1);
    ReadDifFile(p, file);
    /* Precomputations */
    mpz_t* precomputations = (mpz_t*)malloc(B * sizeof(mpz_t));
    mpz_t* cursor = precomputations;
		
    for(int i = 0; i < B; i++, cursor++){
		mpz_init(*cursor);
		mpz_powm_ui(*cursor, b, i, N);
    }
#if DEBUG >= 1
    printf("# Precomputation of phase 2 done.\n");
#endif
    mpz_powm(b, b, p, N);
    while(mpz_cmp_ui(p, bound2) <= 0){
		mpz_sub_ui(bm1, b, 1);
		mpz_gcd(f, bm1, N);
		if(mpz_cmp_ui(f, 1) > 0){
			status = FACTOR_FOUND;
			break;
		}
		fscanf(file, "%d", &dp);
		d = dp << 1;
		mpz_add_ui(p, p, d);		
		if(d < B){
			mpz_mul(b, b, precomputations[d]);
			mpz_mod(b, b, N);
		}
		else{
			printf("Cramer's rule Failed!\n");
			printf("WRITE A PAPER!!!\n");
			return 1;
			//mpz_powm(b, b, precomputations[d], cof);
		}
    }			
    cursor = precomputations;
    for(int i = 0; i < B; i++, cursor++){
		mpz_clear(*cursor);
    }
    free(precomputations);
    mpz_clear(bm1);
    return status;
}

int PollardPminus1(factor_t* result, int *nf, mpz_t cof, const mpz_t N,
				   long bound1, long bound2, FILE* file){
	mpz_t b, p, factor;
	mpz_inits(b, p, factor, NULL);
	mpz_set_ui(b, 2);
	mpz_set_ui(p, 2);
	int status = PollardPminus1Step1(factor, N, bound1, file, b, p);
	if(status != FACTOR_FOUND)
		status = PollardPminus1Step2(factor, N, bound2, file, b, p);
	if(status == FACTOR_FOUND){
		AddFactor(result + *nf, factor, 1, FACTOR_IS_UNKNOWN);
		(*nf)++;
	}
	mpz_clears(b, p, factor, NULL);
	return status;
}
