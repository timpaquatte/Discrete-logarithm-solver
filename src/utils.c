#include <stdio.h>

#include "gmp.h"
#include "utils.h"

/* Return the maximal k s.t. N is N1^k. */
int IsPerfectPower(mpz_t N1, mpz_t N){
    mpz_t tmp, r;
    int k = 2, e = 1;

    mpz_inits(tmp, r, NULL);
    do{
		mpz_rootrem(tmp, r, N, k);
		if(mpz_sgn(r) == 0){
			mpz_set(N1, tmp);
			e = k;
		}
		k++;
    } while(mpz_cmp_ui(tmp, 1) > 0);
    mpz_clears(tmp, r, NULL);
    return e;
}

void UpdateStatus(factor_t *fact){
    fact->status = PrimeStatus(fact->f);
}

void AddFactor(factor_t *fact, mpz_t f, int e, int status){
    mpz_init_set(fact->f, f);
    fact->e = e;
    fact->status = status;
}

void AddSmallFactor(factor_t *fact, int f, int e, int status){
    mpz_init_set_ui(fact->f, f);
    fact->e = e;
    fact->status = status;
}

void PrintFactor(factor_t fact){
    gmp_printf("[%Zd, %d]", fact.f, fact.e);
}

void PrintFactorization(factor_t *fact, int nf){
    int i;

    printf("[");
    for(i = 0; i < nf; i++){
	if(i > 0)
	    printf(", ");
	PrintFactor(fact[i]);
    }
    printf("]");
}

void factor_clear(factor_t* fact, int numberOfFactors){
    for(int i = 0; i < numberOfFactors; i++, fact++){
	mpz_clear(fact->f);
    }
}

