#include <stdio.h>
#include <assert.h>

#include "gmp.h"
#include "utils.h"
// #include "refine.h"
#include "rho.h"



int PollardRhoSteps(mpz_t f, const mpz_t N,
					void (*fct)(mpz_t, mpz_t, const mpz_t), long *iter,
					long nrOfIterations, mpz_t x, mpz_t y){
    int status = FACTOR_NOT_FOUND;
    long i = 0;
    mpz_t abs;
    mpz_init(abs);
    mpz_init(f);
    for(i;i<nrOfIterations;i++){
        fct(x,x,N);
        fct(y,y,N);
        fct(y,y,N);
        if(mpz_cmp(x,y)>=0){
            mpz_sub(abs,x,y);
        } else {
            mpz_sub(abs,y,x);
        }
        mpz_gcd(f,abs,N);
        if(mpz_cmp_ui(f,1) != 0){
            break;
        }
    }
    mpz_clear(abs);
    if(mpz_cmp(f,N) == 0 || mpz_cmp_ui(f,1) == 0){
        return status;
    }
    status = FACTOR_FOUND;
    return status;
}

int PollardRho(factor_t* result, int *nf, mpz_t cof, const mpz_t N,
			   void (*fct)(mpz_t, mpz_t, const mpz_t),
			   long nrOfIterations){
    int status = FACTOR_NOT_FOUND;
    mpz_t x, y, fact;
	long iter = 0;
    mpz_inits(x, y, fact, NULL);
    mpz_set_ui(x, 1);	
    mpz_set_ui(y, 1);
    status = PollardRhoSteps(fact, N, fct, &iter, nrOfIterations, x, y);
	if(status == FACTOR_FOUND){
		AddFactor(result, fact, 1, FACTOR_IS_UNKNOWN);
		(*nf)++;
		mpz_divexact(cof, N, fact);
	}
    mpz_clears(x, y, fact,NULL);
    return status;
}
