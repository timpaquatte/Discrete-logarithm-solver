#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include "time.h"
#include "limits.h"

#include "factor.h"
#include "trialdiv.h"
#include "pminus1.h"
#include "rho.h"


void goodFunction(mpz_t out, mpz_t in, const mpz_t N){
    mpz_powm_ui(out, in, 2, N);
    mpz_add_ui(out, out, 7);
    mpz_mod(out, out, N);
}


void factorize(factor_t factors[], mpz_t N, int* nbFactors) 
{
    int nf, count;
    clock_t finish, start;
    // factor_t factors[MAX_NUM_OF_FACTORS];
    mpz_t cof;
    FILE* file = fopen("ressources/1e7", "r");
    if(file == NULL) {
        printf("Could not find the file ressources/1e7\n");
        exit(0);
    }

    start = clock();
    mpz_init(cof);

    gmp_printf("N : %Zd\n", N);

    if(mpz_probab_prime_p(N, 50) > 0) {
        printf("N est premier\n\n");
        gmp_printf("%Zd 1", N);
        mpz_set_ui(cof, 1);
    }
    else {
        printf("\n--- TRIAL DIVISION ---\n\n");
        
        trialDivision(factors, &nf, cof, N, pow(10, 7), MAX_NUM_OF_FACTORS, file);
        printf("Factors found : %d\n", nf);
        for(int i = 0; i < nf; i++)
            gmp_printf("%Zd %d ", factors[i].f, factors[i].e);
        if(mpz_probab_prime_p(cof, 50) > 0) {
            gmp_printf("%Zd 1", cof);
			AddFactor(factors + nf, cof, 1, FACTOR_IS_PRIME);
            nf++;
            mpz_set_ui(cof, 1);
        }
    }

    printf("\n");
    count = nf;
    nf = 0;
    mpz_set(N, cof);
    rewind(file);

    if(mpz_cmp_ui(cof, 1) != 0){
        gmp_printf("Trial division not sufficient\nRemainder : %Zd\n\n", cof);
        printf("--- p - 1 ---\n\n");
        long B1 = 10000000, B2 = 1000000000;
        
        // if(argc > 3) { B1 = atoi(args[2]); B2 = atoi(args[3]); }

        while(PollardPminus1(factors + count, &nf, cof, N, B1, B2, file) == FACTOR_FOUND) {
            gmp_printf("%Zd %d ", factors[count].f, factors[count].e);
            mpz_div(N, N, factors[count].f);
            if(mpz_cmp_ui(N, 1) == 0)
                break;
            if(mpz_probab_prime_p(N, 50) > 0) {
                gmp_printf("%Zd 1\n", N);
                mpz_set_ui(N, 1);
                break;
            }
            nf = 0;
            count++;
        }
        nf = count;

        if(mpz_cmp_ui(N, 1) != 0) {
            gmp_printf("p-1 method not sufficient\nRemainder : %Zd\n\n", N);
            printf("--- Rho ---\n\n");
            long n = LONG_MAX;

            // if(argc > 4)
            //     n = atoi(args[4]);
            while(PollardRho(factors + count, &nf, cof, N, goodFunction, n) == FACTOR_FOUND) {
                gmp_printf("%Zd %d ", factors[count].f, factors[count].e);
                mpz_div(N, N, factors[count].f);
                if(mpz_cmp_ui(N, 1) == 0)
                    break;
                if(mpz_probab_prime_p(N, 50) > 0) {
                    gmp_printf("%Zd 1\n", N);
                    mpz_set_ui(N, 1);
                    break;
                }
                nf = 0;
                count++;
            }
        }
    }

    finish = clock();
	printf("---------------------------------------\nTime : %fs\n\n", (double)(finish - start)/CLOCKS_PER_SEC);

    if(mpz_cmp_ui(N, 1) != 0)
        gmp_printf("No more factor found\nRemainder : %Zd\n", N);

    *nbFactors = count;

    mpz_clear(cof);
    fclose(file);
}
