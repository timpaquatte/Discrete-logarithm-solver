#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "solve_dlog.h"
#include "factor.h"
#include "crt.h"


int f(mpz_t res, mpz_t x, const mpz_t p, const mpz_t a, const mpz_t g)
{
	mpz_t third;

	mpz_init(third);

	mpz_div_ui(third, p, 3);
	if(mpz_cmp(x, third) < 0) {
		mpz_mul(res, a, x);
		mpz_mod(res, res, p);
		mpz_clear(third);
		return 1;
	}

	mpz_mul_ui(third, third, 2);
	if(mpz_cmp(x, third) < 0) {
		mpz_mul(res, x, x);
		mpz_mod(res, res, p);
		mpz_clear(third);
		return 2;
	}

	mpz_mul(res, x, g);
	mpz_mod(res, res, p);

	mpz_clear(third);
	return 3;
}

void check_invariant(const mpz_t ux, const mpz_t vx, const mpz_t uy, const mpz_t vy, const mpz_t g, const mpz_t a, const mpz_t N, const mpz_t x, const mpz_t y)
{
	mpz_t left, right, temp;

	mpz_inits(left, right, temp, NULL);

	mpz_powm(left, a, ux, N);
	mpz_powm(temp, g, vx, N);
	mpz_mul(left, left, temp);
	mpz_mod(left, left, N);
	if(mpz_cmp(x, left) != 0) {
		gmp_printf("FAAAAIL: x = %Zd\n", x);
		exit(0);
	}

	mpz_powm(right, a, uy, N);
	mpz_powm(temp, g, vy, N);
	mpz_mul(right, right, temp);
	mpz_mod(right, right, N);
	if(mpz_cmp(y, right) != 0) {
		gmp_printf("FAAAAIL: y = %Zd\n", y);
		exit(0);
	}
	mpz_clears(left, right, temp, NULL);
}

void pohligHellman(mpz_t result, const mpz_t g, const mpz_t a, const mpz_t N, mpz_t factors[], int nbFactors)
{
	mpz_t gi, ai, temp;
	mpz_t* dlogs = malloc(nbFactors * sizeof(mpz_t));
	pthread_t *pool = malloc(nbFactors * sizeof(pthread_t));

	mpz_inits(gi, ai, temp, NULL);
	for(int i = 0; i < nbFactors; i++)
		mpz_init(dlogs[i]);

	for(int i = 0; i < nbFactors; i++) {
		mpz_div(temp, N, factors[i]);
		mpz_powm(gi, g, temp, N);
		mpz_powm(ai, a, temp, N);

		PollardRhoArgs args = {
			.result = dlogs[i],
			.factor = factors[i],
			.fct = &f
		};
		mpz_inits(args.result, args.g, args.a, args.N, args.factor, NULL);
		mpz_set(args.g, gi);
		mpz_set(args.a, ai);
		mpz_set(args.N, N);
		mpz_set(args.factor, factors[i]);
		gmp_printf("N = %Zd\na = %Zd\ng = %Zd\n", args.N, args.a, args.g);
		printf("---------------------\n");

		pthread_create(&pool[i], NULL, &PollardRhoDLog, (void*) &args);
		mpz_set(dlogs[i], args.result);
		mpz_clears(args.result, args.g, args.a, args.N, args.factor, NULL);
	}

	for(int i = 0; i < nbFactors; i++)
		pthread_join(pool[i], NULL);

	CRT(result, dlogs, factors, nbFactors);

	mpz_clears(gi, ai, temp, NULL);
	for(int i = 0; i < nbFactors; i++)
		mpz_clear(dlogs[i]);
	free(dlogs);
}

void* PollardRhoDLog(void *args)
{
	mpz_t x0, x, y, temp, N_;
	mpz_t ux, uy, vx, vy;
	mpz_t diff, gcd;
	int ret, cpt = 0;
	PollardRhoArgs *prArgs = (PollardRhoArgs*) args;

	mpz_inits(x0, x, y, temp, N_, NULL);
	mpz_inits(ux, uy, vx, vy, NULL);
	mpz_inits(diff, gcd, NULL);
	mpz_set(x0, prArgs->a);

	do {
		mpz_set(x, x0);
		mpz_set(y, x0);
		mpz_set_ui(ux, 1);
		mpz_set_ui(uy, 1);
		mpz_set_ui(vx, cpt);
		mpz_set_ui(vy, cpt);
		mpz_sub_ui(N_, prArgs->N, 1);

		do {
			// x = f(x)
			// gmp_printf("Avant: prArgs->fct=%p\n", (prArgs->fct));
			// gmp_printf("N = %Zd\na = %Zd\ng = %Zd\n", prArgs->N, prArgs->a, prArgs->g);
			ret = prArgs->fct(temp, x, prArgs->N, prArgs->a, prArgs->g);
			mpz_set(x, temp);
			// update ux, vx
			switch (ret)
			{
				case 1:
					mpz_add_ui(ux, ux, 1);
					break;
				case 2:
					mpz_mul_ui(ux, ux, 2);
					mpz_mul_ui(vx, vx, 2);
					break;
				case 3:
					mpz_add_ui(vx, vx, 1);
					break;
			}

			// y = f(f(y))
			for(int i = 0; i < 2; i++) {
				ret = prArgs->fct(temp, y, prArgs->N, prArgs->a, prArgs->g);
				mpz_set(y, temp);
				switch (ret)
				{
					case 1:
						mpz_add_ui(uy, uy, 1);
						break;
					case 2:
						mpz_mul_ui(uy, uy, 2);
						mpz_mul_ui(vy, vy, 2);
						break;
					case 3:
						mpz_add_ui(vy, vy, 1);
						break;
				}
			}

			mpz_mod(ux, ux, N_);
			mpz_mod(vx, vx, N_);
			mpz_mod(uy, uy, N_);
			mpz_mod(vy, vy, N_);
			// check_invariant(ux, vx, uy, vy, g, a, N, x, y);
		} while(mpz_cmp(x, y) != 0);

		// Now that we've found ux,uy,vx,vy, we have to be sure that they give a solution
		mpz_sub(diff, ux, uy);
		mpz_mod(diff, diff, prArgs->factor);
		mpz_gcd(gcd, diff, prArgs->factor);

		cpt++;
		mpz_mul(x0, x0, prArgs->g);
	} while(mpz_cmp_ui(gcd, 1) != 0);

	mpz_invert(prArgs->result, diff, prArgs->factor);
	mpz_sub(diff, vy, vx);
	mpz_mod(diff, diff, prArgs->factor);
	mpz_mul(prArgs->result, prArgs->result, diff);
	mpz_mod(prArgs->result, prArgs->result, prArgs->factor);

	mpz_clears(x0, x, y, temp, N_, NULL);
	mpz_clears(ux, uy, vx, vy, NULL);
	mpz_clears(diff, gcd, NULL);

	// Print results
	gmp_printf("----------------------\n\ngi = %Zd\nai = %Zd\nfactor = %Zd\nResult: %Zd\n\n", prArgs->g, prArgs->a, prArgs->factor, prArgs->result);
}

void usage(char *s){
    fprintf(stderr, "Usage: %s g a N\nTo solve g^x = a [N]\n", s);
}

int main(int argc, char *argv[]){
    if(argc == 1){
		usage(argv[0]);
		return 0;
    }

	mpz_t a, g, N, order, result;
    clock_t finish, start;

	mpz_inits(a, g, N, order, result, NULL);
	mpz_set_str(g, argv[1], 10);
	mpz_set_str(a, argv[2], 10);
	mpz_set_str(N, argv[3], 10);
	mpz_sub_ui(order, N, 1);

	printf("------------ FACTORISATION ------------\n\n");
	start = clock();

	int nbFactors = 0;
	factor_t factors[MAX_NUM_OF_FACTORS];
	mpz_t factors_mult[MAX_NUM_OF_FACTORS];
	factorize(factors, order, &nbFactors);

	printf("------------ POHLIG-HELLMAN ------------\n\n");

	for(int i = 0; i < nbFactors; i++) {
		mpz_init(factors_mult[i]);
		mpz_pow_ui(factors_mult[i], factors[i].f, factors[i].e);
	}
	// char* factors_str[] = {"2", "38579489651", "173", "194264244901", "5030919566507", "5481173216993", "7819807277899", "212521", "237169", "241081", "48970736047507", "528529", "423564751", "454756609", "549057842453207", "613813499121307", "688768866241"};
	// for(int i = 0; i < 1; i++)
	// 	mpz_init_set_str(factors[i], factors_str[i], 10);

	printf("Solving g^x = a [N]\n");
	gmp_printf("g = %Zd\na = %Zd\nN = %Zd\n\n", g, a, N);

	pohligHellman(result, g, a, N, factors_mult, nbFactors);

	gmp_printf("\nFinal result: %Zd\n", result);
	finish = clock();
	printf("---------------------------------------\nTime : %fs\n\n", (double)(finish - start)/CLOCKS_PER_SEC);

	mpz_clears(a, g, N, order, result, NULL);
	for(int i = 0; i < nbFactors; i++)
		mpz_clear(factors_mult[i]);
	return 0;
}