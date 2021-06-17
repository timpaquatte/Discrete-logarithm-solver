#ifndef SOLVE_DLOG_H
#define SOLVE_DLOG_H

#include "gmp.h"

int f(mpz_t res, mpz_t x, const mpz_t p, const mpz_t a, const mpz_t g);
void check_invariant(const mpz_t ux, const mpz_t vx, const mpz_t uy, const mpz_t vy, const mpz_t g, const mpz_t a, const mpz_t N, const mpz_t x, const mpz_t y);
void pohligHellman(mpz_t result, const mpz_t g, const mpz_t a, const mpz_t N, mpz_t factors[], int nbFactors);
void* PollardRhoDLog(void *args);
void usage(char *s);


typedef struct {
	mpz_t result;
	const mpz_t g;
	const mpz_t a;
	const mpz_t N;
	const mpz_t factor;
	int (*fct)(mpz_t, mpz_t, const mpz_t, const mpz_t, const mpz_t);

} PollardRhoArgs;

#endif