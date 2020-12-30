int PollardRhoSteps(mpz_t f, const mpz_t N,
					void (*fct)(mpz_t, mpz_t, const mpz_t), long *iter,
					long nrOfIterations, mpz_t x, mpz_t y);
int PollardRho(factor_t* result, int *nf, mpz_t cof, const mpz_t N,
	       void (*f)(mpz_t, mpz_t, const mpz_t),
	       long nrOfIterations);
