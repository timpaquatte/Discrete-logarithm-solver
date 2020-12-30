int PollardPminus1Step1(mpz_t f, const mpz_t N,	long bound1, FILE* file,
			mpz_t b, mpz_t p);
int PollardPminus1Step2(mpz_t f, const mpz_t N, long bound2, FILE* file,
			mpz_t b, mpz_t p);
int PollardPminus1(factor_t* result, int *nf, mpz_t cof, const mpz_t N,
		   long bound1, long bound2, FILE* file);

