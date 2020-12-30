#include <stdlib.h>
#include <stdio.h>

#include "gmp.h"
#include "utils.h"
#include "trialdiv.h"

/* OUTPUT: 1 if factorization finished. */
int trialDivision(factor_t* factors, int *nf, mpz_t cof, const mpz_t N,
				  const long bound, uint length, FILE* file)
{
	char* line = malloc(100*sizeof(char));
	int e, offset = 0;
	mpz_t r, factor;

	mpz_init(r);
	mpz_init_set_ui(factor, 2);
	mpz_set(cof, N);

	
	
	while(mpz_cmp_ui(factor, bound) < 0) {
		e = 0;
		mpz_fdiv_r(r, cof, factor);
		while(mpz_cmp_ui(r, 0) == 0) {
			mpz_div(cof, cof, factor);
			e++;

			mpz_fdiv_r(r, cof, factor);
		}
		if(e > 0) {
			AddFactor(factors + offset, factor, e, FACTOR_IS_PRIME);
			offset++;
			if(offset == length)
				break;
		}

		fgets(line, 100, file);

		int diff = atoi(line);
		if(mpz_cmp_ui(factor, 2) != 0)
			diff *= 2;
		mpz_add_ui(factor, factor, diff);
	}

	mpz_clear(r);
	mpz_clear(factor);
	*nf = offset;
	return offset;
}
