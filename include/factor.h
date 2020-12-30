#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include "time.h"
#include "limits.h"

#include "gmp.h"
#include "utils.h"

#define MAX_NUM_OF_FACTORS 20

void goodFunction(mpz_t out, mpz_t in, const mpz_t N);
void factorize(factor_t output[], mpz_t N, int* nbFactors);
