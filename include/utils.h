typedef unsigned int uint;

typedef struct{
    mpz_t f;
    int e, status;
} factor_t;

typedef struct{
    mpz_t N;
    int e, status;
    int next_method;
    int nargs_method;
    mpz_t *args_method;
} working_t;

#define FACTOR_IS_PRIME          1
#define FACTOR_IS_PROBABLE_PRIME 2
#define FACTOR_IS_COMPOSITE      3
#define FACTOR_IS_UNKNOWN        4

#define FACTOR_NOT_FOUND         0
#define FACTOR_FOUND             1
#define FACTOR_FINISHED          2
#define FACTOR_ERROR             3

#define PrimeStatus(N) \
	(mpz_probab_prime_p((N), 2) ? FACTOR_IS_PROBABLE_PRIME : \
	 FACTOR_IS_COMPOSITE)

extern int IsPerfectPower(mpz_t N1, mpz_t N);
extern void UpdateStatus(factor_t *fact);
extern void AddFactor(factor_t *fact, mpz_t f, int e, int status);
extern void AddSmallFactor(factor_t *fact, int f, int e, int status);
extern void PrintFactor(factor_t fact);
extern void PrintFactorization(factor_t *fact, int nf);
extern void factor_clear(factor_t* fact, int numberOfFactors);
