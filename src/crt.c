#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "gmp.h"

#include "xgcd.h"
#include "crt.h"

/* Given (r0, m0) and (r1, m1), compute n such that
   n mod m0 = r0; n mod m1 = r1.  If no such n exists, then this
   function returns 0. Else returns 1.  The moduli m must all be positive.
*/
int CRT2(mpz_t n, mpz_t r0, mpz_t m0, mpz_t r1, mpz_t m1){
    int status = 1;

    mpz_t e0, e1, u, v, g, prod1, prod2, ppcm;
    mpz_init(e0);
    mpz_init(e1);
    mpz_init(u);
    mpz_init(v);
    mpz_init(g);
    mpz_init(prod1);
    mpz_init(prod2);
    mpz_init(ppcm);

    XGCD(g, u, v, m0, m1);          // u*m0 + v*m1 = g

    // Check if the solution exists
    mpz_fdiv_r(e0, r0, g);
    mpz_fdiv_r(e1, r1, g);
    if(mpz_cmp(e0, e1) != 0)
        status = 0;

    else {
        // Computing the ppcm
        mpz_mul(ppcm, m0, m1);
        mpz_fdiv_q(ppcm, ppcm, g);
    
        mpz_mul(e0, v, m1);         // e0 = v*m1
        mpz_mul(e1, u, m0);         // e1 = u*m0

        mpz_mul(prod1, r0, e0); 
        mpz_mul(prod2, r1, e1);

        mpz_add(n, prod1, prod2);   // n = e0*r0 + e1*r1
        mpz_fdiv_q(n, n, g);        // n /= pgcd(m0, m1)
        mpz_fdiv_r(n, n, ppcm);     // mod ppcm(m0, m1)
    }
    
    mpz_clear(e0);
    mpz_clear(e1);
    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(g);
    mpz_clear(prod1);
    mpz_clear(prod2);
    mpz_clear(ppcm);
    return status;
}


/* Given a list S of pairs (r,m), returns an integer n such that n mod
   m = r for each (r,m) in S.  If no such n exists, then this function
   returns 0. Else returns 1.  The moduli m must all be positive.
*/
int CRT(mpz_t n, mpz_t *r, mpz_t *m, int nb_pairs){
    int status = 1;
    mpz_t ppcm, pgcd, u, v;

    mpz_init(ppcm);
    mpz_init(pgcd);
    mpz_init(u);
    mpz_init(v);

    CRT2(n, r[0], m[0], r[1], m[1]);
    XGCD(pgcd, u, v, m[0], m[1]);
    mpz_mul(ppcm, m[0], m[1]);
    mpz_fdiv_q(ppcm, ppcm, pgcd);

    for(int i = 2; i < nb_pairs; i++) {
        if(CRT2(n, n, ppcm, r[i], m[i]) == 0)
            break;
        XGCD(pgcd, u, v, ppcm, m[i]);
        mpz_mul(ppcm, ppcm, m[i]);
        mpz_fdiv_q(ppcm, ppcm, pgcd);
    }

    mpz_clear(ppcm);
    mpz_clear(pgcd);
    mpz_clear(u);
    mpz_clear(v);
    return status;
}
