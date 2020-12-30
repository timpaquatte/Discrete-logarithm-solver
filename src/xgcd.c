#include <stdio.h>
#include <assert.h>

#include "gmp.h"

#include "xgcd.h"

#define DEBUG 0



/* compute g, u and v s.t. a*u+b*v = g = gcd(a, b) */
void XGCD(mpz_t g, mpz_t u, mpz_t v, mpz_t a, mpz_t b)
{
    mpz_t rp, up, vp, q, gcopy, ucopy, vcopy, prod;

    mpz_init(q);
    mpz_init(gcopy);
    mpz_init(ucopy);
    mpz_init(vcopy);
    mpz_init(prod);
    mpz_init_set(rp, b);
    mpz_init_set_ui(up, 0);
    mpz_init_set_ui(vp, 1);
    mpz_set(g, a);
    mpz_set_ui(u, 1);
    mpz_set_ui(v, 0);

    while(mpz_cmp_ui(rp, 0)) {
        mpz_fdiv_q(q, g, rp);

        mpz_set(gcopy, g);
        mpz_set(ucopy, u);
        mpz_set(vcopy, v);

        mpz_set(g, rp);
        mpz_set(u, up);
        mpz_set(v, vp);

        mpz_mul(prod, q, rp);
        mpz_sub(rp, gcopy, prod);

        mpz_mul(prod, q, up);
        mpz_sub(up, ucopy, prod);   

        mpz_mul(prod, q, vp);
        mpz_sub(vp, vcopy, prod);
    }

    mpz_clear(rp);
    mpz_clear(up);
    mpz_clear(vp);
    mpz_clear(q);
    mpz_clear(gcopy);
    mpz_clear(ucopy);
    mpz_clear(vcopy);
    mpz_clear(prod);
}


/* Solve a*x=b mod m if possible. 
   
*/
int linear_equation_mod(mpz_t x, mpz_t a, mpz_t b, mpz_t m){
    mpz_t u, v, g, gcd, r;

    mpz_init(u);
    mpz_init(v);
    mpz_init(g);
    mpz_init(gcd);
    mpz_init(r);

    XGCD(gcd, u, v, a, m);
    mpz_fdiv_r(r, b, gcd);
    
    if(mpz_cmp_ui(r, 0) != 0)
        return -1;

    mpz_div(a, a, gcd);
    mpz_div(b, b, gcd);
    mpz_div(m, m, gcd);

    XGCD(gcd, u, v, a, m);

    mpz_mul(x, u, b);
    mpz_fdiv_r(x, x, m);

    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(g);
    mpz_clear(gcd);
    mpz_clear(r);
    return 1;
}
 
