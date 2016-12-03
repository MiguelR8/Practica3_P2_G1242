#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void toModM(mpz_t x, mpz_t m) {
	mpz_t q;
	mpz_t r;
	mpz_inits(q, r, NULL);
	
	mpz_tdiv_qr(q, r, x, m);
	
	if (mpz_sgn(x) < 0) {
		//|x| < m
		if (mpz_sgn(q) == 0) {
			//gmp_printf("(1) %Zd %% %Zd = %Zd + %Zd\n", x, m, x, m);
			mpz_add(x, x, m);
		} else {
			mpz_abs(q, q);
			if (mpz_sgn(r) != 0) {
				mpz_add_ui(q, q, 1L);
			}
			//gmp_printf(" (2) %Zd %% %Zd = %Zd + %Zd * %Zd\n", x, m, x, q, m);
			mpz_mul(q, m, q);
			mpz_add(x, x, q);
		}
	} else if (mpz_cmp(x, m) >= 0) {
		//gmp_printf(" (3) %Zd %% %Zd = %Zd - %Zd * %Zd\n", x, m, x, q, m);
		mpz_mul(q, m, q);
		mpz_sub(x, x, q);
	}
	//else x is already % m
	
	mpz_clears(q, r, NULL);
}

void exp_mod(mpz_t rop, mpz_t b, mpz_t e, mpz_t m) {
	mp_bitcnt_t l = mpz_sizeinbase(e, 2);	//longitud en bits
	mpz_set_ui(rop, 1L);
	
	do {
		l--;
		mpz_mul(rop, rop, rop);			//x = x^2
		toModM(rop, m);					//x = x^2 % m
		
		if (mpz_tstbit(e, l) == 1) {
			mpz_mul(rop, rop, b);		//x = x * b
			toModM(rop, m);				//x = x * b % m
		}
	} while(l != 0);
}

int main (int argc, char** argv) {
	
	mpz_t b;
	mpz_t e;
	mpz_t m;
	
	if (argc < 4) {
		printf("Uso: %s base exponente modulo\n", argv[0]);
		return EXIT_FAILURE;
	}
	mpz_inits(b, e, m, NULL);
	
	mpz_set_str(b,	   argv[1], 10);
	mpz_set_str(e, argv[2], 10);
	mpz_set_str(m,    argv[3], 10);
	
	if (mpz_sgn(e) < 0 || mpz_sgn(m) < 1) {
		printf("e debe ser mayor o igual que 0 y m mayor que 1\n");
		mpz_clears(b, e, m, NULL);
		return EXIT_FAILURE;
	}
	
	mpz_t res;
	mpz_init(res);
	
	exp_mod(res, b, e, m);
	
	gmp_printf("%Zd\n", res);
	mpz_clear(res);
	mpz_clears(b, e, m, NULL);
	return EXIT_SUCCESS;
}
