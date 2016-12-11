#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <gmp.h>
#include <time.h>
#include "../tables/primes.c"

#define ITERATION_MAX 32768

typedef enum {DEFINITELY_NOT = 0, PROBABLY, DEFINITELY} PRIME_CERTAINTY;
extern char *optarg;
extern int optind, opterr, optopt;

gmp_randstate_t state;
gmp_randstate_t state2;

void generate_prime_candidate_of_size(mpz_t rop, int size) {
	//genera aleatorio entre 0 y 2^(n-1)
	mpz_urandomb(rop, state, size);
	//asegura que el tamaño es el especificado
	mpz_setbit(rop, (mp_bitcnt_t)size-1);
	mpz_setbit(rop, 0);
}

void binary_factor_n(const mpz_t n, unsigned long* k, mpz_t m) {
	*k = mpz_scan1(n, 0L);		//get number of 0s until first 1
	mpz_tdiv_q_2exp(m, n, *k);	//m = floor(n/(2^k))
}

PRIME_CERTAINTY Miller_Rabin_Test(const mpz_t n, int reps) {
	int i, j;
	int r;
	mpz_t p;
	//Primer test, no es multiplo de los primeros 2000 primos
	mpz_init(p);
	//TODO: possible optimization, binary search in table
	//TODO: possible optimization, iterate from top until table element is smaller
	for (i = 0; i < PRIME_TABLE_SIZE; i++) {
		if (mpz_cmp_ui(n, primes[i]) == 0) {
			mpz_clears(p, NULL);
			return DEFINITELY;
		}
		
		mpz_tdiv_r_ui(p, n, primes[i]);
		if (mpz_sgn(p) == 0) {
			mpz_clears(p, NULL);
			return DEFINITELY_NOT;
		}
	}
	//no funciona
	mpz_t base;
	unsigned long k;
	mpz_t m;
	mpz_t rem;
	mpz_t n1;
	mpz_inits(base, m, n1, rem, NULL);
	
	mpz_sub_ui(n1, n, 1);
	//(n-1) = (2^k)*m
	binary_factor_n(n1, &k, m);
	
	mpz_sub_ui(p, n, 4);
	
	for (i = 0; i < reps; i++) {
		//encuentra base aleatoria
		mpz_urandomm(base, state2, p);
		//base <- [2, p-2]
		mpz_add_ui(base, base, 2);
		
		//rem = a^m % n
		mpz_powm(rem, base, m, n);
		
		//pase inmediato: a^m == 1 mod n
		if (mpz_cmpabs_ui(rem, 2) < 0) {
			continue;
		}
		mpz_add_ui(rem, rem, 1L);	//+1
		mpz_sub(rem, rem, n);		//-n
		//pase inmediato: a^m == -1 mod n
		if (mpz_sgn(rem) == 0) {
			continue;
		}
		//revierte si no pasa
		mpz_add(rem, rem, n);		//+n
		mpz_sub_ui(rem, rem, 1L);	//-1
		
		for (j = 1; j < k; j++) {
			mpz_powm_ui(rem, rem, 2L, n);
			//x != +-1 mod n, x == 1 mod n
			if (mpz_cmp_si(rem, 1L) == 0) {
				mpz_clears(base, m, n1, rem, NULL);
				mpz_clear(p);
				return DEFINITELY_NOT;
			}
			
			//rem == -1 mod n
			//rem + 1 == 0 mod n
			mpz_add_ui(rem, rem, 1L);	//+1
			mpz_sub(rem, rem, n);		//-n
			if (mpz_sgn(rem) == 0) {
				break;
			}
			mpz_add(rem, rem, n);		//+n
			mpz_sub_ui(rem, rem, 1L);	//-1
		}
		//no se encontro ningun -1
		if (j == k) {
			mpz_clears(base, m, n1, rem, NULL);
			mpz_clear(p);
			return DEFINITELY_NOT;
		}
		
	}
	mpz_clears(base, m, n1, rem, NULL);
	mpz_clear(p);
	
	return PROBABLY;
}

int main (int argc,char *argv[]) {
	
	if (argc < 5) {
		printf("Uso: %s {-b bits} {-t sec} [-o fileout]\n", argv[0]);
		return 0;
	}
	
	int c = 0;
	int bitlen = -1;
	int test_iters = -1;
	
	FILE* fout = NULL;
	
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
		   {"b", required_argument, 0, 'b'},
		   {"t", required_argument, 0, 't'},
		   {"o", required_argument, 0, 'o'},
		   {0, 0, 0, 0}
		};
		c = getopt_long(argc, argv, "+b:t:o:",
			long_options, &option_index);
		if (c < 0)
			break;
		
		switch (c) {
			case 'b':
				bitlen = atoi(optarg);
				break;
			case 't':
				test_iters = atoi(optarg);
				break;
			case 'o':
				fout = fopen(optarg, "wb");
				if (fout == NULL) {
					printf("Error al abrir %s para escribir\n", optarg);
					return EXIT_FAILURE;
				}
				break;
			default:
				if (fout != NULL) {
					fclose(fout);
				}
				printf("Uso: %s {-b bits} {-t sec} [-o fileout]\n", argv[0]);
				return EXIT_FAILURE;
		}
	}
	
	if (bitlen < 1 || test_iters < 1) {
		printf("Uso: %s {-b bits} {-t sec} [-o fileout]\n", argv[0]);
		printf("bits y sec deben ser mayores que 1\n");
		return -1;
	}
	
	long t = 1481383636L;//time(NULL);
	//printf("%ld\n", t);
	gmp_randinit_default(state);
	gmp_randseed_ui(state, t);
	gmp_randinit_default(state2);
	gmp_randseed_ui(state2, t);
	
	mpz_t n;
	mpz_init(n);
	
	int i;
	
	FILE* f = fopen("suc.txt", "w");
	if (!f) {
		perror(NULL);
		return EXIT_FAILURE;
	}
	
	
	for (i = 0; i < ITERATION_MAX; i++) {
		generate_prime_candidate_of_size(n, bitlen);
	
		//int ret = Miller_Rabin_Test(n, test_iters);
		int gmpret = mpz_probab_prime_p(n, test_iters);
		if (gmpret > DEFINITELY_NOT) {
			fprintf(f, "%d\n", i);
		}
	}
	
	fclose(f);
	
	/*gmp_printf("Para %Zd los resultados son:\n", n);
	printf("\tManual: ");
	if (ret == DEFINITELY_NOT) {
		printf("Compuesto\n");
	} else if (ret == PROBABLY) {
		printf("Probablemente\n");
	} else {
		printf("Definitivamente\n");
	}
	printf("\tGMP: ");
	if (gmpret == DEFINITELY_NOT) {
		printf("Compuesto\n");
	} else if (gmpret == PROBABLY) {
		printf("Probablemente\n");
	} else {
		printf("Definitivamente\n");
	}*/
	
	mpz_clear(n);
	gmp_randclear(state);
	gmp_randclear(state2);
	return 0;
}
