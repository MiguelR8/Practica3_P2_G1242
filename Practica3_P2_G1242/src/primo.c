#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "../tables/primes.c"

typedef enum {DEFINITELY_NOT = 0, PROBABLY, DEFINITELY} PRIME_CERTAINTY;
extern char *optarg;
extern int optind, opterr, optopt;

gmp_randstate_t state;

mpz_t generate_prime_candidate_of_size(int size) {
	mpz_t n;
	//genera aleatorio entre 0 y 2^(n-1)
	mpz_urandomb(n, state, size);
	//asegura que el tamaño es el especificado
	mpz_setbit(n, (mp_bitcnt_t)size);
	mpz_setbit(n, 0);
	return n;
}

void binary_factor_n()

PRIME_CERTAINTY Miller_Rabin_Test(mpz_t n, int reps) {
	int i;
	mpz_t p;
	//Primer test, no es multiplo de los primeros 2000 primos
	mpz_init(p);
	for (i = 0; i < PRIME_TABLE_SIZE; i++) {
		mpz_set_ui(p, primes[i]);
		mpz_tdiv_r(p, n, p);
		if (mpz_sgn(p) == 0) {
			return DEFINITELY_NOT;
		}
	}
	
	
	mpz_t base;
	mpz_t k;
	mpz_t m;
	mpz_inits(base, k, m, NULL);
	
	mpz_set(p, n);
	//p <- [0, p-4]
	mpz_sub_ui(p, n, 4);
	
	for (i = 0; i < reps; i++) {
		//encuentra base aleatoria
		mpz_urandomm(base, state, p);
		//base <- [2, p-2]
		mpz_add_ui(base, base, 2);
		
		//factoriza p-1
	}
	mpz_clears(base, k, m, NULL);
	mpz_clear(p);
	
	return DEFINITELY;
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
	
	gmp_randinit_default(state);
	
	gmp_randclear(state);
	
	return 0;
}
