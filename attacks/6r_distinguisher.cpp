/*
 Proof-of-concept distinguisher on 6 rounds
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <array>

#include <immintrin.h>

// Getrandom
# if __GLIBC__ > 2 || __GLIBC_MINOR__ > 24

#include <sys/random.h>

#else /* older glibc */

#include <unistd.h>
#include <sys/syscall.h>
#include <errno.h>

int getrandom(void *buf, size_t buflen, unsigned int flags) {
    return syscall(SYS_getrandom, buf, buflen, flags);
}
# endif

/* SCARF implementation */
// #include "scarf.hpp"
#include "scarf_c.h"

/* AES-based PRNG */
#include "aes-rng.h"

/* Reduced-round SCARF */
#define REDUCED_ROUNDS 6
inline int enc_partial(int input, uint64_t key3, uint64_t key2, uint64_t key1, uint64_t key0, uint64_t tweak) {
    uint64_t RK[ROUNDS];
    tweakey_schedule(key3, key2, key1, key0, tweak, &RK[0]);

    int ct = input;

    for (int i = 0; i < REDUCED_ROUNDS - 1; i++)
    {
        ct = R1(ct, RK[i]);
    }

    ct = R2(ct, RK[REDUCED_ROUNDS - 1]);

    return ct;
}

/* Definition of the trail */

#define DELTA_IN   (1<<3)
#define NB_SAMPLES (1ULL<<30)
#define THRESHOLD (1.002/1024)

#define REPEAT_CNT 20

int main(int argc, char* argv[]){

    uint64_t key[4];

    if (argc == 5) {
	key[0] = strtoull(argv[1], NULL, 16);
	key[1] = strtoull(argv[2], NULL, 16);
	key[2] = strtoull(argv[3], NULL, 16);
	key[3] = strtoull(argv[4], NULL, 16);
    } else {
	getrandom(&key, sizeof(key), 0);
    }
    // Clean key
    for (int i=0; i<4; i++) {
	key[i] &= (1ULL<<60)-1;
    }
    
    printf ("Key  : %015llx %015llx %015llx %015llx\n",
	    (unsigned long long) key[0], (unsigned long long) key[1],
	    (unsigned long long) key[2], (unsigned long long) key[3]);


    printf ("Running experiments with 6-round SCARF\n");

    int success = 0;
    int collisions = 0;
    
#pragma omp parallel for reduction(+:success,collisions)
    for (uint64_t r = 0; r < REPEAT_CNT; r++) {
        uint32_t seed[2];
        getrandom(&seed, sizeof(seed), 0);
        struct RNG_state* st = init_aesrand_r(seed[0], seed[1]);

	int colls = 0;
	for (uint64_t i=0; i<NB_SAMPLES; i++) {
	    uint64_t tweak = aesrand_int64_r(st) & ((1ULL<<48)-1);
	    int c0 = enc_partial(0, key[3], key[2], key[1], key[0], tweak);
	    int c1 = enc_partial(0, key[3], key[2], key[1], key[0], tweak^DELTA_IN);
	    if (c0 == c1)
		colls++;
	}

	if (colls > THRESHOLD*NB_SAMPLES) {
	    success++;
	}
	#pragma omp critical
	{
	    printf (colls > THRESHOLD*NB_SAMPLES? ".": "!");
	    fflush(stdout);
	}
	collisions += colls;
    }
    
    printf ("\nSuccess: %i/%i\n", success, REPEAT_CNT);
    printf ("Average collision probability: %f*2^-10\n", 1024.*collisions/NB_SAMPLES/REPEAT_CNT);
    
    printf ("Running experiments with random data\n");

    success = 0;
    
#pragma omp parallel for reduction(+:success)
    for (uint64_t r = 0; r < REPEAT_CNT; r++) {
        uint32_t seed[2];
        getrandom(&seed, sizeof(seed), 0);
        struct RNG_state* st = init_aesrand_r(seed[0], seed[1]);

	int colls = 0;
	for (uint64_t i=0; i<NB_SAMPLES; i++) {
	    int c0 = aesrand_int32_r(st) & ((1ULL<<10)-1);
	    int c1 = aesrand_int32_r(st) & ((1ULL<<10)-1);
	    if (c0 == c1)
		colls++;
	}

	if (colls < THRESHOLD*NB_SAMPLES) {
	    success++;
	}
	#pragma omp critical
	{
	    printf (colls < THRESHOLD*NB_SAMPLES? ".": "!");
	    fflush(stdout);
	}
    }
    
    printf ("\nSuccess: %i/%i\n", success, REPEAT_CNT);


}
