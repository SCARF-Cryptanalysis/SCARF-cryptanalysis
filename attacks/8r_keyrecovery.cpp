/*
 Proof-of-concept key recovery
 Starts by guessing the full K[0], and count collisions to check that
 the approach is valid.  Recovers partial K[1].
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

#include "scarf.hpp"
//#include "scarf_c.h"

/* Some useful SCARF components */

// 2 rounds partial encryption
inline int enc_partial(int input, uint64_t RK) {
    int ct = input;
    ct = R1(ct,  RK     &((1<<30)-1));
    ct = R1(ct, (RK>>30)&((1<<30)-1));
    return ct;
}

int sbox_inv[1 << 5];

int R1_inv(int input, uint64_t key)
{
    int SK1 = (key & 0x3E000000) >> 25;

    int xR = extract_bits<0, 5>(input);
    int xL = extract_bits<5, 10>(input);

    int temp = sbox_inv[xR] ^ SK1;

    xR = xL ^ Gfunc(temp, key);
    xL = temp;

    return (xL << 5) ^ xR;
}

// 2 rounds partial encryption
inline int dec_partial(int input, uint64_t RK) {
    int pt = input;
    pt = R1_inv(pt, (RK>>30)&((1<<30)-1));
    pt = R1_inv(pt,  RK     &((1<<30)-1));
    return pt;
}

// Linear part of G
int Glin(int input, uint64_t key)
{
    int x0 = input & extract_bits<0, 5>(key);
    int x1 = rol<1, N / 2>(input) & extract_bits<5, 10>(key);
    int x2 = rol<2, N / 2>(input) & extract_bits<10, 15>(key);
    int x3 = rol<3, N / 2>(input) & extract_bits<15, 20>(key);
    int x4 = rol<4, N / 2>(input) & extract_bits<20, 25>(key);

    int output = x0 ^ x1 ^ x2 ^ x3 ^ x4;

    return output;
}


uint64_t double_sbox_inv[1 << N];

inline uint64_t Slayer_inv(uint64_t input)
{
    return double_sbox_inv[input & 0x3ff] ^ (double_sbox_inv[(input >> N) & 0x3ff] << N) ^ (double_sbox_inv[(input >> (2 * N)) & 0x3ff] << (2 * N)) ^ (double_sbox_inv[(input >> (3 * N)) & 0x3ff] << (3 * N)) ^ (double_sbox_inv[(input >> (4 * N)) & 0x3ff] << (4 * N)) ^ (double_sbox_inv[(input >> (5 * N)) & 0x3ff] << (5 * N));
}


void print_key(uint64_t k) {
    for (int i=5; i>=0; i--) {
	printf ("%02x", (int)((k>>(5*i))&0x1f));
	if (i) {
	    printf (".");
	}
    }
}

/* Definition of the trail */

#define KEY_GUESS_EQ2     0b1111100000111110000000000LL

#define KEY_MASK_ACTIVE   0b111111111100000000001111111111111110000011111111110000011111LL
#define DELTA_IN          0b110011111000000000001011010010100100000001000110110000010110LL
#define DELTA_OUT                                            0b0100100000010100000000000LL
#define DELTA2                                               0b1000000000001000000000000LL
#define GUESS2                                               0b1111100000111110000000000LL
#define GUESS3            0b000000000011111000001111100000000000000000000000000000011111LL
#define GUESS4            0b000000000001000000000001000000000001000000000001000000000001LL
#define DELTA_FIN         0b01110LL


#define TWEAK_BITS        0b011110111101111011110111101111011110111101111011110111101111LL

// 32-bit tweak mask
#define TWEAK_MASK        (KEY_MASK_ACTIVE&TWEAK_BITS)

/* Definition of tweak space: dimension-29 or dimension-30 with at most 32 active bits */

// tweak space of dimension 29
#define TWEAK_DIM 29
#define TWEAK_BASE        0b011110111100000000000111100111001110000001111011110000000111LL
std::array<uint64_t, 3> TWEAK_PARITIES = {
			  0b000000000000000000000000000000000000000000000000000000001111LL,
			  0b000000000000000000000000000000011110000000000000000000000000LL,
			  0b000000000000000000000000001111000000000000000000000000000000LL,
};
#define THRESHOLD (1.4/1024)
#define THRESHOLD2 (2.5/1024)

// // tweak space of dimension 30
// #define TWEAK_DIM 30
// #define TWEAK_BASE        0b011110111100000000000111100111001110000001111011110000001111LL
// std::array<uint64_t, 3> TWEAK_PARITIES = {
// 			  0b000000000000000000000000000000011110000000000000000000000000LL,
// 			  0b000000000000000000000000001111000000000000000000000000000000LL,
// };
// #define THRESHOLD 1080./884736
// #define THRESHOLD2 (1./1024.+49./32768.)


// // tweak space of dimension 32
// #define TWEAK_DIM 32
// #define TWEAK_BASE        0b011110111100000000000111101111011110000001111011110000001111LL
// std::array<uint64_t, 3> TWEAK_PARITIES = {
// };
// #define THRESHOLD (1./1024.+11./32768.)
// #define THRESHOLD2 (1./1024.+49./32768.)


// parity equations for tweak space
int tweak_is_valid(uint64_t t) {
    if (t & (~TWEAK_MASK))
	return 0;
    for (unsigned i=0; i<TWEAK_PARITIES.size(); i++) {
	if (__builtin_popcountll(t&TWEAK_PARITIES[i]) % 2)
	    return 0;
    }
    return 1;
}

// expansion from TWEAK_DIM bits to 60
uint64_t tweak_expand(uint64_t z) {
    uint64_t t = _pdep_u64(z, TWEAK_BASE);
    for (unsigned i=0; i<TWEAK_PARITIES.size(); i++) {
	if (__builtin_popcountll(t&TWEAK_PARITIES[i]) % 2)
	    t ^= TWEAK_PARITIES[i]&(~TWEAK_BASE);
    }
    assert(tweak_is_valid(t));
    return t;
}

// #define RANDOM_KEY 1
// #define DELTA_KEY (1<<5)

int main(int argc, char* argv[]){

    for (int i=0; i<1<<N; i++) {
	double_sbox_inv[double_sbox[i]] = i;
    }
    for (int i=0; i<1<<5; i++) {
	sbox_inv[sbox[i]] = i;
    }
    
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

    // Guess K0
    uint64_t guess = key[0];

#ifdef DELTA_KEY
    guess ^= DELTA_KEY;
#elif defined(RANDOM_KEY)
    getrandom(&guess , sizeof(guess ), 0);
    guess  &= (1ULL<<60)-1;
#endif
    
    printf ("Guess: %015llx\n", (unsigned long long)guess);

    /*
      Precompute pairs of values for the second layer transitions DELTA_OUT -> DELTA2
      4 values corresponding to pairs val2[i],val2[i]+DELTA_OUT
    */
    uint64_t val2[4];
    {
	int i=0;
	for (uint64_t z=0; z < 1<<10; z++) {
	    uint64_t x = _pdep_u64(z, GUESS2);
	    uint64_t y = x ^ DELTA_OUT;
	    if ((Slayer(x) ^ Slayer(y)) == DELTA2) {
		val2[i++] = x;
	    }
	}
	assert(i == 4);
    }

    uint64_t pairs[1<<10] = {0};
    uint64_t colls[1<<10] = {0};
    // Count pairs and collision for key candidates K[1]
    
    /* Iterate over pairs of tweak.
       For each tweak T, build T' so that the difference after S-Box
       layer is equal to DELTA_IN
     */

    int tw = 0;
    
#pragma omp parallel for reduction(+:tw)
    for (uint64_t z = 0; z < 1ULL<<TWEAK_DIM; z++) {
        uint64_t t  = tweak_expand(z);
        uint64_t t0 = _pext_u64(t, TWEAK_BITS);
        // z  is a 32-bit value (only active tweak nibbles)
        // t  is a 60-bit value (key state)
        // t0 is a 48-bit value (full tweak)
	    
        t = Slayer(t ^ (guess&KEY_MASK_ACTIVE));
        t = t ^ DELTA_IN;
        t = Slayer_inv(t) ^ (guess&KEY_MASK_ACTIVE);

        if (tweak_is_valid(t)) {
            uint64_t t1 = _pext_u64(t,TWEAK_BITS);
            // (t0, t1) is a tweak pair leading to difference DELTA_IN
            tw++;
		
            /* Recover K[1] candidates */
            uint64_t guess2[4];
		
            uint64_t w0 = _pdep_u64(t0,TWEAK_BITS);
            w0 ^= guess;
            w0  = Slayer(w0);
            w0  = Sigma(w0);
            for (int i=0; i<4; i++) {
                guess2[i] = (val2[i]^w0)&GUESS2;
            }

            /* Build pairs with collision after 3 rounds */
            int mypairs = 0;
            int mycolls = 0;

            // Use full codebook
            for (int p0=0; p0 < 1<<10; p0++) {
                // Encrypt 2 rounds
                int x0 = enc_partial(p0, guess^_pdep_u64(t0, TWEAK_BITS));
                // Apply diff
                int x1 = x0 ^ Glin(x0>>5, DELTA_OUT);
                // Decrypt 2 rounds
                int p1 = dec_partial(x1, guess^_pdep_u64(t1, TWEAK_BITS));
                mypairs++;
                int c0 = enc(p0, key[3], key[2], key[1], key[0], t0);
                int c1 = enc(p1, key[3], key[2], key[1], key[0], t1);
                if (c0 == c1) {
                    mycolls++;
                }
            }
            // Increment shared counters
            for (int i=0; i<4; i++) {
                __atomic_add_fetch(&colls[_pext_u64(guess2[i],GUESS2)], mycolls, __ATOMIC_RELAXED);
                __atomic_add_fetch(&pairs[_pext_u64(guess2[i],GUESS2)], mypairs, __ATOMIC_RELAXED);
            }
        }
    }

    printf ("Number of tweak pairs passing first layer: %i\n", tw);
	

    {
        /* Print candidate keys passing threshold */
	int candidates = 0;
	int found = 0;
	
	printf ("Actual K1   : ");
        print_key(key[1] & GUESS2);
        printf (" [%lli/%lli]\n", (long long) colls[_pext_u64(key[1], GUESS2)], (long long) pairs[_pext_u64(key[1], GUESS2)]);
        double max = 0;
        for (int i=0; i < 1<<10; i++) {
            double p = 1.*colls[i]/pairs[i];
            if (p > max)
                max = p;
            if (p > THRESHOLD) {
                printf ("K1 candidate: ");
                print_key(_pdep_u64(i, GUESS2));
                printf(" [%lli/%lli]\n", (long long) colls[i], (long long) pairs[i]);
		candidates++;
		if (_pdep_u64(i, GUESS2) == (key[1] & GUESS2))
		    found = 1;
            }
        }
        printf ("Max ratio: %f*2^-10\n", max*1024);
	printf ("PHASE 1   SUMMARY: %i candidates, key %s\n", candidates, found? "FOUND": "NOT FOUND");
    }
        
    printf("*** Phase 1.5\n");

    int candidates = 0;
    int found = 0;
    
    for (int i=0; i < 1<<10; i++) {
	double p = 1.*colls[i]/pairs[i];
	if (p > THRESHOLD) {

	    printf ("K1 candidate: ");
	    print_key(_pdep_u64(i, GUESS2));
	    printf(":");
		
	    uint64_t guess2 = key[1] & GUESS2;

	    uint64_t (*pairs2)[1<<5] = (uint64_t (*)[1<<5]) calloc(sizeof(uint64_t), 1<<20);
	    uint64_t (*colls2)[1<<5] = (uint64_t (*)[1<<5]) calloc(sizeof(uint64_t), 1<<20);
	    // Count pairs and collision for key candidates K[2] K[3]
    
	    /* Iterate over pairs of tweak.
	       For each tweak T, build T' so that the difference after S-Box
	       layer is equal to DELTA_IN
	    */

	    tw = 0;

#pragma omp parallel for reduction(+:tw)
	    for (uint64_t z = 0; z < 1ULL<<TWEAK_DIM; z++) {
		uint64_t t  = tweak_expand(z);
		uint64_t t0 = _pext_u64(t, TWEAK_BITS);
		// z  is a 32-bit value (only active tweak nibbles)
		// t  is a 60-bit value (key state)
		// t0 is a 48-bit value (full tweak)
	    
		t = Slayer(t ^ (guess&KEY_MASK_ACTIVE));
		t = t ^ DELTA_IN;
		t = Slayer_inv(t) ^ (guess&KEY_MASK_ACTIVE);

		if (tweak_is_valid(t)) {
		    uint64_t t1 = _pext_u64(t,TWEAK_BITS);
		    // (t0, t1) is a tweak pair leading to difference DELTA_IN
		    // Test with K[1] guess
		    uint64_t w0 = _pdep_u64(t0,TWEAK_BITS);
		    uint64_t z0;
		    w0 ^= guess;
		    w0  = Slayer(w0);
		    w0  = Sigma(w0);
		    w0 ^= guess2;
		    z0  = Slayer(w0);
                
		    uint64_t w1 = _pdep_u64(t1,TWEAK_BITS);
		    uint64_t z1;
		    w1 ^= guess;
		    w1  = Slayer(w1);
		    w1  = Sigma(w1);
		    w1 ^= guess2;
		    z1  = Slayer(w1);

		    if ((z0 ^ z1) == DELTA2) {

			tw++;
                    
			/* Build pairs with collision after 3 rounds */
			int mycolls = 0;
                
			// Use full codebook and count collisions
			for (int p0=0; p0 < 1<<10; p0++) {
			    // Encrypt 2 rounds
			    int x0 = enc_partial(p0, guess^_pdep_u64(t0, TWEAK_BITS));
			    // Apply diff
			    int x1 = x0 ^ Glin(x0>>5, DELTA_OUT);
			    // Decrypt 2 rounds
			    int p1 = dec_partial(x1, guess^_pdep_u64(t1, TWEAK_BITS));
			    int c0 = enc(p0, key[3], key[2], key[1], key[0], t0);
			    int c1 = enc(p1, key[3], key[2], key[1], key[0], t1);
			    if (c0 == c1) {
				mycolls++;
			    }
			}
			/* Identify key candidates making the tweak pair valid */
			for (int k2 = 0; k2 < 1<<15; k2++) {
			    for (int k3 = 0; k3 < 1<<5; k3++) {
				z0 = w0 ^ _pdep_u64(k2, GUESS3);
				z1 = w1 ^ _pdep_u64(k2, GUESS3);
				z0 = Slayer(z0);
				z1 = Slayer(z1);
				z0 = z0 ^ _pdep_u64(k3, GUESS4);
				z1 = z1 ^ _pdep_u64(k3, GUESS4);
				z0 = permutation(z0);
				z1 = permutation(z1);
				z0 = Slayer(z0);
				z1 = Slayer(z1);

				if (((z0 ^ z1) & ~DELTA_FIN) == 0) {
				    __atomic_add_fetch(&pairs2[k2][k3], 1<<10, __ATOMIC_RELAXED);
				    __atomic_add_fetch(&colls2[k2][k3], mycolls, __ATOMIC_RELAXED);
				}
			    }
			}
		    }
		}
	    }
        
	    {
		double max = 0;
		double sump = 0;
		double sumc = 0;
		int passed = 0;
		for (int i=0; i < 1<<15; i++) {
		    for (int j=0; j < 1<<5; j++) {
			sump += pairs2[i][j];
			sumc += colls2[i][j];
			double p = 1.*colls2[i][j]/pairs2[i][j];
			if (p > max)
			    max = p;
			if (p > THRESHOLD2) {
			    passed = 1;
			}
		    }
		}
                
                printf (" max ratio: %f*2^-10%s\n", max*1024, passed? " - PASSED": "");

		if (passed) {
		    candidates++;
		    if (_pdep_u64(i, GUESS2) == (key[1] & GUESS2))
			found = 1;
		}
	    }
	}
    }

    printf ("PHASE 1.5 SUMMARY: %i candidates, key %s\n", candidates, found? "FOUND": "NOT FOUND");
}
