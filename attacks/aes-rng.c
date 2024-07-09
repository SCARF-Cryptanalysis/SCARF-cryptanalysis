// AES-based RNG
#include "aes-rng.h"

struct RNG_state* init_aesrand_r(uint32_t seed1, uint32_t seed2) {
  // Internal state: count + t
  struct RNG_state *st = malloc(sizeof(struct RNG_state));
  assert(st);
  st->count = _mm_setr_epi32(seed1, seed2, 0, 0);
  st->k = _mm_setr_epi32(seed1+seed2, seed1*seed2, seed1, seed2);
  st->next = 0;
  return st;
}

uint32_t aesrand_int32_r(struct RNG_state* rng) {
  switch(rng->next++) {
  case 0:
    // Increment counter
    rng->count = _mm_add_epi64(rng->count, _mm_setr_epi32(0,0,1,0));

    // 6 AES rounds
    rng->t = rng->count;
    rng->t = _mm_aesenc_si128(rng->t, rng->k);
    rng->t = _mm_aesenc_si128(rng->t, rng->k);
    rng->t = _mm_aesenc_si128(rng->t, rng->k);
    rng->t = _mm_aesenc_si128(rng->t, rng->k);
    rng->t = _mm_aesenc_si128(rng->t, rng->k);
    rng->t = _mm_aesenc_si128(rng->t, rng->k);

    return _mm_extract_epi32(rng->t,0);
  case 1:
    return _mm_extract_epi32(rng->t,1);
  case 2:
    return _mm_extract_epi32(rng->t,2);
  case 3:
    rng->next = 0;
    return _mm_extract_epi32(rng->t,3);
  default:
    assert(0);
  }
}

uint64_t aesrand_int64_r(struct RNG_state* rng) {
  uint32_t a = aesrand_int32_r(rng);
  uint32_t b = aesrand_int32_r(rng);

  return (((uint64_t) a)<<32) + b;
}
