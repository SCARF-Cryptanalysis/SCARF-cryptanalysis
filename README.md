# Cryptanalysis of Full SCARF

This repository contains supplementary meterial (code) for the paper:
> Cryptanalysis of Full SCARF, Antonio Flórez-Gutiérrez, Eran Lambooij, Gaëtan Leurent, Håvard Raddum, Tyge Tiessen, Michiel Verbauwhede, Eurocrypt 2025

## Attack simulation

The `attack` repository contains code to simulate the attacks presented
in the paper.  To compile the code, just run `make` (note: this requires
a compiler with support for C++20 and openmp).

- `6r_distinguisher` is a proof-of-concept implementation of the
  6-round distinguisher in Section 4 of the paper.
  
  The code counts the number of collisions for ciphertext pairs
  corresponding to the encryption of the zero plaintext with two tweaks
  with a specific input difference.  It first picks a random SCARF key
  and repeats 20 experiments with 2^30 samples each, then runs 20
  experiments with random ciphertexts.
  
  Running the code should take about 4 core-hours.
  
- `8r_keyrecovery` is a proof-of-concept implementation of the
  8-round key-recovery in Section 5 of the paper, as described in
  Section 5.5.

  The code guesses the full K1 and recovers a partial K2.  Following
  Section 5.3, we consider a space of 2^29 tweaks, and we filter a set
  of 2^8.2 tweak pairs compatible with the key.  For each compatible
  tweak pair, we consider 2^10 plaintext pairs.

  Then, the code implements the additional filtering of section 5.5 for
  each K2 candidate.

  The code can be modified to run with a random key by defining
  `RANDOM_KEY`, in order to verify that wrong keys do not pass the
  filtering.

  Running the code should take about 3 core-minutes.
	
	
## Counting Minimal Compatible Pairs

The `lower_bound` directory contains code to count the minimal number of
pairs follow a trail, using the theory of quasidifferentials, as
explained in Section 5.7
