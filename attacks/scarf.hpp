#include <cstdint>
#include <bit>
#include <immintrin.h>

// compile with the -DINTEL option to enable the intel specific optimizations
#ifndef AMD
    #ifndef INTEL
        #define AMD
    #endif
#endif

#define ROUNDS 8
#define N 10
#define TL 60
#define unit64 (uint64_t)1

int sbox[1 << (N / 2)] = {0, 2, 4, 12, 8, 14, 24, 21, 16, 19, 28, 5, 17, 20, 11, 23, 1, 6, 7, 26, 25, 18, 10, 27, 3, 13, 9, 29, 22, 30, 15, 31};
uint64_t double_sbox[1 << N] = {0, 2, 4, 12, 8, 14, 24, 21, 16, 19, 28, 5, 17, 20, 11, 23, 1, 6, 7, 26, 25, 18, 10, 27, 3, 13, 9, 29, 22, 30, 15, 31, 64, 66, 68, 76, 72, 78, 88, 85, 80, 83, 92, 69, 81, 84, 75, 87, 65, 70, 71, 90, 89, 82, 74, 91, 67, 77, 73, 93, 86, 94, 79, 95, 128, 130, 132, 140, 136, 142, 152, 149, 144, 147, 156, 133, 145, 148, 139, 151, 129, 134, 135, 154, 153, 146, 138, 155, 131, 141, 137, 157, 150, 158, 143, 159, 384, 386, 388, 396, 392, 398, 408, 405, 400, 403, 412, 389, 401, 404, 395, 407, 385, 390, 391, 410, 409, 402, 394, 411, 387, 397, 393, 413, 406, 414, 399, 415, 256, 258, 260, 268, 264, 270, 280, 277, 272, 275, 284, 261, 273, 276, 267, 279, 257, 262, 263, 282, 281, 274, 266, 283, 259, 269, 265, 285, 278, 286, 271, 287, 448, 450, 452, 460, 456, 462, 472, 469, 464, 467, 476, 453, 465, 468, 459, 471, 449, 454, 455, 474, 473, 466, 458, 475, 451, 461, 457, 477, 470, 478, 463, 479, 768, 770, 772, 780, 776, 782, 792, 789, 784, 787, 796, 773, 785, 788, 779, 791, 769, 774, 775, 794, 793, 786, 778, 795, 771, 781, 777, 797, 790, 798, 783, 799, 672, 674, 676, 684, 680, 686, 696, 693, 688, 691, 700, 677, 689, 692, 683, 695, 673, 678, 679, 698, 697, 690, 682, 699, 675, 685, 681, 701, 694, 702, 687, 703, 512, 514, 516, 524, 520, 526, 536, 533, 528, 531, 540, 517, 529, 532, 523, 535, 513, 518, 519, 538, 537, 530, 522, 539, 515, 525, 521, 541, 534, 542, 527, 543, 608, 610, 612, 620, 616, 622, 632, 629, 624, 627, 636, 613, 625, 628, 619, 631, 609, 614, 615, 634, 633, 626, 618, 635, 611, 621, 617, 637, 630, 638, 623, 639, 896, 898, 900, 908, 904, 910, 920, 917, 912, 915, 924, 901, 913, 916, 907, 919, 897, 902, 903, 922, 921, 914, 906, 923, 899, 909, 905, 925, 918, 926, 911, 927, 160, 162, 164, 172, 168, 174, 184, 181, 176, 179, 188, 165, 177, 180, 171, 183, 161, 166, 167, 186, 185, 178, 170, 187, 163, 173, 169, 189, 182, 190, 175, 191, 544, 546, 548, 556, 552, 558, 568, 565, 560, 563, 572, 549, 561, 564, 555, 567, 545, 550, 551, 570, 569, 562, 554, 571, 547, 557, 553, 573, 566, 574, 559, 575, 640, 642, 644, 652, 648, 654, 664, 661, 656, 659, 668, 645, 657, 660, 651, 663, 641, 646, 647, 666, 665, 658, 650, 667, 643, 653, 649, 669, 662, 670, 655, 671, 352, 354, 356, 364, 360, 366, 376, 373, 368, 371, 380, 357, 369, 372, 363, 375, 353, 358, 359, 378, 377, 370, 362, 379, 355, 365, 361, 381, 374, 382, 367, 383, 736, 738, 740, 748, 744, 750, 760, 757, 752, 755, 764, 741, 753, 756, 747, 759, 737, 742, 743, 762, 761, 754, 746, 763, 739, 749, 745, 765, 758, 766, 751, 767, 32, 34, 36, 44, 40, 46, 56, 53, 48, 51, 60, 37, 49, 52, 43, 55, 33, 38, 39, 58, 57, 50, 42, 59, 35, 45, 41, 61, 54, 62, 47, 63, 192, 194, 196, 204, 200, 206, 216, 213, 208, 211, 220, 197, 209, 212, 203, 215, 193, 198, 199, 218, 217, 210, 202, 219, 195, 205, 201, 221, 214, 222, 207, 223, 224, 226, 228, 236, 232, 238, 248, 245, 240, 243, 252, 229, 241, 244, 235, 247, 225, 230, 231, 250, 249, 242, 234, 251, 227, 237, 233, 253, 246, 254, 239, 255, 832, 834, 836, 844, 840, 846, 856, 853, 848, 851, 860, 837, 849, 852, 843, 855, 833, 838, 839, 858, 857, 850, 842, 859, 835, 845, 841, 861, 854, 862, 847, 863, 800, 802, 804, 812, 808, 814, 824, 821, 816, 819, 828, 805, 817, 820, 811, 823, 801, 806, 807, 826, 825, 818, 810, 827, 803, 813, 809, 829, 822, 830, 815, 831, 576, 578, 580, 588, 584, 590, 600, 597, 592, 595, 604, 581, 593, 596, 587, 599, 577, 582, 583, 602, 601, 594, 586, 603, 579, 589, 585, 605, 598, 606, 591, 607, 320, 322, 324, 332, 328, 334, 344, 341, 336, 339, 348, 325, 337, 340, 331, 343, 321, 326, 327, 346, 345, 338, 330, 347, 323, 333, 329, 349, 342, 350, 335, 351, 864, 866, 868, 876, 872, 878, 888, 885, 880, 883, 892, 869, 881, 884, 875, 887, 865, 870, 871, 890, 889, 882, 874, 891, 867, 877, 873, 893, 886, 894, 879, 895, 96, 98, 100, 108, 104, 110, 120, 117, 112, 115, 124, 101, 113, 116, 107, 119, 97, 102, 103, 122, 121, 114, 106, 123, 99, 109, 105, 125, 118, 126, 111, 127, 416, 418, 420, 428, 424, 430, 440, 437, 432, 435, 444, 421, 433, 436, 427, 439, 417, 422, 423, 442, 441, 434, 426, 443, 419, 429, 425, 445, 438, 446, 431, 447, 288, 290, 292, 300, 296, 302, 312, 309, 304, 307, 316, 293, 305, 308, 299, 311, 289, 294, 295, 314, 313, 306, 298, 315, 291, 301, 297, 317, 310, 318, 303, 319, 928, 930, 932, 940, 936, 942, 952, 949, 944, 947, 956, 933, 945, 948, 939, 951, 929, 934, 935, 954, 953, 946, 938, 955, 931, 941, 937, 957, 950, 958, 943, 959, 704, 706, 708, 716, 712, 718, 728, 725, 720, 723, 732, 709, 721, 724, 715, 727, 705, 710, 711, 730, 729, 722, 714, 731, 707, 717, 713, 733, 726, 734, 719, 735, 960, 962, 964, 972, 968, 974, 984, 981, 976, 979, 988, 965, 977, 980, 971, 983, 961, 966, 967, 986, 985, 978, 970, 987, 963, 973, 969, 989, 982, 990, 975, 991, 480, 482, 484, 492, 488, 494, 504, 501, 496, 499, 508, 485, 497, 500, 491, 503, 481, 486, 487, 506, 505, 498, 490, 507, 483, 493, 489, 509, 502, 510, 495, 511, 992, 994, 996, 1004, 1000, 1006, 1016, 1013, 1008, 1011, 1020, 997, 1009, 1012, 1003, 1015, 993, 998, 999, 1018, 1017, 1010, 1002, 1019, 995, 1005, 1001, 1021, 1014, 1022, 1007, 1023};

uint64_t permutation(uint64_t input);

int R1(int input, uint64_t key);
int R2(int input, uint64_t key);
int Gfunc(int input, uint64_t key);
uint64_t tweakey_schedule(uint64_t key3, uint64_t key2, uint64_t key1, uint64_t key0, uint64_t tweak, uint64_t *RK);
int enc(int input, uint64_t key3, uint64_t key2, uint64_t key1, uint64_t key0, uint64_t tweak);

template <unsigned int start, unsigned int end>
inline uint64_t extract_bits(uint64_t input);

template <unsigned int shift, unsigned int size>
inline uint64_t rol(uint64_t num)
{
    uint64_t output = (num << shift) ^ (num >> (size - shift));
    output = output & ((unit64 << size) - 1);

    return output;
}

#ifdef AMD
template <uint64_t m, uint64_t shift>
inline uint64_t bit_permute_step(uint64_t x)
{
    uint64_t t;
    t = ((x >> shift) ^ x) & m;
    x = (x ^ t) ^ (t << shift);
    return x;
}

inline uint64_t permutation(uint64_t x)
{
    x = bit_permute_step<0x0000555000555000, 1>(x);
    x = bit_permute_step<0x0000222333111000, 2>(x);
    x = bit_permute_step<0x020f070307010500, 4>(x);
    x = bit_permute_step<0x004c004300430000, 8>(x);
    x = bit_permute_step<0x0000010000000000, 16>(x);
    x = bit_permute_step<0x0000000041b35a80, 32>(x);
    x = bit_permute_step<0x00004c3300005970, 16>(x);
    x = bit_permute_step<0x0068002c006f002c, 8>(x);
    x = bit_permute_step<0x020e0b0a0a0e0a0e, 4>(x);
    return x;
}

#endif

#ifdef INTEL
inline uint64_t permutation(uint64_t x)
{
    return _pdep_u64(x, 0x0084210842108421) | _pdep_u64(x >> 12, 0x0108421084210842) | _pdep_u64(x >> 24, 0x0210842108421084) | _pdep_u64(x >> 36, 0x0421084210842108) | _pdep_u64(x >> 48, 0x0842108421084210);
}
#endif

inline uint64_t Slayer(uint64_t input)
{
    return double_sbox[input & 0x3ff] ^ (double_sbox[(input >> N) & 0x3ff] << N) ^ (double_sbox[(input >> (2 * N)) & 0x3ff] << (2 * N)) ^ (double_sbox[(input >> (3 * N)) & 0x3ff] << (3 * N)) ^ (double_sbox[(input >> (4 * N)) & 0x3ff] << (4 * N)) ^ (double_sbox[(input >> (5 * N)) & 0x3ff] << (5 * N));
}

int enc(int input, uint64_t key3, uint64_t key2, uint64_t key1, uint64_t key0, uint64_t tweak)
{
    uint64_t RK[ROUNDS];
    tweakey_schedule(key3, key2, key1, key0, tweak, &RK[0]);

    int ct = input;

    for (int i = 0; i < ROUNDS - 1; i++)
    {
        ct = R1(ct, RK[i]);
    }

    ct = R2(ct, RK[ROUNDS - 1]);

    return ct;
}

uint64_t Sigma(uint64_t input)
{
    uint64_t t_0 = input;
    uint64_t t_1 = rol<6, TL>(input);
    uint64_t t_2 = rol<12, TL>(input);
    uint64_t t_3 = rol<19, TL>(input);
    uint64_t t_4 = rol<29, TL>(input);
    uint64_t t_5 = rol<43, TL>(input);
    uint64_t t_6 = rol<51, TL>(input);

    return t_0 ^ t_1 ^ t_2 ^ t_3 ^ t_4 ^ t_5 ^ t_6;
}

uint64_t tweakey_schedule(uint64_t key3, uint64_t key2, uint64_t key1, uint64_t key0, uint64_t tweak, uint64_t *RK)
{
    tweak &= 0xFFFFFFFFFFFF;

    #ifdef INTEL
    tweak = _pdep_u64(tweak, 0x07bdef7bdef7bdef) | _pdep_u64(tweak >> 48, 0x0042108421084210);
    #endif

    #ifdef AMD
    tweak = std::byteswap(tweak);
    tweak = ((tweak & 0x00f0000000000000) >> 8) | (tweak & 0x0000000000f00000) | ((tweak & 0x0000000000008000) >> 11) | ((tweak & 0x0000000f00000000) << 27) | ((tweak & 0x000000f000000000) >> 36) | std::rotl(tweak & 0x000000000f002000, 45) | ((tweak & 0x000f000000000000) >> 9) | ((tweak & 0xf000000000000000) >> 26) | ((tweak & 0x00000000000f0000) >> 1) | ((tweak & 0x00000000f0000000) >> 18) | ((tweak & 0x0000f00000000000) << 10) | (((std::rotl(tweak, 25) & 0x0000001e001e001e) * 0x0001000000001111) & 0x001f0843e0000000) | (((tweak & 0x0000000000005007) * 0x0002020000022200) & 0x8020000000084200);
    tweak = std::rotl(tweak, 35);
    #endif

    uint64_t K0 = key0 & 0xFFFFFFFFFFFFFFF;
    uint64_t K1 = key1 & 0xFFFFFFFFFFFFFFF;
    uint64_t K2 = key2 & 0xFFFFFFFFFFFFFFF;
    uint64_t K3 = key3 & 0xFFFFFFFFFFFFFFF;

    uint64_t T1 = tweak ^ K0;

    uint64_t T2 = Slayer(T1);
    uint64_t T3 = Sigma(T2) ^ K1;

    uint64_t T4 = Slayer(T3) ^ K2;
    uint64_t T5 = permutation(T4);
    uint64_t T6 = Slayer(T5);

    uint64_t T7 = Sigma(T6) ^ K3;
    uint64_t T8 = Slayer(T7);

    RK[0] = extract_bits<0, 30>(T1);
    RK[1] = extract_bits<30, 60>(T1);

    RK[2] = extract_bits<0, 30>(T3);
    RK[3] = extract_bits<30, 60>(T3);

    RK[4] = extract_bits<0, 30>(T6);
    RK[5] = extract_bits<30, 60>(T6);

    RK[6] = extract_bits<0, 30>(T8);
    RK[7] = extract_bits<30, 60>(T8);

    return T8;
}

template <unsigned int start, unsigned int end>
inline uint64_t extract_bits(uint64_t input)
{
    constexpr uint64_t mask = (unit64 << end) - (unit64 << start);

    return (input & mask) >> start;
}

int R1(int input, uint64_t key)
{
    int SK1 = (key & 0x3E000000) >> 25;

    int xR = extract_bits<0, 5>(input);
    int xL = extract_bits<5, 10>(input);

    int temp_xL = xL;

    xL = xR ^ Gfunc(temp_xL, key);

    xR = sbox[temp_xL ^ SK1];

    return (xL << 5) ^ xR;
}

int R2(int input, uint64_t key)
{
    int SK1 = (key & 0x3E000000) >> 25;

    int xR = extract_bits<0, 5>(input);
    int xL = extract_bits<5, 10>(input);

    int temp_xL = xL;

    xR = xR ^ Gfunc(temp_xL, key);

    xL = sbox[temp_xL] ^ SK1;

    return (xL << 5) ^ xR;
}

int Gfunc(int input, uint64_t key)
{
    int x0 = input & extract_bits<0, 5>(key);
    int x1 = rol<1, N / 2>(input) & extract_bits<5, 10>(key);
    int x2 = rol<2, N / 2>(input) & extract_bits<10, 15>(key);
    int x3 = rol<3, N / 2>(input) & extract_bits<15, 20>(key);
    int x4 = rol<4, N / 2>(input) & extract_bits<20, 25>(key);

    int output = x0 ^ x1 ^ x2 ^ x3 ^ x4 ^ (rol<1, N / 2>(input) & rol<2, N / 2>(input));

    return output;
}