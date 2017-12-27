#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

#ifndef FFT_COMMON_H
#define FFT_COMMON_H

void p(uint32_t * aa, int aa_s, int radix, bool show_all=false);

// reverse the k-bit representation of t
int rev_bin(int t, int k);

// 弃九法验证计算的正误
int qi9_check(uint32_t * aa, int aa_s, int radix, int show=0);

#endif
