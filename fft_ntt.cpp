#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdint.h>
#include "fft_ntt.h"
#include "bignumber_add_sub_div.h"

template <> long Zp<0>::P = 2113929217;
template <> long Zp<0>::G = 5;
template <> long Zp<1>::P = 2013265921;
template <> long Zp<1>::G = 31;
template <> long Zp<2>::P = 1811939329;
template <> long Zp<2>::G = 13;

FntMul fnt;
FftMul fft;

extern "C" {

    void fnt_init(int max_digit) {
        fnt.init(max_digit * 2);
    }

    void fft_init(int max_digit) {
        fft.init(max_digit * 2);
    }

    // mul
    int do_fft(uint32_t * aa, int aa_s, uint32_t * bb, int bb_s, int radix, uint32_t * cc)
    {
        int cc_s;
        fft.fast_prod(aa, aa_s, bb, bb_s, radix, cc, cc_s);
        return cc_s;
    }
    int do_fft_1(uint32_t * aa, int aa_start, int aa_s, 
                 uint32_t * bb, int bb_start, int bb_s, 
                 int radix, uint32_t * cc)
    {
        int cc_s;
        fft.fast_prod(aa + aa_start, aa_s, bb + bb_start, bb_s, radix, cc, cc_s);
        return cc_s;
    }
    int do_fnt(uint32_t * aa, int aa_s, uint32_t * bb, int bb_s, int radix, uint32_t * cc)
    {
        int cc_s;
        fnt.fast_prod(aa, aa_s, bb, bb_s, radix, cc, cc_s);
        return cc_s;
    }
    int do_fnt_1(uint32_t * aa, int aa_start, int aa_s, 
                 uint32_t * bb, int bb_start, int bb_s, 
                 int radix, uint32_t * cc)
    {
        int cc_s;
        fnt.fast_prod(aa + aa_start, aa_s, bb + bb_start, bb_s, radix, cc, cc_s);
        return cc_s;
    }

    // div
    void big_div_by_1digit(uint32_t * divisor, int divisor_size, 
                           uint32_t div_by, 
                           uint32_t * quot, int quot_size, int * out_size, 
                           uint32_t radix)
    {
        BigIntOp::div_by_1digit(divisor, divisor_size, div_by, quot, quot_size, out_size, radix);
    }
    void big_div_by_1digit_1(uint32_t * divisor, int divisor_start, int divisor_size, 
                             uint32_t div_by, 
                             uint32_t * quot, int quot_size, int * out_size, 
                             uint32_t radix)
    {
        BigIntOp::div_by_1digit(divisor + divisor_start, divisor_size, 
                                div_by, quot, quot_size, out_size, radix);
    }

    // add
    void big_add(uint32_t * bigger,  int bigger_size,
                 uint32_t * smaller, int smaller_size,
                 uint32_t *out,      int*out_size,
                 int big_zero_tail_cnt, int radix)
    {
        BigIntOp::add(bigger, bigger_size, smaller, smaller_size, out, out_size, 
                      big_zero_tail_cnt, radix);
    }
    void big_add_1(uint32_t * bigger,  int bigger_start,  int bigger_size,
                   uint32_t * smaller, int smaller_start, int smaller_size,
                   uint32_t *out,      int*out_size,
                   int big_zero_tail_cnt, int radix)
    {
        BigIntOp::add(bigger + bigger_start, bigger_size, smaller + smaller_start, smaller_size, out, out_size, 
                      big_zero_tail_cnt, radix);
    }


    // sub
    void big_sub(uint32_t * bigger,  int bigger_size,
                 uint32_t * smaller, int smaller_size,
                 uint32_t * out,     int * out_size, int * first_is_big,
                 int big_zero_tail_cnt, int radix)
    {
        BigIntOp::sub(bigger, bigger_size, smaller, smaller_size, out, out_size, 
                      first_is_big, big_zero_tail_cnt, radix);
    }
    void big_sub_1(uint32_t * bigger,  int bigger_start,  int bigger_size,
                   uint32_t * smaller, int smaller_start, int smaller_size,
                   uint32_t * out,     int * out_size,    int * first_is_big,
                   int big_zero_tail_cnt, int radix)
    {
        BigIntOp::sub(bigger + bigger_start, bigger_size, 
                      smaller + smaller_start, smaller_size, 
                      out, out_size, 
                      first_is_big, big_zero_tail_cnt, radix);
    }
}

