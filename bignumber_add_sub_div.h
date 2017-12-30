#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifndef BIGNUMBER_ADD_SUB_DIV_H
#define BIGNUMBER_ADD_SUB_DIV_H

class BigIntOp {
public:

static void div_by_1digit(uint32_t * divisor, int divisor_size, 
                   uint32_t div_by, 
                   uint32_t * quot, int quot_size, int * out_size, 
                   uint32_t radix)
{
    if (quot_size <= 0) {
        *out_size = 0;
        return;
    }
    int k = 0;
    uint32_t remain = 0;
    int j = quot_size - 1;
    for (int i = divisor_size - 1; i >= 0 && j >= 0; --i) {
        uint64_t cur_val = 1L * remain * radix + divisor[i];
        quot[j] = cur_val / div_by;
        remain = cur_val % div_by;
        j--;
        k++; //printf("k=%d\n", k);
    }
    for (int i = j; i >= 0; --i) {
        if (remain == 0) {
            break;
        }
        uint64_t cur_val = 1L * remain * radix;
        quot[j] = cur_val / div_by;
        remain = cur_val % div_by;
        k++; //printf("k=%d\n", k);
        j--;
    }
    if (quot[quot_size-1] == 0) {
        *out_size = k - 1;
    } else {
        *out_size = k;
    }
}

// big_zero_tail_cnt: bigger 在 bigger_size之外，尾部还有多少个零
// bigger/small: 补0后的长度哪个更长
static void add(uint32_t * bigger, int bigger_size,
         uint32_t * smaller, int smaller_size,
         uint32_t *out, int*out_size,
         int big_zero_tail_cnt, int radix)
{
    // 1. big_zero_tail_cnt > 0:
    //    [-----------bigger---------]0000
    //          [--------smaller---------]
    // 2. big_zero_tail_cnt == 0:
    //    [-----------bigger-------------]
    //          [--------smaller---------]
    // 3. big_zero_tail_cnt< 0:
    //    [-----------bigger-------------]
    //          [--------smaller-----]0000

    uint32_t remain = 0;

    // == tail part
    uint32_t * tail = big_zero_tail_cnt >= 0 ? smaller : bigger;
    int tail_size = big_zero_tail_cnt >= 0 ? big_zero_tail_cnt : -big_zero_tail_cnt;
    for (int i = 0; i < tail_size; ++i) {
        out[i] = tail[i];
    }
    // == common part
    int common_part_cnt  = (big_zero_tail_cnt >= 0) ? smaller_size - big_zero_tail_cnt : smaller_size;
    uint32_t * bigger_1  = (big_zero_tail_cnt >= 0) ? bigger : bigger+big_zero_tail_cnt;
    uint32_t * smaller_1 = (big_zero_tail_cnt >= 0) ? smaller + big_zero_tail_cnt : smaller;
    uint32_t * out_1 = out + big_zero_tail_cnt;
    for (int i = 0; i < common_part_cnt; ++i) {
        uint32_t res = bigger_1[i] + smaller_1[i] + remain;
        if (res >= radix) {
            out_1[i] = res - radix;
            remain = 1;
        } else {
            out_1[i] = res;
            remain = 0;
        }
    }
    // == head part
    int head_cnt  = bigger_size - smaller_size + big_zero_tail_cnt;
    uint32_t * bigger_2  = bigger - head_cnt;
    uint32_t * out_2 = (big_zero_tail_cnt > 0) ? out + smaller_size : out + (smaller_size - big_zero_tail_cnt);
    for (int i = 0; i < head_cnt; ++i) {
        if (remain == 0) {
            out_2[i] = bigger_2[i];
            continue;
        }
        uint32_t res = bigger_2[i] + remain;
        if (res >= radix) {
            out_2[i] = res - radix;
            remain = 1;
        }
    }
    // == 
    *out_size = bigger_size + big_zero_tail_cnt;
    if (remain) {
        *out_size += 1;
        out[*out_size - 1] = remain;
    }
}

// big_zero_tail_cnt: bigger 在 bigger_size之外，尾部还有多少个零
// bigger/small: 补0后的长度哪个更长. 作为减法，两者长度一样的时候，未必bigger > small
static void sub(uint32_t * bigger, int bigger_size,
         uint32_t * smaller, int smaller_size,
         uint32_t * out, int * out_size, int * first_is_big,
         int big_zero_tail_cnt, int radix)
{
    *first_is_big = 1;
    int eq_cnt = 0;
    if (bigger_size - smaller_size + big_zero_tail_cnt == 0) {
        int min_size = bigger_size > smaller_size ? smaller_size : bigger_size;
        for (int i = 0; i < min_size; ++i) {
            if (bigger[bigger_size - 1 - i] == smaller[smaller_size - 1 - i]) {
                eq_cnt += 1;
            } else if (bigger[bigger_size - 1 - i] > smaller[smaller_size - 1 - i]) {
                *first_is_big = 1;
                break;
            } else {
                *first_is_big = 0;
                break;
            }
        }
        if (*first_is_big == 0) {
            // switch
            uint32_t * x = bigger;
            uint32_t xx = bigger_size;
            bigger = smaller;
            bigger_size = smaller_size;
            smaller = x;
            smaller_size = xx;
            big_zero_tail_cnt *= -1;
        }
    }

    // bigger is really bigger now

    // 1. big_zero_tail_cnt > 0:
    //    [-----------bigger---------]0000
    //          [--------smaller---------]
    // 2. big_zero_tail_cnt == 0:
    //    [-----------bigger-------------]
    //          [--------smaller---------]
    // 3. big_zero_tail_cnt< 0:
    //    [-----------bigger-------------]
    //          [--------smaller-----]0000

    uint32_t borrow = 0;
    int res_size = 0;
    int k = 0;

    // == tail part
    int tail_size = big_zero_tail_cnt >= 0 ? big_zero_tail_cnt : -big_zero_tail_cnt;
    if (big_zero_tail_cnt > 0) {
        for (int i = 0; i < tail_size; ++i) {
            if (borrow || smaller[i]) {
                out[i] = radix - borrow - smaller[i];
                borrow = 1;
            } else {
                out[i] = 0;
                borrow = 0;
            }
            k++;
            if (out[i]) {
                res_size = k;
            }
        }
    } else {
        for (int i = 0; i < tail_size; ++i) {
            out[i] = bigger[i];
            k++;
            if (out[i]) {
                res_size = k;
            }
        }
    }
    // == common part
    int common_part_cnt  = (big_zero_tail_cnt >= 0) ? smaller_size - big_zero_tail_cnt : smaller_size;
    uint32_t * bigger_1  = (big_zero_tail_cnt >= 0) ? bigger : bigger+big_zero_tail_cnt;
    uint32_t * smaller_1 = (big_zero_tail_cnt >= 0) ? smaller + big_zero_tail_cnt : smaller;
    uint32_t * out_1 = out + big_zero_tail_cnt;
    if (eq_cnt > 0) {
        common_part_cnt -= eq_cnt;
        uint32_t * out_1 = out + common_part_cnt;
        for (int i = 0; i < eq_cnt; ++i) {
            out_1[i] = 0;
        }
    }
    for (int i = 0; i < common_part_cnt; ++i) {
        int32_t res = (int32_t)bigger_1[i] - smaller_1[i] - borrow;
        if (res < 0) {
            out_1[i] = uint32_t(res + radix);
            borrow = 1;
        } else {
            out_1[i] = res;
            borrow = 0;
        }
        k++;
        if (out_1[i]) {
            res_size = k;
        }
    }
    // == head part
    int head_cnt  = bigger_size - smaller_size + big_zero_tail_cnt;
    uint32_t * bigger_2  = bigger - head_cnt;
    uint32_t * out_2 = (big_zero_tail_cnt > 0) ? out + smaller_size : out + (smaller_size - big_zero_tail_cnt);
    for (int i = 0; i < head_cnt; ++i) {
        if (borrow == 0) {
            out_2[i] = bigger_2[i];
        } else {
            int32_t res = (int32_t)bigger_2[i] - borrow;
            if (res < 0) {
                out_2[i] = uint32_t(res + radix);
                borrow = 1;
            }
        }
        k++;
        if (out_2[i]) {
            res_size = k;
        }
    }
    *out_size = res_size;
}
};
//int main()
//{
//    int s = 100;
//    uint32_t aa[] = {1, 1};
//    uint32_t bb[s];
//    uint32_t div_by = 7;
//    int out_size;
//    div_by_1digit(aa, sizeof(aa)/sizeof(uint32_t), div_by, bb, s, &out_size, 10);
//    for (int j = s-1; j >=0; --j) {
//        printf("%u ", bb[j]);
//    }
//    printf("\n%d\n", out_size);
//    return 0;
//}

#endif
