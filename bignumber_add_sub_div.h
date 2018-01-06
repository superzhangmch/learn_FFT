#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#ifndef BIGNUMBER_ADD_SUB_DIV_H
#define BIGNUMBER_ADD_SUB_DIV_H


class BigIntOp {
public:

// divisor: 被除数， divisor_size： 被除数在基数radix下的长度
// div_by: 除数，是基数radix内的一个数字
// quot: 除法结果
// quot_size：希望运算结果是多少位的
// out_size: 实际结果是多少位的。如果能整除，往往小于quot_size，否则会算到quot_size的长度
//           不能整除的时候，小数点应该标在哪里，应该由调用方自己算出
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
    int j = quot_size - 1;

    uint32_t remain = 0;
    int sz = 0;
    if (divisor[divisor_size - 1] < div_by) {
        remain = divisor[divisor_size - 1];
        sz = divisor_size - 2;
    } else {
        remain = 0;
        sz = divisor_size - 1;
    }
    for (int i = sz; i >= 0 && j >= 0; --i) {
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
    if (k > quot_size) {
        k = quot_size;
    }
    *out_size = k;
}

// big_zero_tail_cnt: bigger 在 bigger_size之外，尾部还有多少个零
// bigger/small: 补0后的长度哪个更长
// 各参数含义可看 sub(...) 函数
static void add(uint32_t * bigger, int bigger_size,
         uint32_t * smaller, int smaller_size,
         uint32_t *out, int*out_size,
         int big_zero_tail_cnt, int radix)
{
    // 1. big_zero_tail_cnt >= 0:
    //    [-----------bigger---------]0000
    //          [--------smaller---------]
    // 2. big_zero_tail_cnt >= 0 and has gap
    //    [---bigger---]00gap0000000000000
    //                         [-smaller-]
    // 3. big_zero_tail_cnt< 0:
    //    [-----------bigger-------------]
    //          [--------smaller-----]0000

    uint32_t remain = 0;
    bool has_gap = smaller_size < big_zero_tail_cnt;;

    // == tail part
    uint32_t * tail = 0;
    int tail_size = 0;
    int tail_size_1 = 0;
    if (big_zero_tail_cnt < 0) { // case 3
        tail = bigger;
        tail_size = -big_zero_tail_cnt;
        tail_size_1 = tail_size;
    } else {
        if (has_gap) { // case 2
            tail = smaller;
            tail_size = big_zero_tail_cnt;
            tail_size_1 = smaller_size;
        } else { // case 1
            tail = smaller;
            tail_size = big_zero_tail_cnt;
            tail_size_1 = tail_size;
        }
    }
    for (int i = 0; i < tail_size_1; ++i) {
        out[i] = tail[i];
    }
    for (int i = tail_size_1; i < tail_size; ++i) {
        out[i] = 0;
    }
    // == common part
    int common_part_cnt  = 0;
    uint32_t * bigger_1  = 0;
    uint32_t * smaller_1 = 0;
    uint32_t * out_1 = 0;
    if (big_zero_tail_cnt < 0) {
        common_part_cnt = smaller_size;
        bigger_1 = bigger-big_zero_tail_cnt;
        smaller_1 = smaller;
        out_1 = out -big_zero_tail_cnt;
    } else {
        if (has_gap) {
            common_part_cnt = 0;
        } else {
            common_part_cnt = smaller_size - big_zero_tail_cnt;
            bigger_1 = bigger;
            smaller_1 = smaller + big_zero_tail_cnt;
            out_1 = out + big_zero_tail_cnt;
        }
    }
    //printf("common cnt=%d big_sz=%d sm_sz=%d bztc=%d\n", common_part_cnt, bigger_size, smaller_size, big_zero_tail_cnt);
    for (int i = 0; i < common_part_cnt; ++i) {
        uint32_t res = bigger_1[i] + smaller_1[i] + remain;
        //printf(" - common big:%d small:%d remain:%d\n", bigger_1[i], smaller_1[i], remain);
        if (res >= radix) {
            out_1[i] = res - radix;
            remain = 1;
        } else {
            out_1[i] = res;
            remain = 0;
        }
    }
    // == head part
    int head_cnt = 0;
    uint32_t * bigger_2 = 0;
    uint32_t * out_2 = 0;
    if (big_zero_tail_cnt < 0) {
        head_cnt = bigger_size - smaller_size + big_zero_tail_cnt;
        bigger_2 = bigger + bigger_size - head_cnt;
        out_2 = out + (smaller_size - big_zero_tail_cnt);
    } else {
        if (has_gap) {
            head_cnt = bigger_size;
            bigger_2 = bigger;
            out_2 = out + big_zero_tail_cnt;
        } else {
            head_cnt = bigger_size - smaller_size + big_zero_tail_cnt;
            bigger_2 = bigger + bigger_size - head_cnt;
            out_2 = out + smaller_size;
        }
    }
    //printf("head cnt=%d remain=%d\n", head_cnt, remain);
    for (int i = 0; i < head_cnt; ++i) {
        //printf(" - head: big:%d remain:%d\n", bigger_2[i], remain);
        if (remain == 0) {
            out_2[i] = bigger_2[i];
        } else {
            uint32_t res = bigger_2[i] + remain;
            if (res >= radix) {
                out_2[i] = res - radix;
                remain = 1;
            } else {
                out_2[i] = res;
                remain = 0;
            }
        }
    }
    // == 
    *out_size = bigger_size;
    if (big_zero_tail_cnt > 0) {
        *out_size = bigger_size + big_zero_tail_cnt;
    }
    //printf("outsize=%d remain=%d\n", *out_size, remain);
    if (remain) {
        *out_size += 1;
        out[*out_size - 1] = remain;
    }
    //printf("outsize=%d remain=%d\n", *out_size, remain);
}

// big_zero_tail_cnt: bigger 在 bigger_size之外，尾部还有多少个零.如果是smaller补零，则是负数.
// bigger: 指的是两个数字在末尾该补0补足0后，长度更长的那个.
// 举例说明：10进制下，A = 100 * 10^5, B=2000, 那么A末尾需要补充5个零，才能与B计算
//           补0后A更长，所以bigger=A, bigger_size=3, smaller=B, smaller_size=4, big_zero_tail_cnt=5
//           如果 A=10000000, B=200*10^2,则bigger=A, bigger_size=8, smaller=B, smaller_size=3, big_zero_tail_cnt=-2
// bigger/small: 补0后的长度哪个更长. 作为减法，两者长度一样的时候，未必bigger > small
// out_size: 返回最终结果的长度
// first_is_big: 返回bigger是否真的比smaller大. (bigger比smaller只是长度长些.两者一样长的时候，不一定bigger>smaller)
static void sub(uint32_t * bigger, int bigger_size,
         uint32_t * smaller, int smaller_size,
         uint32_t * out, int * out_size, int * first_is_big,
         int big_zero_tail_cnt, int radix)
{
    //struct timeval tpstart, tpend;
    //gettimeofday(&tpstart,NULL);
    *first_is_big = 1;
    int eq_cnt = 0;
    //printf("abc bigger_size=%d, smaller_size=%d, big_zero_tail_cnt=%d\n", bigger_size, smaller_size, big_zero_tail_cnt);
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
    // 2. big_zero_tail_cnt > 0 and has_gap:
    //    [---bigger---]0gap00000000000000
    //                       [--smaller--]
    // 3. big_zero_tail_cnt< 0:
    //    [-----------bigger-------------]
    //          [--------smaller-----]0000

    uint32_t borrow = 0;
    int res_size = 0;
    int k = 0;
    bool has_gap = smaller_size < big_zero_tail_cnt;;

    // == tail part

    int tail_size = big_zero_tail_cnt >= 0 ? big_zero_tail_cnt : -big_zero_tail_cnt;
    if (big_zero_tail_cnt > 0) {
        int tail_size_1 = has_gap ? smaller_size : tail_size; 
        for (int i = 0; i < tail_size_1; ++i) {
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
        for (int i = tail_size_1; i < tail_size; ++i) {
            if (borrow) {
                out[i] = radix - borrow;
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
    int common_part_cnt  = 0;
    uint32_t * bigger_1  = 0;
    uint32_t * smaller_1 = 0;
    uint32_t * out_1 = 0;
    if (big_zero_tail_cnt < 0) {
        common_part_cnt = smaller_size;
        bigger_1 = bigger-big_zero_tail_cnt;
        smaller_1 = smaller;
        out_1 = out -big_zero_tail_cnt;
    } else {
        if (has_gap) {
            common_part_cnt = 0;
        } else {
            common_part_cnt = smaller_size - big_zero_tail_cnt;
            bigger_1 = bigger;
            smaller_1 = smaller + big_zero_tail_cnt;
            out_1 = out + big_zero_tail_cnt;
        }
    }
    if (eq_cnt > 0) {
        // aa=[xxxx---------]
        // bb=[xxxx------]000 or
        // aa=[xxxx------]000
        // bb=[xxxx---------]
        common_part_cnt -= eq_cnt;
        int abs_bztc = big_zero_tail_cnt>0 ? big_zero_tail_cnt : -big_zero_tail_cnt;
        uint32_t * out_1 = out + common_part_cnt + abs_bztc;
        for (int i = 0; i < eq_cnt; ++i) {
            out_1[i] = 0;
        }
    }
    //printf("common=%d eq_cnt=%d\n", common_part_cnt, eq_cnt);
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
    int head_cnt = 0;
    uint32_t * bigger_2 = 0;
    uint32_t * out_2 = 0;
    if (big_zero_tail_cnt < 0) {
        head_cnt = bigger_size - smaller_size + big_zero_tail_cnt;
        bigger_2 = bigger + bigger_size - head_cnt;
        out_2 = out + (smaller_size - big_zero_tail_cnt);
    } else {
        if (has_gap) {
            head_cnt = bigger_size;
            bigger_2 = bigger;
            out_2 = out + big_zero_tail_cnt;
        } else {
            head_cnt = bigger_size - smaller_size + big_zero_tail_cnt;
            bigger_2 = bigger + bigger_size - head_cnt;
            out_2 = out + smaller_size;
        }   
    } 
    //printf("head_cnt=%d\n", head_cnt);
    for (int i = 0; i < head_cnt; ++i) {
        if (borrow == 0) {
            out_2[i] = bigger_2[i];
        } else {
            int32_t res = (int32_t)bigger_2[i] - borrow;
            if (res < 0) {
                out_2[i] = uint32_t(res + radix);
                borrow = 1;
            } else {
                out_2[i] = res;
                borrow = 0;
            }
        }
        k++;
        if (out_2[i]) {
            res_size = k;
        }
    }
    *out_size = res_size;

    //gettimeofday(&tpend,NULL);
    //double tm = ((tpend.tv_sec-tpstart.tv_sec)*1000000+(tpend.tv_usec-tpstart.tv_usec));
    //printf("====tm=%.4f, sz=%d\n", tm/1000000, *out_size);
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
