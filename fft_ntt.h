#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"
#include "zp.h"
#include "fft_common.h"
#include "uint128t.h"
#include <stdint.h>
#include <sys/time.h>
#include <pthread.h>

#ifndef FFT_NTT_H
#define FFT_NTT_H

typedef Zp<0> Zp0;
typedef Zp<1> Zp1;
typedef Zp<2> Zp2;

// 乘法时，乘数长度值大于多少后开始使用多线程
#define MIN2USE_THREAD 5000

// FFT与FNT有共同之处，所以都继承自这个基类。
// TComplex: 如果是FFT则是复数类，如果是FNT则是Zp数域类
template <typename TComplex>
class ProductFFT {
public:
    int _max_result_len;
    int _max_k;
    TComplex * _omega_fft[32];
    uint32_t * _pi_k[32];

    void prepare_data(int k)
    {
        int n = (int)pow(2, k);
        TComplex * omega_fft = new TComplex[n];
        uint32_t * pi_k = new uint32_t[n];

        // if (k == 5) printf("%d--\n", k);
        // if (k == 5) printf("%d->%.5f+%.5fi\n", i, w_ret[i][0], w_ret[i][1]);
        TComplex::get_w_pow(n, omega_fft);

        for (int i = 0; i < n; ++i) {
            pi_k[i] = rev_bin(i, k);
            // if (k == 5) printf("%d->%d\n", i, pi_k_ret[i]);
        }
        _omega_fft[k] = omega_fft;
        _pi_k[k] = pi_k;
    }

public:
    ProductFFT() {
        _max_result_len = _max_k = 0;
        for (int i = 0; i < 32; ++i) {
            _omega_fft[i] = NULL;
            _pi_k[i] = NULL;
        }
    }
    virtual ~ProductFFT() {
        for (int i = 0; i < 32; ++i) {
            if (_omega_fft[i]) {
                delete [] _omega_fft[i];
                _omega_fft[i] = NULL;
            }
            if (_pi_k[i]) {
                delete [] _pi_k[i];
                _pi_k[i] = NULL;
            }
        }
    }

    void init(int max_len)
    {
        int k = int(ceil(log(1.L * max_len) / log(2.)));
        _max_result_len = int(round(2 * pow(2., k)));
        _max_k = k + 1;

        // 两个数相乘，需要适配到最小的一个2次幂长度，以利用FFT。该长度下，需要对应的FFT参数。
        // 但其实只用给最大长度算好参数就可以了，这里为了方便全算了
        for (int i = 1; i < _max_k + 1; ++i) {
            printf("init %d/%d\n", i, _max_k);
            prepare_data(i);
        }
        printf("init done\n");
    }

    virtual int fft_ntt(int k, int n2, 
                        uint32_t * aa, int aa_length, void * trans_aa, int calc_aa,
                        uint32_t * bb, int bb_length, void * trans_bb, int calc_bb,
                        int aa_eq_bb,
                        int radix, uint32_t * out, int &out_size) = 0;
    virtual void calc_trans(int k, uint32_t * aa, int aa_length, void * trans_aa) = 0;

    /*
     * FFT/FNT运算对外暴露的接口
     *
     */
    // trans_aa: 对aa作正向FFT/FNT变换的时候提供的内存buffer。
    //           如果不提供(令取值NULL)，则内部会创建;这时就不能从外部得到变换结果了。
    //           如果要提供trans_aa, 调用方应该保证有足够空间容纳变换结果.
    // 如果trans_aa!=NULL: calc_aa=0, 表示内部不再做正向FFT/FNT而直接复用trans_aa内的内容当做正向结果，参与整体运算。
    //                     这时为了aa多次参与运算的时候，可以复用正向FFT、FNT结果。
    // trans_bb/calc_bb: 同trans_aa、calc_aa
    // out_size: 运算结果长度
    // trans_aa 内存长度估计：aa/bb会扩展成能容纳(aa_length+bb_length)的最小二次幂,记为n2.
    //                        sizeof(Complex) == 16, sizeof(Zp) == 8
    //                        对于FFT：n2*sizeof(Complex) = 16*n2
    //                        对于FNT：n2*sizeof(Zp)*3 = n2 * 8 * 3
    int fast_prod(uint32_t * aa, int aa_length, void * trans_aa, int calc_aa,
                  uint32_t * bb, int bb_length, void * trans_bb, int calc_bb,
                  int aa_eq_bb,
                  int radix, uint32_t * out, int &out_size)
    {
        // aa or bb == 0
        if (aa_length == 0 || bb_length == 0) {
            out_size = 0;
            return 0;
        }
        // aa or bb is small
        // aa或bb是一位数的时候，直接乘就是了
        if (aa_length == 1 || bb_length == 1) {
            int big_size;
            uint32_t small_val, *big;
            if (aa_length == 1) {
                small_val = aa[0];
                big = bb;
                big_size = bb_length;
            } else {
                small_val = bb[0];
                big = aa;
                big_size = aa_length;
            }
            if (small_val == 0) {
                out_size = 0;
                return 0;
            }
            int out_len = big_size;
            int remain = 0;
            for (int i = 0; i < big_size; ++i) {
                uint64_t res = 1L * big[i] * small_val + remain;
                if (res >= radix) {
                    out[i] = res % radix;
                    remain = res / radix;
                } else {
                    out[i] = res;
                    remain = 0;
                }
            }
            if (remain) {
                out[big_size] = remain;
                out_len = big_size + 1;
            }
            out_size = out_len;
            return 0; 
        }

        // aa or bb both big
        // int max_len = aa_length > bb_length ? aa_length : bb_length;
        // int k = int(ceil(log(max_len) / log(2))) + 1;
        // 结果的长度为刚刚超过两个乘数的长度和的2**k即可
        int k = int(ceil(log(aa_length + bb_length) / log(2)));
        int n2 = int(round(pow(2, k)));
        if (k > _max_k) {
            out_size = 0;
            return -1;
        }
        fft_ntt(k, n2, 
                aa, aa_length, trans_aa, calc_aa,
                bb, bb_length, trans_bb, calc_bb,
                aa_eq_bb, radix, out, out_size);
        //printf("--%d %d %d %d %d\n", k, n2, aa_length, bb_length, max_len);
        return 0;
    }
};

class FftMul: public ProductFFT<Complex> {

public:
    // 只做FFT正向变换，结果放入trans_aa。trans_aa大小由调用方保证足够容纳结果
    // trans_aa大小应不小于n2*sizeof(Complex)
    void calc_trans(int k, uint32_t * aa, int aa_length, void * trans_aa)
    {
        int n2 = int(round(pow(2, k)));
        FFT(k, n2, aa, NULL, aa_length, (Complex*)trans_aa);
    }

    // <<-- multi thread to speed up
    struct thread_do_fft_args_t {
        int k; int n2; uint32_t * aa; int aa_length; Complex *transform_aa;
        void *cls_obj;
    };
    static void * thread_do_fft(void * args) {
        thread_do_fft_args_t & a = *((thread_do_fft_args_t*)args);
        ((FftMul*)a.cls_obj)->FFT(a.k, a.n2, a.aa, NULL, a.aa_length, a.transform_aa);
    }
    // -->>

    // trans_aa/trans_bb: 是否提供内存存放aa/bb变换后的结果。
    // calc_fft_aa: 如果提供了trans_aa，同时calc_fft_aa=1，则变换结果放入trans_aa；
    //              如果提供了trans_aa，同时calc_fft_aa=0，直接使用trans_aa中内容作为正向FFT变换结果，
    //              如果没提供trans_aa，忽略calc_fft_aa
    // aa_eq_bb: 指明两个乘数是否是同一个数。如果是，可以少算一次正向FFT，从而加速运算
    // out_size: 返回结果长度
    // k/n2: n2与k满足n2==2^k, aa/bb会扩展到这么长的长度
    // trans_aa长度应不小于n2*sizeof(Complex)
    int fft_ntt(int k, int n2, 
                uint32_t * aa, int aa_length, void * trans_aa, int calc_fft_aa,
                uint32_t * bb, int bb_length, void * trans_bb, int calc_fft_bb,
                int aa_eq_bb,
                int radix, uint32_t * out, int &out_size)
    {
        if (trans_aa == NULL) {
            calc_fft_aa = 1;
        }
        if (trans_bb == NULL) {
            calc_fft_bb = 1;
        }
        Complex *buf = new Complex[2*n2];
        Complex *transform_aa = trans_aa ? (Complex*)trans_aa: buf;
        Complex *transform_bb = trans_bb ? (Complex*)trans_bb: buf + n2;
        Complex *transform_cc = buf + n2;
        Complex *transform_dd = buf;

        if ((calc_fft_aa && calc_fft_bb && (!aa_eq_bb)) && n2 > MIN2USE_THREAD) {
            // multi thread to speed up
            thread_do_fft_args_t args_0 = {k, n2, aa, aa_length, transform_aa, this};
            thread_do_fft_args_t args_1 = {k, n2, bb, bb_length, transform_bb, this};
            pthread_t thread_0;
            pthread_t thread_1;
            pthread_create(&thread_0, NULL, thread_do_fft, &args_0);
            pthread_create(&thread_1, NULL, thread_do_fft, &args_1);
            pthread_join(thread_0, NULL);
            pthread_join(thread_1, NULL);
        } else {
            if (calc_fft_aa) {
                FFT(k, n2, aa, NULL, aa_length, transform_aa);
            }
            if (calc_fft_bb && (!aa_eq_bb)) {
                FFT(k, n2, bb, NULL, bb_length, transform_bb);
            }
        }

        if (aa_eq_bb) {
            for (int i = 0; i < n2; i++) {
                transform_cc[i] = transform_aa[i] * transform_aa[i];
            }
        } else {
            for (int i = 0; i < n2; i++) {
                transform_cc[i] = transform_aa[i] * transform_bb[i];
            }
        }

        rev_FFT(k, n2, NULL, transform_cc, n2, transform_dd);

        int max_res_len = aa_length + bb_length;
        int remain = 0;
        for (int i = 0; i < max_res_len; ++i) {
            int o = transform_dd[i].to_int() + remain;
            if (o >= radix) {
                out[i] = o % radix;
                remain = o / radix;
            } else {
                out[i] = o;
                remain = 0;
            }
        }

        out_size = out[max_res_len-1] != 0 ? max_res_len : max_res_len - 1;
        delete [] buf;
        return 0;
    }

    void FFT(int k, int n, uint32_t * input, Complex *input1, int input_size, Complex *transform)
    {
        Complex *omega = _omega_fft[k];
        uint32_t * pi_k = _pi_k[k];
        uint32_t * P = input;
        Complex * P1 = input1;

        if (P1 == NULL) {
            for (int t = 0; t < n-1; t+= 2) {
                uint32_t P_pi_k_t = pi_k[t] < input_size ? P[pi_k[t]] : 0;
                uint32_t P_pi_k_t_1 = pi_k[t+1] < input_size ? P[pi_k[t+1]] : 0;
                transform[t]   = P_pi_k_t + P_pi_k_t_1;
                transform[t+1] = P_pi_k_t - P_pi_k_t_1;
            }
        } else {
            for (int t = 0; t < n-1; t+= 2) {
                Complex P_pi_k_t = pi_k[t] < input_size ? P1[pi_k[t]] : 0;
                Complex P_pi_k_t_1 = pi_k[t+1] < input_size ? P1[pi_k[t+1]] : 0;

                transform[t]   = P_pi_k_t + P_pi_k_t_1;
                transform[t+1] = P_pi_k_t - P_pi_k_t_1;
            }
        }

        int m = n/2;
        int num = 2;
        int pow_d_2 = (int)round(pow(2, k-2));
        for (int d = k-2; d > -1; d--) {
            m /= 2;
            num *= 2;
            int Max = (pow_d_2-1)*num + 1;
            pow_d_2 /= 2;

            for (int t = 0; t < Max; t += num) {
                for (int j = 0; j < num / 2; ++j) {
                    Complex xPOdd     = omega[m*j] * transform[t+num/2+j];
                    Complex prevTrans = transform[t+j];
                    transform[t+j]       = prevTrans + xPOdd;
                    transform[t+num/2+j] = prevTrans - xPOdd;

                    // omega[m*j] == - omega[m*(j+num/2)], 故这里可节约一次乘法
                    // 如果如下面这样展开了(NTT种就没法节约一次乘法), 实测运行时间翻倍
                    /*
                    Complex xPOdd     = omega[m*j] * transform[t+num/2+j];
                    Complex xPOdd1    = omega[m*(j+num/2)] * transform[t+num/2+j];
                    Complex prevTrans = transform[t+j];
                    transform[t+j]       = prevTrans + xPOdd;
                    transform[t+num/2+j] = prevTrans + xPOdd1;
                    */
                }
            }
        }
    }

    void rev_FFT(int k, int n, uint32_t * input, Complex *input1, int input_size, Complex *transform)
    {
        FFT(k, n, input, input1, input_size, transform);
        for (int i = n-1; i > n/2-1; i--) {
            Complex a = transform[n - i] / n;
            transform[n-i] = transform[i] / n;
            transform[i] = a;
        }
        transform[0] /= n;
    }
};


template <typename TZp>
class NttMul: public ProductFFT<TZp> {
    /*
     *  NTT 与 FFT 两种变换的框架是一样的。但是由于复根的特殊性可以利用以更快加速。因此这里各实现一下
     *  这里对100万位的乘法实测得：NTT比FFT慢，NTT速度是FFT的1/5
     * */

public:
    void FFT(int k, int n, uint32_t * input, TZp *input1, int input_size, TZp *transform)
    {
        TZp *omega = ProductFFT<TZp>::_omega_fft[k];
        uint32_t * pi_k = ProductFFT<TZp>::_pi_k[k];
        uint32_t * P = input;
        TZp * P1 = input1;

        if (P1 == NULL) {
            for (int t = 0; t < n-1; t+= 2) {
                uint32_t P_pi_k_t = pi_k[t] < input_size ? P[pi_k[t]] : 0;
                uint32_t P_pi_k_t_1 = pi_k[t+1] < input_size ? P[pi_k[t+1]] : 0;
                transform[t]   = P_pi_k_t;
                transform[t+1] = P_pi_k_t_1;
            }
        } else {
            for (int t = 0; t < n-1; t+= 2) {
                TZp P_pi_k_t = pi_k[t] < input_size ? P1[pi_k[t]] : 0;
                TZp P_pi_k_t_1 = pi_k[t+1] < input_size ? P1[pi_k[t+1]] : 0;

                transform[t]   = P_pi_k_t;
                transform[t+1] = P_pi_k_t_1;
            }
        }

        int m = n/1;
        int num = 1;
        int pow_d_2 = (int)round(pow(2, k-1));
        for (int d = k-1; d > -1; d--) {
            m /= 2;
            num *= 2;
            int Max = (pow_d_2-1)*num + 1;
            pow_d_2 /= 2;
            for (int t = 0; t < Max; t += num) {
                for (int j = 0; j < num / 2; ++j) {
                    TZp xPOdd     = omega[m*j] * transform[t+num/2+j];
                    TZp xPOdd1    = omega[m*(j+num/2)] * transform[t+num/2+j];
                    TZp prevTrans = transform[t+j];
                    transform[t+j]       = prevTrans + xPOdd;
                    transform[t+num/2+j] = prevTrans + xPOdd1;
                }
            }
        }
    }

    void rev_FFT(int k, int n, uint32_t * input, TZp *input1, int input_size, TZp *transform)
    {
        FFT(k, n, input, input1, input_size, transform);
        // 必须除以mod P 意义上的n. 这里先求倒数，再作乘积
        TZp n_reciprocal = TZp(1) / TZp(n);
        for (int i = n-1; i > n/2-1; i--) {
            TZp a = transform[n - i] * n_reciprocal;
            transform[n-i] = transform[i] * n_reciprocal;
            transform[i] = a;
        }
        transform[0] *= n_reciprocal;
    }

    // <<-- multi thread to speed up
    struct thread_do_fft_args_t {
        int k; int n2; uint32_t * aa; int aa_length; TZp *transform_aa;
        void *cls_obj;
    };
    static void * thread_do_fft(void * args) {
        thread_do_fft_args_t & a = *((thread_do_fft_args_t*)args);
        ((NttMul<TZp>*)a.cls_obj)->FFT(a.k, a.n2, a.aa, NULL, a.aa_length, a.transform_aa);
    }
    // -->>
    //
    // aa_eq_bb==1:算aa平方; ==0: 算aa*bb
    // calc_ntt_aa==1:算aa的NTT变换，结果存入transform_aa; ==0, 不算aa的NTT变换，直接用transform_aa作为aa的NTT变换
    // calc_ntt_bb: 同calc_ntt_aa
    void do_fast_ntt(int k, int n2, 
                     uint32_t * aa, int aa_length, TZp *transform_aa, int calc_ntt_aa,
                     uint32_t * bb, int bb_length, TZp *transform_bb, int calc_ntt_bb,
                     int aa_eq_bb, TZp *transform_temp, TZp *transform_out)
    {
        if ((calc_ntt_aa && calc_ntt_bb && (!aa_eq_bb)) && n2 > MIN2USE_THREAD) {
            // multi thread to speed up
            thread_do_fft_args_t args_0 = {k, n2, aa, aa_length, transform_aa, this};
            thread_do_fft_args_t args_1 = {k, n2, bb, bb_length, transform_bb, this};
            pthread_t thread_0;
            pthread_t thread_1;
            pthread_create(&thread_0, NULL, thread_do_fft, &args_0);
            pthread_create(&thread_1, NULL, thread_do_fft, &args_1);
            pthread_join(thread_0, NULL);
            pthread_join(thread_1, NULL);
        } else {
            if (calc_ntt_aa) {
                FFT(k, n2, aa, NULL, aa_length, transform_aa);
            }
            if (calc_ntt_bb && (!aa_eq_bb)) {
                FFT(k, n2, bb, NULL, bb_length, transform_bb);
            }
        }
        TZp *transform_cc = transform_temp;
        if (aa_eq_bb) {
            for (int i = 0; i < n2; i++) {
                transform_cc[i] = transform_aa[i] * transform_aa[i];
            }
        } else {
            for (int i = 0; i < n2; i++) {
                transform_cc[i] = transform_aa[i] * transform_bb[i];
            }
        }
        rev_FFT(k, n2, NULL, transform_cc, n2, transform_out);
    }
};


// FNT 没有精度问题，但有取模溢出问题。
// 为了解决溢出问题，可以做几组不同的FNT，再用中国剩余定理把溢出结果找回来
class FntMul: public NttMul<Zp0>, NttMul<Zp1>, NttMul<Zp2> {

    // 下面几个数字用于用中国剩余问题恢复overflow的数字
    // apfloat 和 libmpdec 是分三次FNT，且用的以下数字组合。这里也如此。
    // (实际上开始只用了一次FNT，总是有计算出错, 但不次次错。百思不得解。只好看了apfloat等的实现，才知道原因)
    // 之所以用下面这几个数是因为：
    //   1.从计算上, 这三个数字的选取使得整体计算的复杂度可控: 
    //        下面的P0/P1/P2都是31bit数字，可以说是uint32_t内最大的了(是否真的最大三个，不确定；但至少差不多是),
    //        31bit保证了每个FNT可以在uint64_t内轻松进行。P0*P1*p2是93bit的数字, P0*P1*p2*uint32_t 不会爆128bit，
    //        因此在3次FNT且用中国剩余定理恢复结果的时候，uint128_t(128位整数一般系统不支持，但是容易低成本模
    //        拟)可以轻松搞定。
    //   2. uint32_t 的最大十进制基是10^9，他们保证了这个基下，FNT能支持很长的长度，仍能用中国剩余定理给捞回来。
    //        允许的长度应该使FNT单个结果数字不爆P0*P1*P2。Length * 10^9 * Pi * Pi
    //        另外：如果再搞出个P3, 那么P0*P1*P2*P3当然可以支持更长长度数字的乘积。只是P0*P1*P2已经可以支撑
    //        很长的长度了。
    //
    // // Let P0, P1, P2 = 2113929217, 2013265921, 1811939329
    // // M0 = (P1*P2) * ((P1*P2)^(-1) mod P0)
    // // M1 = (P0*P2) * ((P0*P2)^(-1) mod P1)
    // // M2 = (P0*P1) * ((P0*P1)^(-1) mod P2)
    uint128_t M0; //(0x1d, 0x11e00082ec000093LU); // 0x1d11e00082ec000093
    uint128_t M1; //(0x18eabfd7LU, 0x1b97ff4a91ffff39LU); // 0x18eabfd71b97ff4a91ffff39
    uint128_t M2; //(0xc, 0x75600033e4000036LU); // 0xc75600033e4000036
    // // P012 = P0 * P1 * P2
    uint128_t P012; //(0x18eac000LU, 0xa2d8000162000001LU); // 0x18eac000a2d8000162000001

public:
    FntMul():M0(0x1d, 0x11e00082ec000093LU),
            M1(0x18eabfd7LU, 0x1b97ff4a91ffff39LU),
            M2(0xc, 0x75600033e4000036LU),
            P012(0x18eac000LU, 0xa2d8000162000001LU)
    {
    }
    void init(int max_len)
    {
        NttMul<Zp0>::init(max_len);
        NttMul<Zp1>::init(max_len);
        NttMul<Zp2>::init(max_len);
    }

    void calc_trans(int k, uint32_t * aa, int aa_length, void * trans_aa)
    {
        int n2 = int(round(pow(2, k)));
        Zp0 * buf = (Zp0*)trans_aa;
        Zp0 * trans0 = buf + n2 * 0;
        Zp0 * trans1 = buf + n2 * 1;
        Zp0 * trans2 = buf + n2 * 2;
        NttMul<Zp0>::FFT(k, n2, aa, NULL, aa_length, (Zp0*)trans0);
        NttMul<Zp1>::FFT(k, n2, aa, NULL, aa_length, (Zp1*)trans1);
        NttMul<Zp2>::FFT(k, n2, aa, NULL, aa_length, (Zp2*)trans2);
    }

    // <<-- multi thread to speed up
    struct thread_do_ntt_args_t {
        int k; int n2;
        uint32_t * aa; int aa_length; void * trans_aa; int calc_ntt_aa;
        uint32_t * bb; int bb_length; void * trans_bb; int calc_ntt_bb;
        int aa_eq_bb; void * trans_temp; void * fnt_out;
        void *cls_obj;
        int which_ntt;
    };
    static void * thread_do_ntt(void * args) {
        thread_do_ntt_args_t & a = *((thread_do_ntt_args_t*)args);
        if (a.which_ntt == 0) {
            ((NttMul<Zp0>*)a.cls_obj)->NttMul<Zp0>::do_fast_ntt(a.k,  a.n2, 
                                     a.aa, a.aa_length, (Zp0*)a.trans_aa, a.calc_ntt_aa,
                                     a.bb, a.bb_length, (Zp0*)a.trans_bb, a.calc_ntt_bb,
                                     a.aa_eq_bb, (Zp0*)a.trans_temp, (Zp0*)a.fnt_out);
        } else if (a.which_ntt == 1) {
            ((NttMul<Zp1>*)a.cls_obj)->NttMul<Zp1>::do_fast_ntt(a.k,  a.n2, 
                                     a.aa, a.aa_length, (Zp1*)a.trans_aa, a.calc_ntt_aa,
                                     a.bb, a.bb_length, (Zp1*)a.trans_bb, a.calc_ntt_bb,
                                     a.aa_eq_bb, (Zp1*)a.trans_temp, (Zp1*)a.fnt_out);
        } else if (a.which_ntt == 2) {
            ((NttMul<Zp2>*)a.cls_obj)->NttMul<Zp2>::do_fast_ntt(a.k,  a.n2, 
                                     a.aa, a.aa_length, (Zp2*)a.trans_aa, a.calc_ntt_aa,
                                     a.bb, a.bb_length, (Zp2*)a.trans_bb, a.calc_ntt_bb,
                                     a.aa_eq_bb, (Zp2*)a.trans_temp, (Zp2*)a.fnt_out);
        }
        return NULL;
    }
    // -->
    int fft_ntt(int k, int n2, 
                uint32_t * aa, int aa_length, void * trans_aa, int calc_ntt_aa,
                uint32_t * bb, int bb_length, void * trans_bb, int calc_ntt_bb,
                int aa_eq_bb,
                int radix, uint32_t* out, int &out_size)
    {
        int new_cnt = 3;
        int aa_offset = 0;
        int bb_offset = 0;
        int temp_offset = 0;
        if (trans_aa != NULL && trans_bb != NULL) {
            new_cnt += 3; // for temp
            temp_offset = 3;
        }
        if (trans_aa == NULL) {
            new_cnt += 3;
            calc_ntt_aa = 1;
            aa_offset = 3;
            temp_offset = 3;
        }
        if (trans_bb == NULL) {
            new_cnt += 3;
            calc_ntt_bb = 1;
            bb_offset = (trans_aa == NULL) ? 6 : 3;
            temp_offset = bb_offset;
        }

        Zp0 *buf = new Zp0[new_cnt*n2];

        Zp0 *fnt_out0 = buf + n2*0;
        Zp0 *fnt_out1 = buf + n2*1;
        Zp0 *fnt_out2 = buf + n2*2;

        Zp0 * trans_temp_0 = buf + n2*(temp_offset+0);
        Zp0 * trans_temp_1 = buf + n2*(temp_offset+1);
        Zp0 * trans_temp_2 = buf + n2*(temp_offset+2);

        Zp0 *trans_aa_0, *trans_aa_1, *trans_aa_2;
        Zp0 *trans_bb_0, *trans_bb_1, *trans_bb_2;

        if (trans_aa == NULL) {
            trans_aa_0= buf + n2*(aa_offset+0);
            trans_aa_1= buf + n2*(aa_offset+1);
            trans_aa_2= buf + n2*(aa_offset+2);
        } else {
            Zp0 * buf_1 = (Zp0*)trans_aa;
            trans_aa_0= buf_1 + n2*0;
            trans_aa_1= buf_1 + n2*1;
            trans_aa_2= buf_1 + n2*2;
        }

        if (trans_bb == NULL) {
            trans_bb_0= buf + n2*(bb_offset+0);
            trans_bb_1= buf + n2*(bb_offset+1);
            trans_bb_2= buf + n2*(bb_offset+2);
        } else {
            Zp0 * buf_1 = (Zp0*)trans_bb;
            trans_bb_0= buf_1 + n2*0;
            trans_bb_1= buf_1 + n2*1;
            trans_bb_2= buf_1 + n2*2;
        }

        int use_trd_max = MIN2USE_THREAD;
        if (n2 > use_trd_max) { // use multi thread
            pthread_t thread_0;
            pthread_t thread_1;
            pthread_t thread_2;
            thread_do_ntt_args_t args_0 = {k, n2, 
                                           aa, aa_length, trans_aa_0, calc_ntt_aa,
                                           bb, bb_length, trans_bb_0, calc_ntt_bb,
                                           aa_eq_bb, trans_temp_0, fnt_out0, (NttMul<Zp0>*)this, 0};
            thread_do_ntt_args_t args_1 = {k, n2, 
                                           aa, aa_length, trans_aa_1, calc_ntt_aa,
                                           bb, bb_length, trans_bb_1, calc_ntt_bb,
                                           aa_eq_bb, trans_temp_1, fnt_out1, (NttMul<Zp1>*)this, 1};
            thread_do_ntt_args_t args_2 = {k, n2, 
                                           aa, aa_length, trans_aa_2, calc_ntt_aa,
                                           bb, bb_length, trans_bb_2, calc_ntt_bb,
                                           aa_eq_bb, trans_temp_2, fnt_out2, (NttMul<Zp2>*)this, 2};
            pthread_create(&thread_0, NULL, thread_do_ntt, &args_0);
            pthread_create(&thread_1, NULL, thread_do_ntt, &args_1);
            pthread_create(&thread_2, NULL, thread_do_ntt, &args_2);
            pthread_join(thread_0, NULL);
            pthread_join(thread_1, NULL);
            pthread_join(thread_2, NULL);
        } else {
            NttMul<Zp0>::do_fast_ntt(k, n2, 
                                 aa, aa_length, (Zp0*)trans_aa_0, calc_ntt_aa,
                                 bb, bb_length, (Zp0*)trans_bb_0, calc_ntt_bb,
                                 aa_eq_bb, (Zp0*)trans_temp_0, (Zp0*)fnt_out0);
            NttMul<Zp1>::do_fast_ntt(k, n2,
                                 aa, aa_length, (Zp1*)trans_aa_1, calc_ntt_aa,
                                 bb, bb_length, (Zp1*)trans_bb_1, calc_ntt_bb,
                                 aa_eq_bb, (Zp1*)trans_temp_1, (Zp1*)fnt_out1);
            NttMul<Zp2>::do_fast_ntt(k, n2, 
                                 aa, aa_length, (Zp2*)trans_aa_2, calc_ntt_aa,
                                 bb, bb_length, (Zp2*)trans_bb_2, calc_ntt_bb,
                                 aa_eq_bb, (Zp2*)trans_temp_2, (Zp2*)fnt_out2);
        }

        int max_res_len = aa_length + bb_length;

        if (0) { // 只做一个NTT
            int remain = 0;
            for (int i = 0; i < max_res_len; ++i) {
                int o = fnt_out1[i].to_int() + remain;
                if (o >= radix) {
                    out[i] = o % radix;
                    remain = o / radix;
                } else {
                    out[i] = o;
                    remain = 0;
                }
            }
        } else {
            uint128_t remain(0);
            for (int i = 0; i < max_res_len; ++i) {
                uint128_t d = (M0 * fnt_out0[i].n + M1 * fnt_out1[i].n + M2 * fnt_out2[i].n + remain) % P012;

                //uint128_t out_res = d % radix;
                //out[i] = out_res.LOWER;
                //remain = d / radix;

                remain = d.div_mod(radix, out[i]);
            }
        }
        out_size = out[max_res_len-1] != 0 ? max_res_len : max_res_len - 1;
        delete [] buf;
        return 0;
    }

    int fast_prod(uint32_t * aa, int aa_length, void *trans_aa, int calc_ntt_aa,
                  uint32_t * bb, int bb_length, void *trans_bb, int calc_ntt_bb,
                  int aa_eq_bb,
                  int radix, uint32_t * out, int &out_size) 
    {
        //printf("cpp_mul_time begin: aa_size=%d, bb_size=%d\n", aa_length, bb_length);
        //struct timeval tpstart, tpend;
        //gettimeofday(&tpstart,NULL);
        // FNT 类是菱形继承, 因此要访问最上基类的方法，会不知调用Zp<i>中哪个。调用哪个都可以
        int ret = NttMul<Zp0>::fast_prod(aa, aa_length, trans_aa, calc_ntt_aa,
                                         bb, bb_length, trans_bb, calc_ntt_bb,
                                         aa_eq_bb, radix, out, out_size);
        //gettimeofday(&tpend,NULL);
        //double tm = ((tpend.tv_sec-tpstart.tv_sec)*1000000+(tpend.tv_usec-tpstart.tv_usec));
        //printf("cpp_mul_time end:tm=%.4f, out_size=%d\n", tm/1000000, out_size);
        return ret;
    }
};

#endif
