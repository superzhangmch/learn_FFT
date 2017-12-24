#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include<time.h>

class ProductFFT {

    int _max_result_len;
    int _max_k;
    double (*_omega[32])[2];
    int * _pi_k[32];

public:
    ProductFFT() {
        _max_result_len = _max_k = 0;
        for (int i = 0; i < 32; ++i) {
            _omega[i] = NULL;
            _pi_k[i] = NULL;
        }
    }
    ~ProductFFT() {
        for (int i = 0; i < 32; ++i) {
            if (_omega[i]) {
                delete [] _omega[i];
                _omega[i] = NULL;
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

        for (int i = 1; i < _max_k + 1; ++i) {
            int n = (int)pow(2, i);
            _omega[i] = new double[n][2];
            _pi_k[i] = new int[n];
            gen_omega_pi_k(i, _omega[i], _pi_k[i]);
        }
    }


    // reverse the k-bit representation of t
    int rev_bin(int t, int k)
    {
        int ret = 0;
        for (int i = 0; i < k; ++i) {
            int a = (1<<i) & t;
            if (a != 0) {
                ret |= (1<<(k - i-1));
            }
        }
        return ret;
    }

    void gen_omega_pi_k(int k, double (*w_ret)[2], int * pi_k_ret)
    {
        int n = (int)pow(2, k);
        double x, y;

        // if (k == 5) printf("%d--\n", k);
        for (int i = 0; i < n; ++i) {
            w_ret[i][0] = cos(2.L* i * 3.141592653589793238462643383279502884197L / n);
            w_ret[i][1] = sin(2.L* i * 3.141592653589793238462643383279502884197L / n);
            // if (k == 5) printf("%d->%.5f+%.5fi\n", i, w_ret[i][0], w_ret[i][1]);
        }

        for (int i = 0; i < n; ++i) {
            pi_k_ret[i] = rev_bin(i, k);
            // if (k == 5) printf("%d->%d\n", i, pi_k_ret[i]);
        }
    }

    void FFT(int k, int n, int * input, double (*input1)[2], int input_size, double (*transform)[2])
    {
        double (*omega)[2] = _omega[k];
        int * pi_k = _pi_k[k];
        int * P = input;
        double (* P1)[2] = input1;

        // n = 2 ** k
        // printf("nn %d %d\n", k, n);
        if (P1 == NULL) {
            for (int t = 0; t < n-1; t+= 2) {
                double P_pi_k_t = pi_k[t] < input_size ? P[pi_k[t]] : 0.;
                double P_pi_k_t_1 = pi_k[t+1] < input_size ? P[pi_k[t+1]] : 0.;
                transform[t][0]   = P_pi_k_t + P_pi_k_t_1;
                transform[t][1]   = 0.;
                transform[t+1][0] = P_pi_k_t - P_pi_k_t_1;
                transform[t+1][1] = 0.;
            }
        } else {
            for (int t = 0; t < n-1; t+= 2) {
                double P_pi_k_t_0 = pi_k[t] < input_size ? P1[pi_k[t]][0] : 0.;
                double P_pi_k_t_1 = pi_k[t] < input_size ? P1[pi_k[t]][1] : 0.;

                double P_pi_k_t_1_0 = pi_k[t+1] < input_size ? P1[pi_k[t+1]][0] : 0.;
                double P_pi_k_t_1_1 = pi_k[t+1] < input_size ? P1[pi_k[t+1]][1] : 0.;

                transform[t][0] = P_pi_k_t_0 + P_pi_k_t_1_0;
                transform[t][1] = P_pi_k_t_1 + P_pi_k_t_1_1;
                transform[t+1][0] = P_pi_k_t_0 - P_pi_k_t_1_0;
                transform[t+1][1] = P_pi_k_t_1 - P_pi_k_t_1_1;
            }
        }
        //if (input)
        //for (int i = 0; i< input_size; i++) {
        //    printf("%d ", input[i]);
        //}
        //printf("\n");
        //for (int i = 0; i< n; i++) {
        //    printf("%d ", pi_k[i]);
        //}
        //printf("\n");
        //for (int i = 0; i< n; i++) {
        //    printf("%d->%f %f\n", i, transform[i][0], transform[i][1]);
        //}
        //printf("\n");
        int m = n/2;
        int num = 2;
        int pow_d_2 = (int)round(pow(2, k-2));
        for (int d = k-2; d > -1; d--) {
            m /= 2;
            num *= 2;
            int Max = (pow_d_2-1)*num + 1;
            pow_d_2 /= 2;
            //printf("xx %d\n", Max);
            for (int t = 0; t < Max; t += num) {
                for (int j = 0; j < num / 2; ++j) {
                    double a = omega[m*j][0];
                    double b = omega[m*j][1];
                    double c = transform[t+num/2+j][0];
                    double d = transform[t+num/2+j][1];

                    double xPOdd[2];
                    xPOdd[0] = a*c-b*d;
                    xPOdd[1] = b*c+a*d;

                    double prevTrans[2];
                    prevTrans[0] = transform[t+j][0];
                    prevTrans[1] = transform[t+j][1];

                    transform[t+j][0] = prevTrans[0]+xPOdd[0];
                    transform[t+j][1] = prevTrans[1]+xPOdd[1];

                    transform[t+num/2+j][0] = prevTrans[0]-xPOdd[0];
                    transform[t+num/2+j][1] = prevTrans[1]-xPOdd[1];
                }
            }
        }
        //for (int i = 0; i< n; i++) {
        //    printf("%d->%f %f\n", i, transform[i][0], transform[i][1]);
        //}
    }

    void rev_FFT(int k, int n, int * input, double (*input1)[2], int input_size, double (*transform)[2])
    {
        FFT(k, n, input, input1, input_size, transform);
        //for (int i = 0; i< n; i++) {
        //    printf("---%d->%f %f\n", i, transform[i][0], transform[i][1]);
        //}
        for (int i = n-1; i > n/2-1; i--) {
            double a = transform[n - i][0] / n;
            double b = transform[n - i][1] / n;

            transform[n-i][0] = transform[i][0] / n;
            transform[n-i][1] = transform[i][1] / n;

            transform[i][0] = a;
            transform[i][1] = b;
        }
        transform[0][0] /= n;
        transform[0][1] /= n;
        //for (int i = 0; i< n; i++) {
        //    printf("---%d->%f %f\n", i, transform[i][0], transform[i][1]);
        //}
    }

    int fast_prod(int * aa, int aa_length, int * bb, int bb_length, int radix, int * out, int &out_size) 
    {
        // aa or bb == 0
        if (aa_length == 0 || bb_length == 0) {
            out_size = 0;
            return 0;
        }
        // aa or bb is small
        if (aa_length == 1 || bb_length == 1) {
            int small_val, big_size, *big;
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
                int res = big[i] * small_val + remain;
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
        int max_len = aa_length > bb_length ? aa_length : bb_length;
        int k = int(ceil(log(max_len) / log(2))) + 1;
        int n2 = int(round(pow(2, k)));
        if (k > _max_k) {
            out_size = 0;
            return -1;
        }
        //printf("--%d %d %d %d %d\n", k, n2, aa_length, bb_length, max_len);
        double (*buf)[2] = new double[2*n2][2];
        double (*transform_aa)[2] = buf;
        double (*transform_bb)[2] = buf + n2;
        double (*transform_cc)[2] = buf + n2;
        double (*transform_dd)[2] = buf;
        FFT(k, n2, aa, NULL, aa_length, transform_aa);
        //printf("xxxx\n\n");
        FFT(k, n2, bb, NULL, bb_length, transform_bb);

        for (int i = 0; i < n2; i++) {
            double a = transform_aa[i][0];
            double b = transform_aa[i][1];
            double c = transform_bb[i][0];
            double d = transform_bb[i][1];
            transform_cc[i][0] = a*c-b*d;
            transform_cc[i][1] = b*c+a*d;
        }

        rev_FFT(k, n2, NULL, transform_cc, n2, transform_dd);

        int max_res_len = aa_length + bb_length;
        int remain = 0;
        for (int i = 0; i < max_res_len; ++i) {
            int o = int(round(transform_dd[i][0])) + remain;
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
};

int qi9_check(int * aa, int aa_s, int radix)
{
    int qi9_a=0;
    int show = 0;
    if (show)printf("\n");
    for (int i = aa_s - 1; i>=0; --i) {
        qi9_a += aa[i];
        qi9_a = qi9_a % 9;
        if(show)printf("%d", aa[i]);
    }
    if(show)printf("\n");
    return qi9_a;
}
int main()
{
    ProductFFT fft;
    int max_digit = 2000000;
    fft.init(max_digit * 2);
    srand(time(NULL));

    int *aa = new int[max_digit];
    int *bb = new int[max_digit];
    int *cc = new int[max_digit * 2];
    int aa_s, bb_s, cc_s = 0;
    int radix = 10; // for radix=100, max_digit = 200w not ok

    struct timeval tpstart, tpend;
    aa_s = bb_s = max_digit/2;

    for (int k = 0; k < 100; ++k) {
        for (int i = 0; i< aa_s; ++i) {
            aa[i] = rand() % radix;
            bb[i] = rand() % radix;
            if (i == aa_s -1) {
                if (aa[i] == 0) aa[i]++;
                if (bb[i] == 0) bb[i]++;
            }
        }


        gettimeofday(&tpstart,NULL);
        fft.fast_prod(aa, aa_s, bb, bb_s, radix, cc, cc_s); 
        gettimeofday(&tpend,NULL);

        double tm = ((tpend.tv_sec-tpstart.tv_sec)*1000000+(tpend.tv_usec-tpstart.tv_usec));
        int qi9_a=0, qi9_b=0, qi9_c = 0;
        qi9_a = qi9_check(aa, aa_s, radix);
        qi9_b = qi9_check(bb, bb_s, radix);
        qi9_c = qi9_check(cc, cc_s, radix);
        printf("time=%.4f, check=%d*%d=%d ~ %d\n", tm/1000000, qi9_a, qi9_b, (qi9_b * qi9_a) % 9, qi9_c);
        if ((qi9_b * qi9_a) % 9 != qi9_c)
        {
            printf("not match\n");
            return -1;
        }
        for (int i = cc_s - 1; i >= 0; --i) {
            // printf("%d ", cc[i]);
        }
        //printf("\n%d\n", cc_s);
    }
    return 0;
}
