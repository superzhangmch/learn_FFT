#encoding:utf8
import math
import sys

class ProductFFT(object):
    """
    快速傅里叶变换与快速数论变换作乘法
    """

    def __init__(self, max_len, use_fft=True):
        self.use_fft = use_fft
        self.p = 1004535809
        self.g = 3
        # for primitive root of prime: see http://blog.miskcoo.com/2014/07/fft-prime-table
        k = int(math.ceil(math.log(max_len) / math.log(2)))
        self._max_result_len = 2 * 2 **k
        self._max_k = k + 1
        self._omega = {}
        self._omega_ntt = {}
        self._omega_reciprocal = {}
        self._pi_k = {}

        print "init begin"
        for i in xrange(1, self._max_k + 1, 1):
            self.prepare_data(i)
        print "init ok"

    def fast_m(self, a, idx, p):
        """
        calc: (a ** idx) % p
        """
        if idx < 10:
            return a**idx % p
        elif idx % 2 == 0:
            return (self.fast_m(a, idx/2, p)**2) % p
        else:
            return (self.fast_m(a, idx/2, p)**2 * a) % p

    def prepare_data(self, k):
        n = 2 ** k
        self._omega[k] = []
        self._omega_ntt[k] = []
        self._omega_reciprocal[k] = 0

        if self.use_fft:
            w = math.cos(2*math.pi / n) + 1j * math.sin(2*math.pi / n)
            self._omega[k] = [w ** i for i in xrange(n)]
        else: # NTT
            p = self.p
            g = self.g
            skip = p/n
            self._omega_ntt[k] = [int(self.fast_m(g, i*skip, p)) for i in xrange(n)]
            # 计算 n 在 mod p 下的倒数
            # 根据欧拉定理，p为素数，则任意n有 n^(p-1) = 1 mod p, 于是 n * (n^(p-2)) == 1的n倒数为n^(p-2)
            self._omega_reciprocal[k] = self.fast_m(n, p-2, p)

        def rev_bin(t, k):
            """
            对 k-bits 表示的二进制整数t, 进行二进制逆序，染回逆序后的整数
            """
            b = bin(t)[2:]
            b1 = ("0"*(k - len(b))) + b
            b2 = b1[::-1]
            t1 = int("0b"+b2, 2)
            #print t, b, b1, b2, t1
            return t1
        self._pi_k[k] = [rev_bin(t, k) for t in xrange(n)]

    def FFT(self, omega, omega_ntt, pi_k, P, k):
        """ O(n*log(n)). 从 <<computer algorithms: introduction to design and analysis>> 改来，兼容FFT与NTT
         FFT/rev_FFT就是个矩阵运算，所以可以按矩阵计算的方式进行, 但复杂度为O(n^2)
         另外numpy.fft.fft/numpy.fft.ifft 是numpy的FFT实现
        """
        if not self.use_fft: # FTT
            omega = omega_ntt
        transform = [0] * len(P)

        n = 2 ** k
        for t in xrange(0, n-1, 2):
            transform[t]   = P[pi_k[t]]
            transform[t+1] = P[pi_k[t+1]]
        m = n/1
        num = 1
        for d in xrange(k-1, -1, -1):
            m /= 2
            num *= 2
            for t in xrange(0, (2**d-1)*num + 1, num):
                for j in xrange(num/2):

                    xPOdd = transform[t+j+num/2]
                    prevTrans= transform[t+j]
                    # NOTE: for FFT, omega[m*(j+num/2)] == -omega[m*j]
                    transform[t+j] = prevTrans + omega[m*j] * xPOdd
                    transform[t+j+num/2] = prevTrans + omega[m*(j+num/2)] * xPOdd

                    if not self.use_fft: # FFT
                        transform[t+j] %= self.p
                        transform[t+j+num/2] %= self.p

        return transform
    
    def FFT_only(self, omega, pi_k, P, k): 
        """ O(n*log(n)). from <<computer algorithms: introduction to design and analysis>>"""
        transform = [0.] * len(P)
        n = 2 ** k
        for t in xrange(0, n-1, 2): 
            transform[t]   = P[pi_k[t]]  + P[pi_k[t+1]]
            transform[t+1] = P[pi_k[t+1]] - P[pi_k[t+1]]
        m = n/2
        num = 2
        for d in xrange(k-1, -1, -1):
            m /= 2
            num *= 2
            for t in xrange(0, (2**d-1)*num + 1, num):
                for j in xrange(num/2):

                    xPOdd = transform[t+j+num/2] * omega[m*j]
                    prevTrans= transform[t+j]

                    transform[t+j] = prevTrans + xPOdd
                    transform[t+j+num/2] = prevTrans - xPOdd

        return transform
    
    def rev_FFT(self, omega, omega_ntt, pi_k, P, k):
        """ from <<computer algorithms: introduction to design and analysis>>"""
        n = 2 ** k
        transform = self.FFT(omega, omega_ntt, pi_k, P, k)
        if self.use_fft:
            for i in xrange(n-1, n/2 - 1, -1):
                transform[i], transform[n-i] = transform[n - i] / n, transform[i] / n
            transform[0] /= n
        else:
            # rev-NTT 需要mod p方式除以n, 也就是乘以n的mod p 倒数
            reciprocal = self._omega_reciprocal[k]
            for i in xrange(n-1, n/2 - 1, -1):
                transform[i], transform[n-i] = (transform[n - i] * reciprocal) % self.p, (transform[i] *reciprocal) % self.p
            transform[0] = (transform[0] *reciprocal) % self.p
        return transform

    def fast_prod(self, aa, bb, result, radix):
        # ===================
        # aa or bb == 0
        # ===================
        if aa.length == 0 or bb.length == 0:
            return result.set_val([], 0, 0)

        # ===================
        # aa or bb is small
        # ===================
        if aa.length == 1 or bb.length == 1:
            if aa.length == 1:
                small_val, big = aa.val[0], bb
            else:
                small_val, big = bb.val[0], aa
            out = [0] * (big.length + 1)
            out_len = big.length
            remain = 0
            for i in xrange(big.length):
                res = big.val[i] * small_val + remain
                #print i, '->', big.val[i], small_val, remain, res
                if res >= radix:
                    out[i] = res % radix
                    remain = res / radix
                else:
                    out[i] = res
                    remain = 0
            #print ""
            if remain:
                #print 'xxxxxxxxxxxxxxxxxxxxx', small_val
                out[big.length] = remain
                out_len = big.length + 1
            is_neg = False if aa.is_neg == bb.is_neg else True
            result.set_val(out, out_len, aa.exp_idx + bb.exp_idx, is_neg=is_neg)
            #print 'outoutout', out, out_len
            assert out[out_len-1] != 0
            return 

        # =========================
        # both aa and bb are big
        # =========================
        def get_len(len_aa, len_bb):
            k =  int(math.ceil(math.log(max(len_aa, len_bb)) / math.log(2)))
            return 2 * 2 **k, k + 1
        n2, k = get_len(aa.length, bb.length)

        aa.extend_len(n2)
        bb.extend_len(n2)

        w_arr, w_ntt_arr, pi_k_arr = self._omega[k], self._omega_ntt[k], self._pi_k[k]

        # === FFT
        aa2 = self.FFT(w_arr, w_ntt_arr, pi_k_arr, aa.val, k)
        bb2 = self.FFT(w_arr, w_ntt_arr, pi_k_arr, bb.val, k)

        # ==== point-wise multiplication
        cc = [0] * n2
        if self.use_fft:
            for i in xrange(n2):
                cc[i] = aa2[i]*bb2[i]
        else: # NTT
            for i in xrange(n2):
                cc[i] = (aa2[i]*bb2[i]) % self.p

        # ==== rev FFT
        out = self.rev_FFT(w_arr, w_ntt_arr, pi_k_arr, cc, k)

        # ==== 
        max_res_len = aa.length + bb.length
        out = [int(round(o.real)) for o in out]
        out1 = [0] * n2
        remain = 0
        zero_cnt = 0
        not_count = 0
        for i in xrange(len(out)):
            o = out[i]
            o += remain
            if o >= radix:
                out1[i] = o % radix
                remain = o / radix
            else:
                out1[i] = o
                remain = 0

            if (not not_count) and out1[i] == 0:
                zero_cnt += 1
            if out1[i] != 0:
                not_count = 1

            if i >= max_res_len:
                break
        res_len = max_res_len if out1[max_res_len - 1] != 0 else max_res_len - 1
        is_neg = False if aa.is_neg == bb.is_neg else True
        result.set_val(out1, res_len, aa.exp_idx + bb.exp_idx, is_neg=is_neg)
        assert out1[res_len-1] != 0
        #if not_count == 0:
        #    result.set_val(out1, res_len, aa.exp_idx + bb.exp_idx, is_neg=is_neg)
        #else:
        #    result.set_val(out1[zero_cnt:], res_len-zero_cnt, aa.exp_idx + bb.exp_idx+zero_cnt, is_neg=is_neg)
 
