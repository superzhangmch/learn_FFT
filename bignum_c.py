#encoding: utf8
import math
import sys
import time
import random
import ctypes
from ctypes import cdll, c_int, byref, pointer

"""
大数计算：包含加减乘除开平方等
C/C++ 实现背后的基本+-*运算，其他运算由python来组合生成
"""

class BigNum(object):

    max_digit = 100
    extra_digit = 10
    radix = 10
    bignum_so = None
    use_fft = False

    @staticmethod
    def init(max_digit, radix, use_fft, so_file="./fft_ntt.so"):
        BigNum.max_digit = max_digit
        BigNum.radix = radix
        BigNum.use_fft = use_fft
        
        BigNum.bignum_so = cdll.LoadLibrary(so_file)
        fft_max_len = max_digit * 2 + BigNum.extra_digit * 2
        if use_fft:
            BigNum.bignum_so.fft_init(fft_max_len)
        else:
            BigNum.bignum_so.fnt_init(fft_max_len)

    def dec2other_radix(self, num, radix):
        out = []
        while num:
            remain = num % radix
            num = num / radix
            out.append(int(remain))
        return out

    def __init__(self, val, length=None, exp_idx=0, is_neg=False, start_from=0, is_rev=True):
        radix = self.radix
        self.debug = 0
        self.precision = 0
        self.start_from = 0

        if type(val) == type(1) or type(val) == long:
            if val == 0:
                self.val = c_int(0)
                self.length = 0
            else:
                v = self.dec2other_radix(abs(val), radix)
                self.length = len(v)
                self.val = (c_int*len(v))(*v)
            self.exp_idx = exp_idx
            if val == 0:
                self.is_neg = True
            else:
                self.is_neg = True if val / abs(val) == -1 else False
                if not self.is_neg and is_neg:
                    self.is_neg = is_neg
            return
        elif isinstance(val, ctypes.Array):
            self.val = val
            self.length = length
            self.exp_idx = exp_idx
            self.is_neg = is_neg
            self.start_from = start_from
            return

        if type(val) == type(""):
            val = [int(i) for i in val]

        if length is None:
            length = len(val)
        if not is_rev:
            val = val[::-1]

        self.set_val(val, length, exp_idx, is_neg, start_from)

    def set_val(self, val, length, exp_idx, is_neg=False, start_from=0):
        assert len(val) >= length
        self.val = (c_int * length)(*val)
        self.start_from = start_from
        self.length = length
        self.exp_idx = exp_idx
        self.is_neg = is_neg

        max_len = length
        for i in xrange(start_from + max_len-1, start_from-1, -1):
            if val[i] == 0:
                self.length -= 1
            else:
                break

    def cut_len(self, size):
        exp_add = (self.length - size)
        if exp_add > 0:
            old_len = self.length
            self.length = size
            self.exp_idx += exp_add
            self.start_from += exp_add

    def get_normal_val(self):
        aa = list(self.val) + []
        aa = aa[self.start_from:self.length+self.start_from]
        aa.reverse()
        return aa

    def get_all_numbers(self):
        c = int(round(math.log(self.radix) / math.log(10)))
        if self.radix == 10:
            res = "".join([str(i) for i in self.get_normal_val()])
        else:
            prefix = "0" * c
            res = "".join([(prefix + str(i))[-c:] for i in self.get_normal_val()])
        return res

    def get_sum_mod9(self, val=None):
        if not val:
            val = list(self.val)
            return sum(val[self.start_from: self.length+self.start_from]) % 9
        else:
            return sum(val) % 9

    def get_string_val(self, show_all=False, show_real_precision=True, num_only=False):
        if self.length == 0:
            return "0"
        if num_only == True:
            show_all = True
        c = int(round(math.log(self.radix) / math.log(10)))
        if self.precision == 0:
            precision = self.length
        else:
            precision = self.precision
        arr_val = self.get_normal_val()
        if show_real_precision:
            arr_precision = arr_val[:precision]
        else:
            arr_precision = arr_val
        if self.radix == 10:
            res = "".join([str(i) for i in arr_precision])
        else:
            prefix = "0" * c
            res = "".join([(prefix + str(i))[-c:] for i in arr_precision])
        if show_all:
            pass
        elif len(res)> 600:
            res = res[:500] + "..." + res[-100:]
        if self.radix == 10:
            digit_num = self.exp_idx + self.length
        elif arr_val:
            digit_num = len(str(arr_val[0])) + c*(self.exp_idx + self.length - 1)
        else:
            digit_num = 0
        sign = "-" if self.is_neg else ""
        if num_only:
            return sign + "0."+res.lstrip("0") + "E" + str(digit_num) 
        res = sign + ""+res.lstrip("0") + "*E" + str(digit_num) 
        res +=  "(len=%d|start:%d|length=%d|exp_idx:%d|precision:%d|sum_mod9:%d)" % (digit_num, self.start_from, self.length, self.exp_idx, precision, self.get_sum_mod9(arr_val))
        return res

    def __str__(self):
        return self.get_string_val(False)

    # =====================

    def factorial(self, n):
        assert n >= 0
        if n == 0 or n == 1:
            return 1
            
    def pow(self, n):
        assert n >= 0
        if n == 0:
            return BigNum(1)
        elif n == 1:
            return self
        pow_n_half_pow2 = self.pow(n / 2).pow2()
        if n % 2 == 0:
            return pow_n_half_pow2
        else:
            return pow_n_half_pow2 * self

    def pow2(self):
        if self.length:
            assert self.val[self.start_from + self.length - 1] != 0

        if self.length == 0:
            return BigNum(0)
        elif self.length == 1:
            res_c_int = (c_int * (self.length + 1))()
        else:
            k =  int(math.ceil(math.log(self.length) / math.log(2)))
            cc_s = 2 **(k+1)
            res_c_int = (c_int * cc_s)()

        if self.use_fft:
            func = self.bignum_so.do_fft_pow2
        else:
            func = self.bignum_so.do_fnt_pow2

        #tm = time.time()
        res_digit = func(byref(self.val), self.start_from, self.length,
                    self.radix, byref(res_c_int))
        #print "tm_mul %.4f" % (time.time() - tm)
        exp_idx = self.exp_idx * 2
        res = BigNum(res_c_int, length=res_digit, exp_idx=exp_idx, is_neg=False)
        
        res.cut_len(self.max_digit + self.extra_digit)
        if res.length:
            assert res.val[res.start_from + res.length - 1] != 0
        return res

    def __mul__(self, num1):
        if self.length:
            assert self.val[self.start_from + self.length - 1] != 0
        if num1.length:
            assert num1.val[num1.start_from + num1.length - 1] != 0

        if self.length == 0 or num1.length == 0:
            return BigNum(0)
        elif self.length == 1 or num1.length == 1:
            res_c_int = (c_int * (max(self.length, num1.length) + 1))()
        else:
            k =  int(math.ceil(math.log(max(self.length, num1.length)) / math.log(2)))
            cc_s = 2 **(k+1)
            res_c_int = (c_int * cc_s)()

        if self.use_fft:
            func = self.bignum_so.do_fft_1
        else:
            func = self.bignum_so.do_fnt_1

        #tm = time.time()
        res_digit = func(byref(self.val), self.start_from, self.length,
                    byref(num1.val), num1.start_from, num1.length, 
                    self.radix, byref(res_c_int))
        #print "tm_mul %.4f" % (time.time() - tm)
        is_neg = False if self.is_neg == num1.is_neg else True
        exp_idx = self.exp_idx + num1.exp_idx
        res = BigNum(res_c_int, length=res_digit, exp_idx=exp_idx, is_neg=is_neg)
        
        res.cut_len(self.max_digit + self.extra_digit)
        if res.length:
            assert res.val[res.start_from + res.length - 1] != 0
        return res


    def __add__(self, num1):
        if self.length:
            assert self.val[self.start_from + self.length - 1] != 0
        if num1.length:
            assert num1.val[num1.start_from + num1.length - 1] != 0

        #tm = time.time()
        if self.is_neg != num1.is_neg:
            # sub: (+) + (-) || (-) + (+)
            res = self.do_sub(num1)
        # add: (+) + (+) || (-) + (-)
        else:
            res = self.do_add(num1)
        #print "tm_add %.4f" % (time.time() - tm)
        if res.length:
            #print res, res.start_from + res.length - 1
            assert res.val[res.start_from + res.length - 1] != 0
        res.cut_len(self.max_digit + self.extra_digit)
        return res

    def __sub__(self, num1):
        if self.length:
            assert self.val[self.start_from + self.length - 1] != 0
        if num1.length:
            assert num1.val[num1.start_from + num1.length - 1] != 0

        #tm = time.time()
        if self.is_neg != num1.is_neg:
            # add: (+) - (-) => (+) || (-) - (+) => (-)
            res = self.do_add(num1)
        else: # sub: (-) - (-) || (+)-(+)
            res = self.do_sub(num1)
        #print "tm_add %.4f" % (time.time() - tm)
        if res.length:
            assert res.val[res.start_from + res.length - 1] != 0
        res.cut_len(self.max_digit + self.extra_digit)
        return res

    def do_add(self, num1):
        """ same sign, result.sign = self.sign """
        if self.length == 0:
            if self.is_neg == num1.is_neg:
                return num1
            else:
                res = (c_int*num1.length)()
                self.bignum_so.big_copy(byref(num1.val), num1.start_from, num1.length, byref(res))
                return BigNum(res, num1.length, num1.exp_idx, is_neg=self.is_neg)
        if num1.length == 0:
            return self

        if self.exp_idx + self.length > num1.exp_idx + num1.length:
            big, small = self, num1
        else:
            big, small = num1, self
        big_padding_zero_num = big.exp_idx - small.exp_idx

        res = (c_int*(1 + big.length+abs(big_padding_zero_num)))()
        res_size = c_int()
        self.bignum_so.big_add_1(byref(big.val), big.start_from, big.length, 
                                 byref(small.val), small.start_from, small.length, 
                                 byref(res), byref(res_size), big_padding_zero_num, self.radix)
        exp_idx = min(self.exp_idx, num1.exp_idx)
        return BigNum(res, res_size.value, exp_idx, is_neg=self.is_neg)

    def do_sub(self, num1):
        """ self.sign != num1.sign: result.sign = self.sign if (|self| > |num1| else (!self.sign)) """
        if self.length == 0:
            res_neg = not self.is_neg
            if res_neg == num1.is_neg:
                return num1
            else:
                res = (c_int*num1.length)()
                self.bignum_so.big_copy(byref(num1.val), num1.start_from, num1.length, byref(res))
                return BigNum(res, num1.length, num1.exp_idx, is_neg=not self.is_neg)
        if num1.length == 0:
            return self

        if self.exp_idx + self.length > num1.exp_idx + num1.length:
            big, small = self, num1
        else:
            big, small = num1, self
        big_padding_zero_num = big.exp_idx - small.exp_idx

        res = (c_int*(big.length+abs(big_padding_zero_num)))()
        res_size = c_int()
        big_gt_small = c_int()
        self.bignum_so.big_sub_1(byref(big.val), big.start_from, big.length, 
                                 byref(small.val), small.start_from, small.length, 
                                 byref(res), byref(res_size), byref(big_gt_small),
                                 big_padding_zero_num, self.radix)

        exp_idx = min(self.exp_idx, num1.exp_idx)
        if big == self:
            if big_gt_small: # |self| >= |num1|
                res_neg = self.is_neg
            else:            # |self| < |num1|
                res_neg = not self.is_neg
        else: # big == num1
            if big_gt_small: # |self| < |num1|
                res_neg = not self.is_neg
            else:            # |self| >= |num1|
                res_neg = self.is_neg
        return BigNum(res, res_size.value, exp_idx, is_neg=res_neg)

    def __div__(self, num):
        if self.length:
            assert self.val[self.start_from + self.length - 1] != 0
        if num.length:
            assert num.val[num.start_from + num.length - 1] != 0
        assert num.length > 0

        if self.length == 0:
            return BigNum(0)

        if num.length == 1:
            div_by = num.val[0]
            out_size = c_int(0)
            req_out_size = self.max_digit + self.extra_digit
            res = (c_int * (self.max_digit + self.extra_digit + 1))()
            #tm = time.time()
            self.bignum_so.big_div_by_1digit_1(byref(self.val), self.start_from, self.length, div_by, 
                             byref(res), req_out_size, byref(out_size), 
                             self.radix)
            #print "tm_div %.4f" % (time.time() - tm)
            res_exp_idx = self.length if (self.val[self.start_from + self.length - 1] >= div_by) else self.length - 1
            exp_idx = self.exp_idx - (out_size.value - res_exp_idx)
            is_neg = True if self.is_neg != num.is_neg else False
            ret_val = BigNum(res, start_from=req_out_size-out_size.value, 
                             length=out_size.value, exp_idx=exp_idx, is_neg=is_neg)
            ret_val.cut_len(self.max_digit + self.extra_digit)
            return ret_val

        # 除法转化为乘法。先求倒数。
        res = self * num.reciprocal('div_inv')
        if res.length:
            assert res.val[res.start_from + res.length - 1] != 0
        res.cut_len(self.max_digit + self.extra_digit)
        return res

    def reciprocal(self, name=""):
        b = self
        assert b.length > 1

        # print b.val, b.length, 11111111111111
        if b.val[b.start_from + b.length - 1] == 1:
            val = [1]
        else:
            val = [b.val[b.start_from + b.length - 1]]
        init_expidx = -b.length - b.exp_idx - len(val)
        # 上面是优选牛顿迭代的初始值。对待求解数，取其最高位作估计，作为初始值。
        # 注意初始值需要是一个分母是radix的整幂次的数，这样才能规避复杂除法。
        
        def f(x):
            '''
            欲求1/b, 令h(x) =1/x-b. 用牛顿法求出x≈1/b.
            '''
            return x*(BigNum(2) - b*x)
        aa = BigNum(val, 1, init_expidx)
        return self.Newton_method(f, aa, name='inv' if not name else name)

    def sqrt(self):
        b = self
        if b.length == 0:
            return BigNum([], 0, 0)
        assert b.is_neg == False

        # print b.val, b.length, 11111111111111222222222222222
        #----
        first_num = b.val[b.start_from + b.length - 1]
        second_num = 0 if b.length == 1 else b.val[b.start_from + b.length - 2]
        if (b.length + b.exp_idx) % 2 == 0:
            init_expidx = -(b.length + b.exp_idx - 2) / 2
            val = first_num * self.radix + second_num
        else:
            init_expidx = -(b.length + b.exp_idx - 1) / 2
            val = first_num
        if val != 1:
            val = int(math.floor(self.radix * 1/math.sqrt(val)))
            init_expidx -= 1
        val = [val]
        # 上面是优选牛顿迭代的初始值。对待求解数，取其最高位作开方估计，作为初始值。

        #----
        def f(x):
            '''
            欲求sqrt(y), 令h(x) =1/x/x-y. 用牛顿法求出x≈1/sqrt(y).
            不用h(x) = x*x - y 是为了规避迭代公式中的除法
            迭代公式，以计算1/sqrt(.), 直接用牛顿法求sqrt(.)会导致迭代公式里有大数除法。如此则规避了除法。
            详细参https://www.guokr.com/blog/444081/
            '''
            # BigNum([self.radix/2], 1, -1) == 1/2 = 0.5, 为了把除法转为乘法
            b1 = b
            to_be = x.length * 2 + self.extra_digit
            if to_be < b.length:
                exp_idx = b.exp_idx + (b.length - to_be)
                #b1 = BigNum(b.val[b.length-to_be: b.length], to_be, exp_idx)
                b1 = BigNum(b.val, start_from=b.start_from+b.length-to_be, length=to_be, exp_idx=exp_idx)
            ret = x + x * (BigNum(1) - b1*x.pow2()) * BigNum([self.radix/2], 1, -1)
            return ret
        #----
        aa = BigNum(val, 1, init_expidx)
        # print aa, val, init_expidx
        #print aa, aa.length, aa.exp_idx, aa.val
        
        x = self.Newton_method(f, aa, name='sqrt')

        # 先求 1/sqrt(.)
        ret = self * x
        ret.cut_len(self.max_digit + self.extra_digit)

        return ret

    def Newton_method(self, func, init_val, max_loop=100, name=''):
        '''
        init_val似乎初始值。需要先预估一个比较好的初始值。否则有可能导致不收敛。
        '''
        tm = time.time()
        f = func
        x = init_val
        last_2 = [-1, -1]
        has_first_2 = False
        digit_num = 0
        # print "init=", init_val.length, init_val.val, init_val.exp_idx

        last_precise = -1
        last_x = None
        last_realprecise = -1
        tm_start = time.time()
        for i in xrange(max_loop):
            last_x = x
            tm1 = time.time()
            x = f(x)
            x_len = x.length
            if not has_first_2:
                # 再有两个正确数字之前，按小精度计算
                if x_len >= 2:
                    x_len_1 = x.length + x.start_from
                    if last_2[0] == x.val[x_len_1-1] and last_2[1] == x.val[x_len_1-2]:
                        # 前两个数字和上一轮一样的时候，认为这两个数字是对的了
                        has_first_2 = True
                        x.cut_len(2)
                        digit_num = 2
                        last_precise = 2
                        last_realprecise = 2
                    else:
                        last_2[0], last_2[1] = x.val[x_len_1-1], x.val[x_len_1-2]
                        if x.length > 4: 
                            # 再有两个正确数字之前，按小精度计算，所以强制作精度截取。
                            # 这样保证可以很短的时间找到两个有效的数字
                            x.cut_len(4)
            else:
                digit_num *= 2
                #牛顿法是倍增收敛，并不是严格保证正确有效数字倍增，而是基本如此，
                #可能会差或多那么一两位。所以要跟踪精度，需要判断到底倍增后是差一两位还是多一两位
                found_cnt = -1
                max_val = last_x.length - last_precise
                max_val = max_val if max_val >= 3 else 3
                for j in xrange(max_val):
                    if last_x.length - (j + last_precise) < 0 or x.length - (j + last_precise) < 0:
                        continue
                    if last_x.val[last_x.start_from + last_x.length - (j + last_precise)] == \
                                     x.val[x.start_from + x.length - (j + last_precise)]:
                        found_cnt += 1
                    else:
                        break
                assert found_cnt != -1
                last_precise += found_cnt
                last_realprecise = last_precise
                last_precise = last_precise * 2 - 2
                x.cut_len(digit_num)
                x.precision = last_precise
            #print "after, %d->%.4f %.4f" % (i, time.time() - tm1, time.time() - tm_start), name, "len =", x.length
                
            #print "round=%d precise=%d exp_idx=%d length=%d num=%s" % (i, digit_num, x.exp_idx, x.length, str(x)[:200])
            #print "round=%d precise=%d exp_idx=%d length=%d num=%s" % (i, digit_num, x.exp_idx, x.length, x.get_string_val(True, False)), "|", last_realprecise, last_x.get_string_val(True, True)[:last_realprecise+2], last_x.get_string_val(True, False)[2+last_realprecise:]
            #if digit_num >= self.max_digit * 2:
            #print x.val[x.length-20:x.length]
            if x.precision - self.extra_digit >= self.max_digit:
                # print "vvvvvv", digit_num, self.max_digit, last_2
                x.cut_len(self.max_digit + self.extra_digit)
                break
        x.precision = self.max_digit
        #print "nd tm=%.4f name=%s" % (time.time() - tm, name)
        #print "ret", x
        return x

def test_mul():
    BigNum.init(20000, 10000000, use_fft=False)
    #BigNum.init(100000, 10, use_fft=True)
    aa = BigNum(1)
    for i in xrange(1000):
        aa = aa * BigNum(i+1)
    print aa.get_string_val(True)

def test_div_by_one():
    BigNum.init(200, 100, use_fft=False)
    aa = BigNum(17)
    print BigNum(1000) / aa

def test_add():
    from bignum_py import BigNumber
    BigNum.init(100000, 10, use_fft=True)
    #aa = BigNum(0)
    #for i in xrange(1000):
    #    aa = aa + BigNum(i+1, is_neg=True)
    #print aa.get_string_val(True)

    for i in xrange(100):
        aa = random.randint(0, 10000000000000)
        bb = random.randint(0, 10000000000000)
        aa_exp_idx = random.randint(-20, 20)
        bb_exp_idx = random.randint(-20, 20)
        aa_sign = True if random.randint(0, 1) == 1 else False
        bb_sign = True if random.randint(0, 1) == 1 else False

        aaa = BigNum(aa, exp_idx=aa_exp_idx, is_neg=aa_sign)
        bbb = BigNum(bb, exp_idx=bb_exp_idx, is_neg=bb_sign)

        aa1 = BigNumber(aa, exp_idx=aa_exp_idx, is_neg=aa_sign)
        bb1 = BigNumber(bb, exp_idx=bb_exp_idx, is_neg=bb_sign)

        if str(aaa+bbb).split("(")[0] != str(aa1+bb1).split("(")[0]:
            print "~:", str(aaa-bbb).split("(")[0], str(aa1-bb1).split("(")[0]
            print "bignum:", aaa, bbb
            print "bignumber", aa1, bb1
            1/0
            break

if __name__ == "__main__":
    BigNum.init(max_digit=120000, radix=1000000000, use_fft=False)
    tm = time.time()
    aa = BigNum(2)
    print aa.sqrt()
    #print BigNum(2).pow(1230000)
    #print BigNum(33).reciprocal()

    #test_div_by_one()
    #test_add()
    print "tm", time.time() - tm
    sys.exit(0)

    tm = time.time()
    aa = BigNum(2).sqrt()
    print aa
    print time.time() - tm
    xx =  aa.get_string_val(True)
    yy = open("2").read()
    c = 0
    for i in xrange(min(len(xx), len(yy))):
        if xx[i] == yy[i]:
            c += 1
        else:
            break
    print "cnt", c

    tm = time.time()
    aa = aa.sqrt()
    print time.time() - tm
    bb = aa * aa
    xx = bb.get_string_val(True)
    yy = open("2").read()
    c = 0
    for i in xrange(min(len(xx), len(yy))):
        if xx[i] == yy[i]:
            c += 1
        else:
            break
    print "cnt", c
    sys.exit(0)
    a = 981345343
    b = 314159265

    a = random.randint(100000000000000000000000000000, 1000000000000000000000000000000-1)
    b = random.randint(100000000000000000000000000000, 1000000000000000000000000000000-1)
    aa = BigNum(a)
    bb = BigNum(b)
    #a = [350720055,551038198,282600109,233668843,58361539]
    #b = [402348077,174833859,717914086,34944397,546510800]
    #a.reverse()
    #b.reverse()
    #aa = BigNum(a, 5, 0)
    #bb = BigNum(b, 5, 0)
    rr = aa * bb
    print a*b, "should be me."
    print rr, rr.val
    print aa, aa.val
    print bb, bb.val
    assert (aa.get_sum_mod9() * bb.get_sum_mod9()) % 9 == rr.get_sum_mod9()
