#encoding: utf8
import math
import sys
import time
import random
from fft_product import ProductFFT


class BigNumber(object):

    max_digit = 100
    radix = 10
    fft_prod = None

    @staticmethod
    def init(max_digit, radix, use_fft, use_so=True):
        BigNumber.max_digit = max_digit
        BigNumber.radix = radix
        BigNumber.fft_prod = ProductFFT(max_digit * 8, use_fft, use_so)

    def dec2other(self, num, radix):
        out = []
        while num:
            remain = num % radix
            num = num / radix
            out.append(int(remain))
        return out

    def __init__(self, val, length=None, exp_idx=0, is_neg=False, is_rev=True):
        radix = self.radix
        self.debug = 0
        self.precision = 0
        if type(val) == type(1) or type(val) == long:
            if val == 0:
                self.val = []
            else:
                self.val = self.dec2other(abs(val), radix)
            self.length = len(self.val)
            self.exp_idx = exp_idx
            if val == 0:
                self.is_neg = True
            else:
                self.is_neg = True if val / abs(val) == -1 else False
            return

        if type(val) == type(""):
            val = [int(i) for i in val]

        if length is None:
            length = len(val)
        if is_rev:
            self.set_val(val, length, exp_idx, is_neg)
        else:
            self.set_normal_val(val, length, exp_idx, is_neg)

    def set_val(self, val, length, exp_idx, is_neg=False):
        assert len(val) >= length
        self.val = val
        self.length = length
        self.exp_idx = exp_idx
        self.is_neg = is_neg

        max_len = length
        for i in xrange(max_len-1, -1, -1):
            if self.val[i] == 0:
                self.length -= 1
            else:
                break

    def extend_len(self, size):
        if len(self.val) < size:
            self.val += [0] * (size - len(self.val))

    def set_normal_val(self, val, length, exp_idx=0, is_neg=False):
        self.length = length
        aa = val + []
        aa.reverse()
        self.val = aa
        self.exp_idx = exp_idx
        self.is_neg = is_neg

        max_len = length
        for i in xrange(max_len-1, -1, -1):
            if self.val[i] == 0:
                self.length -= 1
            else:
                break

    def get_normal_val(self):
        aa = self.val + []
        aa = aa[:self.length]
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

    def get_sum_mod9(self):
        return sum(self.val[:self.length]) % 9

    def get_string_val(self, show_all=False, show_real_precision=True):
        if self.length == 0:
            return "0"
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
        res = res.lstrip("0") + "*E" + str(c*(self.exp_idx + self.length)) + "("+str(self.length*c)+"|"+str(precision)+")"
        res = ("-" if self.is_neg else "") + res + "sum_mod9=%d" % (self.get_sum_mod9())
        return res

    def __str__(self):
        return self.get_string_val(False)

    def adjust_precision(self, precision):
        if self.length <= precision:
            return
        else:
            self.exp_idx += (self.length - precision)
            self.length = precision
            self.val = self.val[-precision:]

    def adjust_exp_idx(self, new_exp_idx):
        if new_exp_idx == self.exp_idx:
            return
        if new_exp_idx < self.exp_idx:
            adjust_val = self.exp_idx - new_exp_idx
            self.val = ([0] * adjust_val) + self.val
            self.exp_idx = new_exp_idx
            self.length += adjust_val
            return
        else:
            adjust_val = new_exp_idx - self.exp_idx
            assert adjust_val <= self.length
            all_z = True
            for i in xrange(adjust_val):
                if self.val[i] != 0:
                    all_z = False
                    break
            assert all_z == True
            self.length -= adjust_val
            self.val = self.val[adjust_val:]
            self.exp_idx = new_exp_idx
        
    # =====================
    def __mul__(self, num1):
        if self.length:
            assert self.val[self.length - 1] != 0
        if num1.length:
            assert num1.val[num1.length - 1] != 0
        res = BigNumber(0)
        self.fft_prod.fast_prod(self, num1, res, self.radix)
        if res.length:
            assert res.val[res.length - 1] != 0
        return res

    def __add__(self, num1):
        if self.length:
            assert self.val[self.length - 1] != 0
        if num1.length:
            assert num1.val[num1.length - 1] != 0

        if self.is_neg != num1.is_neg:
            # sub: (+) + (-) || (-) + (+)
            res = self.do_sub(num1)
        # add: (+) + (+) || (-) + (-)
        else:
            res = self.do_add(num1)
        if res.length:
            assert res.val[res.length - 1] != 0
        return res

    def __sub__(self, num1):
        if self.length:
            assert self.val[self.length - 1] != 0
        if num1.length:
            assert num1.val[num1.length - 1] != 0

        if self.is_neg != num1.is_neg:
            # add: (+) - (-) => (+) || (-) - (+) => (-)
            res = self.do_add(num1)
        else: # sub: (-) - (-) || (+)-(+)
            res = self.do_sub(num1)
        if res.length:
            assert res.val[res.length - 1] != 0
        return res

    def do_add(self, num1):
        """ same sign, result.sign = self.sign """
        if self.length == 0:
            if self.is_neg == num1.is_neg:
                return num1
            else:
                val = num1.val+[]
                return BigNumber(val, num1.length, num1.exp_idx, is_neg=self.is_neg)
        if num1.length == 0:
            return self

        if self.exp_idx <= num1.exp_idx:
            a, b = self, num1
        else:
            a, b = num1, self
        # now a.exp_idx < b.exp_idx. b should be padded with zero at the begining

        exp_idx_diff = b.exp_idx - a.exp_idx
        if a.length >= b.length + exp_idx_diff:
            a_is_longer = True
        else:
            a_is_longer = False

        use_exp_idx = a.exp_idx
        min_length = min(a.length, b.length + exp_idx_diff)
        res_length = max(a.length, b.length + exp_idx_diff)
        radix = self.radix

        sum_arr = [0] * (1 + res_length)
        length_a = a.length
        a = a.val
        b = b.val
        #print "xx", exp_idx_diff, min_length, res_length
        #print "xx, a=", a
        #print "xx, b=", b
        remain = 0
        for i in xrange(res_length):
            if i < exp_idx_diff:
                if i < length_a:
                    sum_arr[i] = a[i]
            elif i < min_length:
                s = a[i] + b[i - exp_idx_diff] + remain
                if s >= radix:
                    sum_arr[i] = s % radix
                    remain = s / radix
                else:
                    sum_arr[i] = s
                    remain = 0
                #print "s=%d, remain=%d, radix=%d, res=%d" %( s, remain, radix, sum_arr[i])
            else:
                if a_is_longer:
                    s = a[i]
                else:
                    s = b[i - exp_idx_diff]
                if remain:
                    s += remain
                    sum_arr[i] = s % radix
                    remain = s / radix
                else:
                    sum_arr[i] = s
        if remain:
            sum_arr[res_length] = remain
            res_length += 1
        return BigNumber(sum_arr, res_length, use_exp_idx, is_neg=self.is_neg)

    def do_sub(self, num1):
        """ self.sign != num1.sign, result.sign = self.sign if (|self| > |num1| else (!self.sign)) """
        if self.length == 0:
            res_neg = not self.is_neg
            if res_neg == num1.is_neg:
                return num1
            else:
                val = num1.val+[]
                return BigNumber(val, num1.length, num1.exp_idx, is_neg=not self.is_neg)
        if num1.length == 0:
            return self

        is_a_self = False
        if self.exp_idx <= num1.exp_idx:
            a, b = self, num1
            is_a_self = True
        else:
            a, b = num1, self
        # now a.exp_idx < b.exp_idx. b should be padded with zero at the begining

        exp_idx_diff = b.exp_idx - a.exp_idx
        is_a_sub_b = False
        if a.length > b.length + exp_idx_diff:
            is_a_sub_b = True
        elif a.length < b.length + exp_idx_diff:
            is_a_sub_b = False
        else:
            all_equal = True
            for i in xrange(a.length - 1, exp_idx_diff - 1, -1):
                if a.val[i] > b.val[i - exp_idx_diff]:
                    is_a_sub_b = True
                    all_equal = False
                    break
                elif a.val[i] < b.val[i - exp_idx_diff]:
                    all_equal = False
                    break
            if all_equal:
                is_a_sub_b = True

        use_exp_idx = a.exp_idx
        min_length = min(a.length, b.length + exp_idx_diff)
        res_length = max(a.length, b.length + exp_idx_diff)
        radix = self.radix

        sub_arr = [0] * res_length
        length_a = a.length
        a = a.val
        b = b.val
        #print "xx", exp_idx_diff, min_length, res_length
        #print "xx, a=", a
        #print "xx, b=", b
        #print "is_a_sub_b", is_a_sub_b
        borrow = 0
        for i in xrange(res_length):
            if i < exp_idx_diff:
                if is_a_sub_b:
                    sub_arr[i] = a[i]
                else:
                    if i < length_a:
                        s = - a[i] - borrow
                        if s < 0:
                            s += radix
                            borrow = 1
                        else:
                            borrow = 0
                        sub_arr[i] = s
                    else:
                        s = - 0 - borrow
                        if s < 0:
                            s += radix
                            borrow = 1
                        else:
                            borrow = 0
                        sub_arr[i] = s
            elif i < min_length:
                if is_a_sub_b:
                    s = a[i] - b[i - exp_idx_diff] - borrow
                else:
                    s = b[i - exp_idx_diff] - a[i] - borrow
                if s < 0:
                    sub_arr[i] = s + radix
                    borrow = 1
                else:
                    sub_arr[i] = s
                    borrow = 0
                #print "s=%d, remain=%d, radix=%d, res=%d" %( s, remain, radix, sum_arr[i])
            else:
                if is_a_sub_b:
                    s = a[i] - borrow
                else:
                    s = b[i - exp_idx_diff] - borrow
                if s < 0:
                    sub_arr[i] = s + radix
                    borrow = 1
                else:
                    sub_arr[i] = s
                    borrow = 0

        for i in xrange(res_length - 1, -1, -1):
            if sub_arr[i] == 0:
                res_length -= 1
            else:
                break

        is_self_bigger = None
        if is_a_self:
            is_self_bigger = True if is_a_sub_b else False
        else:
            is_self_bigger = False if is_a_sub_b else True
        if is_self_bigger:
            res_neg = self.is_neg
        else:
            res_neg = not self.is_neg
        return BigNumber(sub_arr, res_length, use_exp_idx, is_neg=res_neg)

    def __div__(self, num):
        if self.length:
            assert self.val[self.length - 1] != 0
        if num.length:
            assert num.val[num.length - 1] != 0
        assert num.length > 0

        if self.length == 0:
            return BigNumber([], 0, 0)

        if num.length == 1:
            div_by = num.val[0]
            if self.length > self.max_digit:
                arr = self.val
                range_max = self.length
                res_exp_idx = self.exp_idx - num.exp_idx
            else:
                arr = ([0] * (self.max_digit + 1 - self.length)) + self.val
                res_exp_idx = self.exp_idx - num.exp_idx - (self.max_digit + 1 - self.length)
                range_max = self.max_digit + 1
            res_len = range_max if arr[range_max-1] >= div_by else range_max - 1
            out = [0] * (res_len + 1)
            remain = 0
            for i in xrange(range_max - 1, -1, -1):
                cur_val = remain * self.radix + arr[i]
                out[i] = cur_val / div_by
                remain = cur_val % div_by
            is_neg = True if self.is_neg != num.is_neg else False
            res = BigNumber(out, res_len, res_exp_idx, is_neg=is_neg)
            return res
        # 除法转化为乘法。先求倒数。
        res = self * num.reciprocal('div_inv')
        if res.length:
            assert res.val[res.length - 1] != 0
        return res

    def reciprocal(self, name=""):
        b = self
        assert b.length != 0
        if b.length == 1:
            b.adjust_exp_idx(b.exp_idx - 2)

        # print b.val, b.length, 11111111111111
        if b.val[b.length - 1] == 1:
            val = [1]
        else:
            val = [b.val[b.length - 1]]
        init_expidx = -b.length - b.exp_idx - len(val)
        # 上面是优选牛顿迭代的初始值。对待求解数，取其最高位作估计，作为初始值。注意初始值需要是一个分母是radix的整幂次的数，这样才能规避复杂除法。
        
        def f(x):
            '''
            欲求1/y, 令h(x) =1/x-y. 用牛顿法求出x≈1/y.
            '''
            return x*(BigNumber(2) - b*x)
        aa = BigNumber(val, 1, init_expidx)
        return self.Newton_method(f, aa, name='inv' if not name else name)

    def sqrt(self):
        b = self
        if b.length == 0:
            return BigNumber([], 0, 0)
        assert b.is_neg == False

        # print b.val, b.length, 11111111111111222222222222222
        #----
        first_num = b.val[b.length - 1]
        second_num = 0 if b.length == 1 else b.val[b.length - 2]
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
            # BigNumber([self.radix/2], 1, -1) == 1/2 = 0.5, 为了把除法转为乘法
            ret = x + x * (BigNumber(1) - b*x*x) * BigNumber([self.radix/2], 1, -1)
            return ret
        #----
        aa = BigNumber(val, 1, init_expidx)
        # print aa, val, init_expidx
        #print aa, aa.length, aa.exp_idx, aa.val
        
        x = self.Newton_method(f, aa, name='sqrt')
        
        # 先求 1/sqrt(.), 然后取倒数得解
        return x.reciprocal('sqrt_inv')

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
        for i in xrange(max_loop):
            last_x = x
            tm1 = time.time()
            x = f(x)
            x_len = x.length
            if not has_first_2:
                # 再有两个正确数字之前，按小精度计算
                if x_len >= 2:
                    if last_2[0] == x.val[x_len-1] and last_2[1] == x.val[x_len-2]:
                        # 前两个数字和上一轮一样的时候，认为这两个数字是对的了
                        has_first_2 = True
                        x.val = x.val[x_len-2:x_len]
                        x.exp_idx += (x_len - 2)
                        x.length = 2
                        digit_num = 2
                        last_precise = 2
                        last_realprecise = 2
                    else:
                        last_2[0], last_2[1] = x.val[x_len-1], x.val[x_len-2]
                        if x.length > 4: 
                            # 再有两个正确数字之前，按小精度计算，所以强制作精度截取。这样保证可以很短的时间找到两个有效的数字
                            x.val = x.val[x_len-4:x_len]
                            x.exp_idx += (x_len - 4)
                            x.length = 4
            else:
                digit_num *= 2
                #牛顿法是倍增收敛，并不是严格保证正确有效数字倍增，而是基本如此，可能会差或多那么一两位。所以要跟踪精度，需要判断到底倍增后是差一两位还是多一两位
                found_cnt = -1
                max_val = last_x.length - last_precise
                max_val = max_val if max_val >= 3 else 3
                for j in xrange(max_val):
                    if last_x.length - (j + last_precise) < 0 or x.length - (j + last_precise) < 0:
                        continue
                    if last_x.val[last_x.length - (j + last_precise)] == x.val[x.length - (j + last_precise)]:
                        found_cnt += 1
                    else:
                        break
                assert found_cnt != -1
                last_precise += found_cnt
                last_realprecise = last_precise
                last_precise = last_precise * 2 - 2
                x.precision = last_precise
                x.val = x.val[x_len-digit_num:x_len]
                x.exp_idx += (x_len - digit_num)
                x.length = digit_num
                
            #print "round=%d precise=%d exp_idx=%d length=%d num=%s" % (i, digit_num, x.exp_idx, x.length, str(x)[:200])
            #print "round=%d precise=%d exp_idx=%d length=%d num=%s" % (i, digit_num, x.exp_idx, x.length, x.get_string_val(True, False)), "|", last_realprecise, last_x.get_string_val(True, True)[:last_realprecise+2], last_x.get_string_val(True, False)[2+last_realprecise:]
            #if digit_num >= self.max_digit * 2:
            if x.precision >= self.max_digit:
                # print "vvvvvv", digit_num, self.max_digit, last_2
                break
        x.precision = self.max_digit
        #print "nd tm=%.4f name=%s" % (time.time() - tm, name)
        return x

if __name__ == "__main__":
    #BigNumber.init(100 * 4, 1000000000, False)
    BigNumber.init(100 * 4, 10, use_fft=True, use_so=True)

    print BigNumber(2).sqrt()
    sys.exit(0)
    a = 981345343
    b = 314159265

    a = random.randint(100000000000000000000000000000, 1000000000000000000000000000000-1)
    b = random.randint(100000000000000000000000000000, 1000000000000000000000000000000-1)
    aa = BigNumber(a)
    bb = BigNumber(b)
    #a = [350720055,551038198,282600109,233668843,58361539]
    #b = [402348077,174833859,717914086,34944397,546510800]
    #a.reverse()
    #b.reverse()
    #aa = BigNumber(a, 5, 0)
    #bb = BigNumber(b, 5, 0)
    rr = aa * bb
    print a*b, "should be me."
    print rr, rr.val
    print aa, aa.val
    print bb, bb.val
    assert (aa.get_sum_mod9() * bb.get_sum_mod9()) % 9 == rr.get_sum_mod9()
