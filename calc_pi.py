#encoding: utf8

import math
import sys
import time
import random
from bignum_c import BigNum
#from bignum_py import BigNumber as BigNum
from pi import pi, match_num


def calc_pi_use_chudnovsky():
    """
    A = 13591409; B = 545140134; C = 640320
    1/pi = 12 * sum((-1)^n*[(6n)!/(n!)^3/(3n)!]*(A+Bn)/(C^(3n+3/2)), n = 0..infinity)
    while (6n)!/(n!)^3/(3n)!=product((12*k+2)*(12*k+6)*(12*k+10)/(k+1)^3, k=0..n-1)
    pi = 426880 * sqrt(10005) * (M*L/X)^(-1)

    https://en.wikipedia.org/wiki/Chudnovsky_algorithm
    http://www.numberworld.org/misc_runs/pi-5t/details.html
    https://bellard.org/pi/pi2700e9/pipcrecord.pdf

    http://blog.jobbole.com/67766/: PI的验证。
        一种方法是两种方法互相验证；
        另一种是按16进制计算，再用BBP作片段验证，最后转化为10进制.
    """
    def test6n_n3_3n():
        def calc_6n_n3_3n(n):
            aa = math.factorial(6*n) / math.factorial(3*n) / math.factorial(n) ** 3
            bb = 1
            for k in xrange(0, n, 1):
                bb *= (12*k+2)*(12*k+6)*(12*k+10)
            for k in xrange(0, n, 1):
                bb /= (k+1)**3
            return aa, bb
        for n in xrange(10):
            aa, bb = calc_6n_n3_3n(n)
            print n, len(str(aa)), aa
            assert aa == bb
        print "-----"
    #test6n_n3_3n()

    #1329sec for 1M digits
    #340sec for 0.5M digits
    calc_base10_digit = 1000000
    # 进制基数，最大支持10^9
    radix = 1000000000
    radix_len = int(math.log(radix) / math.log(10))
    max_digit = calc_base10_digit / radix_len + 20

    BigNum.init(max_digit, radix, False)

    round_cnt = calc_base10_digit / 14 # 每次循环产生14位
    A = BigNum(13591409)
    B = BigNum(545140134)
    C = BigNum(640320)
    L = A
    SUM = A
    CC = BigNum(1)
    M = BigNum(1)
    cur_digit = max_digit
    tm = time.time()
    for k in xrange(round_cnt):
        if k % 100 == 0:
            print "=============== %d/%d tm=%.4f cur_digit=%d" % (k, round_cnt, time.time()-tm, cur_digit)
            #print M.length, M.exp_idx
            tm = time.time()

        k_1 = BigNum(k+1)
        #M = M * BigNum(12*k+2)*BigNum(12*k+6)*BigNum(12*k+10) / k_1 / k_1 / k_1 /C/C/C
        # 全部拆解成一个基数内的乘除法
        M = M.mul(BigNum(12*k+2))* BigNum(12*k+6) * BigNum(12*k+10)/ k_1 / k_1 / k_1 / C / C / C
        M.cut_len(cur_digit)

        #L = L + B
        #if k % 2 == 1:
        #    SUM = SUM + L * M
        #else:
        #    SUM = SUM - L * M

        # 全部拆解成一个基数内的乘除法
        if k % 2 == 1:
            SUM = SUM + M*B*BigNum(k+1) + M*A
        else:
            SUM = SUM - M*B*BigNum(k+1) - M*A

        # 每次循环推进14位，所以M的精度可以少14位
        cur_digit = int((calc_base10_digit - 14 * k) / radix_len) + 10

        #PI = BigNum(426880) * BigNum(10005).sqrt() / (SUM)
        #aa = PI.get_all_numbers().lstrip("0")
        #match_num(pi, aa)

    print cur_digit
    PI = BigNum(426880) * BigNum(10005).sqrt() / (SUM)
    print "--------------"
    print PI
    aa = PI.get_all_numbers().lstrip("0")
    match_num(pi, aa)

def calc_e():
    '''
    自然对数: e= sum(1/n!, n = 1 .. N)
    注意:公式中除以的数字都比较小，应该直接做这种除法，不应该用到FFT的乘除法。
    '''
    max_digit = 10000
    radix = 1000000000
    BigNum.init(max_digit, radix, True)

    sum1 = BigNum(0)
    x = BigNum(1)
    for i in xrange(10000):
        sum1 = sum1 + x
        x = x / BigNum(i + 1)
        if i % 1000 == 0:
            print i
    print sum1

def calc_pi_use_machin():
    '''
    Machin公式计算PI.
    Pi = 16*sum(1/(2*i+1)/5^(2*i+1)*(-1)^i)-4*sum(1/(2*i+1)/239^(2*i+1)*(-1)^i) 
    注意:公式中除以的数字都比较小，应该直接做这种除法，不应该用到FFT的乘除法。
    about 1 million, 2434 seconds
    '''
    calc_digits = 100000
    round_cnt = int(calc_digits/1.4) + 10

    radix = 1000000000

    max_digit = calc_digits / int(math.log(radix) / math.log(10)) + 20
    print "ditig", max_digit
    BigNum.init(max_digit, radix, False)

    x = BigNum(5)
    y = BigNum(239)

    sum1 = BigNum(0)
    sum2 = BigNum(0)

    num5 = BigNum(5)
    num239 = BigNum(239)
    num5pow2 = BigNum(5).pow2()
    num239pow2 = BigNum(239).pow2()

    cur_digit_x = max_digit
    cur_digit_y = max_digit
    print "loop_cnt", round_cnt
    for i in xrange(round_cnt):
        tm = time.time()
        a = 2*i+1
        if radix > 239*239:
            #x = x / num5pow2
            #y = y / num239pow2
            if cur_digit_x > 0:
                x = x.div(num5pow2, cur_digit_x)
            if cur_digit_y > 0:
                y = y.div(num239pow2, cur_digit_y)
        else:
            #x = x / num5 / num5
            #y = y / num239 / num239
            if cur_digit_x > 0:
                x = x.div(num5, cur_digit_x).div(num5, cur_digit_x)
            if cur_digit_y > 0:
                y = y.div(num239, cur_digit_y).div(num239, cur_digit_y)
        xx = BigNum(a, is_neg=True if i % 2 == 1 else False)
        xx = BigNum(a, is_neg=True if i % 2 == 1 else False)
        #sum1 = sum1 + x / xx
        #sum2 = sum2 + y / xx
        if cur_digit_x > 0:
            sum1 = sum1 + x.div(xx, cur_digit_x)
        if cur_digit_y > 0:
            sum2 = sum2 + y.div(xx, cur_digit_y)

        def calc_cur_len(x, cur_len):
            x_sum = x.length + x.exp_idx
            if x_sum >= 0:
                return cur_len
            pre_zero_cnt = abs(x_sum)
            if pre_zero_cnt > max_digit + 10:
                return 0
            ret_val = max_digit - pre_zero_cnt + 10
            if ret_val > x.length + 10:
                return cur_len
            return ret_val
        cur_digit_x = calc_cur_len(x, cur_digit_x)
        cur_digit_y = calc_cur_len(y, cur_digit_y)
                
        if i % 1000 == 0:
            #print sum1.length + sum1.exp_idx, "sum_len=", sum1.length, "sum_exp=", sum1.exp_idx, "x_exp=", x.exp_idx,"x_len=",x.length, "x_exp_len=", (x.length + x.exp_idx) , "y_exp=", y.exp_idx, "y_len=", y.length
            #print "y_exp=", y.exp_idx, "y_len=", y.length
            print "i=%d/r_cnt=%d/digits=%d -> %.4f" % (i, round_cnt, max_digit, time.time() - tm)
            #print cur_digit_x, cur_digit_y
            #print sum2
            #print y
    print "y_exp=", y.exp_idx, "y_len=", y.length
    PI = BigNum(16) * sum1 - BigNum(4) * sum2
    aa = PI.get_all_numbers().strip("0")
    #print aa
    match_num(pi, aa)
    print "pi=", PI

def calc_pi_use_agm():
    '''
    AGM 算术几何平均法. 需要第一轮循环各个数字的精度就达到所需要的精度，整个计算过程中精度不变
    '''
    def radix_len(radix):
        return int(round(math.log(radix)/ math.log(10)))
    radix10_digit = 1000000

    # FFT, 10进制下没法算到100万位，有溢出。更大进制更不行了
    # FNT 则无此限制
    radix = 1000000000
    use_fft = False if radix > 10 else True
    max_digit = radix10_digit / radix_len(radix)
    BigNum.init(max_digit, radix, use_fft)

    a = BigNum(1)
    b = (BigNum(1) / BigNum(2)).sqrt()
    t = BigNum(1) / BigNum(4)
    p = BigNum(1)

    round_cnt = int(math.ceil(math.log(radix10_digit) / math.log(2)))

    print "init ok, now to calc %d rounds:" % (round_cnt)

    for i in xrange(round_cnt):
        print "begin %d" % i
        tm = time.time()
        a1 = (a + b) / BigNum(2)
        b = (a*b).sqrt()
        diff = a - a1
        #t = t - p * diff * diff
        # 如果多线程计算乘法，则自乘也未必节约多少时间
        t = t - p * diff.pow2()
        p = BigNum(2)*p
        a = a1
        print "%d -> %.4f" % (i, time.time() - tm)

    #PI = (a + b) * (a + b) / BigNum(4) / t
    PI = (a + b).pow2() / BigNum(4) / t
    aa = PI.get_all_numbers()
    match_num(pi, aa)
    print PI

if __name__ == "__main__":
    time_start = time.time()
    calc_pi_use_agm()
    #calc_pi_use_chudnovsky()
    #calc_pi_use_machin()
    #calc_e()
    print time.time() - time_start
    sys.exit(0)
