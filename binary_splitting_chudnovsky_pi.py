#encoding: utf8
import math
import time

ii = 0
tm = time.time()
# copy from https://www.craig-wood.com/nick/articles/pi-chudnovsky/
# 用gmpy2.mpz算100万位，只需要不到1秒. gmpy2.mpz 用到了GMP库

def pi_chudnovsky_bs(digits, use_mpz):
    """
    Compute int(pi * 10**digits)

    This is done using Chudnovsky's series with binary splitting
    """
    C = 640320
    C3_OVER_24 = C**3 // 24
    def bs(a, b):
        global ii, tm
        """
        Computes the terms for binary splitting the Chudnovsky infinite series

        a(a) = +/- (13591409 + 545140134*a)
        p(a) = (6*a-5)*(2*a-1)*(6*a-1)
        b(a) = 1
        q(a) = a*a*a*C3_OVER_24

        returns P(a,b), Q(a,b) and T(a,b)
        """
        if b - a == 1:
            ii += 1
            if ii % 1000 == 0:
                print ii, time.time() - tm
                #tm = time.time()

            # Directly compute P(a,a+1), Q(a,a+1) and T(a,a+1)
            if a == 0:
                Pab = Qab = mpz(1)
            else:
                Pab = mpz((6*a-5)*(2*a-1)*(6*a-1))
                Qab = mpz(a*a*a*C3_OVER_24)
            Tab = Pab * (13591409 + 545140134*a) # a(a) * p(a)
            if a & 1:
                Tab = -Tab
        else:
            #ii += 1
            #if ii % 100 == 0:
            #    print ii

            # Recursively compute P(a,b), Q(a,b) and T(a,b)
            # m is the midpoint of a and b
            m = (a + b) // 2
            # Recursively calculate P(a,m), Q(a,m) and T(a,m)
            Pam, Qam, Tam = bs(a, m)
            # Recursively calculate P(m,b), Q(m,b) and T(m,b)
            Pmb, Qmb, Tmb = bs(m, b)
            # Now combine
            Pab = Pam * Pmb
            Qab = Qam * Qmb
            Tab = Qmb * Tam + Pam * Tmb
        return Pab, Qab, Tab
    # how many terms to compute
    DIGITS_PER_TERM = math.log10(C3_OVER_24/6/2/6)
    N = int(digits/DIGITS_PER_TERM + 1)
    # Calclate P(0,N) and Q(0,N)
    P, Q, T = bs(0, N)
    print "tm=1", time.time() - tm
    #one_squared = mpz(1) #0)**(2*digits)
    one_squared = mpz(10)**(2*digits)
    print "tm=2", time.time() - tm
    if use_mpz:
        sqrtC = gmpy2.sqrt(10005*one_squared)
    else:
        sqrtC = (10005*one_squared).sqrt()
    print "tm=3", time.time() - tm
    ret = (Q*426880*sqrtC) 
    print "tm=4", time.time() - tm
    ret = ret / T
    print "tm=5", time.time() - tm
    return ret

if __name__ == "__main__":

    calc_digit = 1000000
    use_mpz = True
    use_mpz = False
    if use_mpz:
        # 百万位不到一秒
        from gmpy2 import mpz
        import gmpy2
        PI = pi_chudnovsky_bs(digits=calc_digit, use_mpz=use_mpz)
        print PI
    else:
        # 百万位45秒，不过已经比AGM快了
        from bignum_c import BigNum as mpz 
        from pi import pi, match_num
    
        max_digit = calc_digit / 9 + 10
        mpz.init(max_digit=max_digit, radix=1000000000, use_fft=False)
        print "begin to calc"
    
        PI = pi_chudnovsky_bs(digits=calc_digit, use_mpz=use_mpz)
        aa = PI.get_all_numbers()
        match_num(pi, aa)
        print PI
