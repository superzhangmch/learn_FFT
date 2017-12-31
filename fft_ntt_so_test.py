from ctypes import cdll, c_int, byref, c_float
import sys
import time
import math
import random

libc = cdll.LoadLibrary("./fft_ntt.so")

def gen_rand(sz, radix=10):
    max_val = radix-1
    aa = [random.randint(0, max_val) for _ in xrange(sz)]
    if aa[-1] == 0:
        aa[-1] = random.randint(1, max_val)
    return aa
def a2rs(arr, radix=10):
    num = int(round(math.log(radix)/math.log(10)))
    return "".join([("0"*num+str(a))[-num:] for a in arr][::-1])
    
def init_fft_fnt(is_fft, max_digit_num):
    if is_fft:
        libc.fft_init(max_digit_num)
    else:
        libc.fnt_init(max_digit_num)
    return libc

def do_fft_fnt(is_fft, aa, aa_s, bb, bb_s, radix):
    k =  int(math.ceil(math.log(max(aa_s, bb_s)) / math.log(2)))
    cc_s = 2 **(k+1)
    AA = (c_int * aa_s)(*aa[:aa_s])
    BB = (c_int * bb_s)(*bb[:bb_s])
    CC = (c_int * cc_s)()
    if is_fft:
        cc_s = libc.do_fft(byref(AA), c_int(aa_s), byref(BB), c_int(bb_s), radix, byref(CC))
    else:
        cc_s = libc.do_fnt(byref(AA), c_int(aa_s), byref(BB), c_int(bb_s), radix, byref(CC))
    return list(CC), cc_s

# void big_div_by_1digit(uint32_t * divisor, int divisor_size, 
#                uint32_t div_by, 
#                uint32_t * quot, int quot_size, int * out_size, 
#                uint32_t radix)
def div_by_1digit(aa, aa_s, div_by, out_size, radix):
    AA = (c_int*aa_s)(*aa[:aa_s])
    BB = (c_int*out_size)()
    real_out_size = c_int(0)
    libc.big_div_by_1digit(byref(AA), aa_s, div_by, byref(BB), out_size, byref(real_out_size), radix)
    return list(BB), real_out_size.value
#aa = [1,2,3,4,1]
#bb = 7
#cc,cc_s = div_by_1digit(aa, len(aa), bb, 100, 10)
#print type(cc[0])
#print "".join([str(c) for c in cc][::-1]), cc_s

# void big_add(uint32_t * bigger, int bigger_size,
#      uint32_t * smaller, int smaller_size,
#      uint32_t *out, int*out_size,
#      int big_zero_tail_cnt, int radix)
def add(bigger, big_size, smaller, small_size, big_zero_tail_cnt, radix):
    if big_zero_tail_cnt > 0:
        if big_size + big_zero_tail_cnt < small_size:
            bigger, big_size, smaller, small_size = smaller, small_size, bigger, big_size
            big_zero_tail_cnt *= -1
    elif big_zero_tail_cnt <= 0:
        if small_size + abs(big_zero_tail_cnt) > big_size:
            bigger, big_size, smaller, small_size = smaller, small_size, bigger, big_size
            big_zero_tail_cnt *= -1

    AA = (c_int*big_size)(*bigger[:big_size])
    BB = (c_int*small_size)(*smaller[:small_size])
    CC = (c_int* (1+max(big_size, small_size) + abs(big_zero_tail_cnt)))()
    out_size = c_int(0)
    libc.big_add(byref(AA), big_size, byref(BB), small_size, byref(CC), byref(out_size), big_zero_tail_cnt, radix)
    return list(CC), out_size.value

# void big_sub(uint32_t * bigger, int bigger_size,
#      uint32_t * smaller, int smaller_size,
#      uint32_t * out, int * out_size, int * first_is_big,
#      int big_zero_tail_cnt, int radix)
def sub(bigger, big_size, smaller, small_size, big_zero_tail_cnt, radix):
    if big_zero_tail_cnt > 0:
        if big_size + big_zero_tail_cnt < small_size:
            bigger, big_size, smaller, small_size = smaller, small_size, bigger, big_size
            big_zero_tail_cnt *= -1
    elif big_zero_tail_cnt <= 0:
        if small_size + abs(big_zero_tail_cnt) > big_size:
            bigger, big_size, smaller, small_size = smaller, small_size, bigger, big_size
            big_zero_tail_cnt *= -1

    AA = (c_int*big_size)(*bigger[:big_size])
    BB = (c_int*small_size)(*smaller[:small_size])
    CC = (c_int* (max(big_size, small_size) + abs(big_zero_tail_cnt)))()
    out_size = c_int(0)
    first_is_big = c_int(0)
    libc.big_sub(byref(AA), big_size, byref(BB), small_size, byref(CC), byref(out_size), byref(first_is_big), big_zero_tail_cnt, radix)
    return list(CC), out_size.value, first_is_big.value

def test_add_sub():
    for i in xrange(10000):
        aa_s = random.randint(1, 100)
        bb_s = random.randint(1, 100)
        big_zero_tail_cnt = random.randint(-50, 50)
    
        aa = gen_rand(aa_s, 10)
        bb = gen_rand(bb_s, 10)
    
        #aa = [1]
        #bb = [6, 0, 3, 0, 2, 0, 9, 7, 1, 9, 1, 2, 6, 4, 2, 9, 1, 9, 4]
        #aa_s = len(aa)
        #bb_s = len(bb)
        #big_zero_tail_cnt = 30
    
        aa_prex = [] if (big_zero_tail_cnt<=0) else [0] * big_zero_tail_cnt
        bb_prex = [] if (big_zero_tail_cnt>0) else [0] * abs(big_zero_tail_cnt)

        cc,cc_s, is_neg = sub(aa, aa_s, bb, bb_s, big_zero_tail_cnt, 10)

        print aa_s, bb_s, big_zero_tail_cnt
        res = a2rs(cc[:cc_s]) 
    
        if int(res) != abs(int(a2rs(aa_prex+aa))-int(a2rs(bb_prex+bb))):
            print 'xxxx aa:', aa_prex + aa
            print 'xxxx bb:', bb_prex + bb
            print "xxxx n=", big_zero_tail_cnt
            print "xxxx res", res
            print "xxxx aa", a2rs(aa_prex+aa)
            print "xxxx bb", a2rs(bb_prex+bb)
            print "xxxx rr", cc, cc_s
            print "xxxx real", int(a2rs(aa_prex+aa))-int(a2rs(bb_prex+bb))
            break
        #break
if __name__ == "__main__":
    test_add_sub()
    sys.exit(0)

    aa_s = 1000000
    bb_s = 1000000
    radix = 1000000000
    aa = []
    bb = []
    for i in xrange(aa_s):
        aa.append(random.randint(0, radix - 1))
    for i in xrange(bb_s):
        bb.append(random.randint(0, radix - 1))
    cc = [0] * 1000000
    tm = time.time()
    for i in xrange(1000000):
        cc[i] = aa[i] * bb[i]
    #init_fft_fnt(False, 1000000)
    #cc, cc_s = do_fft_fnt(False, aa, len(aa), bb, len(bb), radix)
    #print len(cc)
    #print ((sum(aa) % 9) * (sum(bb) % 9)) % 9
    #print sum(cc[:cc_s]) % 9
    print time.time() - tm
