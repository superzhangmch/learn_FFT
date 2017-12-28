from ctypes import cdll, c_int, byref, c_float
import time
import math
import random

libc = cdll.LoadLibrary("./fft_ntt.so")

def init_fft_fnt(is_fft, max_digit_num):
    if is_fft:
        libc.fft_init(max_digit_num)
    else:
        libc.fnt_init(max_digit_num)

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

if __name__ == "__main__":
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
