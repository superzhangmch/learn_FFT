#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

class uint128_t {
// partly from https://github.com/calccrypto/uint128_t
public:
    uint64_t UPPER;
    uint64_t LOWER;
    uint128_t(uint64_t up, uint64_t low):UPPER(up), LOWER(low) {}
    uint128_t(uint64_t low): UPPER(0), LOWER(low) {}

private:
    // left -> right, start from 0
    void bit_or(int s)
    {
        if (s < 64) {
            UPPER |= (0x8000000000000000 >> s);
        } else {
            LOWER |= (0x8000000000000000 >> (s - 64));
        }
    }

    int get_len()
    {
        int out = 0;
        if (UPPER){
            out = 64;
            uint64_t up = UPPER;
            while (up){
                up >>= 1;
                out++;
            }
        }
        else{
            uint64_t low = LOWER;
            while (low){
                low >>= 1;
                out++;
            }
        }
        return out;
    }

    void print_byte(int b, char*buf) {
        char a1 = b >> 4;
        char a2 = b & 0x0f;
        char s[16] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
        buf[0] = s[a1];
        buf[1] = s[a2];
    }

public:

    void p(const char * str=0) {
        char buf[50]="";
        for (int i = 0; i < 8; ++i) {
            int j = 7-i;
            print_byte((UPPER & (0xffL<<(8*j)))>>(8*j), buf+i*2);
            print_byte((LOWER & (0xffL<<(8*j)))>>(8*j), buf+18+i*2);
        }

        buf[16]= '\0';
        buf[32+2]= '\0';

        if (str) {
            printf("%s: ", str);
        }
        printf("0x%s,0x%s/0x%s%s\n", buf, buf+18, buf, buf+18);
    }

    uint128_t operator+(uint128_t right) {
        return uint128_t(UPPER + right.UPPER + ((LOWER + right.LOWER) < LOWER), LOWER + right.LOWER);
    }

    uint128_t operator-(uint128_t right) {
        return uint128_t(UPPER - right.UPPER - ((LOWER - right.LOWER) > LOWER), LOWER - right.LOWER);
    }

    bool operator==(uint128_t right) {
        return (LOWER == right.LOWER) && (UPPER == right.UPPER);
    }

    bool operator <(uint128_t right) {
        return (UPPER < right.UPPER) || (UPPER == right.UPPER && LOWER < right.LOWER);
    }
    bool operator >=(uint128_t right) {
        return (UPPER > right.UPPER) || (UPPER == right.UPPER && LOWER >= right.LOWER);
    }

    uint128_t operator>>(int s) {
        uint128_t d = *this;
        if (s == 0) {
            return *this;
        }
        if (s >= 128) {
            return uint128_t(0, 0);
        } else if (s >= 64) {
            return uint128_t(0, d.UPPER >> (s - 64));
        } else {
            uint64_t up = d.UPPER >> s;
            uint64_t low = (d.LOWER >> s) | (d.UPPER << (64 - s));
            return uint128_t(up, low);
        }
    }

    uint128_t operator<<(int s) {
        uint128_t d = *this;
        if (s == 0) {
            return *this;
        }
        if (s >= 128) {
            return uint128_t(0, 0);
        } else if (s >= 64) {
            return uint128_t(d.LOWER << (s - 64), 0);
        } else {
            uint64_t up = (d.UPPER << s) | (d.LOWER >> (64 - s));
            uint64_t low = d.LOWER << s;
            return uint128_t(up, low);
        }
    }

    uint128_t operator %(uint128_t right) {
        if (right == uint128_t(1) || *this == right) {
            return uint128_t(0);
        }
        else if ((*this == uint128_t(0)) || (*this < right)){
            return *this;
        }

        int l_len = this->get_len();
        int r_len = right.get_len();

        // printf("%d %d \n", l_len, r_len);
        uint128_t r = right << (l_len - r_len);
        uint128_t l = *this;
        //printf("%d %d %d\n", l_len, r_len, l_len - r_len);
        //right.p("r0");
        //r.p("r1");
        for (int i = 0; i < l_len - r_len + 1; ++i) {
            if (l >= r) {
                l = l - r;
            }
            r = r >> 1;
        }
        return l;
    }

    uint128_t operator /(uint128_t right) {
        if (right == uint128_t(1)) {
            return *this;
        }
        if (*this == right) {
            return uint128_t(1);
        }
        else if ((*this == uint128_t(0)) || (*this < right)){
            return uint128_t(0);
        }

        int l_len = this->get_len();
        int r_len = right.get_len();

        uint128_t r = right << (l_len - r_len);
        // printf("%d %d, %d\n", l_len, r_len, l_len - r_len);
        uint128_t l = *this;
        uint128_t res(0);
        //l.p("l");
        //r.p("r");
        //right.p("right");
        for (int i = 0; i < l_len - r_len + 1; ++i) {
            if (l >= r) {
                l = l - r;
                res.bit_or(128 - l_len + i + r_len-1);
            }
            r = r >> 1;
        }
        return res;
    }

    uint128_t operator*(uint128_t right) {
        uint128_t a0 = (this->mul(right.LOWER & 0xFFFFFFFF));
        if (right.UPPER == 0 && (right.LOWER >> 32) == 0) {
            return a0;
        }
        uint128_t a1 = (this->mul(right.LOWER >> 32)) << 32;
        if (right.UPPER == 0) {
            return a0 + a1;
        }
        uint128_t a2 = (this->mul(right.UPPER & 0xFFFFFFFF)) << 64;
        if ((right.UPPER >> 32) == 0) {
            return a0 + a1 + a2;
        }
        uint128_t a3 = (this->mul(right.UPPER >> 32)) << (32 * 3);
        return a0 + a1 + a2 + a3;
    }

    uint128_t mul(uint32_t multi_by) {

        // high part -> low part
        uint64_t products[4] = {UPPER >> 32, UPPER & 0xffffffff, LOWER >> 32, LOWER & 0xffffffff};

        //printf("xxxx %u\n", multi_by);
        products[0] *= multi_by;
        products[1] *= multi_by;
        products[2] *= multi_by;
        products[3] *= multi_by;

        uint64_t upper = 0, lower = 0;
        // 1. 
        lower = products[3] & 0xffffffff;
        // 2. 
        products[2] += (products[3] >> 32);
        lower |= (products[2] << 32);
        // 3. 
        products[1] += (products[2] >> 32);
        upper = products[1] & 0xffffffff;
        // 
        products[0] += (products[1] >> 32);
        upper |= (products[0] << 32);
        return uint128_t(upper, lower);
    }
};

/*
int main()
{
    //uint128_t a(0x12345678, 0x987654);
    uint128_t a(0x987654, 0x12345678);
    uint128_t b(0x34323, 0x12345678);
    printf("\n");
    (a/b).p("/");
    (a%b).p("%");
    //printf("\n");
    //(a*b).p("*");
    //(a+b).p("+");
    //(a-b).p("-");
    a.p("a");
    b.p("b");
    return 0;
}
*/
