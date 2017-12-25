#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// calc: (a ** idx) % p
long fast_m(long a, long idx, long p)
{
    if (idx >= 4) {
        long b = fast_m(a, idx/2, p);
        b = (b * b) % p;
        if (idx % 2 == 0) {
            return b;
        } else {
            return (b * a) % p;
        }
    }

    if (idx < 0) {
        return 0;
    } else if (idx == 1) {
        return a;
    } else if (idx == 2) {
        return (a * a) % p;
    } else if (idx == 3) {
        long b = (a * a) % p;
        return (b * a) % p;
    }
}

class Zp {
public:
    static int P;
    static int G;

    long n;

    Zp(long num) {n = num % P;}
    Zp() {}

    // trans to int
    int to_int() {
        return (int)n;
    }
    // + 
    Zp operator+=(Zp a) {
        n += a.n;
        n = n % P;
        return *this;
    }
    Zp operator+(Zp a) {
        return Zp(n + a.n);
    }

    // -
    Zp operator-=(Zp a) {
        n -= a.n;
        n = n % P;
        return *this;
    }
    Zp operator-(Zp a) {
        return Zp(n - a.n);
    }

    // *
    Zp operator*=(Zp a) {
        n *= a.n;
        n = n % P;
        return *this;
    }
    Zp operator*(Zp a) {
        return Zp(n*a.n);
    }

    // = 
    Zp operator=(Zp a) {
        n = a.n;
        return *this;
    }
    Zp operator=(long a) {
        n = a % P;
        return *this;
    }

    // div
    Zp operator/=(Zp a) {
        n *= a.reciprocal().n;
        n = n % P;
        return *this;
    }
    Zp operator/(Zp a) {
        return (*this) * a.reciprocal();
    }

    // inverse or reciprocal
    Zp reciprocal() {
        return Zp(fast_m(n, P-2, P));
    }

    // pow
    Zp pow(int a) {
        a = a % (P - 1);
        return Zp(fast_m(n, a, P));
    }

    static Zp get_w_pow(long skip, long i)
    {
        return Zp(fast_m(G, (i*skip) % (P-1), P));
    }

    static void print(Zp a)
    {
        printf("%lu", a.n);
    }

    static void print(Zp * aa, int aa_s, bool show_all=false)
    {
        for (int i = aa_s - 1; i>=0; --i) {
            if (show_all) {
                printf("%lu ", aa[i].n);
            } else {
                if (i > aa_s - 100 && i > 100) printf("%lu ", aa[i].n);
                if (i == 100) printf("...");
                if (i < 100) printf("%lu ", aa[i].n);
            }
        }
        printf("(%d)\n", aa_s);
    }

};

