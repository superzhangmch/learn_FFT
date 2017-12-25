#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class Complex {
public:
    double real;
    double image;

    Complex(double real, double img):real(real), image(img) {}
    Complex(double real):real(real), image(0.) {}
    Complex() {}

    int to_int() {
        return int(round(real));
    }
    Complex operator+=(Complex a) {
        real += a.real;
        image += a.image;
        return *this;
    }
    Complex operator+(Complex a) {
        return Complex(real+a.real, image+a.image);
    }

    Complex operator-=(Complex a) {
        real -= a.real;
        image -= a.image;
        return *this;
    }
    Complex operator-(Complex a) {
        return Complex(real-a.real, image-a.image);
    }

    Complex operator*=(Complex aa) {
        double a = real, b = image, c = aa.real, d = aa.image;
        real = a*c-b*d;
        image = b*c+a*d;
        return *this;
    }
    Complex operator*(Complex aa) {
        double a = real, b = image, c = aa.real, d = aa.image;
        return Complex(a*c-b*d, b*c+a*d);
    }

    Complex operator=(Complex aa) {
        real = aa.real;
        image = aa.image;
        return *this;
    }
    Complex operator=(double a) {
        real = a;
        image = 0.;
        return *this;
    }

    Complex operator/=(Complex aa) {
        printf("not supported\n");
        return *this;
    }
    Complex operator/=(double a) {
        real /= a;
        image /= a;
        return *this;
    }
    Complex operator/(Complex aa) {
        printf("not supported\n");
        return Complex(real, image);
    }
    Complex operator/(double a) {
        return Complex(real/a, image/a);
    }

    static int get_w_pow(int n, Complex * c)
    {
        for (int i = 0; i < n; ++i) {
            c[i].real = cos(2.L* i * 3.141592653589793238462643383279502884197L / n);
            c[i].image = sin(2.L* i * 3.141592653589793238462643383279502884197L / n);
        }
        return 0;
    }
    static void print(Complex aa)
    {
        double real = aa.real;
        double image = aa.image;
        if ((real < 0.000001 && real > -0.000001) && (image < 0.000001 && image > -0.000001)) {
            printf("0 ");
            return;
        }
        if (real < 0.000001 && real > -0.000001)
        {
            printf("%.3fj ", image);
            return;
        }
        if (image < 0.000001 && image > -0.000001)
        {
            printf("%.3f ", real);
            return;
        }
        printf("%.3f+%.3fj ", real, image);
    }

    static void print(Complex * aa, int aa_s, bool show_all=false)
    {
        for (int i = aa_s - 1; i>=0; --i) {
            if (show_all) {
                print(aa[i]);
            } else {
                if (i > aa_s - 100 && i > 100) print(aa[i]);
                if (i == 100) printf("...");
                if (i < 100) print(aa[i]);
            }
        }
        printf("(%d)\n", aa_s);
    }

};

