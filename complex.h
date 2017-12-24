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

    static Complex get_w_pow(int n, int i)
    {
        double real = cos(2.L* i * 3.141592653589793238462643383279502884197L / n);
        double image = sin(2.L* i * 3.141592653589793238462643383279502884197L / n);
        return Complex(real, image);
    }
    static void print_complex(Complex aa)
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

    static void print_complex_list(Complex * aa, int aa_s, bool show_all=false)
    {
        for (int i = aa_s - 1; i>=0; --i) {
            if (show_all) {
                print_complex(aa[i]);
            } else {
                if (i > aa_s - 100 && i > 100) print_complex(aa[i]);
                if (i == 100) printf("...");
                if (i < 100) print_complex(aa[i]);
            }
        }
        printf("(%d)\n", aa_s);
    }

};

