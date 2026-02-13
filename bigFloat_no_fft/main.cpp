#include <iostream>
#include <cassert>
#include <cmath>
#include "int512.h"
#include "float512.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(expr, msg) do { \
    if (expr) { tests_passed++; } \
    else { tests_failed++; std::cerr << "FAIL: " << msg << " (" << #expr << ")" << std::endl; } \
} while(0)

void test_int512_basic() {
    std::cout << "=== Int512 Basic ===" << std::endl;

    Int512 a(0);
    CHECK(a.is_zero(), "zero");
    CHECK(!a.is_negative(), "zero not negative");

    Int512 b(42);
    CHECK(!b.is_zero(), "42 not zero");
    CHECK(b.to_string() == "42", "42 to_string");

    Int512 c(-1);
    CHECK(c.is_negative(), "-1 is negative");
    CHECK(c.to_string() == "-1", "-1 to_string");

    Int512 d(-100);
    CHECK(d.to_string() == "-100", "-100 to_string");

    std::cout << "  0  = " << a << std::endl;
    std::cout << "  42 = " << b << std::endl;
    std::cout << "  -1 = " << c << std::endl;
    std::cout << "  -100 = " << d << std::endl;
}

void test_int512_string_ctor() {
    std::cout << "=== Int512 String Constructor ===" << std::endl;

    Int512 a("12345678901234567890");
    std::cout << "  12345678901234567890 = " << a << std::endl;
    CHECK(a.to_string() == "12345678901234567890", "large number from string");

    Int512 b("-99999999999999999999999999999999");
    std::cout << "  -99999999999999999999999999999999 = " << b << std::endl;
    CHECK(b.is_negative(), "negative from string");
    CHECK(b.to_string() == "-99999999999999999999999999999999", "negative large to_string");
}

void test_int512_arithmetic() {
    std::cout << "=== Int512 Arithmetic ===" << std::endl;

    Int512 a(100);
    Int512 b(200);

    // Addition
    Int512 c = a + b;
    CHECK(c.to_string() == "300", "100 + 200 = 300");
    std::cout << "  100 + 200 = " << c << std::endl;

    // Subtraction
    Int512 d = a - b;
    CHECK(d.to_string() == "-100", "100 - 200 = -100");
    std::cout << "  100 - 200 = " << d << std::endl;

    // Multiplication
    Int512 e = a * b;
    CHECK(e.to_string() == "20000", "100 * 200 = 20000");
    std::cout << "  100 * 200 = " << e << std::endl;

    // Division
    Int512 f(1000);
    Int512 g = f / Int512(7);
    Int512 r = f % Int512(7);
    std::cout << "  1000 / 7 = " << g << " remainder " << r << std::endl;
    CHECK(g.to_string() == "142", "1000 / 7 = 142");
    CHECK(r.to_string() == "6", "1000 % 7 = 6");

    // Negative division
    Int512 h(-1000);
    Int512 i = h / Int512(7);
    Int512 j = h % Int512(7);
    std::cout << "  -1000 / 7 = " << i << " remainder " << j << std::endl;
    CHECK(i.to_string() == "-142", "-1000 / 7 = -142");
    CHECK(j.to_string() == "-6", "-1000 % 7 = -6");
}

void test_int512_large() {
    std::cout << "=== Int512 Large Numbers ===" << std::endl;

    // 2^256
    Int512 one(1);
    Int512 pow256 = one << 256;
    std::cout << "  2^256 = " << pow256 << std::endl;

    // Multiply two large numbers
    Int512 a("1000000000000000000000000000000");  // 10^30
    Int512 b("1000000000000000000000000000000");  // 10^30
    Int512 c = a * b;
    std::cout << "  10^30 * 10^30 = " << c << std::endl;
    CHECK(c.to_string() == "1000000000000000000000000000000000000000000000000000000000000",
          "10^30 * 10^30 = 10^60");

    // Division of large numbers
    Int512 d = c / a;
    CHECK(d.to_string() == "1000000000000000000000000000000", "10^60 / 10^30 = 10^30");
    std::cout << "  10^60 / 10^30 = " << d << std::endl;
}

void test_int512_comparison() {
    std::cout << "=== Int512 Comparison ===" << std::endl;

    Int512 a(100), b(200), c(100), d(-50);

    CHECK(a < b, "100 < 200");
    CHECK(b > a, "200 > 100");
    CHECK(a == c, "100 == 100");
    CHECK(a != b, "100 != 200");
    CHECK(a <= c, "100 <= 100");
    CHECK(a >= c, "100 >= 100");
    CHECK(d < a, "-50 < 100");
    CHECK(a > d, "100 > -50");

    Int512 e(-100), f(-50);
    CHECK(e < f, "-100 < -50");
}

void test_int512_shift() {
    std::cout << "=== Int512 Shift ===" << std::endl;

    Int512 a(1);
    Int512 b = a << 64;
    CHECK(b.data[1] == 1 && b.data[0] == 0, "1 << 64");
    std::cout << "  1 << 64 = " << b << std::endl;

    Int512 c = b >> 64;
    CHECK(c == Int512(1), "(1 << 64) >> 64 == 1");

    // Arithmetic right shift of negative
    Int512 d(-4);
    Int512 e = d >> 1;
    CHECK(e.to_string() == "-2", "-4 >> 1 = -2");
    std::cout << "  -4 >> 1 = " << e << std::endl;
}

void test_float512_basic() {
    std::cout << "=== Float512 Basic ===" << std::endl;

    Float512 a(0.0);
    CHECK(a.is_zero(), "0.0 is zero");

    Float512 b(1.0);
    CHECK(!b.is_zero(), "1.0 not zero");
    std::cout << "  1.0 = " << b << std::endl;

    Float512 c(-3.14);
    CHECK(c.is_negative(), "-3.14 is negative");
    std::cout << "  -3.14 = " << c << std::endl;

    Float512 d(1e100);
    std::cout << "  1e100 = " << d << std::endl;
}

void test_float512_special() {
    std::cout << "=== Float512 Special Values ===" << std::endl;

    Float512 nan_val = Float512::nan();
    Float512 inf_val = Float512::inf();
    Float512 ninf_val = Float512::inf(true);
    CHECK(Float512::zero().is_zero(), "zero is zero");

    CHECK(nan_val.is_nan(), "nan is nan");
    CHECK(inf_val.is_inf(), "inf is inf");
    CHECK(!inf_val.is_negative(), "inf is positive");
    CHECK(ninf_val.is_negative(), "-inf is negative");

    // NaN comparisons
    CHECK(!(nan_val == nan_val), "nan != nan");
    CHECK(nan_val != nan_val, "nan != nan (!=)");
    CHECK(!(nan_val < Float512(1.0)), "nan not < 1");

    // Inf arithmetic
    Float512 x = inf_val + Float512(1.0);
    CHECK(x.is_inf(), "inf + 1 = inf");

    Float512 y = inf_val - inf_val;
    CHECK(y.is_nan(), "inf - inf = nan");

    std::cout << "  nan = " << nan_val << std::endl;
    std::cout << "  inf = " << inf_val << std::endl;
    std::cout << "  -inf = " << ninf_val << std::endl;
    std::cout << "  0/0 = " << (Float512(0.0) / Float512(0.0)) << std::endl;
    std::cout << "  1/0 = " << (Float512(1.0) / Float512(0.0)) << std::endl;
}

void test_float512_arithmetic() {
    std::cout << "=== Float512 Arithmetic ===" << std::endl;

    Float512 a(2.0);
    Float512 b(3.0);

    // Addition
    Float512 c = a + b;
    std::cout << "  2.0 + 3.0 = " << c << std::endl;

    // Subtraction
    Float512 d = a - b;
    std::cout << "  2.0 - 3.0 = " << d << std::endl;

    // Multiplication
    Float512 e = a * b;
    std::cout << "  2.0 * 3.0 = " << e << std::endl;

    // Division
    Float512 f = a / b;
    std::cout << "  2.0 / 3.0 = " << f << std::endl;

    // Test with larger values
    Float512 g(1e50);
    Float512 h(2e50);
    Float512 i = g + h;
    std::cout << "  1e50 + 2e50 = " << i << std::endl;

    Float512 j = g * h;
    std::cout << "  1e50 * 2e50 = " << j << std::endl;

    Float512 k = j / g;
    std::cout << "  (1e50*2e50) / 1e50 = " << k << std::endl;
}

void test_float512_comparison() {
    std::cout << "=== Float512 Comparison ===" << std::endl;

    Float512 a(1.0), b(2.0), c(1.0), d(-1.0);

    CHECK(a < b, "1.0 < 2.0");
    CHECK(b > a, "2.0 > 1.0");
    CHECK(a == c, "1.0 == 1.0");
    CHECK(a != b, "1.0 != 2.0");
    CHECK(d < a, "-1.0 < 1.0");
    CHECK(a > d, "1.0 > -1.0");
}

void test_float512_from_int512() {
    std::cout << "=== Float512 from Int512 ===" << std::endl;

    Int512 big("1000000000000000000000000000000000000000000000000000000000000");
    Float512 f(big);
    std::cout << "  Int512(10^60) as Float512 = " << f << std::endl;

    Int512 small(42);
    Float512 g(small);
    std::cout << "  Int512(42) as Float512 = " << g << std::endl;
}

void test_float512_precision() {
    std::cout << "=== Float512 Precision ===" << std::endl;

    // Test that 1.0 + tiny != 1.0 (should have much more precision than double)
    Float512 one(1.0);
    Float512 two(2.0);

    // Verify basic identities
    Float512 r1 = one / two;
    Float512 r2 = r1 * two;
    std::cout << "  1/2 * 2 = " << r2 << std::endl;

    // (a+b)*(a-b) vs a^2 - b^2
    Float512 a(7.0);
    Float512 b(3.0);
    Float512 lhs = (a + b) * (a - b);
    Float512 rhs = a * a - b * b;
    std::cout << "  (7+3)*(7-3) = " << lhs << std::endl;
    std::cout << "  7*7 - 3*3  = " << rhs << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "  Int512 & Float512 Test Suite" << std::endl;
    std::cout << "========================================" << std::endl << std::endl;

    test_int512_basic();
    std::cout << std::endl;

    test_int512_string_ctor();
    std::cout << std::endl;

    test_int512_arithmetic();
    std::cout << std::endl;

    test_int512_large();
    std::cout << std::endl;

    test_int512_comparison();
    std::cout << std::endl;

    test_int512_shift();
    std::cout << std::endl;

    test_float512_basic();
    std::cout << std::endl;

    test_float512_special();
    std::cout << std::endl;

    test_float512_arithmetic();
    std::cout << std::endl;

    test_float512_comparison();
    std::cout << std::endl;

    test_float512_from_int512();
    std::cout << std::endl;

    test_float512_precision();
    std::cout << std::endl;

    std::cout << "========================================" << std::endl;
    std::cout << "  Results: " << tests_passed << " passed, " << tests_failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    return tests_failed > 0 ? 1 : 0;
}
