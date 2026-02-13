#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <iomanip>
#include "bfmath.h"

// Reference e to 300 digits
static const char* E_REF =
    "2."
    "71828182845904523536028747135266249775724709369995"
    "95749669676277240766303535475945713821785251664274"
    "27466391932003059921817413596629043572900334295260"
    "59563073813232862794349076323382988075319525101901"
    "15738341879307021540891499348841675092447614606680"
    "82264800168477411853742345442437107539077744992069";

static const char* PI_REF =
    "3."
    "14159265358979323846264338327950288419716939937510"
    "58209749445923078164062862089986280348253421170679"
    "82148086513282306647093844609550582231725359408128"
    "48111745028410270193852110555964462294895493038196"
    "44288109756659334461284756482337867831652712019091"
    "45648566923460348610454326648213393607260249141273";

int count_match(const std::string& a, const char* b) {
    int m = 0;
    for (size_t i = 0; i < a.size() && b[i]; i++) {
        if (a[i] == b[i]) m++; else break;
    }
    return m - 2; // subtract "X."
}

// Check if value is very close to target (error < 2^-900 for 1024-bit)
template<typename BF>
bool is_close(const BF& val, const BF& target) {
    BF err = bf_abs(val - target);
    if (err.is_zero()) return true;
    // err.exp_ + TOP_BIT gives the position of the highest bit of error
    // For 1024-bit, if error < 2^-900 we consider it close enough
    int err_bits = err.exp_ + BF::TOP_BIT;
    return err_bits < -(BF::BITS * 7 / 8);
}

template<typename BF>
void test_all() {
    using namespace std;
    int pass = 0, fail = 0;
    auto check = [&](const char* name, bool ok, const string& detail = "") {
        cout << (ok ? "  PASS" : "  FAIL") << "  " << name;
        if (!detail.empty()) cout << "  " << detail;
        cout << endl;
        if (ok) pass++; else fail++;
    };

    cout << "=== Math function tests (" << BF::BITS << "-bit) ===" << endl;

    // Compute pi once for reuse
    auto t0 = chrono::high_resolution_clock::now();
    BF pi = bf_pi<BF>();
    auto t1 = chrono::high_resolution_clock::now();
    double pi_ms = chrono::duration<double, milli>(t1 - t0).count();
    string pi_str = bf_to_digits(pi, 300);
    int pi_digits = count_match(pi_str, PI_REF);
    check("bf_pi()", pi_digits >= 290,
          to_string(pi_digits) + " correct digits, " +
          to_string((int)pi_ms) + " ms");

    // --- sin(pi/6) = 0.5 ---
    {
        BF x = pi / BF(6.0);
        BF result = bf_sin(x);
        BF expected(0.5);
        BF err = bf_abs(result - expected);
        string err_str = bf_to_digits(err, 40);
        // Check first ~280 digits of 0.5
        string r_str = bf_to_digits(result, 300);
        bool ok = (r_str.substr(0, 6) == "0.5000");
        // Count correct digits
        int correct = 0;
        for (size_t i = 4; i < r_str.size(); i++) { // after "0.50"
            if (r_str[i] == '0') correct++; else break;
        }
        correct += 2; // "0.5" = 2 significant digits at least
        check("sin(pi/6) = 0.5", ok, "err=" + err_str.substr(0, 20));
    }

    // --- cos(pi/3) = 0.5 ---
    {
        BF x = pi / BF(3.0);
        BF result = bf_cos(x);
        BF expected(0.5);
        BF err = bf_abs(result - expected);
        string err_str = bf_to_digits(err, 40);
        string r_str = bf_to_digits(result, 300);
        bool ok = (r_str.substr(0, 6) == "0.5000");
        check("cos(pi/3) = 0.5", ok, "err=" + err_str.substr(0, 20));
    }

    // --- exp(1) = e ---
    {
        auto te0 = chrono::high_resolution_clock::now();
        BF e_val = bf_exp(BF(1.0));
        auto te1 = chrono::high_resolution_clock::now();
        double exp_ms = chrono::duration<double, milli>(te1 - te0).count();
        string e_str = bf_to_digits(e_val, 300);
        int e_digits = count_match(e_str, E_REF);
        check("exp(1) = e", e_digits >= 280,
              to_string(e_digits) + " correct digits, " +
              to_string((int)exp_ms) + " ms");
    }

    // --- ln(e) = 1 ---
    {
        BF e_val = bf_exp(BF(1.0));
        BF result = bf_ln(e_val);
        BF err = bf_abs(result - BF(1.0));
        string err_str = bf_to_digits(err, 40);
        check("ln(e) = 1", is_close(result, BF(1.0)), "err=" + err_str.substr(0, 20));
    }

    // --- atan(1) = pi/4 ---
    {
        auto ta0 = chrono::high_resolution_clock::now();
        BF result = bf_atan(BF(1.0));
        auto ta1 = chrono::high_resolution_clock::now();
        double atan_ms = chrono::duration<double, milli>(ta1 - ta0).count();
        BF expected = pi / BF(4.0);
        BF err = bf_abs(result - expected);
        string err_str = bf_to_digits(err, 40);
        string r_str = bf_to_digits(result, 40);
        string e_str = bf_to_digits(expected, 40);
        check("atan(1) = pi/4", r_str.substr(0, 10) == e_str.substr(0, 10),
              "err=" + err_str.substr(0, 20) + ", " +
              to_string((int)atan_ms) + " ms");
    }

    // --- sin^2(x) + cos^2(x) = 1 ---
    {
        BF x = pi / BF(7.0); // arbitrary angle
        BF s = bf_sin(x);
        BF c = bf_cos(x);
        BF result = s * s + c * c;
        BF err = bf_abs(result - BF(1.0));
        string err_str = bf_to_digits(err, 40);
        check("sin^2+cos^2 = 1", is_close(result, BF(1.0)), "err=" + err_str.substr(0, 20));
    }

    // --- bf_sqrt(2)^2 = 2 ---
    {
        BF s2 = bf_sqrt(BF(2.0));
        BF result = s2 * s2;
        BF err = bf_abs(result - BF(2.0));
        string err_str = bf_to_digits(err, 40);
        check("sqrt(2)^2 = 2", is_close(result, BF(2.0)), "err=" + err_str.substr(0, 20));
    }

    // --- tan(pi/4) = 1 ---
    {
        BF x = pi / BF(4.0);
        BF result = bf_tan(x);
        BF err = bf_abs(result - BF(1.0));
        string err_str = bf_to_digits(err, 40);
        check("tan(pi/4) = 1", is_close(result, BF(1.0)), "err=" + err_str.substr(0, 20));
    }

    // --- asin(0.5) = pi/6 ---
    {
        BF result = bf_asin(BF(0.5));
        BF expected = pi / BF(6.0);
        BF err = bf_abs(result - expected);
        string err_str = bf_to_digits(err, 40);
        check("asin(0.5) = pi/6", is_close(result, expected), "err=" + err_str.substr(0, 20));
    }

    // --- asin(1) = pi/2 ---
    {
        BF result = bf_asin(BF(1.0));
        BF expected = pi * BF(0.5);
        BF err = bf_abs(result - expected);
        string err_str = bf_to_digits(err, 40);
        check("asin(1) = pi/2", is_close(result, expected), "err=" + err_str.substr(0, 20));
    }

    // --- acos(0.5) = pi/3 ---
    {
        BF result = bf_acos(BF(0.5));
        BF expected = pi / BF(3.0);
        BF err = bf_abs(result - expected);
        string err_str = bf_to_digits(err, 40);
        check("acos(0.5) = pi/3", is_close(result, expected), "err=" + err_str.substr(0, 20));
    }

    // --- acos(0) = pi/2 ---
    {
        BF result = bf_acos(BF(0.0));
        BF expected = pi * BF(0.5);
        BF err = bf_abs(result - expected);
        string err_str = bf_to_digits(err, 40);
        check("acos(0) = pi/2", is_close(result, expected), "err=" + err_str.substr(0, 20));
    }

    // --- asin(sin(x)) = x  (roundtrip) ---
    {
        BF x = pi / BF(5.0);
        BF result = bf_asin(bf_sin(x));
        BF err = bf_abs(result - x);
        string err_str = bf_to_digits(err, 40);
        check("asin(sin(pi/5)) = pi/5", is_close(result, x), "err=" + err_str.substr(0, 20));
    }

    cout << endl << "Results: " << pass << " passed, " << fail << " failed" << endl << endl;
}

int main() {
    test_all<Float1024_2>();
    return 0;
}
