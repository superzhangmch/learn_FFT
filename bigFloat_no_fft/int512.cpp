#include "int512.h"
#include <cstring>
#include <algorithm>
#include <stdexcept>

// --- Constructors ---

Int512::Int512() {
    memset(data, 0, sizeof(data));
}

Int512::Int512(int val) : Int512(int64_t(val)) {}

Int512::Int512(int64_t val) {
    if (val >= 0) {
        data[0] = static_cast<uint64_t>(val);
        for (int i = 1; i < 8; i++) data[i] = 0;
    } else {
        // Sign extend: fill with 0xFFFFFFFFFFFFFFFF
        data[0] = static_cast<uint64_t>(val);
        for (int i = 1; i < 8; i++) data[i] = ~uint64_t(0);
    }
}

Int512::Int512(const char* decimal_str) {
    memset(data, 0, sizeof(data));
    if (!decimal_str || *decimal_str == '\0') return;

    bool negative = false;
    const char* p = decimal_str;
    if (*p == '-') { negative = true; p++; }
    else if (*p == '+') { p++; }

    // Parse decimal digits: result = result * 10 + digit
    for (; *p != '\0'; p++) {
        if (*p < '0' || *p > '9')
            throw std::invalid_argument("Int512: invalid character in decimal string");
        int digit = *p - '0';

        // Multiply by 10
        uint64_t carry = 0;
        for (int i = 0; i < 8; i++) {
            __uint128_t prod = (__uint128_t)data[i] * 10 + carry;
            data[i] = (uint64_t)prod;
            carry = (uint64_t)(prod >> 64);
        }

        // Add digit
        uint64_t c = digit;
        for (int i = 0; i < 8 && c; i++) {
            __uint128_t sum = (__uint128_t)data[i] + c;
            data[i] = (uint64_t)sum;
            c = (uint64_t)(sum >> 64);
        }
    }

    if (negative && !is_zero()) {
        // Negate (two's complement)
        *this = -(*this);
    }
}

Int512::Int512(const std::string& decimal_str) : Int512(decimal_str.c_str()) {}

// --- Helpers ---

bool Int512::is_negative() const {
    return (data[7] >> 63) & 1;
}

bool Int512::is_zero() const {
    for (int i = 0; i < 8; i++)
        if (data[i] != 0) return false;
    return true;
}

int Int512::get_bit(int pos) const {
    if (pos < 0 || pos >= 512) return 0;
    return (data[pos / 64] >> (pos % 64)) & 1;
}

void Int512::set_bit(int pos, int val) {
    if (pos < 0 || pos >= 512) return;
    if (val)
        data[pos / 64] |= (uint64_t(1) << (pos % 64));
    else
        data[pos / 64] &= ~(uint64_t(1) << (pos % 64));
}

int Int512::highest_bit() const {
    for (int i = 7; i >= 0; i--) {
        if (data[i] != 0) {
            // Find highest set bit in this limb
            int bit = 63;
            while (bit >= 0 && !((data[i] >> bit) & 1)) bit--;
            return i * 64 + bit;
        }
    }
    return -1;
}

Int512 Int512::abs() const {
    if (is_negative()) return -(*this);
    return *this;
}

// --- Unary ---

Int512 Int512::operator-() const {
    Int512 result = ~(*this);
    // Add 1
    uint64_t carry = 1;
    for (int i = 0; i < 8; i++) {
        __uint128_t sum = (__uint128_t)result.data[i] + carry;
        result.data[i] = (uint64_t)sum;
        carry = (uint64_t)(sum >> 64);
    }
    return result;
}

Int512 Int512::operator~() const {
    Int512 result;
    for (int i = 0; i < 8; i++)
        result.data[i] = ~data[i];
    return result;
}

// --- Addition ---

Int512 Int512::operator+(const Int512& rhs) const {
    Int512 result;
    uint64_t carry = 0;
    for (int i = 0; i < 8; i++) {
        __uint128_t sum = (__uint128_t)data[i] + rhs.data[i] + carry;
        result.data[i] = (uint64_t)sum;
        carry = (uint64_t)(sum >> 64);
    }
    return result;
}

Int512 Int512::operator-(const Int512& rhs) const {
    return *this + (-rhs);
}

// --- Multiplication ---

Int512 Int512::operator*(const Int512& rhs) const {
    // Schoolbook multiplication, keeping only lower 512 bits
    Int512 result;
    for (int i = 0; i < 8; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 8 - i; j++) {
            __uint128_t prod = (__uint128_t)data[i] * rhs.data[j]
                             + result.data[i + j] + carry;
            result.data[i + j] = (uint64_t)prod;
            carry = (uint64_t)(prod >> 64);
        }
    }
    return result;
}

// --- Unsigned comparison ---

int Int512::unsigned_compare(const Int512& a, const Int512& b) {
    for (int i = 7; i >= 0; i--) {
        if (a.data[i] < b.data[i]) return -1;
        if (a.data[i] > b.data[i]) return 1;
    }
    return 0;
}

// --- Unsigned division (binary long division) ---

void Int512::unsigned_divmod(const Int512& dividend, const Int512& divisor,
                             Int512& quotient, Int512& remainder) {
    quotient = Int512(0);
    remainder = Int512(0);

    if (divisor.is_zero())
        throw std::domain_error("Int512: division by zero");

    if (unsigned_compare(dividend, divisor) < 0) {
        remainder = dividend;
        return;
    }

    int dividend_bits = dividend.highest_bit();
    for (int i = dividend_bits; i >= 0; i--) {
        // Left shift remainder by 1
        remainder = remainder << 1;
        // Set lowest bit of remainder to bit i of dividend
        remainder.set_bit(0, dividend.get_bit(i));
        // If remainder >= divisor, subtract
        if (unsigned_compare(remainder, divisor) >= 0) {
            remainder = remainder - divisor;
            quotient.set_bit(i, 1);
        }
    }
}

// --- Signed division ---

Int512 Int512::operator/(const Int512& rhs) const {
    if (rhs.is_zero())
        throw std::domain_error("Int512: division by zero");

    bool result_negative = is_negative() != rhs.is_negative();
    Int512 a = this->abs();
    Int512 b = rhs.abs();

    Int512 q, r;
    unsigned_divmod(a, b, q, r);

    return result_negative ? -q : q;
}

Int512 Int512::operator%(const Int512& rhs) const {
    if (rhs.is_zero())
        throw std::domain_error("Int512: modulo by zero");

    bool dividend_negative = is_negative();
    Int512 a = this->abs();
    Int512 b = rhs.abs();

    Int512 q, r;
    unsigned_divmod(a, b, q, r);

    // Remainder has same sign as dividend (C++ convention)
    return dividend_negative ? -r : r;
}

// --- Compound assignment ---

Int512& Int512::operator+=(const Int512& rhs) { *this = *this + rhs; return *this; }
Int512& Int512::operator-=(const Int512& rhs) { *this = *this - rhs; return *this; }
Int512& Int512::operator*=(const Int512& rhs) { *this = *this * rhs; return *this; }
Int512& Int512::operator/=(const Int512& rhs) { *this = *this / rhs; return *this; }
Int512& Int512::operator%=(const Int512& rhs) { *this = *this % rhs; return *this; }

// --- Shift ---

Int512 Int512::operator<<(int shift) const {
    if (shift <= 0) return *this;
    if (shift >= 512) return Int512(0);

    Int512 result;
    int limb_shift = shift / 64;
    int bit_shift = shift % 64;

    for (int i = 7; i >= 0; i--) {
        int src = i - limb_shift;
        if (src >= 0) {
            result.data[i] = data[src] << bit_shift;
            if (bit_shift > 0 && src - 1 >= 0)
                result.data[i] |= data[src - 1] >> (64 - bit_shift);
        }
    }
    return result;
}

Int512 Int512::operator>>(int shift) const {
    if (shift <= 0) return *this;
    if (shift >= 512) {
        // Arithmetic shift: fill with sign bit
        return is_negative() ? Int512(-1) : Int512(0);
    }

    Int512 result;
    int limb_shift = shift / 64;
    int bit_shift = shift % 64;
    bool negative = is_negative();

    // Fill upper limbs with sign extension
    uint64_t fill = negative ? ~uint64_t(0) : 0;

    for (int i = 0; i < 8; i++) {
        int src = i + limb_shift;
        if (src < 8) {
            result.data[i] = data[src] >> bit_shift;
            if (bit_shift > 0) {
                uint64_t upper = (src + 1 < 8) ? data[src + 1] : fill;
                result.data[i] |= upper << (64 - bit_shift);
            }
        } else {
            result.data[i] = fill;
        }
    }
    return result;
}

// --- Comparison ---

bool Int512::operator==(const Int512& rhs) const {
    for (int i = 0; i < 8; i++)
        if (data[i] != rhs.data[i]) return false;
    return true;
}

bool Int512::operator!=(const Int512& rhs) const {
    return !(*this == rhs);
}

bool Int512::operator<(const Int512& rhs) const {
    bool a_neg = is_negative();
    bool b_neg = rhs.is_negative();

    if (a_neg && !b_neg) return true;   // negative < positive
    if (!a_neg && b_neg) return false;   // positive > negative

    // Same sign: compare from most significant limb
    for (int i = 7; i >= 0; i--) {
        if (data[i] < rhs.data[i]) return true;
        if (data[i] > rhs.data[i]) return false;
    }
    return false; // equal
}

bool Int512::operator<=(const Int512& rhs) const { return !(rhs < *this); }
bool Int512::operator>(const Int512& rhs) const { return rhs < *this; }
bool Int512::operator>=(const Int512& rhs) const { return !(*this < rhs); }

// --- String conversion ---

std::string Int512::to_string() const {
    if (is_zero()) return "0";

    bool negative = is_negative();
    Int512 val = negative ? -(*this) : *this;

    std::string digits;
    Int512 ten(10);

    while (!val.is_zero()) {
        Int512 q, r;
        unsigned_divmod(val, ten, q, r);
        digits.push_back('0' + (int)r.data[0]);
        val = q;
    }

    if (negative) digits.push_back('-');
    std::reverse(digits.begin(), digits.end());
    return digits;
}

std::ostream& operator<<(std::ostream& os, const Int512& val) {
    return os << val.to_string();
}
