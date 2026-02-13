#include <iostream>
#include <string>
#include <cctype>
#include <cstdlib>
#include "bfmath.h"

using BF = Float1024_2;

// ============================================================
// Recursive descent expression parser
//
// Grammar:
//   expr   = term (('+' | '-') term)*
//   term   = unary (('*' | '/') unary)*
//   unary  = '-' unary | power
//   power  = atom ('^' atom)?
//   atom   = NUMBER | FUNC '(' expr ')' | '(' expr ')' | 'pi' | 'e'
// ============================================================

struct Parser {
    const char* s;
    int pos;

    Parser(const char* input) : s(input), pos(0) {}

    void skip_spaces() {
        while (s[pos] == ' ' || s[pos] == '\t') pos++;
    }

    char peek() { skip_spaces(); return s[pos]; }
    char get()  { skip_spaces(); return s[pos++]; }

    bool match(char c) {
        if (peek() == c) { pos++; return true; }
        return false;
    }

    // Parse a decimal number: [0-9]+ ('.' [0-9]*)?
    BF parse_number() {
        skip_spaces();
        int start = pos;
        while (std::isdigit(s[pos])) pos++;
        if (s[pos] == '.') {
            pos++;
            while (std::isdigit(s[pos])) pos++;
        }
        std::string num_str(s + start, s + pos);

        // Convert decimal string to BF
        // Split into integer and fractional parts
        size_t dot = num_str.find('.');
        std::string int_part = (dot == std::string::npos) ? num_str : num_str.substr(0, dot);
        std::string frac_part = (dot == std::string::npos) ? "" : num_str.substr(dot + 1);

        // Integer part
        BF result(0.0);
        BF ten(10.0);
        for (char c : int_part) {
            result = result * ten + BF((double)(c - '0'));
        }

        // Fractional part
        if (!frac_part.empty()) {
            BF frac(0.0);
            BF place(1.0);
            for (char c : frac_part) {
                place = place / ten;
                frac = frac + place * BF((double)(c - '0'));
            }
            result = result + frac;
        }
        return result;
    }

    // Parse an identifier (function name or constant)
    std::string parse_ident() {
        skip_spaces();
        int start = pos;
        while (std::isalpha(s[pos]) || s[pos] == '_') pos++;
        return std::string(s + start, s + pos);
    }

    BF parse_expr() {
        BF left = parse_term();
        while (true) {
            if (match('+'))      left = left + parse_term();
            else if (match('-')) left = left - parse_term();
            else break;
        }
        return left;
    }

    BF parse_term() {
        BF left = parse_unary();
        while (true) {
            if (match('*'))      left = left * parse_unary();
            else if (match('/')) left = left / parse_unary();
            else break;
        }
        return left;
    }

    BF parse_unary() {
        if (match('-')) return -parse_unary();
        return parse_power();
    }

    BF parse_power() {
        BF base = parse_atom();
        if (match('^')) {
            BF exp_val = parse_atom();
            // a^b = exp(b * ln(a))
            return bf_exp(exp_val * bf_ln(base));
        }
        return base;
    }

    BF call_func(const std::string& name, const BF& arg) {
        if (name == "sin")  return bf_sin(arg);
        if (name == "cos")  return bf_cos(arg);
        if (name == "tan")  return bf_tan(arg);
        if (name == "asin") return bf_asin(arg);
        if (name == "acos") return bf_acos(arg);
        if (name == "atan") return bf_atan(arg);
        if (name == "exp")  return bf_exp(arg);
        if (name == "ln")   return bf_ln(arg);
        if (name == "log")  return bf_ln(arg);
        if (name == "sqrt") return bf_sqrt(arg);
        if (name == "abs")  return bf_abs(arg);
        std::cerr << "Unknown function: " << name << std::endl;
        return BF::nan();
    }

    BF parse_atom() {
        skip_spaces();

        // Number
        if (std::isdigit(s[pos]) || (s[pos] == '.' && std::isdigit(s[pos + 1]))) {
            return parse_number();
        }

        // Parenthesized expression
        if (match('(')) {
            BF val = parse_expr();
            if (!match(')')) std::cerr << "Missing ')'" << std::endl;
            return val;
        }

        // Identifier: function or constant
        if (std::isalpha(s[pos])) {
            std::string id = parse_ident();
            if (id == "pi")  return bf_pi<BF>();
            if (id == "e")   return bf_exp(BF(1.0));
            // Function call
            if (!match('(')) {
                std::cerr << "Expected '(' after " << id << std::endl;
                return BF::nan();
            }
            BF arg = parse_expr();
            if (!match(')')) std::cerr << "Missing ')'" << std::endl;
            return call_func(id, arg);
        }

        std::cerr << "Unexpected character: '" << s[pos] << "'" << std::endl;
        pos++;
        return BF::nan();
    }
};

int main() {
    std::cout << "BigFloat Calculator (1024-bit, 2-bit limbs)" << std::endl;
    std::cout << "Functions: sin cos tan asin acos atan exp ln sqrt abs" << std::endl;
    std::cout << "Constants: pi e" << std::endl;
    std::cout << "Operators: + - * / ^ ()" << std::endl;
    std::cout << "Type 'quit' to exit." << std::endl << std::endl;

    // Precompute pi and e
    std::cout << "Initializing... " << std::flush;
    BF pi_val = bf_pi<BF>();
    std::cout << "ready." << std::endl << std::endl;

    std::string line;
    while (true) {
        std::cout << "> " << std::flush;
        if (!std::getline(std::cin, line)) break;
        if (line.empty()) continue;
        if (line == "quit" || line == "exit" || line == "q") break;

        Parser parser(line.c_str());
        BF result = parser.parse_expr();

        if (result.is_nan()) {
            std::cout << "  NaN" << std::endl;
        } else if (result.is_inf()) {
            std::cout << (result.is_negative() ? "  -Inf" : "  Inf") << std::endl;
        } else {
            std::cout << "  " << bf_to_digits(result, 300) << std::endl;
        }
    }

    return 0;
}
