此文件夹下的内容, 都是 AI 生成. 初衷是早就想用短 bit (比如2/4/8) 的 int 类型的 +-*/ 来实现(且不用 FFT)一个长 bit (128/256/etc)的浮点数, 然后用这个长 bit float 来实现较高精度的基本函数, pi 等的计算(权当一玩). 但是真要做, 还是很费时间. 现在有 ai 帮助, 终于了此宿缘. 

而为啥想这样做呢? 单纯为了体验下假如有一个极致的 2bit CPU, 要用软件的方式实现高精度浮点计算, 会是如何(当然, 可行性毫无问题; 所以只是为了体验下, 单纯为了玩).

# BigFloat：任意精度浮点数与数学函数实现

## 数据表示

BigFloat 的设计类似 IEEE 754 浮点数，但尾数宽度可以任意大：

```
值 = (-1)^sign × mantissa × 2^exp
```

```cpp
template<typename T, int LIMBS>
class BigFloat {
    bool     sign_;              // 符号位
    int32_t  exp_;               // 二进制指数
    T        mantissa_[LIMBS];   // 尾数，由 LIMBS 个 T 类型的 limb 组成
    Special  special_;           // 特殊值：ZERO, NORMAL, INF, NAN_
};
```

**归一化**：和 IEEE 浮点数的隐含前导 1 类似，mantissa 始终保持最高位为 1。每次运算后调用 `normalize()` 左移/右移尾数并调整指数。

**Limb 类型的魔法**：通过模板参数 `T`，同一套代码可以用 uint64_t（64位 limb）也可以用 uint2_t（2位 limb）。关键在于 `DoubleWidth<T>` 提供双倍宽度类型来接住进位和乘积：

```
uint2_t  × uint2_t  → uint4_t      (2位 ALU)
uint64_t × uint64_t → __uint128_t  (64位 ALU)
```

所以即使是 2 位的"CPU"，也能完成 1024 位精度的运算——只是需要 512 个 limb，运算量是 64 位版本的约 1000 倍。

---

## 四则运算

### 加法 / 减法

和**手算**浮点加法完全一样：

**第一步：对齐指数**

两个操作数的指数可能不同。将指数较小的那个尾数右移，使两者指数相等：

```
  1.0110 × 2^5
+ 1.1010 × 2^3  →  右移2位  →  0.0110 × 2^5
```

```cpp
int ed = a.exp_ - b.exp_;   // 指数差
m_shift_right(bm, ed);       // 将较小数的尾数右移 ed 位
```

**第二步：尾数加减**

- **同号** → 尾数相加。若产生进位，右移1位，指数+1
- **异号** → 大减小。结果可能前导零，左移归一化

```cpp
if (同号) {
    m_add(a, b, result, carry);         // limb逐个相加，传递进位链
    if (carry) { shift_right(1); exp_++; }
} else {
    m_sub(bigger, smaller, result);     // limb逐个相减，传递借位链
    normalize();                         // 左移消除前导零
}
```

逐 limb 加法的核心——每一步用双倍宽度类型 T2 来捕获进位：

```cpp
for (int i = 0; i < LIMBS; i++) {
    T2 s = (T2)a[i] + (T2)b[i] + (T2)carry;
    r[i] = (T)s;            // 低位存结果
    carry = T(s >> TBITS);   // 高位就是进位
}
```

### 乘法

小学**竖式**乘法，在 limb 层面执行：

```
     a[3] a[2] a[1] a[0]
  ×  b[3] b[2] b[1] b[0]
  ────────────────────────
  → full[0..7]  (双倍宽度的积)
```

```cpp
for (i = 0..LIMBS-1)
    for (j = 0..LIMBS-1)
        full[i+j] += a[i] × b[j]     // T2 双倍宽度接住乘积
```

得到 2×LIMBS 个 limb 的完整乘积后，从最高有效位开始截取 LIMBS 个 limb 作为结果尾数，指数相加。

复杂度 **O(LIMBS²)**。对于 uint2_t 的 512 个 limb，就是 512×512 = 262,144 次微型乘加。

> 这正是后续引入 FFT 的动机——FFT 乘法可以把 O(N²) 降到 O(N log N)。

### 除法

逐位恢复除法（Restoring Division），每次迭代处理 1 个二进制位：

```
对每一位 i（从 2×BITS-1 到 0）:
    余数左移 1 位
    从被除数取下一位放入余数最低位
    若 余数 ≥ 除数:
        余数 -= 除数
        商的第 i 位 = 1
```

```cpp
for (int i = 2 * BITS - 1; i >= 0; i--) {
    rem <<= 1;                    // 余数左移
    rem[0] |= 被除数第i位;         // 取下一位
    if (rem >= div) {
        rem -= div;               // 减去除数
        q[i] = 1;                 // 商位置1
    }
}
```

产生 2×BITS 位精度的商，截取高 BITS 位作为结果。指数相减。

这是最昂贵的运算——每一位都需要 O(LIMBS) 的比较和减法，总复杂度 **O(BITS × LIMBS)**。

---

## 数学函数

所有数学函数都建立在四则运算之上，核心思想是两个：
1. **参数规约（Argument Reduction）**：把输入变小，让级数快速收敛
2. **Taylor 级数**：用递推公式逐项累加，避免计算大数阶乘

### bf_sqrt(x) — 牛顿迭代

```
x_{n+1} = (x_n + S/x_n) / 2
```

初始值用 double 精度的 `std::sqrt` 估算，然后迭代 12 次。牛顿法每次迭代有效精度翻倍（二次收敛），12 次迭代 → 2^12 = 4096 位精度，足够覆盖 1024 位。

### bf_pi() — Chudnovsky 公式

```
1/π = 12 × Σ (-1)^k × (6k)! × (13591409 + 545140134k) / ((3k)! × (k!)^3 × 640320^(3k+3/2))
```

每项贡献约 14.18 个十进制位。1024 位二进制 ≈ 308 位十进制 → 只需约 22 项。

### bf_exp(x) — 参数规约 + Taylor

**问题**：Taylor 级数 `e^x = 1 + x + x²/2! + x³/3! + ...` 当 x 很大时收敛极慢。

**解决**：

第一步——**参数规约**：反复将 x 除以 2，直到 |x| < 1：

```
x → x/2 → x/4 → ... → x/2^n （其中 |x/2^n| < 1）
```

第二步——**Taylor 级数**：对小的 r = x/2^n 计算级数，用递推避免阶乘：

```
term_k = term_{k-1} × r / k
```

第三步——**平方还原**：将结果连续平方 n 次：

```
exp(x) = exp(r)^{2^n} = ((exp(r))²)²)²...
```

### bf_ln(x) — 牛顿迭代

利用 exp 反推 ln，牛顿迭代公式：

```
y_{n+1} = y_n + x × exp(-y_n) - 1
```

初始值用 `std::log` 估算，迭代 14 次达到全精度。每次迭代调用一次 bf_exp，所以 ln 比 exp 慢得多。

### bf_sin(x) / bf_cos(x) — 参数规约 + Taylor

**第一步——规约到 [0, 2π)**：

```
x → x mod 2π    （利用高精度 π）
```

**第二步——象限折叠到 [0, π/2]**：

```
sin: 若 x > π → sin(x) = -sin(x-π)
     若 x > π/2 → sin(x) = sin(π-x)

cos: 若 x > π → cos(x) = cos(2π-x)
     若 x > π/2 → cos(x) = -cos(π-x)
```

**第三步——Taylor 级数**（此时 x ∈ [0, π/2]，级数收敛快）：

```
sin(x) = x - x³/3! + x⁵/5! - ...
    递推: term *= -x² / ((2k)(2k+1))

cos(x) = 1 - x²/2! + x⁴/4! - ...
    递推: term *= -x² / ((2k-1)(2k))
```

迭代次数 ≈ BITS/2（1024位约 500 次），直到 term 小到不影响累加和。

### bf_atan(x) — 三级参数规约 + Taylor

atan 的 Taylor 级数 `x - x³/3 + x⁵/5 - ...` 仅在 |x| < 1 时收敛，且 |x| 越小越快。

**规约1**：若 |x| > 1，用恒等式翻转：

```
atan(x) = π/2 - atan(1/x)    →  现在 |参数| < 1
```

**规约2**：若 |x| > 0.5，用半角公式缩小：

```
atan(x) = 2 × atan(x / (1 + √(1+x²)))    →  参数约缩小一半
```

**Taylor 级数**（|x| < 0.5，收敛快）：

```
atan(x) = x - x³/3 + x⁵/5 - ...
    递推: term *= -x² × (2k-1) / (2k+1)
```

### bf_asin(x) / bf_acos(x) — 归结为 atan

利用恒等式直接复用 atan 和 sqrt：

```
asin(x) = atan(x / √(1-x²))       （|x| < 1）
asin(1) = π/2                       （特殊情况）

acos(x) = π/2 - asin(x)
```

### bf_tan(x)

最简单的一个：

```
tan(x) = sin(x) / cos(x)
```

---

## 各函数的依赖关系

```
           bf_pi (Chudnovsky)
            │
            ├──── bf_sqrt (Newton)
            │
            ▼
bf_exp ◄──── bf_ln (Newton, 调用 exp)
  │
  ▼
bf_sin ◄──── bf_pi (参数规约需要 π)
bf_cos ◄──── bf_pi
  │
  ▼
bf_tan = sin / cos
  │
  ▼
bf_atan ◄── bf_pi + bf_sqrt (规约)
  │
  ▼
bf_asin = atan + sqrt
bf_acos = π/2 - asin
```

## 性能（1024-bit, Float1024_2, 2-bit limbs）-- Mac 上测试

| 函数 | 耗时 | 说明 |
|------|------|------|
| pi | ~100 ms | Chudnovsky 22项 |
| exp(1) | ~75 ms | 规约 + ~500项 Taylor |
| atan(1) | ~254 ms | 三级规约 + Taylor |
| sin, cos | ~200 ms | 含 pi 计算 + Taylor |
| ln(e) | ~1 s | 14次牛顿迭代，每次调用 exp |

用 uint64_t limb（Float1024）则快约 20-50 倍。

(看来 2-bit 也挺快)
