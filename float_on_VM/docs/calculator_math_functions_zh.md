# 科学计算器函数实现说明（面向 int16 VM）

本文面向你的目标：在仅支持 `int16` 指令的 VM 上，实现接近科学计算器精度的 `sin/cos/tan/asin/acos/atan/log/exp/pow/sqrt`。

## 1. 目标与约束

- 硬件模型：16-bit 寄存器 + 小 RAM + ROM 固件（`asm -> bin -> load -> run`）
- VM 原语：整数加减乘除、位运算、分支、循环
- 浮点能力：不依赖宿主语言浮点，采用软件浮点（建议 binary64/Soft-FP64）
- 精度目标：通常函数达到 10-12 位有效数字（范围相关）

## 2. 总体架构

1. 表达式层：词法/语法解析，生成计算序列（可选）
2. 数值层：Soft-FP64（加减乘除、比较、规格化、舍入）
3. 函数层：`sin/cos/log/exp...` 的近似算法
4. 固件层：汇编程序与查表常量（ROM）

关键原则：
- 所有循环和分支最终都应落到 VM 指令上
- 高层函数只组合 Soft-FP64 原语，不直接“偷用”宿主浮点

## 3. 三类常用近似方法

## 3.1 CORDIC

特点：
- 核心只用加减、移位、查表
- 适合没有高性能乘法器的硬件

圆模式旋转（求 sin/cos）典型迭代：
- `x_{i+1} = x_i - d_i * y_i * 2^-i`
- `y_{i+1} = y_i + d_i * x_i * 2^-i`
- `z_{i+1} = z_i - d_i * atan(2^-i)`
- `d_i = sign(z_i)`

注意：
- 需处理增益 `K`（预缩放或结果补偿）
- 先做象限约简，再进入迭代
- 迭代次数决定精度

## 3.2 Minimax 多项式（Remez）

特点：
- 在给定区间上最小化最大误差（Chebyshev 意义下最优）
- 通常比同阶 Taylor 误差更小
  - 它的展开系数和 Taylor 系数, 数值上往往差不多

流程：
1. 选区间（如 `sin` 在 `[-pi/4, pi/4]`）
2. 用 Remez 求系数
3. 用 Horner/Estrin 在 VM 上求值

建议：
- 对 `sin/cos/exp/log/atan` 都可用 minimax
- 系数离线生成，固化到 ROM

## 3.3 Taylor 展开

特点：
- 推导简单、系数直观
- 但“同阶最坏误差”通常不如 minimax

用途：
- 作为基线实现/验证实现正确性
- 在很小区间内可接受

## 4. 每个函数怎么做

## 4.1 sin/cos/tan

推荐主线：
1. `x` 先做范围约简到小区间（非常关键）
2. 在核心区间用 `CORDIC` 或 `minimax` 计算 `sin/cos`
3. `tan = sin/cos`，并处理 `cos` 接近 0 的异常

工程建议：
- 若追求“硬件感”和少乘法：优先 CORDIC
- 若追求速度/ROM 平衡和高精度：优先 minimax

## 4.2 atan/asin/acos

- `atan`：
  - 先变换到小区间（如 `|x|>1` 用 `atan(x)=pi/2-atan(1/x)`）
  - 小区间用 minimax 或级数
- `asin`：`asin(x)=atan(x/sqrt(1-x^2))`
- `acos`：`acos(x)=pi/2-asin(x)`

## 4.3 log/exp

- `log(x)`：
  - 分解 `x = m * 2^e`
  - `log(x)=log(m)+e*ln2`
  - 对 `m` 在窄区间做多项式近似（minimax）
- `exp(x)`：
  - `x = k*ln2 + r`
  - `exp(x)=2^k * exp(r)`
  - `exp(r)` 用 minimax

## 4.4 sqrt/pow

- `sqrt(x)`：牛顿迭代（`g=(g+x/g)/2`）
- `pow(a,b)`：
  - 整数幂：平方-乘法（binary exponentiation）
  - 实数幂：`pow(a,b)=exp(b*log(a))`（要求 `a>0`）

## 5. CORDIC vs Minimax 怎么选

- CORDIC：
  - 优点：移位+加减友好，结构规则
  - 缺点：迭代轮数多，吞吐常不如多项式
- Minimax：
  - 优点：同阶精度高、步数少
  - 缺点：需要离线求系数，系数量化要小心

实务上常见混合：
- `sin/cos` 用 CORDIC（硬件友好）或 minimax（更快）
- `log/exp/atan` 常用 minimax + 变换

## 6. 在 VM 上落地时最容易踩坑

1. 范围约简不够精确（大输入误差会暴涨）
2. 系数量化位数不足
3. Soft-FP64 的舍入/规格化不完整
4. 特殊值漏处理（NaN/Inf/±0/次正规数）
5. 只看平均误差，不看 worst-case/ULP

## 7. 验证标准（建议）

- 分层测试：
  - Soft-FP64 原语：加减乘除、比较、转换
  - 函数核：核心区间误差
  - 端到端：随机输入 + 边界输入
- 指标：
  - 绝对误差、相对误差
  - ULP 误差（更适合浮点）
- 用例：
  - 小区间高精度
  - 大范围随机（验证范围约简）
  - 奇异点附近（如 `tan(pi/2)`、`log(0+)`）

## 8. 推荐实现路线（你这个项目）

1. 先把 Soft-FP64 原语做完整（含特殊值和舍入）
2. `sin/cos` 先上一个稳定版本（CORDIC 或 minimax 二选一）
3. `log/exp` 用“分解 + minimax”
4. `atan -> asin/acos` 复用
5. `pow/sqrt` 最后收口
6. 最后补大规模随机与 ULP 回归

---

如果你要，我可以再补一份“纯实现文档”：
- 每个函数具体到伪代码 + 所需 ROM 常量表格式 + VM 调用约定（输入输出寄存器/内存地址）。
