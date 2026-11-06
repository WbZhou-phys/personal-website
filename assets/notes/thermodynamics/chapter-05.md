---
tags:
  - 热统
  - 热力学
  - 第五章
created: 2026-10-22
---

# 第五章　不可逆过程热力学简介

> **教材**：汪志诚《热力学·统计物理》（第六版）
> **范围**：§5.1–§5.7（教材 pp. 112–133；其中 §5.3–§5.7 为带星号的选读内容）

## 本章地图

| 节 | 标题 | 一句话 |
| --- | --- | --- |
| §5.1 | 局域平衡　熵流密度与局域熵产生率 | 局域平衡假设；$\mathrm{d}S=\mathrm{d}_eS+\mathrm{d}_iS$；$\Theta=\sum_k\boldsymbol{J}_k\cdot\boldsymbol{X}_k$ |
| §5.2 | 线性与非线性过程　昂萨格倒易关系 | 线性唯象律 $J_k=\sum_l L_{kl}X_l$；$L_{kl}=L_{lk}$ |
| *§5.3 | 温差电现象 | 泽贝克、佩尔捷、汤姆孙效应；开尔文第一、第二关系 |
| *§5.4 | 最小熵产生定理 | 线性区非平衡定态是熵产生率极小且稳定的状态 |
| *§5.5 | 化学反应与扩散过程 | 局部化学亲和势 $a_\rho$；反应流与力的非线性关系 |
| *§5.6 | 非平衡系统在非线性区的发展判据 | $\mathrm{d}P/\mathrm{d}t=\mathrm{d}_XP/\mathrm{d}t+\mathrm{d}_JP/\mathrm{d}t$；非线性区符号不定 |
| *§5.7 | 三分子模型与耗散结构的概念 | Brusselator；线性稳定性分析；耗散结构与自组织 |

**主线**：§5.1 建立局域平衡与熵产生率 → §5.2 引入线性唯象律与昂萨格关系 → *§5.3 用温差电作范例 → *§5.4 得到线性区的最小熵产生定理 → *§5.5 转向反应–扩散（线性与非线性）→ *§5.6 给出非线性区的发展判据 → *§5.7 以三分子模型说明耗散结构的产生。

---

## §5.1　局域平衡　熵流密度与局域熵产生率

### 5.1.1 从平衡态热力学到不可逆过程

前面几章讲述了平衡态热力学。对于不可逆过程，从平衡态热力学只能得到非常有限的信息。例如，可以根据热力学函数的不等式判断过程的方向；如果不可逆过程的初态和终态都是平衡态，可以通过初态和终态热力学函数之间的关系求得整个过程的总效应；在不可逆性不是过程本质特征的情形下可以将过程近似地当作可逆过程处理，等等。

我们知道，自然界中存在大量的不可逆过程，例如热传导、扩散等输运过程、化学反应过程，乃至生命过程，等等。**在不可逆过程中系统经历一系列的非平衡态。** 将热力学方法推广到非平衡的情形，对不可逆过程本身进行研究无疑是重要和有意义的。

§1.16 中根据热力学第二定律得到下述不等式：

$$\mathrm{d}S \geqslant \frac{\mathrm{d}Q}{T} \tag{5.1.1}$$

式中等号适用于可逆过程，不等号适用于不可逆过程。可以将上式推广为下述等式：

$$\mathrm{d}S = \mathrm{d}_eS + \mathrm{d}_iS \tag{5.1.2}$$

式中 $\mathrm{d}_eS$ 是由于系统与外界交换物质和能量所引起的系统熵变，**它是可正可负的**；$\mathrm{d}_iS$ 表示系统内部发生的过程引起的**熵产生**，**它不会取负值**。

| 情形 | $\mathrm{d}_eS$ | $\mathrm{d}_iS$ | 结果 |
| --- | --- | --- | --- |
| 内部过程可逆 | — | $0$ | $\mathrm{d}S=\mathrm{d}_eS$ |
| 内部过程不可逆 | — | $>0$ | $\mathrm{d}S>\mathrm{d}_eS$ |
| 孤立系统 | $0$（无交换） | $\geqslant 0$ | $\mathrm{d}S=\mathrm{d}_iS\geqslant 0$，即熵增加原理 |
| 闭系 | $\mathrm{d}_eS=\dfrac{\mathrm{d}Q}{T}$ | $\geqslant 0$ | 回到式 (5.1.1) |
| 开系 | 热量交换 + 物质交换 | $\geqslant 0$ | $\mathrm{d}_eS$ 还包含物质交换的贡献 |

对于闭系，$\mathrm{d}_eS=\dfrac{\mathrm{d}Q}{T}$，就得到式 (5.1.1)，这时 $\mathrm{d}_eS$ 的正负取决于系统是吸热还是放热。对于开系，除热量交换外，系统与外界的物质交换也会引起 $\mathrm{d}_eS$。**为了建立不可逆过程热力学，需要计算各种不可逆过程的 $\mathrm{d}_iS$ 和 $\mathrm{d}_eS$。**

### 5.1.2 局域平衡假设与局域基本方程

在平衡态热力学中，我们得到了热力学基本方程 (4.1.12)。在本章中我们将该式改写为

$$T\,\mathrm{d}S = \mathrm{d}U + p\,\mathrm{d}V - \sum_i\mu_i\,\mathrm{d}N_i \tag{5.1.3}$$

式中 $N_i$ 是 $i$ 组元的分子数，相应地，$\mu_i$ 是一个 $i$ 分子的化学势。上式给出系统在两个相邻平衡态的熵、内能、体积和分子数之差的关系。

对于系统在不可逆过程中所经历的非平衡状态，我们限于讨论下述情形：**整个系统虽然处在非平衡状态，但如果将系统分成若干个小部分，使每一小部分仍然是含有大量粒子的宏观系统，由于各部分之间只通过界面区域的分子发生相互作用，且各小部分的弛豫时间比整个系统的弛豫时间要短得多，可以认为各小部分近似处在局域平衡的状态。** 在这种情形下，每一小部分的温度、压强、化学势、内能、熵、粒子数等就都有确定的意义。我们假设这些局域热力学量的改变仍然满足热力学基本方程式 (5.1.3)。

> **脚注**（教材 p.113）：在本章后面要讨论的问题中，可以略去基本方程中的 $p\,\mathrm{d}V$ 项。如果涉及流体力学问题，该项不能略去。

将全式除以局域体积可以得到联系局域熵密度 $s$、内能密度 $u$ 和粒子数密度 $n_i$ 的方程式：

$$T\,\mathrm{d}s = \mathrm{d}u - \sum_i\mu_i\,\mathrm{d}n_i \tag{5.1.4}$$

对于内能、熵、粒子数等广延量，整个系统的量可以表示为

$$U = \int u\,\mathrm{d}\tau, \qquad S = \int s\,\mathrm{d}\tau, \qquad N_i = \int n_i\,\mathrm{d}\tau \tag{5.1.5}$$

式中 $\mathrm{d}\tau$ 是体积元。**对于强度量（例如温度、化学势），系统不具有统一的数值。**

式 (5.1.4) 对于局域热力学量仍然成立，在不可逆过程热力学中是一个**假设**，其正确性由其推论与实际相符而得到肯定。统计物理学可以分析上述假设的正确性及其适用的限度。

> **脚注**（教材 p.113）：de Groot S R, Mazur P. *Non-equilibrium Thermodynamics*[M]. Amsterdam: North Holland Publishing Company, 1962; Chap 9.

### 5.1.3 熵流密度与局域熵产生率

在局域平衡情形下，可以将局域熵密度的增加率写成如下的形式：

$$\frac{\partial s}{\partial t} = -\nabla\cdot\boldsymbol{J}_s + \Theta \tag{5.1.6}$$

式中 $\boldsymbol{J}_s$ 是单位时间内流过单位截面的熵，称为**熵流密度**；$\Theta$ 是单位时间内单位体积中产生的熵，称为**局域熵产生率**。

根据式 (5.1.5)，整个系统熵的增加率可以表示为

$$\frac{\mathrm{d}S}{\mathrm{d}t} = \frac{\mathrm{d}}{\mathrm{d}t}\int s\,\mathrm{d}\tau = \int\frac{\partial s}{\partial t}\mathrm{d}\tau = \int\left(-\nabla\cdot\boldsymbol{J}_s+\Theta\right)\mathrm{d}\tau$$

利用高斯定理将右方第一项化为面积分，得

$$\frac{\mathrm{d}S}{\mathrm{d}t} = -\oint\boldsymbol{J}_s\cdot\mathrm{d}\boldsymbol{\sigma} + \int\Theta\,\mathrm{d}\tau \tag{5.1.7}$$

上式右方第一项表示单位时间内通过系统表面从外界流入的熵，第二项表示单位时间内系统各体积元的熵产生之和。与式 (5.1.2) 比较，可知

$$\frac{\mathrm{d}_eS}{\mathrm{d}t} = -\oint\boldsymbol{J}_s\cdot\mathrm{d}\boldsymbol{\sigma}, \qquad \frac{\mathrm{d}_iS}{\mathrm{d}t} = \int\Theta\,\mathrm{d}\tau \tag{5.1.8}$$

由于在任何宏观区域中熵产生都是正定的，故有

$$\Theta \geqslant 0$$

式 (5.1.6) 或式 (5.1.7) 只是一种形式的表示，需要对具体的不可逆过程求得熵流密度和局域熵产生率的具体表达式。下面介绍两个例子。

### 5.1.4 例一：单纯的热传导

当物体各处温度不均匀时，物体内部将发生热传导过程。先讨论单纯的热传导过程。考虑物体中一个固定的体积元，在单纯的热传导过程中，体积元中物质内能的增加是热量流入的结果。以 $\boldsymbol{J}_q$ 表示单位时间内流过单位截面的热量，称为**热流密度**，则内能密度的增加率为

$$\frac{\partial u}{\partial t} = -\nabla\cdot\boldsymbol{J}_q \tag{5.1.9}$$

式 (5.1.9) 表达了能量守恒定律。

对于单纯的热传导过程，式 (5.1.4) 简化为

$$T\,\mathrm{d}s = \mathrm{d}u \tag{5.1.10}$$

由上式得局域熵密度的增加率为

$$\frac{\partial s}{\partial t} = \frac{1}{T}\frac{\partial u}{\partial t} \tag{5.1.11}$$

将式 (5.1.9) 代入上式，得

$$\frac{\partial s}{\partial t} = -\frac{1}{T}\nabla\cdot\boldsymbol{J}_q$$

而

$$\frac{1}{T}\nabla\cdot\boldsymbol{J}_q = \nabla\cdot\frac{\boldsymbol{J}_q}{T} - \boldsymbol{J}_q\cdot\nabla\frac{1}{T}$$

故有

$$\frac{\partial s}{\partial t} = -\nabla\cdot\frac{\boldsymbol{J}_q}{T} + \boldsymbol{J}_q\cdot\nabla\frac{1}{T} \tag{5.1.12}$$

上式表明，局域熵密度的增加率可以分为两部分：$-\nabla\cdot\dfrac{\boldsymbol{J}_q}{T}$ 是由于热量流入而引起的局域熵密度的增加率；$\boldsymbol{J}_q\cdot\nabla\dfrac{1}{T}$ 是温度梯度导致的热传导过程所引起的局域熵密度的产生率。与式 (5.1.6) 比较，有

$$\boldsymbol{J}_s = \frac{\boldsymbol{J}_q}{T}, \qquad \Theta = \boldsymbol{J}_q\cdot\nabla\frac{1}{T} \tag{5.1.13}$$

前面说过，温度不均匀性是引起热传导的原因。定义

$$\boldsymbol{X}_q = \nabla\frac{1}{T}$$

称为**热流动力**。局域熵密度的产生率 $\Theta$ 可以表示为热流密度与热流动力的乘积：

$$\Theta = \boldsymbol{J}_q\cdot\boldsymbol{X}_q \tag{5.1.14}$$

假设热传导过程遵从**傅里叶定律**：

$$\boldsymbol{J}_q = -\kappa\nabla T$$

$\kappa$ 称为**导热系数**，则式 (5.1.14) 可以表示为

$$\Theta = \boldsymbol{J}_q\cdot\nabla\frac{1}{T} = -\boldsymbol{J}_q\cdot\frac{\nabla T}{T^2} = \kappa\frac{(\nabla T)^2}{T^2} \geqslant 0$$

由于导热系数 $\kappa$ 恒为正值，故热传导过程中局域熵产生率是正定的。

### 5.1.5 例二：热传导与物质输运

如果系统内部除了温度不均匀外，化学势也不均匀，则除了热传导外，还将有物质的输运。现在讨论同时存在热传导和物质输运时的熵流密度和局域熵产生率。

考虑物体中一个固定的体积元。体积元中粒子数密度 $n$ 的变化率满足**连续性方程**：

$$\frac{\partial n}{\partial t} + \nabla\cdot\boldsymbol{J}_n = 0 \tag{5.1.15}$$

其中 $\boldsymbol{J}_n$ 是**粒子流密度**。式 (5.1.15) 是物质守恒定律的表达式。

类似地，体积元中内能密度 $u$ 的变化率满足连续性方程：

$$\frac{\partial u}{\partial t} + \nabla\cdot\boldsymbol{J}_u = 0 \tag{5.1.16}$$

$\boldsymbol{J}_u$ 称为**内能流密度**。式 (5.1.16) 是能量守恒定律的表达式。

由式 (5.1.4) 可知，当粒子数密度增加 $\mathrm{d}n$ 时，内能密度的增加为 $\mu\,\mathrm{d}n$，$\mu$ 是一个粒子的化学势。因此，当存在粒子流时，内能流密度 $\boldsymbol{J}_u$ 可以表示为

$$\boldsymbol{J}_u = \boldsymbol{J}_q + \mu\boldsymbol{J}_n \tag{5.1.17}$$

即内能流密度是热流密度与粒子流所携带的能流密度之和。

将式 (5.1.17) 代入式 (5.1.16)，得

$$\frac{\partial u}{\partial t} = -\nabla\cdot\boldsymbol{J}_q - \nabla\cdot\left(\mu\boldsymbol{J}_n\right) \tag{5.1.18}$$

由式 (5.1.4) 知，局域熵密度的增加率为

$$\frac{\partial s}{\partial t} = \frac{1}{T}\frac{\partial u}{\partial t} - \frac{\mu}{T}\frac{\partial n}{\partial t} \tag{5.1.19}$$

将式 (5.1.15) 和式 (5.1.18) 代入上式，得

$$\begin{aligned}\frac{\partial s}{\partial t} &= -\frac{1}{T}\nabla\cdot\boldsymbol{J}_q - \frac{1}{T}\nabla\cdot\left(\mu\boldsymbol{J}_n\right) + \frac{\mu}{T}\nabla\cdot\boldsymbol{J}_n \\ &= -\nabla\cdot\left(\frac{\boldsymbol{J}_q}{T}\right) + \boldsymbol{J}_q\cdot\nabla\frac{1}{T} - \frac{\boldsymbol{J}_n}{T}\cdot\nabla\mu \end{aligned} \tag{5.1.20}$$

上式右方第一项是由于热量流入而引起的局域熵密度的增加率，第二项是温度梯度导致的热传导过程所引起的局域熵密度产生率，第三项是化学势梯度导致的物质输运过程所引起的局域熵密度产生率。将上式与式 (5.1.6) 比较，可知

$$\boldsymbol{J}_s = \frac{\boldsymbol{J}_q}{T}, \qquad \Theta = \boldsymbol{J}_q\cdot\nabla\frac{1}{T} - \frac{\boldsymbol{J}_n}{T}\cdot\nabla\mu \tag{5.1.21}$$

前面说过，化学势的不均匀性是引起物质输运的原因。定义

$$\boldsymbol{X}_n = -\frac{1}{T}\nabla\mu$$

称为**粒子流动力**。局域熵密度产生率 $\Theta$ 可以表示为两种流与力的乘积之和：

$$\Theta = \boldsymbol{J}_q\cdot\boldsymbol{X}_q + \boldsymbol{J}_n\cdot\boldsymbol{X}_n \tag{5.1.22}$$

式 (5.1.22) 的形式具有普遍性。**当多个不可逆过程同时存在时，局域熵密度可以表示为各种不可逆过程的流与力的双线性函数**：

$$\Theta = \sum_k\boldsymbol{J}_k\cdot\boldsymbol{X}_k \tag{5.1.23}$$

如前所述，局域熵产生率 $\Theta$ 满足 $\Theta\geqslant0$，其中等号适用于所有的动力与流量均为零的情形。

---

## §5.2　线性与非线性过程　昂萨格倒易关系

### 5.2.1 四种经验输运定律

许多不可逆过程都是因物体内部某种性质不均匀而引起的输运过程。例如，物体中温度不均匀引起能量的输运，称为**热传导过程**；混合物中各组元浓度不均匀引起物质的输运，称为**扩散过程**；流体流动时速度不均匀引起动量的输运，称为**黏性现象**；导体中的电势差引起电荷的输运，称为**导电过程**，等等。对于一系列输运过程，都已建立了经验规律。**当物体中的不均匀性较小，即偏离平衡不远时，这些经验规律都是线性的。**

| 过程 | 经验定律 | 系数 |
| --- | --- | --- |
| 热传导 | 傅里叶（Fourier）定律 | 导热系数 $\kappa$ |
| 扩散 | 菲克（Fick）定律 | 扩散系数 $D$ |
| 导电 | 欧姆（Ohm）定律 | 电导率 $\sigma$ |
| 黏性 | 牛顿黏性定律 | 黏度 $\eta$ |

$$\boldsymbol{J}_q = -\kappa\nabla T \tag{5.2.1}$$

$$\boldsymbol{J}_n = -D\nabla n \tag{5.2.2}$$

$$\boldsymbol{J}_e = \sigma\boldsymbol{E} = -\sigma\nabla V \tag{5.2.3}$$

$$P_{xy} = \eta\frac{\mathrm{d}v}{\mathrm{d}x} \tag{5.2.4}$$

其中：

- $\boldsymbol{J}_q$ 是单位时间内流过单位截面的热量，称为**热流密度**；
- $\boldsymbol{J}_n$ 是混合物中某组元物质在单位时间内流过单位截面的分子数，称为**粒子流密度**；$n$ 是该组元的浓度；
- $\boldsymbol{J}_e$ 是单位时间内流过单位截面的电荷量，称为**电流密度**；$E$ 是电场强度，$V$ 是电势；
- $P_{xy}$ 是黏性应力，它等于单位时间内沿法线方向通过 $x$ 的单位面积所输运的 $y$ 方向的动量。

### 5.2.2 线性唯象律与昂萨格倒易关系

我们把在单位时间内通过单位截面所输运的物理量（分子数、电荷量、动量、能量等）统称为**热力学流**，以 $\boldsymbol{J}$ 表示；把引起物理量输运的物体中某种性质的梯度（浓度梯度、电势梯度、速度梯度、温度梯度等）统称为**热力学力**，以 $\boldsymbol{X}$ 表示。**在各向同性物体中，上述各种输运过程的经验规律都可表述为"流量与动力成正比"**，即

$$\boldsymbol{J} = L\boldsymbol{X} \tag{5.2.5}$$

不过，在许多情形下，往往有几种力与几种流同时存在，这时将出现不同过程的**交叉现象**。例如，当温度梯度和浓度梯度同时存在时，温度梯度与浓度梯度都会引起热流，也都会引起粒子流，所以式 (5.2.5) 应推广为

$$\boldsymbol{J}_k = \sum_l L_{kl}\boldsymbol{X}_l \tag{5.2.6}$$

上式称为**线性唯象律**，系数 $L_{kl}$ 称为**动理系数**。$L_{kl}$ 等于一个单位的第 $l$ 种动力所引起的第 $k$ 种流量。$L_{kl}$ 一般是**局域强度量**（例如温度、化学势）的函数。

统计物理学可以证明，适当选择流量和动力，使局域熵产生率表达为式 (5.1.23) 的形式，即

$$\Theta = \sum_k\boldsymbol{J}_k\cdot\boldsymbol{X}_k \tag{5.2.7}$$

则动理系数满足关系

$$L_{kl} = L_{lk} \tag{5.2.8}$$

上式称为**昂萨格（Onsager）倒易关系**，简称为**昂萨格关系**。**它描述了第 $l$ 种力对第 $k$ 种流与第 $k$ 种力对第 $l$ 种流所产生的线性效应的对称性。这个关系是微观可逆性在宏观规律上的表现**，它不可能根据热力学理论推导出来，在不可逆过程热力学中我们将直接引用这个公式。

> **脚注**（教材 p.116）：关于昂萨格关系的证明可参阅林宗涵. 热力学与统计物理学[M]. 北京：北京大学出版社，2007: §11.7.

### 5.2.3 熵产生率正定对动理系数的限制

现在讨论条件 $\Theta\geqslant0$ 对动理系数的限制。将式 (5.2.6) 代入式 (5.2.7) 可得

$$\Theta = \sum_{kl}L_{kl}X_kX_l \tag{5.2.9}$$

$\Theta\geqslant0$ 意味着上式是**正定二次型**。为简单起见，我们讨论存在两个耦合的不可逆过程的情形。这时，式 (5.2.9) 为

$$\Theta = L_{11}X_1^2 + \left(L_{12}+L_{21}\right)X_1X_2 + L_{22}X_2^2 \tag{5.2.10}$$

根据线性代数，上式是正定二次型的充要条件为

$$L_{11} > 0, \qquad \begin{vmatrix} L_{11} & L_{12} \\ L_{21} & L_{22} \end{vmatrix} > 0 \tag{5.2.11}$$

在上式得到满足的情形下，仅当 $X_1=X_2=0$，即不存在力，因而不存在流时，$\Theta=0$。**式 (5.2.11) 是热力学第二定律对动理系数的限制。** 将昂萨格关系代入，该式简化为

$$L_{11} > 0, \qquad L_{11}L_{22} > L_{12}^2 \tag{5.2.12}$$

### 5.2.4 线性区与非线性区

前面介绍了线性不可逆过程的动理系数及其性质。如前所述，线性过程相应于动力小、系统偏离平衡不远的情形。一般地说，以 $J_k(\{X_l\})=J_k(X_1,\cdots,X_l,\cdots)$ 表示流量 $J_k$ 作为各种动力的函数，如果将 $J_k(\{X_l\})$ 在 $\{X_l\}$ 的零点展开，得

$$J_k\left(\{X_l\}\right) = J_k(0) + \sum_l\left(\frac{\partial J_k}{\partial X_l}\right)_0X_l + \frac{1}{2}\sum_{l,n}\left(\frac{\partial^2 J_k}{\partial X_l\partial X_n}\right)_0X_lX_n + \cdots \tag{5.2.13}$$

当所有的动力都为零时，流量也将为零，因此上式右方第一项为零。定义

$$L_{kl} = \left(\frac{\partial J_k}{\partial X_l}\right)_0, \qquad L_{kln} = \left(\frac{\partial^2 J_k}{\partial X_l\partial X_n}\right)_0, \qquad \cdots \tag{5.2.14}$$

分别为**动理系数**、**二阶动理系数**……它们一般是局域强度量的函数，则式 (5.2.13) 可改写为

$$J_k = \sum_l L_{kl}X_l + \frac{1}{2}\sum_{l,n}L_{kln}X_lX_n + \cdots \tag{5.2.15}$$

当动力较小时，只需保留展开的一阶项，流与力呈线性关系。如果需要保留二阶以上的项，则流与力呈**非线性关系**。实际问题中，在诸如热传导、电导等输运过程中，流与力的关系一般是线性的；**化学反应中，流与力一般呈非线性关系**。处在**线性区**（流与力呈线性关系）与处在**非线性区**（流与力呈非线性关系）的非平衡系统，其行为有质的不同。我们在后面将分别加以介绍。

---

## *§5.3　温差电现象

> 带星号的选读内容。

### 5.3.1 三种温差电效应

将两种不同的金属相连接，并在两接头处保持不同的温度，电路中将在温度梯度和化学势梯度，因而同时产生热流和粒子流（电流），出现 §5.2 中所说的交叉现象。在这种情形下，实验观察到五种效应：**泽贝克效应、佩尔捷效应、汤姆孙效应、焦耳热效应和热传导过程**。焦耳热效应和热传导过程是我们所熟悉的，现在对前三个效应进行介绍。

**（1）泽贝克效应**（Seebeck，1827 年）。如图 5.3.1 所示，由金属 A、B 接成的热电偶在两个接头处保持不同的温度 $T$ 和 $T+\mathrm{d}T$，实验发现，在电容器中将有电势差：

$$\mathrm{d}V = \varepsilon_{AB}\,\mathrm{d}T \tag{5.3.1}$$

式中 $\varepsilon_{AB}$ 是这两种金属的**温差电动势系数**。我们约定这样选择 $\varepsilon_{AB}$ 的符号：**如果在高温端电动势驱使电流由金属 A 流向金属 B，$\varepsilon_{AB}$ 为正。** $\varepsilon_{AB}$ 取决于两种金属 A、B 的性质，并与温度有关。

**（2）佩尔捷效应**（Peltier，1834 年）。如图 5.3.2 所示，将两种金属 A、B 相连接，并保持其温度恒定不变。实验发现，当有电流通过电路时，在一个接头处有热量放出，在另一个接头处则吸收热量；如果将电流反向，则原来吸热的一端变为放热，原来放热的一端变为吸热。以 $\boldsymbol{J}_e$ 表示由 A 流往 B 的电流密度，单位时间内，单位截面的导线在接头处所吸收的热量为

$$\boldsymbol{J}_{q\pi} = \pi_{AB}\boldsymbol{J}_e \tag{5.3.2}$$

式中 $\boldsymbol{J}_{q\pi}$ 称为**佩尔捷热流密度**，$\pi_{AB}$ 是这两种金属的**佩尔捷系数**，取决于两种金属的性质，并与温度有关。

**（3）汤姆孙效应**（Thomson，1854 年）。当电流通过具有温度梯度的均匀导体时，导体除了放出焦耳热外，还要放出或吸收另外的热量，称为**汤姆孙热**。以 $\boldsymbol{J}_e$ 表示电流密度，在单位时间内，单位体积的导体放出的汤姆孙热 $q_T$ 为

$$q_T = -\tau\boldsymbol{J}_e\cdot\nabla T \tag{5.3.3}$$

$\tau$ 称为**汤姆孙系数**。

**佩尔捷效应和汤姆孙效应都与电流密度成正比，当电流反向时，吸热效应便变为放热效应，所以这两个效应是可逆的。但是，当电路中存在电流和热流时，其中发生的焦耳热效应和热传导过程是不可逆的。因此应该用不可逆过程热力学理论全面地研究整个温差电现象。**

### 5.3.2 电子流与热流的线性唯象律

以 $-e$ 表示电子的电荷量，电流密度 $\boldsymbol{J}_e$ 与电子流 $\boldsymbol{J}_n$ 的关系为 $\boldsymbol{J}_e=-e\boldsymbol{J}_n$。当电路同时存在电子流与热流时，流与力应表示为式 (5.2.6) 的线性律的形式，其中流与力根据式 (5.1.21) 确定。设电子流与热流都平行于 $x$ 轴，略去指标 $x$ 不写，即有

$$\begin{cases} \boldsymbol{J}_n = -L_{11}\dfrac{1}{T}\nabla\mu + L_{12}\nabla\dfrac{1}{T} \\[3mm] \boldsymbol{J}_q = -L_{12}\dfrac{1}{T}\nabla\mu + L_{22}\nabla\dfrac{1}{T} \end{cases} \tag{5.3.4}$$

其中已应用昂萨格关系令 $L_{21}=L_{12}$。

当存在电场时，电子在迁移时除携带通常的化学势外，还携带电能，因此上式中的化学势是**电化学势**，它包括两项：

$$\mu = \mu_c + \mu_e \tag{5.3.5}$$

其中 $\mu_e=-eV$ 是电子的静电势能，$V$ 是电势；$\mu_c$ 是通常的化学势，是温度和电子浓度的函数。

在根据式 (5.3.4) 进行分析之前，要将其中的动理系数 $L_{11}$、$L_{12}$、$L_{22}$ 换为实验测得的经验常量。

### 5.3.3 动理系数与经验常量的关系

**（1）电导率。** $\sigma$ 是在温度均匀的条件下，单位电场强度在导体中产生的电流密度，即有

$$\boldsymbol{J}_e = \sigma\boldsymbol{E}$$

在导体性质和温度都是均匀的情形下，$\nabla\mu_c=0$。因此 $\boldsymbol{E}=-\nabla V=\dfrac{1}{e}\nabla\mu_c=\dfrac{1}{e}\nabla\mu$。所以

$$\sigma = -\frac{eJ_n}{\frac{1}{e}\nabla\mu} \tag{5.3.6}$$

令式 (5.3.4) 第一式的 $\nabla T=0$，与式 (5.3.6) 比较得

$$\sigma = \frac{e^2L_{11}}{T} \tag{5.3.7}$$

**（2）导热系数。** $\kappa$ 是在不存在电流的条件下，单位温度梯度所产生的热流密度，即有

$$\boldsymbol{J}_q = -\kappa\nabla T$$

令式 (5.3.4) 第一式的 $\boldsymbol{J}_n=0$，然后将式 (5.3.4) 的两式联立，消去 $\nabla\mu$，与上式比较得

$$\kappa = \frac{L_{11}L_{22}-L_{12}^2}{T^2L_{11}} \tag{5.3.8}$$

**（3）温差电动势系数。** 现在求动理系数与温差电动势系数的关系。温差电动势是热电偶中不存在电流时所产生的电势差。如图 5.3.3 所示，以 $T_1$、$T_2$（$T_2>T_1$）分别表示导体 A、B 两端接头处的温度。在导体 B 中温度为 $T'$ 处接一电压表，电压表的电阻很高使其无电流通过，但热阻为零。

令式 (5.3.4) 第一式的 $\boldsymbol{J}_n=0$，可得

$$\nabla\mu = -\frac{L_{12}}{TL_{11}}\nabla T \tag{5.3.9}$$

上式对导体 A、B 都成立，因此有

$$\begin{cases} \mu_2-\mu_1 = \displaystyle\int_{T_1}^{T_2}\left(-\frac{L_{12}^A}{TL_{11}^A}\right)\mathrm{d}T \\[3mm] \mu_2-\mu'_R = \displaystyle\int_{T'}^{T_2}\left(-\frac{L_{12}^B}{TL_{11}^B}\right)\mathrm{d}T \\[3mm] \mu'_1-\mu_1 = \displaystyle\int_{T_1}^{T'}\left(-\frac{L_{12}^B}{TL_{11}^B}\right)\mathrm{d}T \end{cases} \tag{5.3.10}$$

消去 $\mu_1$ 和 $\mu_2$，可得

$$\mu'_R-\mu'_L = \int_{T_1}^{T_2}\left(\frac{L_{12}^B}{TL_{11}^B}-\frac{L_{12}^A}{TL_{11}^A}\right)\mathrm{d}T \tag{5.3.11}$$

L、R 两端的 $\mu_c$ 相等，故有

$$V = \frac{1}{-e}\left(\mu'_R-\mu'_L\right) = \int_{T_1}^{T_2}\left(\frac{-L_{12}^B}{eTL_{11}^B}-\frac{-L_{12}^A}{eTL_{11}^A}\right)\mathrm{d}T$$

由此可得温差电动势系数为

$$\varepsilon_{AB} = \frac{-L_{12}^B}{eTL_{11}^B} - \frac{-L_{12}^A}{eTL_{11}^A} \tag{5.3.12}$$

定义导体 A、B 的**绝对温差电动势系数** $\varepsilon$ 为

$$\varepsilon_A = \frac{-L_{12}^A}{eTL_{11}^A}, \qquad \varepsilon_B = \frac{-L_{12}^B}{eTL_{11}^B} \tag{5.3.13}$$

即有

$$\varepsilon_{AB} = \varepsilon_B - \varepsilon_A \tag{5.3.14}$$

由此可知，当 $\varepsilon_B>\varepsilon_A$ 时，在高温端电流密度 $J_e$ 由金属 A 流向金属 B。

### 5.3.4 开尔文第一、第二关系

上面求出了电导率 $\sigma$、导热系数 $\kappa$ 和绝对温差电动势系数 $\varepsilon$ 与动理系数的关系。式 (5.3.4) 含有三个动理系数，可以将式 (5.3.4) 中的动理系数换为经验常量而将式 (5.3.4) 改写为

$$\begin{cases} \boldsymbol{J}_n = -\left(\dfrac{T\sigma}{e^2}\right)\dfrac{1}{T}\nabla\mu - \dfrac{T^2\sigma\varepsilon}{e}\nabla\dfrac{1}{T} \\[3mm] \boldsymbol{J}_q = \left(\dfrac{T^2\sigma\varepsilon}{e}\right)\dfrac{1}{T}\nabla\mu + \left(T^3\sigma\varepsilon^2+\kappa T^2\right)\nabla\dfrac{1}{T} \end{cases} \tag{5.3.15}$$

将式 (5.3.15) 的两式联立，消去 $\nabla\mu$，可得

$$\boldsymbol{J}_q = T\varepsilon\left(-e\boldsymbol{J}_n\right) - \kappa\nabla T$$

或

$$\boldsymbol{J}_q = T\varepsilon\boldsymbol{J}_e - \kappa\nabla T \tag{5.3.16}$$

注意 $\boldsymbol{J}_s=\dfrac{\boldsymbol{J}_q}{T}$，故由上式得

$$\boldsymbol{J}_s = \varepsilon\boldsymbol{J}_e - \frac{1}{T}\kappa\nabla T \tag{5.3.17}$$

式 (5.3.17) 右方第二项是热传导引起的熵流，第一项是电流携带的熵流。由此可知，**绝对温差电动势系数 $\varepsilon$ 是单位电流密度所携带的熵流密度。**

同理，式 (5.3.16) 右方第二项是热传导的热流密度，第一项是伴随着电流密度 $\boldsymbol{J}_e$ 的热流密度，即**佩尔捷热流密度**：

$$\boldsymbol{J}_{q\pi} = T\varepsilon\boldsymbol{J}_e$$

在接头处两金属的 $\varepsilon$ 不相等，使 $\boldsymbol{J}_{q\pi}$ 在接头处有突变而出现吸热或放热：

$$\boldsymbol{J}_{q\pi} = T\left(\varepsilon_B-\varepsilon_A\right)\boldsymbol{J}_e \tag{5.3.18}$$

与式 (5.3.2) 比较可知，**佩尔捷系数 $\pi_{AB}$ 与绝对温差电动势系数之间存在关系**：

$$\pi_{AB} = T\left(\varepsilon_B-\varepsilon_A\right) \tag{5.3.19}$$

上式称为**开尔文第二关系**。

现在讨论汤姆孙效应。当定常电流流过具有温度梯度的均匀导体时，导体内能密度的增加率为

$$\frac{\partial u}{\partial t} = -\nabla\cdot\boldsymbol{J}_u = -\nabla\cdot\left(\boldsymbol{J}_q+\mu\boldsymbol{J}_n\right) = -\nabla\cdot\boldsymbol{J}_q - \boldsymbol{J}_n\cdot\nabla\mu \tag{5.3.20}$$

其中用到了定常电流条件 $\nabla\cdot\boldsymbol{J}_n=0$。将式 (5.3.16) 和式 (5.3.15) 的第一式代入式 (5.3.20)，可得

$$\begin{aligned}\frac{\partial u}{\partial t} &= -\nabla\cdot\left(T\varepsilon\boldsymbol{J}_e-\kappa\nabla T\right) + \varepsilon\boldsymbol{J}_e\cdot\nabla T + \frac{J_e^2}{\sigma} \\ &= -T\boldsymbol{J}_e\cdot\nabla\varepsilon - \nabla\cdot\left(-\kappa\nabla T\right) + \frac{1}{\sigma}J_e^2 \end{aligned} \tag{5.3.21}$$

式 (5.3.21) 右方第三项是焦耳热，第二项代表热传导过程所聚集的热，第一项就是汤姆孙热。注意绝对温差电动势系数 $\varepsilon$ 是温度的函数：

$$\nabla\varepsilon = \frac{\mathrm{d}\varepsilon}{\mathrm{d}T}\nabla T$$

因此，式 (5.3.21) 右方第一项可表示为

$$-T\frac{\mathrm{d}\varepsilon}{\mathrm{d}T}\boldsymbol{J}_e\cdot\nabla T \tag{5.3.22}$$

与式 (5.3.3) 比较，得

$$\tau = T\frac{\mathrm{d}\varepsilon}{\mathrm{d}T} \tag{5.3.23}$$

**上式给出了汤姆孙系数 $\tau$ 与绝对温差电动势系数 $\varepsilon$ 之间的关系。**

将式 (5.3.19) 对 $T$ 求导数，得

$$\frac{\mathrm{d}\pi_{AB}}{\mathrm{d}T} = \left(\varepsilon_B-\varepsilon_A\right) + T\frac{\mathrm{d}}{\mathrm{d}T}\left(\varepsilon_B-\varepsilon_A\right)$$

再将式 (5.3.23) 代入上式，即得

$$\frac{\mathrm{d}\pi_{AB}}{\mathrm{d}T} + \left(\tau_A-\tau_B\right) = \varepsilon_B-\varepsilon_A \tag{5.3.24}$$

上式称为**开尔文第一关系**。

### 5.3.5 线性不可逆过程热力学的解题程序

通过以上讨论可以看出，用线性不可逆过程热力学处理问题的一般程序为：

1. 写出线性唯象律，利用昂萨格关系减少其中出现的动理系数的数目；
2. 分析一些物理效应（其数目等于线性唯象律所含动理系数的数目），求出经验常量与动理系数的关系，从而将线性唯象律用经验常量表出；
3. 进一步分析其他物理效应，即可找出经验常量间的关系。

---

## *§5.4　最小熵产生定理

> 带星号的选读内容。

### 5.4.1 非平衡定态

上节以温差电效应为例，介绍了用线性不可逆过程热力学处理问题的方法。在流与力呈线性关系的情形下，如果外界施加某种恒定的动力，系统将处在某种**定常（不随时间变化）的非平衡态**，称为**非平衡定态**。

我们先讨论一个最简单的例子。考虑处在两块面积很大的金属平面板之间的液体薄层。外界以定常的速率均匀地向下板供给热量，从上板吸取热量，使下板保持恒温 $T_2$，上板保持恒温 $T_1$（$T_2>T_1$）。在 $\Delta T=T_2-T_1$ 不大的情形下，液体内部将建立起定常的温度分布，其中存在定常的温度梯度和热流。这就是一种非平衡定态。

下面我们证明，**这种非平衡定态是一种熵产生率最小的状态**。根据式 (5.1.13)，在单纯的热传导过程中，局域熵密度产生率为

$$\Theta = \boldsymbol{J}_q\cdot\nabla\frac{1}{T} \tag{5.4.1}$$

在热流密度与热流动力呈线性关系的情形下，有

$$\boldsymbol{J}_q = L_{qq}\nabla\frac{1}{T} \tag{5.4.2}$$

系数 $L_{qq}=\left(\dfrac{\partial J_q}{\partial X_q}\right)_0$ 是动理系数。整个系统的熵产生率为

$$P = \int\Theta\,\mathrm{d}\tau = \int\boldsymbol{J}_q\cdot\nabla\frac{1}{T}\mathrm{d}\tau = \int L_{qq}\left(\nabla\frac{1}{T}\right)^2\mathrm{d}\tau \tag{5.4.3}$$

### 5.4.2 最小熵产生定理的证明

将上式对时间求导数，在 $L_{qq}$ 不随时间变化的情形下，有

$$\begin{aligned}\frac{\mathrm{d}P}{\mathrm{d}t} &= 2\int L_{qq}\nabla\left(\frac{1}{T}\right)\cdot\nabla\left(\frac{\partial}{\partial t}\frac{1}{T}\right)\mathrm{d}\tau \\ &= 2\int\boldsymbol{J}_q\cdot\nabla\left(\frac{\partial}{\partial t}\frac{1}{T}\right)\mathrm{d}\tau \\ &= 2\int\nabla\cdot\left(\boldsymbol{J}_q\frac{\partial}{\partial t}\frac{1}{T}\right)\mathrm{d}\tau - 2\int\left(\frac{\partial}{\partial t}\frac{1}{T}\right)\nabla\cdot\boldsymbol{J}_q\,\mathrm{d}\tau\end{aligned}$$

上式右方第一项可换为面积积分

$$2\oint\frac{\partial}{\partial t}\left(\frac{1}{T}\right)\boldsymbol{J}_q\cdot\mathrm{d}\boldsymbol{\sigma}$$

在边界温度不随时间变化的情形下面积分为零，故有

$$\frac{\mathrm{d}P}{\mathrm{d}t} = -2\int\left(\frac{\partial}{\partial t}\frac{1}{T}\right)\nabla\cdot\boldsymbol{J}_q\,\mathrm{d}\tau \tag{5.4.4}$$

在体积变化可以忽略时，单位体积的内能可表示为 $\mathrm{d}u=c_V\mathrm{d}T$，其中 $c_V$ 是单位体积的定容热容，因此

$$\frac{\mathrm{d}u}{\mathrm{d}t} = c_V\frac{\partial T}{\partial t}$$

与式 (5.1.9) 比较得

$$c_V\frac{\partial T}{\partial t} = -\nabla\cdot\boldsymbol{J}_q \tag{5.4.5}$$

代入式 (5.4.4) 得

$$\frac{\mathrm{d}P}{\mathrm{d}t} = -2\int\frac{c_V}{T^2}\left(\frac{\partial T}{\partial t}\right)^2\mathrm{d}\tau$$

由于被积函数非负，故有

$$\frac{\mathrm{d}P}{\mathrm{d}t} \leqslant 0 \qquad \text{或} \qquad \frac{\mathrm{d}\Theta}{\mathrm{d}t} \leqslant 0 \tag{5.4.6}$$

上式表明，**如果系统的温度分布随时间变化，即 $\dfrac{\partial T}{\partial t}\neq0$，其中发生的（线性）热传导过程将使系统的熵产生率随时间而减小，直到熵产生率达到极小值，系统处在具有定常温度分布的非平衡定态为止。这就是最小熵产生定理。**

不仅如此，根据最小熵产生定理，系统处在非平衡定态时，如果由于某种外界扰动或内部涨落使系统离开了这一状态，只要扰动或涨落不大，未破坏流与力的线性关系，**系统会回到熵产生率最小的非平衡定态**。这就是说，**在流与力呈线性关系的范围内，这种具有最小熵产生率的非平衡定态是稳定的。**

### 5.4.3 两个耦合过程的情形

前面就单纯的热传导过程介绍了最小熵产生定理，现在讨论存在两个耦合的不可逆过程的情形。根据式 (5.2.9)，局域熵产生率为

$$\Theta = L_{11}X_1^2 + 2L_{12}X_1X_2 + L_{22}X_2^2 \tag{5.4.7}$$

如果对力未加约束，最小熵产生要求（假设动理系数是常量）：

$$\begin{cases} \left(\dfrac{\partial\Theta}{\partial X_1}\right)_{X_2} = 2L_{11}X_1+2L_{12}X_2 = 2J_1 = 0 \\[3mm] \left(\dfrac{\partial\Theta}{\partial X_2}\right)_{X_1} = 2L_{22}X_2+2L_{21}X_1 = 2J_2 = 0 \end{cases} \tag{5.4.8}$$

在上式得到满足即 $J_1=J_2=0$ 的情形下，如果 $L_{11}L_{22}-L_{12}^2\neq0$，由线性代数知，式 (5.4.8) 只有平庸解 $X_1=X_2=0$。这就是说，**在对动力未加约束的情形下，熵产生最小的状态是动力和流量均为零的平衡态，其熵产生率为零。**

如果对动力加以约束，例如令 $X_1$ 为常量，最小熵产生条件要求：

$$\left(\frac{\partial\Theta}{\partial X_2}\right)_{X_1} = 2L_{22}X_2+2L_{21}X_1 = 2J_2 = 0 \tag{5.4.9}$$

由上式得 $X_2=-\dfrac{L_{21}}{L_{22}}X_1$，由此可得

$$J_1 = L_{11}X_1+L_{12}X_2 = \left(L_{11}-\frac{L_{12}L_{21}}{L_{22}}\right)X_1 \tag{5.4.10}$$

及

$$\Theta = L_{11}X_1^2+2L_{12}X_1X_2+L_{22}X_2^2 = \left(L_{11}-\frac{L_{12}L_{21}}{L_{22}}\right)X_1^2 \tag{5.4.11}$$

在这种情形下，系统处在具有定常的 $X_1$、$X_2$ 和定常的 $J_1$、$J_2$（$J_2=0$）的非平衡定态。**在 $X_1$ 是常量的约束条件下，这状态的熵产生率最小。**

> **脚注**（教材 p.123）：根据式 (5.2.12)，热力学第二定律要求 $L_{11}>0$、$L_{11}L_{22}-L_{12}^2>0$，满足 $L_{11}L_{22}-L_{12}^2\neq0$ 的要求。

前面我们通过两个例子介绍了最小熵产生定理。在习题 5.4 和 §5.6 中，我们会看到其他例子。更普遍的讨论请参阅其他书籍。

> **脚注**（教材 p.124）：De Groot S R, Mazur P. *Non-equilibrium Thermodynamics*[M]. Amsterdam: North Holland Publishing Company, 1962; Chap 5.

---

## *§5.5　化学反应与扩散过程

> 带星号的选读内容。

### 5.5.1 反应速率与连续性方程

考虑恒温恒压的理想流体（理想气体或理想溶液），其中含有 $k$ 个组元且同时存在 $r$ 个化学反应。由于压强和温度保持恒定，系统内部不存在因压强或温度不均匀引起的物质流动或热传导，但存在由于组元浓度在空间分布不均匀引起的扩散过程。我们假设整个系统虽然处在非平衡状态，但将系统分为若干小部分，各小部分处在局域平衡。本节将在局域平衡假设下导出同时存在化学反应和扩散过程时的熵流密度和熵产生率。**我们将会看到，化学反应中流与力之间的关系往往是非线性的。**

系统中某体积元内，化学反应

$$\mathrm{A}+\mathrm{X}_i \xrightarrow{k_1} \mathrm{Y}+\mathrm{B} \tag{5.5.1}$$

的反应速率 $\omega_1$ 显然与体积元内分子 A 与分子 $\mathrm{X}_i$ 发生碰撞的频率成正比，因而与其中反应物 A 和 $\mathrm{X}_i$ 的分子数密度 $n_A$ 和 $n_i$ 成正比，即

$$\omega_1 = k_1n_An_i \tag{5.5.2}$$

式中 $k_1$ 是比例系数。一般来说，在非平衡系统中，$n_A$、$n_i$ 和 $\omega_1$ 可以是时间和坐标的函数。

同理，体积元内，化学反应

$$\mathrm{A}+2\mathrm{X}_i \xrightarrow{k_2} \mathrm{Z}+\mathrm{B}+\mathrm{C} \tag{5.5.3}$$

的反应速率 $\omega_2$ 可表示为

$$\omega_2 = k_2n_An_in_i = k_2n_An_i^2 \tag{5.5.4}$$

式中 $k_2$ 是比例系数。

两个反应同时发生时，体积元内 $\mathrm{X}_i$ 的分子数密度 $n_i$ 的变化率为

$$\left(\frac{\partial n_i}{\partial t}\right)_{\mathrm{ch}} = -k_1n_An_i - 2k_2n_An_i^2 = -\omega_1 - 2\omega_2 \tag{5.5.5}$$

另一方面，由扩散引起分子数密度 $n_i$ 的变化率可表示为

$$\left(\frac{\partial n_i}{\partial t}\right)_{\mathrm{di}} = -\nabla\cdot\boldsymbol{J}_i \tag{5.5.6}$$

$\boldsymbol{J}_i$ 是 $\mathrm{X}_i$ 的粒子流密度。将上两式相加，得体积元中分子数密度 $n_i$ 的变化率为

$$\frac{\partial n_i}{\partial t} = -\nabla\cdot\boldsymbol{J}_i - \omega_1 - 2\omega_2 \tag{5.5.7}$$

一般地，当体积元中存在 $r$ 个化学反应时，组元 $i$ 的分子数密度 $n_i$ 的变化率可表示为

$$\frac{\partial n_i}{\partial t} = -\nabla\cdot\boldsymbol{J}_i + \sum_{\rho=1}^{r}\nu_{i\rho}\omega_\rho \tag{5.5.8}$$

式中 $\omega_\rho$ 是第 $\rho$ 个化学反应的反应速率，$\nu_{i\rho}$ 是在第 $\rho$ 个反应方程中组元 $i$ 的系数。**当组元 $i$ 在反应方程 $\rho$ 中是生成物时，$\nu_{i\rho}$ 为正；是反应物时，$\nu_{i\rho}$ 为负。**

### 5.5.2 熵流密度与熵产生率

对于反应–扩散过程，由式 (5.1.4) 知，局域熵密度的变化率为

$$\begin{aligned}\frac{\partial s}{\partial t} &= -\sum_i\frac{\mu_i}{T}\frac{\partial n_i}{\partial t} \\ &= \sum_i\frac{\mu_i}{T}\nabla\cdot\boldsymbol{J}_i - \sum_i\sum_\rho\frac{\mu_i}{T}\nu_{i\rho}\omega_\rho \\ &= \nabla\cdot\left(\sum_i\frac{\mu_i}{T}\boldsymbol{J}_i\right) - \sum_i\boldsymbol{J}_i\cdot\nabla\frac{\mu_i}{T} + \sum_\rho\left(-\sum_i\frac{\mu_i}{T}\nu_{i\rho}\right)\omega_\rho \end{aligned} \tag{5.5.9}$$

引入反应 $\rho$ 的**局部化学亲和势**

$$a_\rho = -\sum_{i=1}^{k}\mu_i\nu_{i\rho} \tag{5.5.10}$$

可将式 (5.5.9) 表示为

$$\frac{\partial s}{\partial t} = \nabla\cdot\left(\sum_i\frac{\mu_i}{T}\boldsymbol{J}_i\right) - \sum_i\boldsymbol{J}_i\cdot\nabla\frac{\mu_i}{T} + \sum_\rho\frac{a_\rho}{T}\omega_\rho \tag{5.5.11}$$

与式 (5.1.6) 比较知，熵流密度为

$$\boldsymbol{J}_s = -\sum_i\frac{\mu_i}{T}\boldsymbol{J}_i \tag{5.5.12}$$

上式给出粒子流所携带的熵流。局域熵产生率为

$$\Theta = -\sum_i\boldsymbol{J}_i\cdot\nabla\frac{\mu_i}{T} + \sum_\rho\frac{a_\rho}{T}\omega_\rho \tag{5.5.13}$$

两项分别表示**扩散过程**和**化学反应过程**的局域熵产生率。将上式与式 (5.1.23) 比较知，反应扩散过程的流量 $J$ 与动力 $X$ 为

$$\begin{cases} \boldsymbol{J}_i^{\mathrm{di}} = \boldsymbol{J}_i, & \boldsymbol{X}_i^{\mathrm{di}} = -\nabla\dfrac{\mu_i}{T} \\[3mm] J_\rho^{\mathrm{ch}} = \omega_\rho, & X_\rho^{\mathrm{ch}} = \dfrac{a_\rho}{T} \end{cases} \tag{5.5.14}$$

### 5.5.3 局部化学亲和势

现在讨论式 (5.5.10) 引入的局部化学亲和势的意义。假设只存在反应 $\rho$，当 $a_\rho=0$ 时，式 (5.5.10) 与式 (4.5.5) 相当，意味着体积元内反应 $\rho$ 达到**局部平衡**。局域平衡下，理想流体局域化学势的函数形式为

$$\mu_i = g_i(T,p) + kT\ln x_i \tag{5.5.15}$$

> **脚注**（教材 p.125）：参阅式 (4.6.6) 和式 (4.6.17)，其中 $k=\dfrac{R}{N_A}$，$N_A$ 是阿伏伽德罗常量。

因此，局部化学亲和势 $a_\rho$ 可表示为

$$a_\rho = -\sum_i\nu_{i\rho}g_i(T,p) - kT\ln\prod_i x_i^{\nu_{i\rho}} \tag{5.5.16}$$

定义反应 $\rho$ 的**局部平衡常量** $K_\rho(T,p)$ 为

$$\ln K_\rho(T,p) = -\frac{1}{kT}\sum_i\nu_{i\rho}g_i(T,p) \tag{5.5.17}$$

则 $a_\rho$ 可表示为

$$a_\rho = kT\ln\frac{K_\rho(T,p)}{\prod_i x_i^{\nu_{i\rho}}} \tag{5.5.18}$$

显然，当

$$K_\rho(T,p) = \prod_i x_i^{\nu_{i\rho}} \tag{5.5.19}$$

时，$a_\rho=0$，化学反应 $\rho$ 达到局部平衡。式 (5.5.19) 与理想气体或理想溶液中化学反应的质量作用律具有相同的形式。

### 5.5.4 亲和势与反应速率的关系

现在讨论化学亲和势与反应速率的关系。我们考虑下述反应：

$$\mathrm{A} \underset{k_-}{\overset{k_+}{\rightleftharpoons}} \mathrm{B} \tag{5.5.20}$$

与式 (5.5.2) 或式 (5.5.4) 相似，反应正向进行的速率 $\omega_+$ 为

$$\omega_+ = k_+x_A \tag{5.5.21}$$

逆向进行的速率 $\omega_-$ 为

$$\omega_- = k_-x_B \tag{5.5.22}$$

式中 $k_+$、$k_-$ 是比例系数。净反应速率 $\omega_\rho$ 为

$$\omega_\rho = \omega_+-\omega_- = k_+x_A-k_-x_B = k_+x_A\left(1-\frac{k_-x_B}{k_+x_A}\right) \tag{5.5.23}$$

利用反应达到局部平衡时 $\omega_\rho=0$ 和 $K_\rho(T,p)=\dfrac{x_B}{x_A}$ 的条件，可由上式得出

$$\frac{k_-}{k_+} = \frac{1}{K_\rho(T,p)}$$

代回式 (5.5.23)，得

$$\omega_\rho = k_+x_A\left[1-\frac{1}{K_\rho(T,p)}\frac{x_B}{x_A}\right] = k_+x_A\left(1-e^{-\frac{a_\rho}{kT}}\right) \tag{5.5.24}$$

最后一步用到了式 (5.5.18)。通过这个例子可以引出化学亲和势 $a_\rho$ 与反应速率 $\omega_\rho$ 的一般性关系。如前所述，$a_\rho=0$ 或 $\omega_\rho=0$ 相当于局域平衡的情形。如果

$$a_\rho \ll kT \tag{5.5.25}$$

式 (5.5.24) 可近似为

$$\omega_\rho = k_+x_A\frac{a_\rho}{kT} = \omega_+\frac{a_\rho}{kT} \tag{5.5.26}$$

这时反应速率 $\omega_\rho$ 与化学亲和势 $a_\rho$ 呈线性关系，亦即化学反应的流 $J_\rho^{\mathrm{ch}}$ 与力 $X_\rho^{\mathrm{ch}}$ 呈线性关系。

在相反的情形下，如果 $a_\rho\gg kT$，由式 (5.5.24) 得

$$\omega_\rho = k_+x_A = \omega_+$$

反应将单向进行。由此可知，**化学反应中流与力的线性关系只在 $a_\rho$ 很小即非常接近化学平衡时成立。实际的化学反应过程大都不满足式 (5.5.15) 而属于非线性情形。**

后面我们将就反应扩散过程介绍非线性不可逆过程热力学。该理论是普里高津（Prigogine）和他的合作者在二十世纪六七十年代发展起来的。非线性不可逆过程热力学也被用于流体力学、激光和生物等领域的研究。

> **脚注**（教材 p.127）：① Nicolis G, Prigogine I. *Self-organization in Nonequilibrium Systems*[M]. New York: Wiley-Interscience, 1977.
> ② Glansdorff P, Prigogine I. *Thermodynamic Theory of Structure, Stability and Fluctuations*[M]. New York: Wiley-Interscience, 1971.
> ③ De Hemptinne X. Thermodynamics of Laser Systems[J]. Phys. Rept., 1985:1.

---

## *§5.6　非平衡系统在非线性区的发展判据

> 带星号的选读内容。

### 5.6.1 熵产生率变化率的分解

最小熵产生定理指出，对于处在线性区的非平衡系统，如果保持某种恒定的动力，系统将稳定地处在一个非平衡定态。现在讨论处在非线性区的非平衡系统的行为。

根据式 (5.1.23)，系统熵产生率随时间的变化率可表示为

$$\begin{aligned}\frac{\mathrm{d}P}{\mathrm{d}t} &= \int\frac{\mathrm{d}\Theta}{\mathrm{d}t}\mathrm{d}\tau = \int\sum_k J_k\frac{\mathrm{d}X_k}{\mathrm{d}t}\mathrm{d}\tau + \int\sum_k\frac{\mathrm{d}J_k}{\mathrm{d}t}X_k\,\mathrm{d}\tau \\ &= \frac{\mathrm{d}_XP}{\mathrm{d}t} + \frac{\mathrm{d}_JP}{\mathrm{d}t} \end{aligned} \tag{5.6.1}$$

两项分别表示**力与流随时间变化引起的系统熵产生率的变化率**。

### 5.6.2 力随时间变化引起的熵产生率变化率

我们先就恒温恒压下反应扩散过程的情形计算由于力随时间变化引起的熵产生率的变化率 $\dfrac{\mathrm{d}_XP}{\mathrm{d}t}$。根据式 (5.5.13)，得

$$\frac{\mathrm{d}_XP}{\mathrm{d}t} = \frac{1}{T}\int\mathrm{d}\tau\left(-\sum_i\boldsymbol{J}_i\cdot\nabla\frac{\partial\mu_i}{\partial t} + \sum_\rho\omega_\rho\frac{\partial a_\rho}{\partial t}\right) \tag{5.6.2}$$

对右方第一项进行分部积分并将第二项中的 $a_\rho$ 用式 (5.5.10) 代入，得

$$\frac{\mathrm{d}_XP}{\mathrm{d}t} = -\frac{1}{T}\int\nabla\cdot\left(\sum_i\boldsymbol{J}_i\frac{\partial\mu_i}{\partial t}\right)\mathrm{d}\tau + \frac{1}{T}\int\left(\sum_i\frac{\partial\mu_i}{\partial t}\nabla\cdot\boldsymbol{J}_i - \sum_{i,\rho}\omega_\rho\nu_{i\rho}\frac{\partial\mu_i}{\partial t}\right)\mathrm{d}\tau$$

> **⚠ 教材笔误提示**：教材 p.127 此式括号内的第二项印作 $-\dfrac{1}{T}\sum_{i,\rho}\omega_\rho\nu_{i\rho}\dfrac{\partial\mu_i}{\partial t}$，多出一个 $1/T$。由式 (5.5.10) 求导得 $\dfrac{\partial a_\rho}{\partial t}=-\sum_i\nu_{i\rho}\dfrac{\partial\mu_i}{\partial t}$，代入式 (5.6.2) 后外层 $1/T$ 已经用掉，故括号内不应再有 $1/T$；只有取上式（无内层 $1/T$）代回连续性方程 (5.5.8)，才能得到教材中正确的式 (5.6.3)。

利用高斯定理将右方第一项化为系统边界上的面积分，如果边界条件不随时间变化，此项为零。在恒温恒压条件下，$\dfrac{\partial\mu_i}{\partial t}=\sum_j\dfrac{\partial\mu_i}{\partial n_j}\dfrac{\partial n_j}{\partial t}$，故上式可表示为

$$\frac{\mathrm{d}_XP}{\mathrm{d}t} = -\frac{1}{T}\int\left[\sum_{i,j}\frac{\partial\mu_i}{\partial n_j}\frac{\partial n_j}{\partial t}\left(-\nabla\cdot\boldsymbol{J}_i + \sum_\rho\nu_{i\rho}\omega_\rho\right)\right]\mathrm{d}\tau$$

利用反应扩散过程的连续性方程式 (5.5.8)，可将上式表示为

$$\frac{\mathrm{d}_XP}{\mathrm{d}t} = -\frac{1}{T}\int\sum_{i,j}\frac{\partial\mu_i}{\partial n_j}\frac{\partial n_j}{\partial t}\frac{\partial n_i}{\partial t}\mathrm{d}\tau \tag{5.6.3}$$

> **脚注**（教材 p.127）：式中的 $\dfrac{\partial\mu_i}{\partial n_j}$ 不是完全独立的，要满足 $\sum_j n_j\dfrac{\partial\mu_i}{\partial n_j}=0$。

现在讨论上式被积函数的符号。由于系统各小部分处在局域平衡，在恒温恒压条件下，局域吉布斯函数密度 $g$ 应具有极小值，即它的一级微分 $\delta g=0$，二级微分 $\delta^2g\geqslant0$。因此根据式 (4.1.11)，在恒温恒压条件下，有

$$\delta g = \sum_i\mu_i\delta n_i = 0$$

及

$$\delta^2g = \sum_{i,j}\frac{\partial\mu_i}{\partial n_j}\delta n_i\delta n_j \geqslant 0$$

由此可知，式 (5.6.3) 的被积函数不为负，故

$$\frac{\mathrm{d}_XP}{\mathrm{d}t} \leqslant 0 \tag{5.6.4}$$

上式意味着，**力随时间变化将导致系统的熵产生率减小。这一结论对于处在线性区和非线性区的情形都适用。**

### 5.6.3 线性区与非线性区的区别

对于处在线性区的情形，流与力存在式 (5.2.6) 的关系：

$$\boldsymbol{J}_k = \sum_l L_{kl}\boldsymbol{X}_l$$

在动理系数为常量的情形下，有

$$\begin{aligned}\frac{\mathrm{d}_JP}{\mathrm{d}t} &= \int\sum_{k,l}L_{kl}\frac{\mathrm{d}X_l}{\mathrm{d}t}X_k\,\mathrm{d}\tau \\ &= \int\sum_{k,l}L_{lk}X_k\frac{\mathrm{d}X_l}{\mathrm{d}t}\mathrm{d}\tau = \int\sum_l J_l\frac{\mathrm{d}X_l}{\mathrm{d}t}\mathrm{d}\tau \\ &= \frac{\mathrm{d}_XP}{\mathrm{d}t} \end{aligned} \tag{5.6.5}$$

其中第二步用到了昂萨格关系。因此

$$\frac{\mathrm{d}P}{\mathrm{d}t} = 2\frac{\mathrm{d}_XP}{\mathrm{d}t} \leqslant 0 \tag{5.6.6}$$

**式 (5.6.6) 就是 §5.4 中讲述的最小熵产生定理。** 这意味着，对于处在线性区的反应扩散过程，系统的定态对于外界的扰动或内部的涨落是稳定的。

对于处在非线性区的非平衡系统，由于 $J$ 与 $X$ 的非线性关系，$\dfrac{\mathrm{d}_JP}{\mathrm{d}t}$ 的符号是不定的，因而 $\dfrac{\mathrm{d}P}{\mathrm{d}t}$ 的符号也不定。**这意味着，存在这样的可能性，当发生扰动或涨落时，系统原来所处的定态会变得不稳定而演化到另一个新的定态，即发生非平衡相变。** 在 §5.7 中我们将以三分子模型为例，介绍非平衡相变和耗散结构的概念。

---

## *§5.7　三分子模型与耗散结构的概念

> 带星号的选读内容。

### 5.7.1 三分子模型（Brusselator）

实验观察到一些化学反应会出现空间、时间或空间–时间的有序图案。**别洛乌索夫–扎博京斯基（Belousov–Zhabotinskii）反应**（简称为 BZ 反应）就是一个著名的例子。将溴酸钾（$\mathrm{KBrO_3}$）、丙二酸［$\mathrm{CH_2(COOH)_2}$］、硫酸铈［$\mathrm{Ce_2(SO_4)_3}$］、硫酸（$\mathrm{H_2SO_4}$）及几滴亚铁灵试剂加以混合并搅拌，溶液的颜色会在红色和蓝色之间振荡。振荡的周期是分钟的量级，整个现象的寿命是小时量级。颜色的振荡反映了 $\mathrm{Br^-}$ 浓度和 $\mathrm{Ce^{4+}/Ce^{3+}}$ 相对浓度的周期性变化。$\mathrm{Ce^{3+}}$ 过多呈红色，$\mathrm{Ce^{4+}}$ 过多呈蓝色。如果反应在未被搅动的浅盘内进行，实验观察到 $\mathrm{Br^-}$ 浓度和 $\mathrm{Ce^{4+}/Ce^{3+}}$ 相对浓度的行波。

BZ 反应的行为可以用俄勒冈模型（Oregonator）描述。不过较为简单的三分子模型也可显示类似的定性行为。本节对三分子模型加以介绍。**三分子模型是由普里高津和勒菲弗（Lefever）提出来的，常被称为布鲁塞尔模型（Brusselator）。** 它包含下述四步化学反应：

$$\begin{cases} \mathrm{A} \xrightarrow{k_1} \mathrm{X} \\[1mm] \mathrm{B}+\mathrm{X} \xrightarrow{k_2} \mathrm{Y}+\mathrm{D} \\[1mm] 2\mathrm{X}+\mathrm{Y} \xrightarrow{k_3} 3\mathrm{X} \\[1mm] \mathrm{X} \xrightarrow{k_4} \mathrm{E} \end{cases} \tag{5.7.1}$$

在反应中不断供给反应物 A 和 B，使其浓度保持恒定，并不断将生成物 D 和 E 排除。于是反应单向进行，系统处在**远离平衡**的状态。在 A、B、D、E 的浓度保持恒定的情形下，只有 X 和 Y 的浓度 $n_x$ 和 $n_y$ 随时间变化，其变化率为

$$\begin{cases} \dfrac{\mathrm{d}n_x}{\mathrm{d}t'} = k_1n_A - \left(k_2n_B+k_4\right)n_x + k_3n_x^2n_y + D'_1\nabla^2n_x \\[3mm] \dfrac{\mathrm{d}n_y}{\mathrm{d}t'} = k_2n_Bn_x - k_3n_x^2n_y + D'_2\nabla^2n_y \end{cases} \tag{5.7.2}$$

> **脚注**（教材 p.129）：式 (5.7.2) 两式右方末项是扩散项。由于 $\dfrac{\partial n}{\partial t}=-\nabla\cdot\boldsymbol{J}$，而 $\boldsymbol{J}=-D\nabla n$，因此由扩散引起的浓度变化率为 $\left(\dfrac{\partial n}{\partial t}\right)_{\mathrm{di}}=D\nabla^2n$。

作变量代换：

$$\begin{cases} t = k_4t', \qquad X = \left(\dfrac{k_3}{k_4}\right)^{\frac{1}{2}}n_x, \qquad Y = \left(\dfrac{k_3}{k_4}\right)^{\frac{1}{2}}n_y \\[3mm] A = \left(\dfrac{k_1^2k_3}{k_4^3}\right)^{\frac{1}{2}}n_A, \qquad B = \left(\dfrac{k_2}{k_4}\right)n_B, \qquad D_i = \dfrac{D'_i}{k_4} \end{cases} \tag{5.7.3}$$

可将式 (5.7.2) 表示为

$$\begin{cases} \dfrac{\mathrm{d}X}{\mathrm{d}t} = A-\left(B+1\right)X+X^2Y+D_1\nabla^2X \\[3mm] \dfrac{\mathrm{d}Y}{\mathrm{d}t} = BX-X^2Y+D_2\nabla^2Y \end{cases} \tag{5.7.4}$$

### 5.7.2 定态解与线性化

容易验明，方程 (5.7.4) 有下述均匀的定常解：

$$X_0 = A, \qquad Y_0 = \frac{B}{A} \tag{5.7.5}$$

下面分析在什么情形下式 (5.7.5) 的定常解会失稳。我们考虑两类边界条件：

**（一）$X$ 和 $Y$ 在边界上是常量。**

此时有

$$X = A, \qquad Y = \frac{B}{A}$$

可以令边界上存在 $X$ 和 $Y$ 的流来实现上述边界条件。

**（二）$X$ 和 $Y$ 在边界上不存在法线方向的流量。**

此时有

$$\boldsymbol{e}_n\cdot\nabla X = \boldsymbol{e}_n\cdot\nabla Y = 0$$

其中 $\boldsymbol{e}_n$ 是边界法线方向的单位矢量。

假设由于扰动或涨落，定常解 $X_0$ 和 $Y_0$ 发生偏离，使

$$\begin{cases} X(r,t) = A+\alpha(r,t) \\[1mm] Y(r,t) = \dfrac{B}{A}+\beta(r,t) \end{cases} \tag{5.7.6}$$

其中 $\alpha$ 和 $\beta$ 是一阶小量。将上式代入式 (5.7.4)，保留 $\alpha$ 和 $\beta$ 的线性项，可得

$$\begin{cases} \dfrac{\mathrm{d}\alpha}{\mathrm{d}t} = \left(B-1+D_1\nabla^2\right)\alpha + A^2\beta \\[3mm] \dfrac{\mathrm{d}\beta}{\mathrm{d}t} = -B\alpha + \left(-A^2+D_2\nabla^2\right)\beta \end{cases} \tag{5.7.7}$$

上式是 $\alpha$ 和 $\beta$ 的线性联立方程，可以将函数 $\alpha$ 和 $\beta$ 作傅里叶展开而讨论其傅里叶分量的行为。**一个傅里叶分量称为一个模。**

为简明起见，我们考虑一维问题。设容器的长度为 $L$。对于边界条件（一），模 $n$ 可表示为

$$\begin{cases} \alpha_n(z,t) = \alpha_n e^{\omega_nt}\sin k_nz \\[1mm] \beta_n(z,t) = \beta_n e^{\omega_nt}\sin k_nz \end{cases} \tag{5.7.8a}$$

其中

$$k_n = \frac{n\pi}{L}, \qquad n=1,2,\cdots$$

对于边界条件（二），则有

$$\begin{cases} \alpha_n(z,t) = \alpha_n e^{\omega_nt}\cos k_nz \\[1mm] \beta_n(z,t) = \beta_n e^{\omega_nt}\cos k_nz \end{cases} \tag{5.7.8b}$$

其中

$$k_n = \frac{n\pi}{L}, \qquad n=0,1,2,\cdots$$

### 5.7.3 线性稳定性分析

将式 (5.7.8a) 或式 (5.7.8b) 代入式 (5.7.7)，可得

$$\begin{cases} \omega_n\alpha_n = \left(B-1-D_1k_n^2\right)\alpha_n + A^2\beta_n \\[1mm] \omega_n\beta_n = -B\alpha_n + \left(-A^2-D_2k_n^2\right)\beta_n \end{cases} \tag{5.7.9}$$

由线性代数知道，上式具有非零解的条件是

$$\begin{vmatrix} B-1-D_1k_n^2-\omega_n & A^2 \\ -B & -A^2-D_2k_n^2-\omega_n \end{vmatrix} = 0$$

整理得

$$\omega_n^2-\left(C_1-C_2\right)\omega_n + A^2B-C_1C_2 = 0 \tag{5.7.10}$$

其中

$$\begin{cases} C_1 = B-1-k_n^2D_1 \\ C_2 = A^2+k_n^2D_2 \end{cases} \tag{5.7.11}$$

二次方程 (5.7.10) 的根是

$$\omega_n^{\pm} = \frac{1}{2}\left[C_1-C_2\pm\sqrt{\left(C_1+C_2\right)^2-4A^2B}\right] \tag{5.7.12}$$

根 $\omega_n^{\pm}$ 可以是复数或实数。复根的实部或实根可以为正或为负，取决于参数 $D_1$、$D_2$、$A$、$B$ 和 $L$ 的取值。

| $\omega_n$ 的情况 | 扰动的时间行为 |
| --- | --- |
| 复数 | $\alpha_n(z,t)$、$\beta_n(z,t)$ 以 $\omega_n$ 的虚部为频率随时间振动 |
| 实部为负 | $\alpha_n$、$\beta_n$ 随时间衰减（定态稳定） |
| 实部为正 | $\alpha_n$、$\beta_n$ 随时间增强（定态失稳） |

如果 $\omega_n$ 的实部为正，$\alpha_n(z,t)$ 和 $\beta_n(z,t)$ 将随时间增强，使均匀定常解 $X_0$、$Y_0$ 失稳。**上述分析称为线性稳定性分析。它可以告诉我们在参数取何值时原来的定常解会失稳，但不能告知新态的具体形式。**

### 5.7.4 失稳条件与耗散结构

下面分别讨论 $\omega_n$ 为实根或复根的两种情形。由式 (5.7.12) 知，$\omega_n$ 为实根的条件是

$$\left(C_1+C_2\right)^2-4A^2B > 0 \tag{5.7.13}$$

将上式改写为

$$\left(C_1-C_2\right)^2+4\left(C_1C_2-A^2B\right) > 0 \tag{5.7.14}$$

比较式 (5.7.12) 和式 (5.7.14) 可知，当

$$C_1C_2-A^2B > 0 \tag{5.7.15}$$

时，$\omega_n^+$ 为正实根，上式相当于

$$B > B_n = 1+\frac{D_1}{D_2}A^2+\frac{A^2}{D_2k_n^2}+D_1k_n^2 \tag{5.7.16}$$

图 5.7.1 画出了在给定 $L$、$A$、$D_1$ 和 $D_2$ 时，$B_n$ 随 $n$ 变化的曲线。当 $B>B_n$ 时，模 $n$ 将失稳。由图中曲线可以看出，模 $n=2$ 失稳所需的 $B$ 值最小，其 $B$ 值以 $B_c$ 表示。当一个或多个模失稳时，均匀的定常解 $X_0$、$Y_0$ 将失稳。如前所述，由线性稳定性分析不能获知新态的具体形式。图 5.7.2 画出了由数值计算得到的 X 的浓度随空间变化的定常解。参量取值为 $B=4.6$、$A=2$、$L=1$、$D_1=1.6\times10^{-3}$、$D_2=6.0\times10^{-3}$，边界条件为边界上 $X=A$、$Y=\dfrac{B}{A}$。图形显示，在上述参量取值下，系统处在具有空间结构的定常态。

由式 (5.7.12) 知，$\omega_n$ 为复数的条件为

$$\left(C_1+C_2\right)^2-4A^2B < 0 \tag{5.7.17}$$

将式 (5.7.11) 的 $C_1$ 和 $C_2$ 代入，并令

$$\delta_n = 1+k_n^2\left(D_1-D_2\right) = 1+\frac{n^2\pi^2}{L^2}\left(D_1-D_2\right) \tag{5.7.18}$$

可将式 (5.7.17) 化为

$$\left(B+A^2-\delta_n\right)^2-4A^2B < 0$$

整理得

$$B^2-2\left(A^2+\delta_n\right)B+\left(A^2-\delta_n\right)^2 < 0 \tag{5.7.19}$$

令上式左方 $B$ 的二次式等于零，即

$$B^2-2\left(A^2+\delta_n\right)B+\left(A^2-\delta_n\right)^2 = 0$$

它的两个根为

$$B_{\pm} = A^2+\delta_n\pm2A\sqrt{\delta_n} = \left(A\pm\sqrt{\delta_n}\right)^2$$

则式 (5.7.19) 可表示为

$$\left(B-B_+\right)\left(B-B_-\right) < 0$$

这就是说，要使式 (5.7.17) 成立，$B$ 必然在下述区间之内：

$$\left(A-\sqrt{\delta_n}\right)^2 < B < \left(A+\sqrt{\delta_n}\right)^2 \tag{5.7.20}$$

其中隐含着另一个必要条件：

$$\delta_n > 0 \tag{5.7.21}$$

由式 (5.7.18) 知，$\delta_n>0$ 意味着

$$\frac{n^2\pi^2}{L^2}\left(D_2-D_1\right) < 1 \tag{5.7.22}$$

在满足式 (5.7.20) 和式 (5.7.21) 的情形下，模 $n$ 将随时间而空间呈周期性变化。对于边界条件（二），$n$ 可为零，因而可存在 X 和 Y 的浓度在空间均匀分布但随时间振荡的解，称为**化学钟**。

由式 (5.7.12) 知，复根 $\omega_n$ 具有正实部的条件为 $C_1-C_2>0$，或

$$B > B_n = 1+A^2+\frac{n^2\pi^2}{L^2}\left(D_1+D_2\right) \tag{5.7.23}$$

对于给定的 $L$、$A$、$D_1$ 和 $D_2$，可以根据上式画出 $B_n$ 随 $n$ 变化的曲线，如图 5.7.3 所示。当 $B>B_n$ 时，模 $n$ 将失稳。图 5.7.4 画出了由数值计算得到的在不同时刻 X 的浓度随坐标的分布，所取参量值为 $L=1$、$A=2$、$B=5.45$、$D_1=8\times10^{-3}$、$D_2=4\times10^{-3}$，边界条件为边界上 $X=A$、$Y=\dfrac{B}{A}$。图形显示，这些解类似于弦的驻波，由于非线性而呈现更为复杂的结构。

对于二维的情形，计算机模拟显示，在适当的参量取值下存在非线性行波解。

### 5.7.5 耗散结构

**上述空间、时间或空间–时间的有序结构是远离平衡的开放系统在耗散过程中通过涨落的增强和自组织而实现的。普里高津把它们称为耗散结构。**

我们知道，熵是混乱度的量度，有序状态的熵低于无序状态的熵。耗散过程存在熵产生 $\mathrm{d}_iS>0$，因此耗散过程中有序状态只能**在开放系统形成**。**开放系统与外界不断地交换物质和能量形成负熵流，导致 $\mathrm{d}_eS<0$，使 $\mathrm{d}S=\mathrm{d}_eS+\mathrm{d}_iS<0$。** 这样，耗散结构理论为远离平衡的开放系统中有序状态的形成提供了解释，深化了人们对热力学第二定律的认识。

除了反应扩散过程外，自然界中存在大量其他的耗散结构或自组织现象。

| 例子 | 现象 |
| --- | --- |
| **贝纳尔（Bénard）对流**（1900 年） | 处在温差为 $\Delta T$ 的两极间的液体薄层：$\Delta T$ 小于临界值 $\Delta T_c$ 时液体静止、只有热传导；$\Delta T>\Delta T_c$ 时液体突然开始有组织地宏观流动，由许多对流元胞组成。$\Delta T=\Delta T_c$ 时，液体由均匀无序的状态突然转变为空间有序的状态 |
| **激光** | 抽运未达到阈值时，原子主要通过自发辐射发光，光场无序；抽运超过阈值后，原子主要通过受激辐射发光，产生有序的相干光场。光场在阈值处突然由无序状态转变为有序状态，是非平衡相变的又一个例子 |
| **生命现象** | 从细胞到生命机体都是开放系统，结构极其有序；生命过程还表现出随时间作周期性变化的振荡行为（生命节奏）；生命机体不断地进行着新陈代谢、细胞分裂、形成新的机体……这一切都是持续进行的**自组织过程** |

热力学和进化论是 19 世纪自然科学的两大成就。根据热力学，孤立系统总是朝着熵增加的方向变化，因而变得越来越无序；而根据进化论，自然界的物种总是朝着从简单到复杂、从低级到高级的方向进化，因而结构变得越来越复杂。二者表象上的不同曾经使人们认为生命现象与非生命现象遵循不同的规律。**耗散结构理论指出，远离平衡的开放系统可以形成有序的状态，打开了从物理科学通向生命科学的窗口。**

---

## 本章小结

### 一、公式速查

**熵产生与局域平衡**

| 内容 | 公式 |
| --- | --- |
| 熵的分解 | $\mathrm{d}S=\mathrm{d}_eS+\mathrm{d}_iS$ (5.1.2)，$\mathrm{d}_iS\geqslant0$ |
| 局域基本方程 | $T\mathrm{d}S=\mathrm{d}U+p\mathrm{d}V-\sum_i\mu_i\mathrm{d}N_i$ (5.1.3)；$T\mathrm{d}s=\mathrm{d}u-\sum_i\mu_i\mathrm{d}n_i$ (5.1.4) |
| 局域熵平衡 | $\dfrac{\partial s}{\partial t}=-\nabla\cdot\boldsymbol{J}_s+\Theta$ (5.1.6)，$\Theta\geqslant0$ |
| 总熵平衡 | $\dfrac{\mathrm{d}S}{\mathrm{d}t}=-\oint\boldsymbol{J}_s\cdot\mathrm{d}\boldsymbol{\sigma}+\int\Theta\mathrm{d}\tau$ (5.1.7) |
| 熵流与熵产生 | $\dfrac{\mathrm{d}_eS}{\mathrm{d}t}=-\oint\boldsymbol{J}_s\cdot\mathrm{d}\boldsymbol{\sigma}$，$\dfrac{\mathrm{d}_iS}{\mathrm{d}t}=\int\Theta\mathrm{d}\tau$ (5.1.8) |
| 普遍形式 | $\Theta=\sum_k\boldsymbol{J}_k\cdot\boldsymbol{X}_k$ (5.1.23) |

**两个基本例子**

| 过程 | 流 | 力 | 熵产生率 |
| --- | --- | --- | --- |
| 单纯热传导 | $\boldsymbol{J}_q$ | $\boldsymbol{X}_q=\nabla\dfrac{1}{T}$ | $\Theta=\boldsymbol{J}_q\cdot\nabla\dfrac{1}{T}=\kappa\dfrac{(\nabla T)^2}{T^2}\geqslant0$ (5.1.14) |
| 热传导 + 扩散 | $\boldsymbol{J}_q$，$\boldsymbol{J}_n$ | $\boldsymbol{X}_q=\nabla\dfrac{1}{T}$，$\boldsymbol{X}_n=-\dfrac{1}{T}\nabla\mu$ | $\Theta=\boldsymbol{J}_q\cdot\nabla\dfrac{1}{T}-\dfrac{\boldsymbol{J}_n}{T}\cdot\nabla\mu$ (5.1.21)(5.1.22) |

**线性唯象律与昂萨格关系**

| 内容 | 公式 |
| --- | --- |
| 四种经验定律 | $\boldsymbol{J}_q=-\kappa\nabla T$ (5.2.1)；$\boldsymbol{J}_n=-D\nabla n$ (5.2.2)；$\boldsymbol{J}_e=\sigma\boldsymbol{E}$ (5.2.3)；$P_{xy}=\eta\dfrac{\mathrm{d}v}{\mathrm{d}x}$ (5.2.4) |
| 线性唯象律 | $\boldsymbol{J}_k=\sum_l L_{kl}\boldsymbol{X}_l$ (5.2.6) |
| 昂萨格倒易关系 | $L_{kl}=L_{lk}$ (5.2.8) |
| 熵产生率 | $\Theta=\sum_{kl}L_{kl}X_kX_l$ (5.2.9) |
| 正定条件 | $L_{11}>0$，$L_{11}L_{22}>L_{12}^2$ (5.2.12) |
| 泰勒展开 | $J_k=\sum_l L_{kl}X_l+\dfrac{1}{2}\sum_{l,n}L_{kln}X_lX_n+\cdots$ (5.2.15) |

**温差电（*§5.3）**

| 内容 | 公式 |
| --- | --- |
| 泽贝克效应 | $\mathrm{d}V=\varepsilon_{AB}\mathrm{d}T$ (5.3.1) |
| 佩尔捷效应 | $\boldsymbol{J}_{q\pi}=\pi_{AB}\boldsymbol{J}_e$ (5.3.2) |
| 汤姆孙效应 | $q_T=-\tau\boldsymbol{J}_e\cdot\nabla T$ (5.3.3) |
| 线性唯象律 | (5.3.4) |
| 电化学势 | $\mu=\mu_c+\mu_e$，$\mu_e=-eV$ (5.3.5) |
| 经验常量 | $\sigma=\dfrac{e^2L_{11}}{T}$ (5.3.7)；$\kappa=\dfrac{L_{11}L_{22}-L_{12}^2}{T^2L_{11}}$ (5.3.8) |
| 绝对温差电动势系数 | $\varepsilon_A=\dfrac{-L_{12}^A}{eTL_{11}^A}$，$\varepsilon_{AB}=\varepsilon_B-\varepsilon_A$ (5.3.13)(5.3.14) |
| 开尔文第二关系 | $\pi_{AB}=T\left(\varepsilon_B-\varepsilon_A\right)$ (5.3.19) |
| 开尔文第一关系 | $\dfrac{\mathrm{d}\pi_{AB}}{\mathrm{d}T}+\left(\tau_A-\tau_B\right)=\varepsilon_B-\varepsilon_A$ (5.3.24)；$\tau=T\dfrac{\mathrm{d}\varepsilon}{\mathrm{d}T}$ (5.3.23) |

**最小熵产生定理（*§5.4）**

| 内容 | 公式 |
| --- | --- |
| 熵产生率 | $P=\int L_{qq}\left(\nabla\dfrac{1}{T}\right)^2\mathrm{d}\tau$ (5.4.3) |
| 定理 | $\dfrac{\mathrm{d}P}{\mathrm{d}t}\leqslant0$ (5.4.6) |
| 耦合情形 | $\Theta=L_{11}X_1^2+2L_{12}X_1X_2+L_{22}X_2^2$ (5.4.7) |
| 有约束时 | $X_2=-\dfrac{L_{21}}{L_{22}}X_1$，$\Theta=\left(L_{11}-\dfrac{L_{12}L_{21}}{L_{22}}\right)X_1^2$ (5.4.10)(5.4.11) |

**化学反应与扩散（*§5.5）**

| 内容 | 公式 |
| --- | --- |
| 连续性方程 | $\dfrac{\partial n_i}{\partial t}=-\nabla\cdot\boldsymbol{J}_i+\sum_\rho\nu_{i\rho}\omega_\rho$ (5.5.8) |
| 熵流密度 | $\boldsymbol{J}_s=-\sum_i\dfrac{\mu_i}{T}\boldsymbol{J}_i$ (5.5.12) |
| 熵产生率 | $\Theta=-\sum_i\boldsymbol{J}_i\cdot\nabla\dfrac{\mu_i}{T}+\sum_\rho\dfrac{a_\rho}{T}\omega_\rho$ (5.5.13) |
| 局部化学亲和势 | $a_\rho=-\sum_i\mu_i\nu_{i\rho}$ (5.5.10) |
| 流与力 | $\boldsymbol{J}_i^{\mathrm{di}}=\boldsymbol{J}_i$，$\boldsymbol{X}_i^{\mathrm{di}}=-\nabla\dfrac{\mu_i}{T}$；$J_\rho^{\mathrm{ch}}=\omega_\rho$，$X_\rho^{\mathrm{ch}}=\dfrac{a_\rho}{T}$ (5.5.14) |
| 局部平衡常量 | $\ln K_\rho=-\dfrac{1}{kT}\sum_i\nu_{i\rho}g_i$ (5.5.17)；$a_\rho=kT\ln\dfrac{K_\rho}{\prod_i x_i^{\nu_{i\rho}}}$ (5.5.18) |
| 速率与亲和势 | $\omega_\rho=k_+x_A\left(1-e^{-a_\rho/(kT)}\right)$ (5.5.24)；$a_\rho\ll kT$ 时 $\omega_\rho=\omega_+\dfrac{a_\rho}{kT}$ (5.5.26) |

**非线性区发展判据（*§5.6）**

| 内容 | 公式 |
| --- | --- |
| 分解 | $\dfrac{\mathrm{d}P}{\mathrm{d}t}=\dfrac{\mathrm{d}_XP}{\mathrm{d}t}+\dfrac{\mathrm{d}_JP}{\mathrm{d}t}$ (5.6.1) |
| 力的贡献 | $\dfrac{\mathrm{d}_XP}{\mathrm{d}t}=-\dfrac{1}{T}\displaystyle\int\sum_{i,j}\dfrac{\partial\mu_i}{\partial n_j}\dfrac{\partial n_j}{\partial t}\dfrac{\partial n_i}{\partial t}\mathrm{d}\tau\leqslant0$ (5.6.3)(5.6.4) |
| 线性区 | $\dfrac{\mathrm{d}_JP}{\mathrm{d}t}=\dfrac{\mathrm{d}_XP}{\mathrm{d}t}$，$\dfrac{\mathrm{d}P}{\mathrm{d}t}=2\dfrac{\mathrm{d}_XP}{\mathrm{d}t}\leqslant0$ (5.6.5)(5.6.6) |
| 非线性区 | $\dfrac{\mathrm{d}_JP}{\mathrm{d}t}$ 符号不定，$\dfrac{\mathrm{d}P}{\mathrm{d}t}$ 符号不定 → 可发生非平衡相变 |

**三分子模型（*§5.7）**

| 内容 | 公式 |
| --- | --- |
| 反应步骤 | (5.7.1) |
| 动力学方程 | $\dfrac{\mathrm{d}X}{\mathrm{d}t}=A-(B+1)X+X^2Y+D_1\nabla^2X$；$\dfrac{\mathrm{d}Y}{\mathrm{d}t}=BX-X^2Y+D_2\nabla^2Y$ (5.7.4) |
| 均匀定态 | $X_0=A$，$Y_0=\dfrac{B}{A}$ (5.7.5) |
| 线性化 | (5.7.7) |
| 色散关系 | $\omega_n^2-(C_1-C_2)\omega_n+A^2B-C_1C_2=0$ (5.7.10)；$C_1=B-1-k_n^2D_1$，$C_2=A^2+k_n^2D_2$ (5.7.11) |
| 实根失稳 | $B>B_n=1+\dfrac{D_1}{D_2}A^2+\dfrac{A^2}{D_2k_n^2}+D_1k_n^2$ (5.7.16) |
| 复根失稳（化学钟） | $B>B_n=1+A^2+\dfrac{n^2\pi^2}{L^2}\left(D_1+D_2\right)$ (5.7.23) |

### 二、概念脉络

```
§1.16 dS ≥ đQ/T（第二定律的数学表述）
        │
        ▼
§5.1 局域平衡假设：把小部分当作平衡态，仍用 T ds = du − Σ μ_i dn_i
        │
        ├──► dS = d_eS + d_iS，d_iS ≥ 0
        ├──► ∂s/∂t = −∇·J_s + Θ，Θ ≥ 0
        └──► Θ = Σ J_k·X_k
        │
        ├──► 例：热传导  → X_q = ∇(1/T)，J_q = −κ∇T
        └──► 例：热传导 + 扩散 → X_n = −(1/T)∇μ
        │
        ▼
§5.2 线性唯象律 J_k = Σ L_kl X_l
     昂萨格倒易关系 L_kl = L_lk
     第二定律 → Θ 正定 → L_11 > 0，L_11L_22 > L_12²
        │
        ├──► *§5.3 温差电：σ、κ、ε_AB；开尔文第一、第二关系
        │
        ├──► *§5.4 最小熵产生定理（线性区）dP/dt ≤ 0，非平衡定态稳定
        │
        ▼
*§5.5 反应–扩散
     a_ρ = −Σ μ_i ν_iρ；J_s = −Σ(μ_i/T)J_i
        ├──► a_ρ ≪ kT：ω_ρ ∝ a_ρ（线性区）
        └──► a_ρ ≫ kT：ω_ρ = ω_+（单向，非线性）
        │
        ▼
*§5.6 非线性区：dP/dt = d_XP/dt + d_JP/dt
     d_XP/dt ≤ 0 恒成立；d_JP/dt 符号不定
        │
        ▼
*§5.7 三分子模型（Brusselator）
     dX/dt = A − (B+1)X + X²Y + D_1∇²X
     dY/dt = BX − X²Y + D_2∇²Y
     线性稳定性分析 → B > B_n 时均匀定态失稳
        ├──► 实根失稳 → 空间结构（图 5.7.2）
        └──► 复根失稳 → 时间振荡 / 化学钟
        │
        ▼
     耗散结构：远离平衡的开放系统，d_eS < 0、d_iS > 0，dS < 0
     例：BZ 反应、Bénard 对流、激光、生命自组织
```

### 三、易错点

1. **$\mathrm{d}_eS$ 与 $\mathrm{d}_iS$ 的性质不同**：$\mathrm{d}_eS$ 是交换项，可正可负；$\mathrm{d}_iS$ 是产生项，恒 $\geqslant0$。孤立系 $\mathrm{d}_eS=0$，才有 $\mathrm{d}S=\mathrm{d}_iS\geqslant0$。
2. **式 (5.1.3) 的符号**：$T\mathrm{d}S=\mathrm{d}U+p\,\mathrm{d}V-\sum_i\mu_i\mathrm{d}N_i$，与 $dU=TdS-pdV+\sum\mu_i dN_i$ 是同一式的移项；不要把 $p\,\mathrm{d}V$ 的符号抄错。
3. **局域平衡是假设，不是定理**：它要求小部分的弛豫时间远短于整个系统的弛豫时间；其正确性由推论与实验相符来检验。
4. **略去 $p\,\mathrm{d}V$ 是有条件的**：只有在不涉及流体力学问题时才能略去；涉及流体运动时必须保留。
5. **$\Theta$ 与 $P$ 不同**：$\Theta$ 是单位时间单位体积的熵产生率，$P=\int\Theta\,\mathrm{d}\tau$ 才是整个系统的熵产生率；$\mathrm{d}_iS/\mathrm{d}t=P$。
6. **$\boldsymbol{J}_s$ 与 $\boldsymbol{J}_q$ 不是一回事**：单纯热传导时 $\boldsymbol{J}_s=\boldsymbol{J}_q/T$；有物质流时还要加上 $-\sum_i(\mu_i/T)\boldsymbol{J}_i$（式 (5.5.12)）。
7. **热力学力不是"力"**：热流动力是 $\boldsymbol{X}_q=\nabla(1/T)$ 而不是 $\nabla T$；粒子流动力是 $\boldsymbol{X}_n=-(1/T)\nabla\mu$ 而不是 $\nabla\mu$。换一种"流–力"配对，$\Theta=\sum JX$ 与昂萨格关系都要随之改写。
8. **昂萨格关系只在"按式 (5.1.23) 配对"时成立**：若把流或力作线性组合，动理系数矩阵要做相应的相似变换，$L_{kl}=L_{lk}$ 不再直接适用于新系数。
9. **$\Theta\geqslant0$ 给出的是矩阵正定**（$L_{11}>0$，$L_{11}L_{22}>L_{12}^2$），不是每个 $L_{kl}$ 都为正；交叉系数完全可以为负。
10. **最小熵产生定理只适用于线性区**（动理系数为常量）。非线性区 $\mathrm{d}_JP/\mathrm{d}t$ 符号不定，$\mathrm{d}P/\mathrm{d}t$ 可正可负，这正是耗散结构可能出现的前提。
11. **温差电的符号约定**：高温端电动势驱使电流由 A 流向 B 时 $\varepsilon_{AB}$ 为正；$\varepsilon_{AB}=\varepsilon_B-\varepsilon_A$ 是两个"绝对温差电动势系数"之差。
12. **佩尔捷热只出现在接头处**：均匀导体内部没有佩尔捷热；接头处因两侧 $\varepsilon$ 不同，$\boldsymbol{J}_{q\pi}=T\varepsilon\boldsymbol{J}_e$ 发生突变，才表现为吸热或放热。
13. **汤姆孙系数 $\tau=T\dfrac{\mathrm{d}\varepsilon}{\mathrm{d}T}$ 可正可负**，取决于 $\mathrm{d}\varepsilon/\mathrm{d}T$ 的符号，不要与佩尔捷系数 $\pi_{AB}$ 混淆。
14. **局部化学亲和势 $a_\rho=-\sum_i\mu_i\nu_{i\rho}$ 与 §4.8 的亲和势 $A=-\Delta G$ 是同一思想**：$a_\rho=0$ 就是局部平衡，式 (5.5.19) 就是质量作用律。
15. **反应速率与亲和势的线性关系只在 $a_\rho\ll kT$ 时成立**；$a_\rho\gg kT$ 时 $\omega_\rho=\omega_+$，反应单向进行，不能再用线性唯象律。
16. **三分子模型中 $B$ 是控制参量**：$B$ 超过临界值 $B_c$ 时均匀定态失稳；线性稳定性分析只能给出失稳条件，**给不出新态的具体形式**，必须解非线性方程。
17. **边界条件决定模的形式**：定值边界取 $\sin k_nz$（$n=1,2,\cdots$），零流量边界取 $\cos k_nz$（$n=0,1,2,\cdots$）；只有零流量边界允许 $n=0$，对应空间均匀、时间振荡的"化学钟"。
18. **耗散结构要求开放**：$\mathrm{d}_iS>0$ 恒成立，只有靠外界提供的负熵流 $\mathrm{d}_eS<0$ 才能维持 $\mathrm{d}S<0$ 的有序状态；孤立系中不可能出现耗散结构。

### 四、名词中英对照

| 中文 | 英文 | 中文 | 英文 |
| --- | --- | --- | --- |
| 不可逆过程热力学 | thermodynamics of irreversible processes | 局域平衡 | local equilibrium |
| 熵产生 | entropy production | 熵流密度 | entropy flux density |
| 局域熵产生率 | local entropy production rate | 热力学流 | thermodynamic flux |
| 热力学力 | thermodynamic force | 线性唯象律 | linear phenomenological law |
| 动理系数 | kinetic coefficient | 昂萨格倒易关系 | Onsager reciprocal relations |
| 傅里叶定律 | Fourier's law | 菲克定律 | Fick's law |
| 欧姆定律 | Ohm's law | 牛顿黏性定律 | Newton's law of viscosity |
| 导热系数 | thermal conductivity | 扩散系数 | diffusion coefficient |
| 电导率 | electrical conductivity | 黏度 | viscosity |
| 热流密度 | heat flux density | 粒子流密度 | particle flux density |
| 内能流密度 | internal energy flux density | 电流密度 | electric current density |
| 温差电现象 | thermoelectric phenomena | 泽贝克效应 | Seebeck effect |
| 佩尔捷效应 | Peltier effect | 汤姆孙效应 | Thomson effect |
| 温差电动势系数 | thermoelectric power | 佩尔捷系数 | Peltier coefficient |
| 汤姆孙系数 | Thomson coefficient | 绝对温差电动势系数 | absolute thermoelectric power |
| 焦耳热 | Joule heat | 开尔文第一关系 | Kelvin's first relation |
| 开尔文第二关系 | Kelvin's second relation | 最小熵产生定理 | minimum entropy production theorem |
| 非平衡定态 | nonequilibrium steady state | 局部化学亲和势 | local chemical affinity |
| 局部平衡常量 | local equilibrium constant | 反应速率 | reaction rate |
| 扩散过程 | diffusion process | 反应–扩散过程 | reaction–diffusion process |
| 线性稳定性分析 | linear stability analysis | 控制参量 | control parameter |
| 三分子模型 | trimolecular model | 布鲁塞尔模型 | Brusselator |
| 化学钟 | chemical clock | 耗散结构 | dissipative structure |
| 自组织 | self-organization | 非平衡相变 | nonequilibrium phase transition |
| 贝纳尔对流 | Bénard convection | BZ 反应 | Belousov–Zhabotinskii reaction |
