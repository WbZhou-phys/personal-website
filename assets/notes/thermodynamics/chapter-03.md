---
tags:
  - 热统
  - 热力学
  - 第三章
created: 2026-10-08
---

# 第三章　单元系的相变

> **教材**：汪志诚《热力学·统计物理》（第六版）
> **范围**：§3.1–§3.9（教材 pp. 63–86；其中 §3.6、§3.8 为带星号的选读内容）

## 本章地图

| 节     | 标题            | 一句话                                                                 |
| ----- | ------------- | ------------------------------------------------------------------- |
| §3.1  | 热动平衡判据        | 熵判据、自由能判据、吉布斯函数判据；均匀系的平衡条件与稳定性条件                                    |
| §3.2  | 开系的热力学基本方程    | 引入化学势 $\mu$，把基本方程推广到物质量可变的情形；巨热力势 $J$                               |
| §3.3  | 单元系的复相平衡条件    | 由熵判据导出 $T^\alpha=T^\beta$、$p^\alpha=p^\beta$、$\mu^\alpha=\mu^\beta$ |
| §3.4  | 单元复相系的平衡性质    | 相图与三相点；克拉珀龙方程、蒸气压方程                                                 |
| §3.5  | 临界点和气液两相的转变   | 杠杆规则、麦克斯韦等面积法则、临界系数、对应态定律                                           |
| *§3.6 | 液滴的形成（选读）     | 曲面分界面的力学平衡 $p^\alpha=p^\beta+2\sigma/r$；过饱和蒸气与过热液体                  |
| §3.7  | 相变的分类         | 一级相变与连续相变；埃伦菲斯特方程                                                   |
| *§3.8 | 临界现象和临界指数（选读） | $\alpha$、$\beta$、$\gamma$、$\gamma'$、$\delta$ 的实验规律与普适性              |
| §3.9  | 朗道连续相变理论      | 序参量、对称破缺、朗道自由能；平均场临界指数                                              |

**主线**：§3.1 建立判据 → §3.2 把基本方程推广到开系（引入 $\mu$）→ §3.3–§3.4 导出复相平衡条件并给出相图上的定量关系 → §3.5 在 $p$-$V$ 图上讨论临界点 → *§3.6 计入表面相 → §3.7 给相变分类 → *§3.8–§3.9 讨论连续相变的临界现象与朗道理论。

---

## §3.1　热动平衡判据

### 3.1.1 熵判据

本章和第四章讨论相变和化学变化问题。作为基础，本节讲述在热力学中如何判定一个系统的平衡状态。

§1.16 中根据热力学第二定律证明了熵增加原理。熵增加原理指出，**孤立系统的熵永不减少。孤立系统中发生的趋向平衡过程，必朝着熵增加的方向进行。如果孤立系统已经达到了熵为极大的状态，就不可能再发生任何热力学意义上的变化，系统就达到了平衡态。** 我们可以利用熵函数这一性质来判定孤立系统的平衡态，这称为**熵判据**。

为了判定孤立系统的某一状态是否为平衡态，可以设想系统围绕该状态发生各种可能的**虚变动**。所谓虚变动是理论上假想的、满足外加约束条件的各种可能的变动。在应用数学方法求各种可能的虚变动所引起的熵变时，外加约束条件（现在是孤立系条件）需要用函数形式表示。孤立系与其他物体既没有热量的交换，也没有功的交换；如果只有体积变化功，**孤立系条件相当于体积不变和内能不变**。

| 虚变动引起的熵变 | 状态 |
| --- | --- |
| $\Delta S<0$ | 熵具有极大值，是**稳定平衡态** |
| $\Delta S=0$ | **中性平衡态** |

因此，**孤立系统处在稳定平衡状态的充分必要条件为**

$$\Delta S < 0 \tag{3.1.1}$$

将 $S$ 作泰勒展开，准确到二级，有

$$\Delta S = \delta S + \frac{1}{2}\delta^2 S \tag{3.1.2}$$

根据数学上熟知的结果：当熵函数的一级微分 $\delta S=0$ 时，熵函数有极值；当一级微分 $\delta S=0$、二级微分 $\delta^2 S<0$ 时，熵函数有极大值。**由 $\delta S=0$ 可以得到平衡条件，由 $\delta^2 S<0$ 可以得到平衡的稳定性条件。**

如果熵函数的极大值不止一个，则其中**最大的极大值**相应于稳定平衡，其他较小的极大值相应于**亚稳平衡**。亚稳平衡是这样一种平衡，它对于无穷小的变动是稳定的，对于有限大的变动则是不稳定的；如果发生较大的涨落或者通过某种触发作用，系统就可能由亚稳平衡状态过渡到更加稳定（熵更大）的平衡状态。

> **脚注**（教材 p.63）：如果 $\delta S=0$，$\delta^2 S=0$，为了判定状态是否是稳定的平衡态，要考察熵的高级微分。在 §3.5 中讨论液气临界点将遇到这种情形。

**熵判据是基本的平衡判据，适用于孤立系统。** 如果把参与变化的全部物体都包括在系统之内，原则上可以对各种热动平衡问题作出回答。不过在实际上，对于某些经常遇到的物理条件，引入其他判据是更为方便的。

### 3.1.2 自由能判据与吉布斯函数判据

§1.18 中讲过，在等温等体条件下，系统的自由能永不增加；在等温等压条件下，系统的吉布斯函数永不增加。因此在上述条件下，系统中发生的趋向平衡过程，必朝着自由能或吉布斯函数减小的方向进行。**如果系统达到了自由能或吉布斯函数极小的状态，就不可能再发生任何热力学意义上的变化，系统就达到了平衡态。** 根据自由能或吉布斯函数的上述性质，对等温等体系统或等温等压系统进行判断，称为**自由能判据**和**吉布斯函数判据**。

| 判据 | 适用条件 | 稳定平衡的充要条件 | 泰勒展开（准确到二级） | 平衡条件 / 稳定性条件 |
| --- | --- | --- | --- | --- |
| 熵判据 | 孤立系 | $\Delta S<0$ (3.1.1) | $\Delta S=\delta S+\dfrac{1}{2}\delta^2 S$ (3.1.2) | $\delta S=0$ / $\delta^2 S<0$ |
| 自由能判据 | 等温等体 | $\Delta F>0$ (3.1.3) | $\Delta F=\delta F+\dfrac{1}{2}\delta^2 F$ (3.1.4) | $\delta F=0$ / $\delta^2 F>0$ |
| 吉布斯函数判据 | 等温等压 | $\Delta G>0$ (3.1.5) | $\Delta G=\delta G+\dfrac{1}{2}\delta^2 G$ (3.1.6) | $\delta G=0$ / $\delta^2 G>0$ |

通过类似的分析可以知道，等温等体系统处在稳定平衡状态的充分必要条件为

$$\Delta F > 0 \tag{3.1.3}$$

将 $F$ 作泰勒展开，准确到二级，有

$$\Delta F = \delta F + \frac{1}{2}\delta^2 F \tag{3.1.4}$$

由 $\delta F=0$ 和 $\delta^2 F>0$ 可以确定平衡条件和平衡的稳定性条件。

等温等压系统处在稳定平衡状态的充分必要条件为

$$\Delta G > 0 \tag{3.1.5}$$

将 $G$ 作泰勒展开，准确到二级，有

$$\Delta G = \delta G + \frac{1}{2}\delta^2 G \tag{3.1.6}$$

由 $\delta G=0$ 和 $\delta^2 G>0$ 可以确定平衡条件和平衡的稳定性条件。

类似地，对于等温等体系统和等温等压系统，也可能出现**亚稳平衡**或**中性平衡**等情况。

除了熵、自由能和吉布斯函数判据以外，还可以根据其他热力学函数的性质（参阅习题 3.1）进行判断。例如，根据在熵和体积不变的条件下系统的内能永不增加的性质，可以得到**内能判据**。

### 3.1.3 均匀系统的热动平衡条件与平衡的稳定性条件

作为热动平衡判据的应用，我们讨论均匀系统的热动平衡条件和平衡的稳定性条件。

设有一个孤立的均匀系统，考虑系统中任意一个小部分（图 3.1.1）。这部分虽小，但仍含有大量的微观粒子，可以看作一个宏观系统。我们把这个小部分称作**子系统**，而把系统的其他部分看作子系统的**介质**。以不带下标的量表示子系统的热力学量，带有下标 0 的量表示介质的热力学量，例如 $T$、$p$ 和 $T_0$、$p_0$ 分别表示子系统和介质的温度、压强。

设想子系统发生一个虚变动，其内能和体积的变化分别为 $\delta U$ 和 $\delta V$。由于整个系统是孤立的，介质的内能和体积应有相应的变化 $\delta U_0$ 和 $\delta V_0$，使

$$\begin{cases} \delta U + \delta U_0 = 0 \\ \delta V + \delta V_0 = 0 \end{cases} \tag{3.1.7}$$

熵是广延量，虚变动引起整个系统的熵变 $\Delta\widetilde{S}$ 为 $\Delta\widetilde{S}=\Delta S+\Delta S_0$。将 $S$ 和 $S_0$ 作泰勒展开，准确到二级，有

$$\Delta S = \delta S + \frac{1}{2}\delta^2 S, \qquad \Delta S_0 = \delta S_0 + \frac{1}{2}\delta^2 S_0$$

在稳定的平衡状态下，整个孤立系统的熵应取极大值。熵函数的极值要求

$$\delta\widetilde{S} = \delta S + \delta S_0 = 0 \tag{3.1.8}$$

根据热力学基本方程，有

$$\delta S = \frac{\delta U + p\,\delta V}{T}, \qquad \delta S_0 = \frac{\delta U_0 + p_0\,\delta V_0}{T_0}$$

将以上两式代入式 (3.1.8)，并考虑到式 (3.1.7)，可得

$$\delta\widetilde{S} = \delta U\left(\frac{1}{T} - \frac{1}{T_0}\right) + \delta V\left(\frac{p}{T} - \frac{p_0}{T_0}\right) = 0$$

因为在虚变动中 $\delta U$ 和 $\delta V$ 可以独立地改变，$\delta\widetilde{S}=0$ 要求

$$T = T_0, \qquad p = p_0 \tag{3.1.9}$$

上式表明，**达到平衡时子系统与介质具有相同的温度和压强**。如前所述，子系统是整个系统中任意的一个小部分，这意味着**达到平衡时整个系统的温度和压强是均匀的**。

如果熵函数的二级微分是负的，即

$$\delta^2\widetilde{S} = \delta^2 S + \delta^2 S_0 < 0 \tag{3.1.10}$$

则熵函数将具有极大值。由于介质比子系统大得多（物质的量 $n_0 \gg n$），当子系统发生变动而内能和体积有 $\delta U$ 和 $\delta V$ 的变化时，$|\delta^2 S_0| \ll |\delta^2 S|$。因此可以忽略 $\delta^2 S_0$，将式 (3.1.10) 近似为

$$\delta^2\widetilde{S} \approx \delta^2 S < 0 \tag{3.1.11}$$

> **脚注**（教材 p.65）：这一点可以这样理解——$S$、$U$、$V$ 都是广延量（与物质的量成正比），其一阶偏导数 $\dfrac{\partial S}{\partial U}$、$\dfrac{\partial S}{\partial V}$ 是强度量（与物质的量无关），二阶偏导数 $\dfrac{\partial^2 S}{\partial U^2}$、$\dfrac{\partial^2 S}{\partial U\partial V}$、$\dfrac{\partial^2 S}{\partial V^2}$ 与物质的量成反比。既然 $n_0 \gg n$，即有 $|\delta^2 S_0| \ll |\delta^2 S|$。$\delta^2 S$ 的具体表达式参阅习题 3.2，它的各项分别与 $C_V$ 或 $V$ 成反比。

根据泰勒展开公式，式 (3.1.11) 为

$$\delta^2 S = \left[\left(\frac{\partial^2 S}{\partial U^2}\right)(\delta U)^2 + 2\frac{\partial^2 S}{\partial U\partial V}\delta U\,\delta V + \left(\frac{\partial^2 S}{\partial V^2}\right)(\delta V)^2\right] < 0 \tag{3.1.12}$$

选 $T$、$V$ 为独立变量，通过导数变换可以将式 (3.1.12) 的二次型化为平方和，而有（证明见本节末）

$$\delta^2 S = -\frac{C_V}{T^2}(\delta T)^2 + \frac{1}{T}\left(\frac{\partial p}{\partial V}\right)_T (\delta V)^2 < 0 \tag{3.1.13}$$

如果要求 $\delta^2 S$ 对于各种可能的虚变动都小于零，应有

$$C_V > 0, \qquad \left(\frac{\partial p}{\partial V}\right)_T < 0 \tag{3.1.14}$$

式 (3.1.14) 是**平衡的稳定性条件**。

- 假如子系统的温度由于涨落或某种外界影响而略高于介质，热量将从子系统传递到介质；根据平衡稳定性条件 $C_V>0$，热量的传递将使子系统的温度降低，从而恢复平衡。
- 假如子系统的体积由于某种原因发生收缩，根据平衡稳定性条件 $\left(\dfrac{\partial p}{\partial V}\right)_T<0$，子系统的压强将增高而略高于介质的压强，于是子系统膨胀而恢复平衡。

这就是说，**如果平衡稳定性条件得到满足，当系统对平衡发生某种偏离时，系统中将会自发产生相应的过程，以恢复系统的平衡**。

平衡的稳定性条件既适用于均匀系统的任何部分，显然也适用于整个均匀系统。

### 3.1.4 式 (3.1.13) 的证明

现在证明式 (3.1.13)。将 $\delta^2 S$ 改写为

$$\delta^2 S = \left[\frac{\partial}{\partial U}\left(\frac{\partial S}{\partial U}\right)\delta U + \frac{\partial}{\partial V}\left(\frac{\partial S}{\partial U}\right)\delta V\right]\delta U + \left[\frac{\partial}{\partial U}\left(\frac{\partial S}{\partial V}\right)\delta U + \frac{\partial}{\partial V}\left(\frac{\partial S}{\partial V}\right)\delta V\right]\delta V \tag{3.1.15}$$

但由热力学基本方程

$$T\,\mathrm{d}S = \mathrm{d}U + p\,\mathrm{d}V$$

可得

$$\left(\frac{\partial S}{\partial U}\right)_V = \frac{1}{T}, \qquad \left(\frac{\partial S}{\partial V}\right)_U = \frac{p}{T}$$

代入式 (3.1.15)，可将式 (3.1.12) 表达为

$$\begin{aligned}\delta^2 S &= \left[\frac{\partial}{\partial U}\left(\frac{1}{T}\right)\delta U + \frac{\partial}{\partial V}\left(\frac{1}{T}\right)\delta V\right]\delta U + \left[\frac{\partial}{\partial U}\left(\frac{p}{T}\right)\delta U + \frac{\partial}{\partial V}\left(\frac{p}{T}\right)\delta V\right]\delta V \\ &= \left[\delta\left(\frac{1}{T}\right)\delta U + \delta\left(\frac{p}{T}\right)\delta V\right] < 0 \end{aligned} \tag{3.1.16}$$

以 $T$、$V$ 为自变量，有

$$\delta U = \left(\frac{\partial U}{\partial T}\right)_V \delta T + \left(\frac{\partial U}{\partial V}\right)_T \delta V = C_V\,\delta T + \left[T\left(\frac{\partial p}{\partial T}\right)_V - p\right]\delta V$$

$$\delta\left(\frac{1}{T}\right) = \frac{\partial}{\partial T}\left(\frac{1}{T}\right)\delta T + \left(\frac{\partial}{\partial V}\frac{1}{T}\right)\delta V = -\frac{1}{T^2}\delta T$$

$$\delta\left(\frac{p}{T}\right) = \left(\frac{\partial}{\partial T}\frac{p}{T}\right)_V \delta T + \left(\frac{\partial}{\partial V}\frac{p}{T}\right)_T \delta V = \frac{1}{T^2}\left[T\left(\frac{\partial p}{\partial T}\right)_V - p\right]\delta T + \frac{1}{T}\left(\frac{\partial p}{\partial V}\right)_T \delta V$$

将以上各式代入式 (3.1.16) 即得式 (3.1.13)。

---

## §3.2　开系的热力学基本方程

### 3.2.1 开系的吉布斯函数与化学势

在讲述如何判定系统的平衡状态之后，我们将在本章中讨论单元系的相变问题。

**单元系**指化学上纯的物质系统，它只含一种化学组分（一个组元）。§1.1 中讲过，如果一个系统不是均匀的，但可以分为若干个均匀的部分，该系统称为**复相系**。例如，水和水蒸气共存构成一个单元两相系，水为一个相，水蒸气为另一个相；冰、水和水蒸气共存构成一个单元三相系，冰、水和水蒸气各为一个相。

§1.1 中还讲过，在热力学中需要用四类参量来描述一个均匀系统的平衡状态，均匀系统的热力学函数可以表达为这四类参量的函数。对于复相系中的每一个相，也需要用四类参量来描述它的平衡态，各相的热力学函数可以表达为各自参量的函数。但是，这里有两点很重要的区别：

1. **以前所讨论的均匀系都是闭系**，物质的量是不变的。现在物质可以由一相变到另一相，一个相的质量或物质的量是可变的，是一个**开系**。
2. **整个复相系要处于平衡，必须满足一定的平衡条件**，各相的状态参量不完全是独立的变量。

本节先讨论开系的热力学方程，复相系的平衡条件将在下节讨论。

先考虑吉布斯函数。根据式 (2.1.10)，吉布斯函数的全微分为

$$\mathrm{d}G = -S\,\mathrm{d}T + V\,\mathrm{d}p \tag{3.2.1}$$

上式适用于物质的量不发生变化的情况。它给出在系统的两个邻近的平衡态，其吉布斯函数之差与温度、压强之差的关系。**吉布斯函数是一个广延量。当物质的量发生变化时，吉布斯函数显然也将发生变化。** 所以对于开系，式 (3.2.1) 应推广为

$$\mathrm{d}G = -S\,\mathrm{d}T + V\,\mathrm{d}p + \mu\,\mathrm{d}n \tag{3.2.2}$$

式中右方第三项代表由于物质的量改变了 $\mathrm{d}n$ 所引起的吉布斯函数的改变，其中

$$\mu = \left(\frac{\partial G}{\partial n}\right)_{T,p} \tag{3.2.3}$$

称为**化学势**。**它等于在温度和压强保持不变的条件下，增加 1 mol 物质时吉布斯函数的增量。**

由于吉布斯函数是广延量，系统的吉布斯函数等于物质的量 $n$ 与摩尔吉布斯函数 $G_m(T,p)$ 之积：

$$G(T,p,n) = n\,G_m(T,p) \tag{3.2.4}$$

因此

$$\mu = \left(\frac{\partial G}{\partial n}\right)_{T,p} = G_m \tag{3.2.5}$$

这就是说，**化学势 $\mu$ 等于摩尔吉布斯函数**。这个结论适用于单元系。对于含有多种化学组分的系统，其化学势将在第四章中讨论。

由式 (3.2.2) 可知，$G$ 是以 $T$、$p$、$n$ 为独立变量的特性函数。如果已知 $G(T,p,n)$，其他热力学量可以通过下列偏导数分别求得：

$$S = -\left(\frac{\partial G}{\partial T}\right)_{p,n}, \qquad V = \left(\frac{\partial G}{\partial p}\right)_{T,n}, \qquad \mu = \left(\frac{\partial G}{\partial n}\right)_{T,p} \tag{3.2.6}$$

### 3.2.2 开系的热力学基本方程

根据 $U = G + TS - pV$ 及式 (3.2.2)，容易求得内能的全微分为

$$\mathrm{d}U = T\,\mathrm{d}S - p\,\mathrm{d}V + \mu\,\mathrm{d}n \tag{3.2.7}$$

**式 (3.2.7) 就是开系的热力学基本方程。它是式 (1.14.6) 的推广。** 由式 (3.2.7) 可知，$U$ 是以 $S$、$V$、$n$ 为独立变量的特性函数。同理可以求得焓和自由能的全微分：

$$\mathrm{d}H = T\,\mathrm{d}S + V\,\mathrm{d}p + \mu\,\mathrm{d}n \tag{3.2.8}$$

（$H$ 是以 $S$、$p$、$n$ 为独立变量的特性函数。）

$$\mathrm{d}F = -S\,\mathrm{d}T - p\,\mathrm{d}V + \mu\,\mathrm{d}n \tag{3.2.9}$$

（$F$ 是以 $T$、$V$、$n$ 为独立变量的特性函数。）

### 3.2.3 巨热力势

定义一个热力学函数

$$J = F - \mu n \tag{3.2.10}$$

称为**巨热力势**。它的全微分为

$$\mathrm{d}J = -S\,\mathrm{d}T - p\,\mathrm{d}V - n\,\mathrm{d}\mu \tag{3.2.11}$$

$J$ 是以 $T$、$V$、$\mu$ 为独立变量的特性函数。如果已知 $J(T,V,\mu)$，其他热力学量可以通过下列偏导数分别求得：

$$S = -\left(\frac{\partial J}{\partial T}\right)_{V,\mu}, \qquad p = -\left(\frac{\partial J}{\partial V}\right)_{T,\mu}, \qquad n = -\left(\frac{\partial J}{\partial \mu}\right)_{T,V} \tag{3.2.12}$$

由式 (3.2.10) 知，巨热力势 $J$ 也可表为

$$J = F - G = -pV \tag{3.2.13}$$

第八章的玻色（费米）统计理论和第九章的巨正则系综理论就是用统计物理方法求出巨热力势 $J$ 作为 $T$、$V$、$\mu$ 的函数，进而求其他热力学函数的。

---

## §3.3　单元系的复相平衡条件

现在讨论单元复相系达到平衡所要满足的条件。

考虑一个单元两相系，这个单元两相系构成一个孤立系统。我们用指标 $\alpha$ 和 $\beta$ 表示两个相，用 $U^\alpha$、$V^\alpha$、$n^\alpha$ 和 $U^\beta$、$V^\beta$、$n^\beta$ 分别表示 $\alpha$ 相和 $\beta$ 相的内能、体积、物质的量。整个系统既然是孤立系统，它的总内能、总体积、总物质的量应是固定的，即

$$\begin{cases} U^\alpha + U^\beta = \text{常量} \\ V^\alpha + V^\beta = \text{常量} \\ n^\alpha + n^\beta = \text{常量} \end{cases} \tag{3.3.1}$$

设想系统发生一个虚变动。在虚变动中 $\alpha$ 相和 $\beta$ 相的内能、体积、物质的量分别发生虚变动 $\delta U^\alpha$、$\delta V^\alpha$、$\delta n^\alpha$ 和 $\delta U^\beta$、$\delta V^\beta$、$\delta n^\beta$。孤立系条件要求

$$\begin{cases} \delta U^\alpha + \delta U^\beta = 0 \\ \delta V^\alpha + \delta V^\beta = 0 \\ \delta n^\alpha + \delta n^\beta = 0 \end{cases} \tag{3.3.2}$$

由式 (3.2.7) 知，两相的熵变分别为

$$\begin{cases} \delta S^\alpha = \dfrac{\delta U^\alpha + p^\alpha \delta V^\alpha - \mu^\alpha \delta n^\alpha}{T^\alpha} \\[3mm] \delta S^\beta = \dfrac{\delta U^\beta + p^\beta \delta V^\beta - \mu^\beta \delta n^\beta}{T^\beta} \end{cases} \tag{3.3.3}$$

根据熵的广延性质，整个系统的熵变是

$$\delta S = \delta S^\alpha + \delta S^\beta = \delta U^\alpha\left(\frac{1}{T^\alpha} - \frac{1}{T^\beta}\right) + \delta V^\alpha\left(\frac{p^\alpha}{T^\alpha} - \frac{p^\beta}{T^\beta}\right) - \delta n^\alpha\left(\frac{\mu^\alpha}{T^\alpha} - \frac{\mu^\beta}{T^\beta}\right) \tag{3.3.4}$$

其中应用了式 (3.3.2)。**整个系统达到平衡时，总熵有极大值**，必有

$$\delta S = 0$$

因为式 (3.3.4) 中的 $\delta U^\alpha$、$\delta V^\alpha$、$\delta n^\alpha$ 是可以独立改变的，$\delta S=0$ 要求

$$\begin{cases} \dfrac{1}{T^\alpha} - \dfrac{1}{T^\beta} = 0 \\[3mm] \dfrac{p^\alpha}{T^\alpha} - \dfrac{p^\beta}{T^\beta} = 0 \\[3mm] \dfrac{\mu^\alpha}{T^\alpha} - \dfrac{\mu^\beta}{T^\beta} = 0 \end{cases} \tag{3.3.5}$$

即

$$\begin{cases} T^\alpha = T^\beta & (\text{热平衡条件}) \\ p^\alpha = p^\beta & (\text{力学平衡条件}) \\ \mu^\alpha = \mu^\beta & (\text{相变平衡条件}) \end{cases} \tag{3.3.6}$$

式 (3.3.6) 指出，**整个系统达到平衡时，两相的温度、压强和化学势必须分别相等。这就是单元复相系达到平衡所要满足的平衡条件。**

如果平衡条件未能满足，复相系将发生变化，变化是朝着熵增加的方向进行的。

| 未满足的条件 | 变化方向 |
| --- | --- |
| 热平衡条件 | 朝 $\delta U^\alpha\left(\dfrac{1}{T^\alpha}-\dfrac{1}{T^\beta}\right)>0$ 进行。例如 $T^\alpha>T^\beta$ 时朝 $\delta U^\alpha<0$ 进行，即能量由高温相传递给低温相 |
| 力学平衡条件（热平衡已满足） | 朝 $\delta V^\alpha\left(\dfrac{p^\alpha}{T^\alpha}-\dfrac{p^\beta}{T^\beta}\right)>0$ 进行。例如 $p^\alpha>p^\beta$ 时朝 $\delta V^\alpha>0$ 进行，即压强大的相膨胀、压强小的相被压缩 |
| 相变平衡条件（热平衡已满足） | 朝 $-\delta n^\alpha\left(\dfrac{\mu^\alpha}{T^\alpha}-\dfrac{\mu^\beta}{T^\beta}\right)>0$ 进行。例如 $\mu^\alpha>\mu^\beta$ 时朝 $\delta n^\alpha<0$ 进行，即物质由化学势高的相转移到化学势低的相。**这是 $\mu$ 被称为化学势的原因** |

---

## §3.4　单元复相系的平衡性质

### 3.4.1 单元系的相图

实验指出，在不同的温度和压强范围，一个单元系可以分别处在气相、液相或固相。有些物质的固相还可以具有不同的晶格结构，不同的晶格结构也是不同的相。**用温度和压强作为直角坐标可以画出单元系的相图。**

图 3.4.1 是单元系相图的示意图。三条曲线将图分为三个区域，分别是固相、液相和气相单相存在时的温度和压强范围。在各自的区域内，温度和压强可以独立改变。

| 曲线      | 分隔的两相         | 说明                                                        |
| ------- | ------------- | --------------------------------------------------------- |
| **汽化线** | 液相 / 气相       | 其上液、气两相平衡共存；有一终点 $C$（**临界点**），温度高于 $C$ 点温度时液相即不存在，汽化线也不存在 |
| **熔化线** | 固相 / 液相       | 固、液两相平衡曲线                                                 |
| **升华线** | 固相 / 气相       | 固、气两相平衡曲线                                                 |
| **三相点** | 汽化线、熔化线、升华线交点 | 固、液、气三相平衡共存，其温度和压强是确定的                                    |

例如，水的临界温度是 647.05 K，临界压强是 22.09×10⁶ Pa；水的三相点温度为 273.16 K，压强为 610.9 Pa。

图 3.4.2 是水的相图，其中左方的图是高压下冰的相图，画出高压下八种不同的冰（没有画上气相的原因是压强的单位太大，把气相挤到 $t$ 轴去了）；右方的图用不同的压强单位画出气、液两相的相图。

**由一相到另一相的转变过程**（以液-气两相的转变为例）：如图 3.4.3 所示，设系统开始处在由点 1 所代表的气相，压强为 $p$，温度为 $T$。如果维持温度不变，缓慢地增加外界的压强，系统的体积将被压缩，压强则相应增大以维持其与外界的平衡。这样，系统的状态沿直线 1→2 变化，直到与汽化线相交于点 2，这时开始有液体凝结，并放出热量（相变潜热）。在点 2，气、液两相平衡共存。如果系统放出的热量不断被外界吸收，物质将不断地由气相转变为液相，而保持其温度和压强不变；直到系统全部转变为液相后，如果仍保持温度不变而增加外界的压强，系统的压强将相应地增大，其状态沿直线 2→3 变化。

现在根据热力学理论对单元系的相图加以解释。§3.1 中说过，**在一定的温度和压强下，系统的平衡状态是其化学势最小的状态**，各相的化学势是其温度和压强的确定的函数。如果在某一温度和压强范围内，$\alpha$ 相的化学势 $\mu^\alpha(T,p)$ 较其他相的化学势低，系统将以 $\alpha$ 相单独存在。这个温度和压强范围就是 $\alpha$ 相的**单相区域**。在这个区域内温度和压强是独立的状态参量。

### 3.4.2 两相平衡与三相平衡

单元系两相平衡共存时，必须满足 §3.3 中所讲的热平衡条件、力学平衡条件和相变平衡条件：

$$\begin{cases} T^\alpha = T^\beta = T \\ p^\alpha = p^\beta = p \\ \mu^\alpha(T,p) = \mu^\beta(T,p) \end{cases} \tag{3.4.1}$$

**式 (3.4.1) 给出两相平衡共存时压强与温度的关系，就是两相平衡曲线的方程式。在平衡曲线上，温度和压强两个参量中只有一个可以独立改变。** 由于在平衡曲线上两相的化学势相等，两相以任意比例共存，整个系统的吉布斯函数都是相等的。**这两相平衡就是 §3.1 中所说的中性平衡的例子。** 当系统缓慢地从外界吸收或放出热量时，物质将由一相转变到另一相而始终保持在平衡态，称为**平衡相变**。

单元系三相共存时，三相的温度、压强和化学势都必须相等，即

$$\begin{cases} T^\alpha = T^\beta = T^\gamma = T \\ p^\alpha = p^\beta = p^\gamma = p \\ \mu^\alpha(T,p) = \mu^\beta(T,p) = \mu^\gamma(T,p) \end{cases} \tag{3.4.2}$$

三相点的温度和压强由式 (3.4.2) 确定。

### 3.4.3 克拉珀龙方程

如果已知两相的化学势的表达式，由式 (3.4.1) 即可确定相图的两相平衡曲线。**由于缺乏化学势的全部知识，实际上相图上的平衡曲线是由实验直接测定的。** 不过，根据热力学理论可以求出两相平衡曲线的斜率。

设 $(T,p)$ 和 $(T+\mathrm{d}T, p+\mathrm{d}p)$ 是两相平衡曲线上邻近的两点，如图 3.4.4 所示。在这两点上，两相的化学势都相等：

$$\mu^\alpha(T,p) = \mu^\beta(T,p), \qquad \mu^\alpha(T+\mathrm{d}T, p+\mathrm{d}p) = \mu^\beta(T+\mathrm{d}T, p+\mathrm{d}p)$$

两式相减，得

$$\mathrm{d}\mu^\alpha = \mathrm{d}\mu^\beta \tag{3.4.3}$$

式 (3.4.3) 表示，**当沿着平衡曲线由 $(T,p)$ 变到 $(T+\mathrm{d}T, p+\mathrm{d}p)$ 时，两相的化学势的变化相等**。

化学势的全微分为

$$\mathrm{d}\mu = -S_m\,\mathrm{d}T + V_m\,\mathrm{d}p$$

其中 $S_m$ 和 $V_m$ 分别是摩尔熵和摩尔体积。代入式 (3.4.3) 得

$$-S_m^\alpha\,\mathrm{d}T + V_m^\alpha\,\mathrm{d}p = -S_m^\beta\,\mathrm{d}T + V_m^\beta\,\mathrm{d}p$$

或

$$\frac{\mathrm{d}p}{\mathrm{d}T} = \frac{S_m^\beta - S_m^\alpha}{V_m^\beta - V_m^\alpha} \tag{3.4.4}$$

以 $L$ 表示 1 mol 物质由 $\alpha$ 相转变到 $\beta$ 相时所吸收的**相变潜热**，因为相变时物质的温度不变，由式 (1.14.3) 得

$$L = T\left(S_m^\beta - S_m^\alpha\right) \tag{3.4.5}$$

代入式 (3.4.4) 得

$$\frac{\mathrm{d}p}{\mathrm{d}T} = \frac{L}{T\left(V_m^\beta - V_m^\alpha\right)} \tag{3.4.6}$$

式 (3.4.6) 称为**克拉珀龙（Clapeyron）方程**。**它给出两相平衡曲线的斜率。克拉珀龙方程与实验结果符合得很好，为热力学的正确性提供了一个直接的实验验证。**

**例：冰的熔点随压强的变化。** 在 1 atm 下冰的熔点为 273.15 K，此时冰的熔化热为 $L = 3.35\times10^5\ \mathrm{J\cdot kg^{-1}}$，冰的比体积为 $v^\alpha = 1.090\times10^{-3}\ \mathrm{m^3\cdot kg^{-1}}$，水的比体积为 $v^\beta = 1.000\,13\times10^{-3}\ \mathrm{m^3\cdot kg^{-1}}$，代入式 (3.4.6) 可算得

$$\frac{\mathrm{d}T}{\mathrm{d}p} = -\frac{273.15\ \mathrm{K}\times 0.090\,57\times10^{-3}\ \mathrm{m^3\cdot kg^{-1}}}{3.35\times10^5\ \mathrm{J\cdot kg^{-1}}} = -0.738\times10^{-7}\ \mathrm{K\cdot Pa^{-1}} = -0.007\,48\ \mathrm{K\cdot atm^{-1}}$$

这个结果与实验观测值 $\dfrac{\mathrm{d}T}{\mathrm{d}p} = -0.0075\ \mathrm{K\cdot atm^{-1}}$ 相符合。

**例：水的沸点随压强的变化。** 在 1 atm 下水的沸点为 373.15 K。此时水的汽化热为 $L = 2.257\times10^6\ \mathrm{J\cdot kg^{-1}}$，水的比体积为 $v^\alpha = 1.043\times10^{-3}\ \mathrm{m^3\cdot kg^{-1}}$，水蒸气的比体积为 $v^\beta = 1\,673\times10^{-3}\ \mathrm{m^3\cdot kg^{-1}}$。代入式 (3.4.6) 可算得

$$\frac{\mathrm{d}p}{\mathrm{d}T} = \frac{2.257\times10^6\ \mathrm{J\cdot kg^{-1}}}{373.15\ \mathrm{K}\times 1\,672\times10^{-3}\ \mathrm{m^3\cdot kg^{-1}}} = 3.62\times10^3\ \mathrm{Pa\cdot K^{-1}} = 0.0357\ \mathrm{atm\cdot K^{-1}}$$

这个结果与实验观测值 $\dfrac{\mathrm{d}p}{\mathrm{d}T} = 0.0356\ \mathrm{atm\cdot K^{-1}}$ 相符合。

**平衡曲线的斜率**：当物质发生熔化、蒸发或升华时，通常比体积增大，且相变潜热是正的（混乱度增加，因而比熵增加）。因此平衡曲线的斜率 $\dfrac{\mathrm{d}p}{\mathrm{d}T}$ 通常是正的。不过在某些情形下，熔化线具有负的斜率。例如冰熔化时比体积变小，因而熔化线的斜率 $\dfrac{\mathrm{d}p}{\mathrm{d}T}$ 是负的。$^3$He 熔化时比体积增大，但在 0.3 K 以下固相的比熵大于液相，也使熔化线具有负的斜率，如图 3.4.5 所示。

### 3.4.4 蒸气压方程

由克拉珀龙方程可以推导**蒸气压方程**。与凝聚相（液相或固相）达到平衡的蒸气称为**饱和蒸气**。由于两相平衡时压强与温度间存在一定的关系，饱和蒸气的压强是温度的函数。**描述饱和蒸气压与温度的关系的方程称为蒸气压方程。**

以 $\alpha$ 相表示凝聚相，$\beta$ 相表示气相。凝聚相的摩尔体积远小于气相的摩尔体积。如果在式 (3.4.6) 中略去 $V_m^\alpha$，并把气相看作理想气体，$pV_m^\beta = RT$，式 (3.4.6) 可简化为

$$\frac{1}{p}\frac{\mathrm{d}p}{\mathrm{d}T} = \frac{L}{RT^2} \tag{3.4.7}$$

如果更进一步近似地认为相变潜热与温度无关（这个近似是十分粗糙的），就可以将上式积分得

$$\ln p = -\frac{L}{RT} + A \tag{3.4.8}$$

式 (3.4.8) 是**蒸气压方程的近似表达式**。关于蒸气压方程的进一步讨论可参阅本章习题。

---

## §3.5　临界点和气液两相的转变

### 3.5.1 等温线与杠杆规则

§3.4 中用温度和压强为坐标画出了单元系的相图。图 3.4.1 中的汽化线是液、气两相的平衡曲线，汽化线终止于临界点。本节中我们再用 $p$-$V$ 图的等温线分析液、气两相的转变，可以更清楚地显示出其中的某些特性。

图 3.5.1 是安德鲁斯（Andrews）于 1869 年得到的二氧化碳在高温下的等温线。在临界温度 31.1 °C 以上，等温线的形状与玻意耳定律给出的双曲线近似，是气相的等温线。在临界温度以下，等温线包括三段：

- 左边的一段几乎与 $p$ 轴平行（其等温压缩系数很小），代表液相；
- 右边的一段代表气相；
- 中间的一段是与 $v$ 轴平行的直线。§3.4 中说过，两相共存时，在一定的温度下压强是一定的，因此这一段代表液、气共存的状态。

对于单位质量的物质，这段直线左端的横坐标就是液相的比体积 $v_l$，右端的横坐标是气相的比体积 $v_g$。直线中比体积为 $v$ 的一点，相应的液相比例 $x$ 和气相比例 $1-x$ 由下式给出：

$$v = x v_l + (1-x) v_g \tag{3.5.1}$$

**等温线中的水平段随温度的升高而缩短**，说明液、气两相的比体积随温度升高而接近。当温度达到某一极限温度时，水平段的左右端重合。这时两相的比体积相等，两相的其他差别也不再存在，物质处在液气不分的状态。这一极限温度就是**临界温度** $T_c$，相应的压强是**临界压强** $p_c$。在温度为 $T_c$ 的等温线上，压强小于 $p_c$ 时物质处在气相；压强等于或高于 $p_c$ 时物质处在液气不分的状态。**当温度高于 $T_c$ 时，无论处在多大的压强下，物质都处于气态，液态不可能存在。** 由于有了临界点，我们可以像图 3.4.3 中的 4→5 线那样绕过临界点，使气相连续地转变为液相而不必经过气、液两相共存的阶段。

由上面的讨论可知，**临界等温线在临界点的切线是水平的**，即 $\left(\dfrac{\partial p}{\partial V}\right)_T = 0$。

### 3.5.2 范德瓦耳斯等温线与麦克斯韦等面积法则

范德瓦耳斯在 1873 年用他的方程统一地描述液、气两态，并据此讨论了液、气两相转变和临界点问题。对于 1 mol 物质，范德瓦耳斯方程是

$$\left(p + \frac{a}{V_m^2}\right)\left(V_m - b\right) = RT \tag{3.5.2}$$

图 3.5.2 画出了范德瓦耳斯方程的等温线。比较图 3.5.1 和图 3.5.2 可以看出，范德瓦耳斯气体的等温线与实际观测到的等温线很像。不过在温度低于 $T_c$ 时，范德瓦耳斯气体的等温线在 $p_1<p<p_2$ 的范围，对应于一个 $p$ 值有三个可能的 $V_m$ 值，如图 3.5.3 所示。当然，在 $V_{m1}<V_m<V_{m2}$ 的范围内，$\left(\dfrac{\partial p}{\partial V_m}\right)_T>0$。**由于不满足平衡稳定性条件的要求，这些状态是不可能作为均匀系而实现的。** 现在我们根据吉布斯函数最小的要求，讨论在 $p_1<p<p_2$ 的范围内，在给定的 $T$、$p$ 下，什么状态是稳定的平衡状态。

化学势的全微分是

$$\mathrm{d}\mu = -S_m\,\mathrm{d}T + V_m\,\mathrm{d}p \tag{3.5.3}$$

由此可知，等温线上压强为 $p$ 与压强为 $p_0$ 的两个状态的化学势之差为

$$\mu - \mu_0 = \int_{p_0}^{p} V_m\,\mathrm{d}p \tag{3.5.4}$$

积分 (3.5.4) 等于等温线与 $p$ 轴之间由 $p_0$ 到 $p$ 的面积。由图 3.5.3 可以看出，当式 (3.5.4) 的积分下限固定为 $Q$ 点的压强 $p_0$ 而沿等温线积分时，积分的数值由 $Q$ 点出发后增加，到 $N$ 点后减少，到 $J$ 点后又再次增加。因此在温度保持为恒定时，$\mu$ 随 $p$ 的改变如图 3.5.4 所示。

图 3.5.4 中 $Q$、$K$、$B$、$N$、$J$、$A$、$M$、$R$ 各点分别与图 3.5.3 中各点相应，其中 $A$、$B$ 两点在图 3.5.4 中重合，说明 $A$、$B$ 两态的 $\mu$ 值相等。从图 3.5.4 可以看出，在 $p_1<p<p_2$ 的范围内，对应于一个 $p$ 值，$\mu$ 有三个可能的值。这与图 3.5.3 中在 $p_1<p<p_2$ 的范围内，对应于一个 $p$ 值有三个可能的 $V_m$ 值是相应的。**根据吉布斯函数判据，在给定的 $T$、$p$ 下，稳定平衡态的吉布斯函数最小**，因此线段 $QKBAMR$ 上各点代表系统的稳定平衡状态。

在 $B$ 点物质全部处在气态，在 $A$ 点物质全部处在液态。$B$ 点和 $A$ 点的 $\mu$ 值相等，这正是在等温线的温度和 $A$、$B$ 两点的压强下气、液两相的相变平衡条件。$A$、$B$ 两态既然满足条件

$$\mu_A = \mu_B$$

由图 3.5.3 可以看出，这相当于积分

$$\int_{BNDJA} V_m\,\mathrm{d}p = 0$$

或

$$\text{面积}(BND) = \text{面积}(DJA) \tag{3.5.5}$$

这就是说，$A$、$B$ 两点在图 3.5.3 中的位置可以由条件 (3.5.5) 确定，**称为麦克斯韦等面积法则**。根据等面积法则，将范德瓦耳斯气体等温线中的 $BNDJA$ 段换为直线 $BA$ 就与图 3.5.1 中的实测等温线相符了。

如前所述，线段 $JDN$ 上的状态不满足平衡稳定性的要求，物质不可能作为均匀系存在而必将发生相变；线段 $BN$ 和 $AJ$ 上的状态满足平衡稳定性的要求，由于其化学势高于两相平衡的化学势，它们可以作为**亚稳态**单相存在，分别对应于**过饱和蒸气**和**过热液体**。我们将在 §3.6 中讨论在什么条件下可能出现过饱和蒸气和过热液体。

### 3.5.3 临界点与临界系数

在等温线的极大点 $N$，有 $\left(\dfrac{\partial p}{\partial V_m}\right)_T = 0$，$\left(\dfrac{\partial^2 p}{\partial V_m^2}\right)_T < 0$；在极小点 $J$，有 $\left(\dfrac{\partial p}{\partial V_m}\right)_T = 0$，$\left(\dfrac{\partial^2 p}{\partial V_m^2}\right)_T > 0$。随着温度的升高，极大点与极小点逐渐靠近，达到临界温度 $T_c$ 时，这两点重合，而形成**拐点**。因此临界点的温度 $T_c$ 和压强 $p_c$ 满足方程：

$$\left(\frac{\partial p}{\partial V_m}\right)_T = 0, \qquad \left(\frac{\partial^2 p}{\partial V_m^2}\right)_T = 0 \tag{3.5.6}$$

将范德瓦耳斯方程代入，可得

$$\left(\frac{\partial p}{\partial V_m}\right)_T = -\frac{RT}{(V_m-b)^2} + \frac{2a}{V_m^3} = 0$$

$$\left(\frac{\partial^2 p}{\partial V_m^2}\right)_T = \frac{2RT}{(V_m-b)^3} - \frac{6a}{V_m^4} = 0$$

再加上范德瓦耳斯方程本身共三个方程，可以把临界点的温度 $T_c$、压强 $p_c$ 和体积 $V_{mc}$ 定出来。结果为

$$T_c = \frac{8a}{27Rb}, \qquad p_c = \frac{a}{27b^2}, \qquad V_{mc} = 3b$$

$T_c$、$p_c$ 和 $V_{mc}$ 之间存在以下关系：

$$\frac{RT_c}{p_c V_{mc}} = \frac{8}{3} = 2.667 \tag{3.5.7}$$

这个量纲一的比值叫做**临界系数**。**根据范德瓦耳斯方程，临界系数对各种气（液）体应相同。** 实测结果中，临界系数的数值如下：

| 气体 | He | H₂ | Ne | Ar | O₂ | CO₂ | H₂O |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $RT_c/(p_cV_{mc})$ | 3.28 | 3.27 | 3.43 | 3.42 | 3.42 | 3.65 | 4.37 |

### 3.5.4 对应态定律

引进新的变量：

$$t^* = \frac{T}{T_c}, \qquad p^* = \frac{p}{p_c}, \qquad v^* = \frac{V_m}{V_{mc}}$$

分别称为**对比温度**、**对比压强**和**对比体积**，可将范德瓦耳斯方程化为

$$\left(p^* + \frac{3}{v^{*2}}\right)\left(v^* - \frac{1}{3}\right) = \frac{8}{3}t^* \tag{3.5.8}$$

上式称为**范德瓦耳斯对比方程**。在这个方程中，不含与具体物质性质有关的常量。也就是说，**如果采用对比变量，范德瓦耳斯方程是普适的。这个结果称为对应态定律。**

### 3.5.5 临界态的平衡稳定性

前面根据范德瓦耳斯方程讨论了液气两相及其临界点，得到了临界点所满足的方程：

$$\left(\frac{\partial p}{\partial V_m}\right)_T = 0, \qquad \left(\frac{\partial^2 p}{\partial V_m^2}\right)_T = 0$$

实际上，根据热力学的一般论据就可以证明，液气流体系统**临界态的平衡稳定条件**（与物态方程的具体表达式无关）为

$$\left(\frac{\partial p}{\partial V_m}\right)_T = 0, \qquad \left(\frac{\partial^2 p}{\partial V_m^2}\right)_T = 0, \qquad \left(\frac{\partial^3 p}{\partial V_m^3}\right)_T < 0 \tag{3.5.9}$$

如前所述，液气两相平衡时，两相具有相同的温度和压强。在接近临界点的等温线上，两相平衡时其摩尔体积是接近的。以 $V_m$ 和 $V_m+\Delta V_m$ 分别表示液相和气相的摩尔体积，两相压强相等可以表示为

$$p(V_m, T) = p(V_m + \Delta V_m, T) \tag{3.5.10}$$

将上式右方按 $\Delta V_m$ 作泰勒展开，准确到二级，有

$$p(V_m+\Delta V_m, T) = p(V_m, T) + \left(\frac{\partial p}{\partial V_m}\right)_T \Delta V_m + \frac{1}{2}\left(\frac{\partial^2 p}{\partial V_m^2}\right)_T (\Delta V_m)^2$$

代入式 (3.5.10)，全式除以 $\Delta V_m$，可得

$$\left(\frac{\partial p}{\partial V_m}\right)_T + \frac{\Delta V_m}{2}\left(\frac{\partial^2 p}{\partial V_m^2}\right)_T = 0 \tag{3.5.11}$$

当 $T$ 趋于临界温度 $T_c$ 时，$\Delta V_m$ 趋于零，由此可知，在临界点有

$$\left(\frac{\partial p}{\partial V_m}\right)_T = 0 \tag{3.5.12}$$

上式指出，**处在临界点的液气流体系统不满足平衡稳定条件 (3.1.14) 的第二式**。将上式代入式 (3.1.13) 可知，如果将系统约束在临界温度 $T_c$ 上（即 $\delta T=0$），熵函数的二级微分为零，即 $\delta^2 S_m=0$。

为了分析液气流体系统临界态的平衡稳定性，将 $\Delta S_m$ 作泰勒展开，得

$$\Delta S_m = \delta S_m + \frac{1}{2!}\delta^2 S_m + \frac{1}{3!}\delta^3 S_m + \frac{1}{4!}\delta^4 S_m + \cdots$$

根据式 (3.1.13)，有

$$\delta^2 S_m = -\frac{C_{V,m}}{T^2}(\delta T)^2 + \frac{1}{T}\left(\frac{\partial p}{\partial V_m}\right)_T (\delta V_m)^2$$

由

$$\delta^3 S_m = \left(\frac{\partial}{\partial T}\delta^2 S_m\right)\delta T + \left(\frac{\partial}{\partial V_m}\delta^2 S_m\right)\delta V_m$$

可知，将系统约束在临界温度 $T_c$ 的情形下有

$$\delta^3 S_m = \frac{1}{T}\left(\frac{\partial^2 p}{\partial V_m^2}\right)_T (\delta V_m)^3 \tag{3.5.13}$$

如果 $\left(\dfrac{\partial^2 p}{\partial V_m^2}\right)_T \neq 0$，系统的平衡稳定性将取决于在虚变动中体积是膨胀还是缩小，这显然是不合理的。由此可以推断

$$\left(\frac{\partial^2 p}{\partial V_m^2}\right)_T = 0 \tag{3.5.14}$$

在 $\delta^2 S_m=0$ 和 $\delta^3 S_m=0$ 的情形下，我们考察 $\delta^4 S_m$。类似地约束在临界温度 $T_c$ 上：

$$\delta^4 S_m = \frac{1}{T}\left(\frac{\partial^3 p}{\partial V_m^3}\right)_T (\delta V_m)^4$$

平衡稳定性要求 $\delta^4 S_m<0$，由此可知

$$\left(\frac{\partial^3 p}{\partial V_m^3}\right)_T < 0 \tag{3.5.15}$$

这就证明了式 (3.5.9)。注意，式 (3.5.12) 是由临界态的特性决定的，它不满足平衡稳定条件的第二式；式 (3.5.14) 和式 (3.5.15) 则是临界态具有平衡稳定性的要求。

---

## *§3.6　液滴的形成

> 带星号的选读内容。

### 3.6.1 表面相与三相的平衡条件

前面讨论两相平衡时没有考虑表面相的影响，适用于分界面为平面的情形。本节以液滴在蒸气中的形成为例，讨论表面相对相变过程的影响。

我们先讨论在考虑表面相以后系统在达到平衡时所要满足的平衡条件。设液滴为 $\alpha$ 相，蒸气为 $\beta$ 相，表面为 $\gamma$ 相。根据 §3.2 和 §2.5 中的讨论，三相的热力学基本方程分别为

$$\begin{cases} \mathrm{d}U^\alpha = T^\alpha\,\mathrm{d}S^\alpha - p^\alpha\,\mathrm{d}V^\alpha + \mu^\alpha\,\mathrm{d}n^\alpha \\ \mathrm{d}U^\beta = T^\beta\,\mathrm{d}S^\beta - p^\beta\,\mathrm{d}V^\beta + \mu^\beta\,\mathrm{d}n^\beta \\ \mathrm{d}U^\gamma = T^\gamma\,\mathrm{d}S^\gamma + \sigma\,\mathrm{d}A \end{cases} \tag{3.6.1}$$

在热力学中我们把表面理想化为几何面，因此表面相的物质的量 $n^\gamma=0$，在基本方程中不含 $\mathrm{d}n^\gamma$ 的项。

系统的热平衡条件为三相的温度相等，即

$$T^\alpha = T^\beta = T^\gamma \tag{3.6.2}$$

假定热平衡条件已经满足，温度保持不变，我们用自由能判据推求系统的力学平衡条件和相变平衡条件。设想在温度和总体积保持不变的条件下，系统发生一个虚变动。在这个虚变动中，三相的物质的量、体积（或面积）分别有 $\delta n^\alpha$、$\delta V^\alpha$；$\delta n^\beta$、$\delta V^\beta$；$\delta A$ 的变化。由于在虚变动中系统的总物质的量和总体积保持不变，应有

$$\begin{cases} \delta n^\alpha + \delta n^\beta = 0 \\ \delta V^\alpha + \delta V^\beta = 0 \end{cases} \tag{3.6.3}$$

在这个虚变动中，三相自由能的变化分别为

$$\begin{cases} \delta F^\alpha = -p^\alpha\,\delta V^\alpha + \mu^\alpha\,\delta n^\alpha \\ \delta F^\beta = -p^\beta\,\delta V^\beta + \mu^\beta\,\delta n^\beta \\ \delta F^\gamma = \sigma\,\delta A \end{cases} \tag{3.6.4}$$

在三相温度相等的情形下，整个系统的自由能是三相的自由能之和，因此整个系统的自由能变化是

$$\delta F = \delta F^\alpha + \delta F^\beta + \delta F^\gamma = -\left(p^\alpha - p^\beta\right)\delta V^\alpha + \sigma\,\delta A + \left(\mu^\alpha - \mu^\beta\right)\delta n^\alpha \tag{3.6.5}$$

其中用了式 (3.6.3)。如果假设液滴是球形的，半径为 $r$，即有

$$V^\alpha = \frac{4\pi}{3}r^3, \qquad A = 4\pi r^2$$

$$\delta V^\alpha = 4\pi r^2\,\delta r, \qquad \delta A = 8\pi r\,\delta r$$

于是式 (3.6.5) 可简化为

$$\delta F = -\left(p^\alpha - p^\beta - \frac{2\sigma}{r}\right)\delta V^\alpha + \left(\mu^\alpha - \mu^\beta\right)\delta n^\alpha$$

根据自由能判据，在温度和总体积不变的条件下，平衡态的自由能最小，必有 $\delta F=0$。因为 $\delta V^\alpha$ 和 $\delta n^\alpha$ 是独立的，所以有

$$p^\alpha = p^\beta + \frac{2\sigma}{r} \tag{3.6.6}$$

$$\mu^\alpha = \mu^\beta \tag{3.6.7}$$

- 式 (3.6.6) 是**力学平衡条件**。它指出，**由于表面张力有使液滴收缩的趋势，液滴的压强必须大于蒸气的压强才能维持力学平衡。** 当 $r\to\infty$（相当于分界面为平面）时，式 (3.6.6) 给出 $p^\alpha=p^\beta$。这就是说，**当分界面为平面时，力学平衡条件是两相的压强相等。**
- 式 (3.6.7) 指出，**相变平衡条件仍然是两相的化学势相等**。但是必须注意，在式 (3.6.7) 两方的化学势中，压强 $p^\alpha$ 和 $p^\beta$ 满足式 (3.6.6)，其数值是不同的。假如式 (3.6.7) 不能满足，物质将由化学势高的相转变到化学势低的相。

### 3.6.2 曲界面上的平衡蒸气压

现在讨论液滴的形成问题。我们先讨论气液两相平衡时分界面为曲面的蒸气压强与分界面为平面的饱和蒸气压的关系。

当液面为平面时，力学平衡条件是两相的压强相等。以 $p$ 表示这时两相的压强，相变平衡条件为

$$\mu^\alpha(T,p) = \mu^\beta(T,p) \tag{3.6.8}$$

上式确定了**饱和蒸气压与温度的关系**。

在液面为曲面的情形下，设气液两相平衡时蒸气的压强为 $p'$。由式 (3.6.6) 知，这时液滴的压强为 $p' + \dfrac{2\sigma}{r}$。相变平衡条件式 (3.6.7) 应为

$$\mu^\alpha\left(p' + \frac{2\sigma}{r}, T\right) = \mu^\beta(p', T) \tag{3.6.9}$$

上式给出分界面为曲面时的平衡蒸气压强 $p'$ 与温度 $T$ 及曲面半径 $r$ 的关系。

现在讨论 $p'$ 与 $p$ 的关系。当压强改变时，液体的性质改变很小，我们可以将液滴的化学势 $\mu^\alpha$ 按压强展开，只取线性项，有

$$\begin{aligned}\mu^\alpha\left(p' + \frac{2\sigma}{r}, T\right) &= \mu^\alpha(p,T) + \left(p' - p + \frac{2\sigma}{r}\right)\frac{\partial \mu^\alpha}{\partial p} \\ &= \mu^\alpha(p,T) + \left(p' - p + \frac{2\sigma}{r}\right)v^\alpha \end{aligned} \tag{3.6.10}$$

如果把蒸气看成理想气体，根据式 (2.4.15)，蒸气的化学势为

$$\mu^\beta(p,T) = RT(\varphi + \ln p) \tag{3.6.11}$$

其中 $\varphi$ 是温度的函数。由式 (3.6.11) 得

$$\mu^\beta(p',T) = \mu^\beta(p,T) + RT\ln\frac{p'}{p} \tag{3.6.12}$$

将式 (3.6.10) 和式 (3.6.12) 代入式 (3.6.9)，考虑到式 (3.6.8) 可得

$$\left(p' - p + \frac{2\sigma}{r}\right)v^\alpha = RT\ln\frac{p'}{p} \tag{3.6.13}$$

在实际问题中，通常有 $p'-p \ll \dfrac{2\sigma}{r}$。在这种情形下，式 (3.6.13) 可近似为

$$\ln\frac{p'}{p} = \frac{2\sigma v^\alpha}{RTr} \tag{3.6.14}$$

**数值估算（水滴）**：在 $T=291\ \mathrm{K}$ 时，水的表面张力系数 $\sigma = 0.073\ \mathrm{N\cdot m^{-1}}$，$v^\alpha = 18.016\times10^{-6}\ \mathrm{m^3\cdot mol^{-1}}$，代入式 (3.6.14) 得

$$\ln\frac{p'}{p} = \frac{1.087\times10^{-9}}{r}$$

或

$$\lg\frac{p'}{p} = \frac{4.72\times10^{-10}}{r}$$

当 $r=10^{-7}\ \mathrm{m}$ 时，$\dfrac{p'}{p}=1.011$，但 $p = 2.042\times10^{-3}\ \mathrm{Pa}$，可见 $p'-p \ll \dfrac{2\sigma}{r}$，说明略去式 (3.6.13) 中的 $p'-p$ 是容许的。

> **⚠ 教材笔误提示**：教材 p.79 此处印为 $p=2.042\times10^{-3}$ Pa。291 K（18 °C）时水的饱和蒸气压约为 $2.042\times10^{3}$ Pa（2.042 kPa），故该处指数应为笔误。无论取何值，$p'-p\approx0.011p$ 都远小于 $\dfrac{2\sigma}{r}=1.46\times10^{6}\ \mathrm{Pa}$，结论不受影响。

根据式 (3.6.14)，可计算在各种不同半径下水滴与蒸气达到平衡时所需的蒸气压。当 $r=10^{-7}\ \mathrm{m}$ 时，$\dfrac{p'}{p}=1.011$；当 $r=10^{-8}\ \mathrm{m}$ 时，$\dfrac{p'}{p}=1.115$；当 $r=10^{-9}\ \mathrm{m}$ 时，$\dfrac{p'}{p}=2.966$。由此可知，**水滴越小，与水滴达到平衡所需的蒸气压强就越高**。

### 3.6.3 临界半径与过饱和蒸气

在一定的温度和蒸气压强 $p'$ 下，与蒸气达到平衡的液滴半径 $r_c$ 为

$$r_c = \frac{2\sigma v^\alpha}{RT\ln\frac{p'}{p}} \tag{3.6.15}$$

$r_c$ 称为**临界半径**。由式 (3.6.10) 可以看出：

- 对于 $r>r_c$ 的液滴，有 $\mu^\alpha<\mu^\beta$，因而液滴将继续凝结而增大；
- 对于 $r<r_c$ 的液滴，有 $\mu^\alpha>\mu^\beta$，因而液滴将汽化而消失。

在蒸气中，液体的凝结是通过先形成微小液滴然后逐渐生长的方式发生的。**如果在蒸气中不存在凝结核（例如灰尘或带电微粒等），由涨落而形成的液滴往往过小，不能增大，因此在非常干净的蒸气中，蒸气的压强可以超过饱和蒸气压而不凝结，形成过饱和蒸气。**

### 3.6.4 气泡与过热液体

液体中的气泡可以同样考虑。如果仍然令 $\alpha$ 相表示液相，$\beta$ 相表示气相，则在式 (3.6.6) 和式 (3.6.14) 中要将 $r$ 换成 $-r$。将 $r$ 换成 $-r$ 后，由式 (3.6.6) 得

$$p^\beta = p^\alpha + \frac{2\sigma}{r} \tag{3.6.16}$$

式 (3.6.16) 指出，气泡内蒸气的压强必须大于液体的压强才能维持力学平衡。将 $r$ 换成 $-r$ 后，由式 (3.6.14) 得

$$\ln\frac{p}{p'} = \frac{2\sigma v^\alpha}{RTr} \tag{3.6.17}$$

式 (3.6.17) 指出，为满足相变平衡条件，气泡内的压强必须小于同温度的饱和蒸气压。

根据式 (3.6.16) 和式 (3.6.17) 可以说明液体沸腾前的**过热现象**。液体沸腾时，液体内部有大量的气泡形成，使气液分界面大大增加，于是整个液体剧烈汽化。在一般情形下，液体中溶有空气。以这些已有的空气泡作核而形成的气泡具有足够大的半径，接近于分界面为平面的情形，只要气泡中的蒸气压等于液体的压强，即发生沸腾。**如果液体中没有现存的空气泡作核**，由涨落而形成的气泡半径很小，当达到正常沸点的温度即饱和蒸气压 $p$ 等于液体的压强 $p^\alpha$ 时，力学平衡条件 (3.6.16) 要求气泡内的蒸气压强 $p^\beta$ 大于液体的压强即大于饱和蒸气压 $p$，而相变平衡条件 (3.6.17) 又要求气泡内的蒸气压强 $p'$ 小于饱和蒸气压 $p$。因此在正常的沸点温度，式 (3.6.16) 和式 (3.6.17) 不可能同时满足。除非液体的温度高于正常的沸点，使相应的饱和蒸气压 $p$ 大于液体的压强 $p^\alpha$，式 (3.6.16) 和式 (3.6.17) 才可能同时满足。这就是形成**过热液体**的原因。

通过以上的讨论可以知道，**在新相生成时，表面相起着重要的作用。**

---

## §3.7　相变的分类

### 3.7.1 一级相变

如前所述，固相、液相和气相之间的转变（包括固相不同晶格结构之间的转变）存在相变潜热 $L = T\left(s^{(2)}-s^{(1)}\right)$ 和体积突变 $\Delta v = v^{(2)}-v^{(1)}$，而且可能出现亚稳态。自然界还存在另一类相变，在转变时既无潜热又无体积突变。例如液-气通过临界点的转变、铁磁顺磁的转变、合金的有序无序转变、液 He Ⅰ 和 He Ⅱ 的转变、零磁场下金属超导状态和正常状态的转变，等等。

1933 年，埃伦菲斯特（Ehrenfest）试图对相变进行分类。因为 $s = -\dfrac{\partial \mu}{\partial T}$，$v = \dfrac{\partial \mu}{\partial p}$，埃伦菲斯特将前述第一类相变概括为：**在相变点，两相的化学势连续，但化学势的一级偏导数存在突变**，即

$$\begin{cases} \mu^{(1)}(T,p) = \mu^{(2)}(T,p) \\[2mm] \dfrac{\partial \mu^{(1)}}{\partial T} \neq \dfrac{\partial \mu^{(2)}}{\partial T}, \qquad \dfrac{\partial \mu^{(1)}}{\partial p} \neq \dfrac{\partial \mu^{(2)}}{\partial p} \end{cases} \tag{3.7.1}$$

称这一类相变为**一级相变**。图 3.7.1 形象地表达了上述特征和亚稳态存在的可能性。概括地说，**在一级相变中，两相有其各自的非奇异的化学势函数，相变点是两相化学势函数的交点**。

在一级相变中，两相的化学势相等，两相可以平衡共存；但是两相化学势的一级导数不等，转变时有潜热和比体积突变。**在相变点的两侧，化学势较低的相是稳定相，化学势较高的相可以作为亚稳态存在。**

### 3.7.2 二级相变与埃伦菲斯特方程

根据埃伦菲斯特的分类，**如果在相变点两相的化学势和化学势的一级偏导数连续，但化学势的二级偏导数存在突变，称为二级相变**。因为

$$\begin{cases} c_p = T\left(\dfrac{\partial s}{\partial T}\right)_p = -T\dfrac{\partial^2 \mu}{\partial T^2} \\[3mm] \alpha = \dfrac{1}{v}\left(\dfrac{\partial v}{\partial T}\right)_p = \dfrac{1}{v}\dfrac{\partial^2 \mu}{\partial T\partial p} \\[3mm] \kappa_T = -\dfrac{1}{v}\left(\dfrac{\partial v}{\partial p}\right)_T = -\dfrac{1}{v}\dfrac{\partial^2 \mu}{\partial p^2} \end{cases} \tag{3.7.2}$$

所以**二级相变没有相变潜热和比体积突变，但是定压比热、定压膨胀系数和等温压缩系数存在突变**。埃伦菲斯特根据二级相变在邻近的相变点 $(T,p)$ 和 $(T+\mathrm{d}T, p+\mathrm{d}p)$ 两相的比熵和比体积变化相等，即 $\mathrm{d}s^{(1)}=\mathrm{d}s^{(2)}$ 和 $\mathrm{d}v^{(1)}=\mathrm{d}v^{(2)}$ 的条件，导出了二级相变压强随温度变化的斜率公式：

$$\begin{cases} \dfrac{\mathrm{d}p}{\mathrm{d}T} = \dfrac{\alpha^{(2)} - \alpha^{(1)}}{\kappa_T^{(2)} - \kappa_T^{(1)}} \\[4mm] \dfrac{\mathrm{d}p}{\mathrm{d}T} = \dfrac{c_p^{(2)} - c_p^{(1)}}{Tv\left(\alpha^{(2)} - \alpha^{(1)}\right)} \end{cases} \tag{3.7.3}$$

称为**埃伦菲斯特方程**。

根据埃伦菲斯特的分类，如果在相变点两相的化学势和化学势的一级、二级……直到 $n-1$ 级的偏导数连续，但化学势的 $n$ 级偏导数存在突变，则称为 **$n$ 级相变**。

**补充证明

在二级相变中，两相的化学势及其一级偏导数连续，所以在相变点两相的比熵和比体积分别相等：

$$
s^{(1)}=s^{(2)},\qquad v^{(1)}=v^{(2)}.
$$

沿相平衡曲线取邻近的两个相变点 $(T,p)$ 与 $(T+\mathrm{d}T,p+\mathrm{d}p)$。两相在第二个点仍满足相同的相平衡条件，而 $v$、$s$ 在相变点连续，因此一、二两相从第一个相变点到第二个相变点时，比体积和比熵的变化分别相等：

$$
\mathrm{d}v^{(1)}=\mathrm{d}v^{(2)},\qquad
\mathrm{d}s^{(1)}=\mathrm{d}s^{(2)}.
$$

**先导出第一个斜率公式。** 由

$$
\alpha=\frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_p,
\qquad
\kappa_T=-\frac{1}{v}\left(\frac{\partial v}{\partial p}\right)_T
$$

可得两相的比体积全微分：

$$
\begin{aligned}
\mathrm{d}v^{(i)}
&=\left(\frac{\partial v^{(i)}}{\partial T}\right)_p\mathrm{d}T
+\left(\frac{\partial v^{(i)}}{\partial p}\right)_T\mathrm{d}p\\
&=\alpha^{(i)}v^{(i)}\,\mathrm{d}T
-\kappa_T^{(i)}v^{(i)}\,\mathrm{d}p,
\qquad i=1,2.
\end{aligned}
$$

在相变点 $v^{(1)}=v^{(2)}=v$，代入 $\mathrm{d}v^{(1)}=\mathrm{d}v^{(2)}$ 并约去公共因子 $v$，得到

$$
\alpha^{(1)}\mathrm{d}T-\kappa_T^{(1)}\mathrm{d}p
=
\alpha^{(2)}\mathrm{d}T-\kappa_T^{(2)}\mathrm{d}p.
$$

移项即得

$$
\frac{\mathrm{d}p}{\mathrm{d}T}
=
\frac{\alpha^{(2)}-\alpha^{(1)}}
{\kappa_T^{(2)}-\kappa_T^{(1)}}.
$$

**再导出第二个斜率公式。** 由

$$
c_p=T\left(\frac{\partial s}{\partial T}\right)_p,
\qquad
\left(\frac{\partial s}{\partial p}\right)_T
=-\left(\frac{\partial v}{\partial T}\right)_p
=-\alpha v
$$

可得两相的比熵全微分：

$$
\begin{aligned}
\mathrm{d}s^{(i)}
&=\left(\frac{\partial s^{(i)}}{\partial T}\right)_p\mathrm{d}T
+\left(\frac{\partial s^{(i)}}{\partial p}\right)_T\mathrm{d}p\\
&=\frac{c_p^{(i)}}{T}\mathrm{d}T
-\alpha^{(i)}v^{(i)}\,\mathrm{d}p,
\qquad i=1,2.
\end{aligned}
$$

代入 $\mathrm{d}s^{(1)}=\mathrm{d}s^{(2)}$，并再次利用相变点处 $v^{(1)}=v^{(2)}=v$，有

$$
\frac{c_p^{(1)}}{T}\mathrm{d}T-v\alpha^{(1)}\mathrm{d}p
=
\frac{c_p^{(2)}}{T}\mathrm{d}T-v\alpha^{(2)}\mathrm{d}p.
$$

移项后得到

$$
\frac{\mathrm{d}p}{\mathrm{d}T}
=
\frac{c_p^{(2)}-c_p^{(1)}}
{Tv\left(\alpha^{(2)}-\alpha^{(1)}\right)}.
$$

以上两个斜率公式合起来就是式 (3.7.3)，即**埃伦菲斯特方程**。证明只用到了二级相变中 $\mu$ 及一阶偏导数的连续性、相变点处 $v^{(1)}=v^{(2)}$，以及沿相平衡曲线 $\mathrm{d}s^{(1)}=\mathrm{d}s^{(2)}$、$\mathrm{d}v^{(1)}=\mathrm{d}v^{(2)}$。

### 3.7.3 连续相变

埃伦菲斯特的分类适用于突变为有限的情形。后来发现，在前面提到的第二类相变中，**热容、等温压缩系数、磁化率等在趋近相变点时往往趋于无穷**。现在人们习惯上只把相变区分为**一级相变**和**连续相变**两类，把非一级相变统称为**连续相变**。

在连续相变的相变点，**两相的化学势和化学势的一级偏导数连续。因此在温度、压强为 $T$、$p$ 的相变点，两相不仅 $\mu$、$s$、$v$ 相等，而且 $u$、$h$、$f$ 也相等**。这就是说，在相变点两相的热力学状态是相同的，**不存在两相共存现象**。但是**相变点是热力学函数的奇异性点**，与化学势二级偏导数相应的热容、等温压缩系数等在相变点出现跃变或无穷尖峰。同时，在连续相变相变点的两侧，各存在热力学函数的一个分支与物质系统的一个相对应，**不存在亚稳态现象**。

图 3.7.2 就热容在相变点出现无穷尖峰的情形示意地画出了 $\mu$、$s = -\left(\dfrac{\partial \mu}{\partial T}\right)_p$ 和 $c_p = T\left(\dfrac{\partial s}{\partial T}\right)_p$ 随温度变化的典型特征。

**表 3.7.1　连续相变的例子**

| 相变 | 序参量 | 例子 | $T_c/\mathrm{K}$ |
| --- | --- | --- | --- |
| 液气 | $\rho_l - \rho_g$ | $\mathrm{H_2O}$ | 647.05 |
| 铁磁 | 磁化强度 | Fe | 1 044.0 |
| 反铁磁 | 子晶格磁化强度 | $\mathrm{FeF_2}$ | 78.26 |
| 超流 | He 原子的量子概率幅度 | $^4\mathrm{He}$ | 1.8 ~ 2.1 |
| 超导 | 电子对的量子概率幅度 | Pb | 7.19 |
| 二元合金 | 次晶格中某组元的密度 | Cu-Zn | 739 |

---

## *§3.8　临界现象和临界指数

> 带星号的选读内容。

连续相变的相变点也称**临界点**。**临界现象指物质在连续相变临界点邻域的行为。** 如前所述，在连续相变临界点的邻域，与化学势二阶导数相应的热容、等温压缩系数、磁化率等出现跃变或无穷尖峰。**人们用幂函数表述一些热力学量在临界点邻域的特性，其幂次（负幂次）称为临界指数。**

### 3.8.1 液气系统的临界现象与临界指数

先介绍液气流体系统。图 3.5.2 以体积和压强为坐标画出了流体系统的等温线。改以密度和压强为坐标，等温线将如图 3.8.1 所示。$\rho_c$ 表示物质在临界点的密度，两侧的虚线称为**共存线**，分别表示两相平衡时气相和液相的密度 $\rho_g$ 和 $\rho_l$。以

$$t = \frac{T - T_c}{T_c}$$

表示温度与临界温度差的对比值，称为**对比温度**或**约化温度**。人们发现，各种液气流体系统在临界点的邻域存在以下几个共同的实验规律：

（1）当 $t\to-0$ 时，$\rho_l$ 与 $\rho_g$ 之差随 $-t$ 的变化遵从如下的规律：

$$\rho_l - \rho_g \propto (-t)^\beta, \qquad t\to-0 \tag{3.8.1}$$

临界指数 $\beta$ 的实验值在 0.32~0.35。如前所述，$t>0$ 时物质处在液气不分的状态，$\rho_l-\rho_g$ 为零。

（2）当 $t\to\pm0$ 时，物质的等温压缩系数 $\kappa_T = -\dfrac{1}{v}\left(\dfrac{\partial v}{\partial p}\right)_T = \dfrac{1}{\rho}\left(\dfrac{\partial \rho}{\partial p}\right)_T$ 是发散的。这意味着在临界点的邻域，偶然的压强涨落将导致显著的密度涨落。$\kappa_T$ 随 $t$ 的变化规律为

$$\begin{cases} \kappa_T \propto t^{-\gamma}, & t\to+0 \\ \kappa_T \propto (-t)^{-\gamma'}, & t\to-0 \end{cases} \tag{3.8.2}$$

式中在 $t>0$ 时，沿临界等容线 $\rho=\rho_c$ 趋于临界点；在 $t<0$ 时，沿两相共存线趋于临界点。临界指数 $\gamma$ 的实验值在 1.2~1.3，$\gamma'$ 的实验值在 1.1~1.2。

（3）在临界等温线 $t=0$，压强与临界压强之差 $p-p_c$ 和密度与临界密度之差 $\rho-\rho_c$ 在临界点的邻域遵从以下规律：

$$p - p_c \propto \pm\left|\rho - \rho_c\right|^\delta \tag{3.8.3}$$

临界指数 $\delta$ 的实验值在 4.6~5.0。

（4）当 $t\to\pm0$ 时，物质的定容比热是发散的。这意味着，在临界点的邻域系统达到热平衡非常困难；为了使整个系统保持在恒定的温度（例如在 $10^{-2}$ K 甚至在 $10^{-1}$ K 的范围），往往需要很长的时间，并不断进行搅拌。$c_V$ 随 $t$ 的变化规律为

$$\begin{cases} c_V \propto t^{-\alpha}, & t\to+0 \\ c_V \propto (-t)^{-\alpha'}, & t\to-0 \end{cases} \tag{3.8.4}$$

式中 $t\to\pm0$ 沿临界等容线即 $\rho=\rho_c$ 趋于临界点。临界指数 $\alpha$ 和 $\alpha'$ 的实验值在 0.1~0.2。

### 3.8.2 铁磁系统的临界现象与临界指数

现在介绍铁磁系统。以 $T_c$ 表示铁磁顺磁转变的临界温度（也称为**居里温度**）。在 $T_c$ 以下，物质处在铁磁状态；在 $T_c$ 以上，处在顺磁状态。铁磁物质的特征是在外磁场为零时，物质的磁化强度不为零，称为**自发磁化强度**。自发磁化强度 $\mathcal{M}$ 是温度 $T$ 的函数，$\mathcal{M}(T)$ 随温度的升高而减小。当温度达到临界温度 $T_c$ 时，自发磁化强度为零，物质转变为顺磁状态。顺磁状态没有自发磁化，但在外磁场作用下可以发生磁化。

人们发现各种铁磁系统在临界点的邻域存在以下几个共同的实验规律：

（1）当 $t\to-0$ 时，自发磁化强度随 $-t$ 的变化遵从以下规律：

$$\mathcal{M} \propto (-t)^\beta, \qquad t\to-0 \tag{3.8.1'}$$

临界指数 $\beta$ 的实验值在 0.30~0.36。如前所述，在临界温度以上，$\mathcal{M}=0$。

（2）当 $t\to\pm0$ 时，各种铁磁物质的零场磁化率 $\chi = \left(\dfrac{\partial \mathcal{M}}{\partial \mathcal{H}}\right)_T$ 是发散的，$\chi$ 随 $t$ 的变化规律为

$$\begin{cases} \chi \propto t^{-\gamma}, & t\to+0 \\ \chi \propto (-t)^{-\gamma'}, & t\to-0 \end{cases} \tag{3.8.2'}$$

临界指数 $\gamma$ 的实验值在 1.2~1.4，$\gamma'$ 的实验值在 1.0~1.2。

（3）在 $t=0$ 和十分弱的磁场下，磁化强度 $\mathcal{M}$ 与外加磁场 $\mathcal{H}$ 的关系为

$$\mathcal{M} \propto \mathcal{H}^{1/\delta} \tag{3.8.3'}$$

临界指数 $\delta$ 的实验值在 4.2~4.8。

（4）当 $t\to\pm0$ 时，铁磁物质的零场比热容 $c_{\mathcal{H}}(\mathcal{H}=0)$ 遵从以下规律：

$$\begin{cases} c_{\mathcal{H}} \propto t^{-\alpha}, & t>0 \\ c_{\mathcal{H}} \propto (-t)^{-\alpha'}, & t<0 \end{cases} \tag{3.8.4'}$$

临界指数 $\alpha$ 和 $\alpha'$ 的实验值在 0.0~0.2。

### 3.8.3 临界现象的普适性

我们看到，**不仅各种流体系统和各种铁磁系统在临界点的邻域遵从相同的规律，而且，如果将液气密度差 $\rho_l-\rho_g$ 比作磁化强度 $\mathcal{M}$，压强 $p$ 比作磁场强度 $\mathcal{H}$，等温压缩系数 $\kappa_T$ 比作磁化率 $\chi$，上述两个系统的物理特性虽然很不相同，但在临界点邻域的行为却有极大的相似性，不仅变化规律相同，临界指数也大致相等。这一事实显示临界现象具有某种普适性。**

| 液气系统 | 铁磁系统 |
| --- | --- |
| 密度差 $\rho_l - \rho_g$ | 磁化强度 $\mathcal{M}$ |
| 压强 $p$ | 磁场强度 $\mathcal{H}$ |
| 等温压缩系数 $\kappa_T$ | 磁化率 $\chi$ |

---

## §3.9　朗道连续相变理论

### 3.9.1 序参量与对称破缺

1937 年，朗道（Landau）试图对连续相变提供一个统一的描述。他提出了**序参量**的概念，认为连续相变的特征是物质有序程度的改变及与之相伴随的物质对称性质的变化。我们以三维各向同性的铁磁体为例加以说明。

铁磁物质的原子具有固有的磁矩。两个相邻原子在其磁矩平行时相互作用能量较低。绝对零度下物质处在能量最低的状态，所有原子的磁矩取向都相同，是磁完全有序的状态。温度升高时，热运动有减弱有序取向的趋势。不过只要温度不太高，仍有为数较多的原子磁矩沿某一取向。这就是铁磁物质存在**自发磁化强度** $\mathcal{M}$，且 $\mathcal{M}$ 随温度升高而减小的原因。我们可以用 $\mathcal{M}(T)$ 作为**序参量**来描述铁磁物质的有序程度。

- 当温度升高达到临界温度 $T_c$ 时，自发磁化强度减小为零，物质转变为顺磁体。**在临界温度以上自发磁化强度为零时，系统是各向同性的，或者说系统具有空间转动的对称性。**
- **在临界温度以下自发磁化强度非零时，自发磁化强度向量 $\mathcal{M}$ 处在空间某个方向而破坏了系统对于空间转动的对称性。**

通常在**临界温度以下的相，对称性较低，有序度较高，序参量非零；临界温度以上的相，对称性较高，有序度较低，序参量为零**。随着温度的降低，**序参量在临界点连续地从零变到非零**。

下面我们着重讨论单轴各向异性的铁磁体。单轴铁磁体具有一个容易磁化的晶轴。原子磁矩的取向只能平行或反平行于这个晶轴，因此它的序参量 $\mathcal{M}(T)$ 也只能沿这个轴。$\mathcal{M}$ 取正值对应于磁矩向上，$\mathcal{M}$ 取负值对应于磁矩向下，由偶然的因素决定。**对于自发磁化强度为零的顺磁状态，上和下两个方向是对称或等价的。对于自发磁化强度非零的铁磁状态，自发磁化强度只能在上、下之中取一个方向。当系统由顺磁状态转变为铁磁状态时，系统的对称性就突然降低，发生了对称破缺。** 图 3.9.1 示意地画出了单轴铁磁体的 $\mathcal{M}(T)$ 随温度的变化。

### 3.9.2 朗道自由能与自发磁化强度

朗道认为，在临界点 $T_c$ 附近，序参量 $\mathcal{M}$ 是一个小量，可以将自由能 $F(T,\mathcal{M})$ 在 $T_c$ 附近按 $\mathcal{M}$ 的幂展开，称为**朗道自由能**：

$$F(T,\mathcal{M}) = F_0(T) + \frac{1}{2}a(T)\mathcal{M}^2 + \frac{1}{4}b(T)\mathcal{M}^4 + \cdots \tag{3.9.1}$$

其中 $F_0(T)$ 是 $\mathcal{M}=0$ 时的自由能。**由于系统对于变换 $\mathcal{M}\to-\mathcal{M}$ 是对称的，展开式不含 $\mathcal{M}$ 的奇次幂。**

先讨论不存在外磁场的情形。在此情形下，自由能的自变量是温度 $T$ 和体积 $V$。假设体积的变化可以忽略，$V$ 是常量，并不失普遍性地令 $V$ 等于单位体积。这样状态参量就只有温度 $T$。式中展开系数 $a(T)$ 和 $b(T)$ 是温度 $T$ 的函数。注意，我们将自由能写作 $F(T,\mathcal{M})$，其中的 $\mathcal{M}$ 并不是自变量。我们将利用在 $T$、$V$ 不变的条件下稳定平衡自由能最小的条件确定 $\mathcal{M}$ 为温度 $T$ 的函数。

稳定平衡态下 $F$ 具有极小值，应有

$$\frac{\partial F}{\partial \mathcal{M}} = \mathcal{M}\left(a + b\mathcal{M}^2\right) = 0 \tag{3.9.2}$$

$$\frac{\partial^2 F}{\partial \mathcal{M}^2} = a + 3b\mathcal{M}^2 > 0 \tag{3.9.3}$$

式 (3.9.2) 有三个解：

$$\mathcal{M} = 0, \qquad \mathcal{M} = \pm\sqrt{-\frac{a}{b}} \tag{3.9.4}$$

$\mathcal{M}=0$ 代表无序态，相应于 $T>T_c$ 的温度范围。将 $\mathcal{M}=0$ 代入式 (3.9.3) 可知，在 $T>T_c$ 时，$a>0$。非零解 $\mathcal{M}=\pm\sqrt{-\dfrac{a}{b}}$ 代表有序态，相应于 $T<T_c$ 的温度范围。将 $\mathcal{M}=\pm\sqrt{-\dfrac{a}{b}}$ 代入式 (3.9.3) 可知，在 $T<T_c$ 时，$a<0$。序参量在 $T_c$ 处连续地由零转变到非零，所以在 $T=T_c$ 处应有 $a=0$。**在临界点的邻域可以简单地假设**

$$a = a_0\left(\frac{T-T_c}{T_c}\right) = a_0 t, \qquad a_0 > 0 \tag{3.9.5}$$

和

$$b(T) = b\ (\text{常量}) \tag{3.9.6}$$

其中 $t$ 为对比温度。因为式 (3.9.4) 给出的 $\mathcal{M}=\pm\left(-\dfrac{a}{b}\right)^{\frac{1}{2}}$ 应是实数，而在 $T<T_c$ 时，$a<0$，故常量 $b>0$。

将式 (3.9.5) 和式 (3.9.6) 代入式 (3.9.4) 可知，在临界点的邻域，单轴铁磁体的自发磁化强度 $\mathcal{M}$ 为

$$\begin{aligned} &\mathcal{M} = 0, && t>0 \\ &\mathcal{M} = \pm\left(\frac{a_0}{b}\right)^{1/2}(-t)^{1/2}, && t<0 \end{aligned} \tag{3.9.7}$$

式 (3.9.7) 所给出的 $\mathcal{M}$ 对 $t$ 的依赖关系与式 (3.8.1′) 相同，**临界指数 $\beta=\dfrac{1}{2}$**。

由此可知，**在 $T_c$ 的两侧热力学函数各存在一个分支**。对于 $T>T_c$ 的无序相（$a>0$，$b>0$），朗道自由能的极小值在 $\mathcal{M}=0$ 处，即无序相的自由能为

$$F = F_0, \qquad T>T_c \tag{3.9.8a}$$

对于 $T<T_c$ 的有序相（$a<0$，$b>0$），朗道自由能的极小值在 $\mathcal{M}=\pm\left(-\dfrac{a}{b}\right)^{\frac{1}{2}}$ 处，即有序相的自由能为

$$F = F_0 - \frac{a^2}{4b}, \qquad T<T_c \tag{3.9.8b}$$

在临界点 $T=T_c$（$a=0$，$b>0$），自由能的两个分支重合。图 3.9.2 画出了式 (3.9.1) 朗道自由能在 $T>T_c$［图 (a)］和 $T<T_c$［图 (b)］随 $\mathcal{M}$ 的变化曲线。

### 3.9.3 临界指数的导出

现在讨论铁磁体的零场比热容。利用公式 $c = -T\dfrac{\partial^2 F}{\partial T^2}$ 可得，在 $T=T_c$ 处，两相的比热容之差为

$$c(t\to-0) - c(t\to+0) = \frac{a_0^2}{2bT_c} \tag{3.9.9}$$

上式表明，**有序相的比热容大于无序相的比热容，且在 $t=0$ 处比热容的突变是有限的。由此可知，临界指数 $\alpha=\alpha'=0$。**

存在外磁场的情形下，如果体积的变化可以忽略，铁磁体自由能的自变量是温度 $T$ 和磁化强度 $\mathcal{M}$。根据式 (1.18.3) 和式 (2.7.3)，自由能的全微分为

$$\mathrm{d}F = -S\,\mathrm{d}T + \mu_0\mathcal{H}\,\mathrm{d}\mathcal{M} \tag{3.9.10}$$

因此

$$\mu_0\mathcal{H} = \left(\frac{\partial F}{\partial \mathcal{M}}\right)_T = a\mathcal{M} + b\mathcal{M}^3 \tag{3.9.11}$$

式中第二步利用了朗道自由能表达式 (3.9.1)。我们假设外磁场十分弱，式 (3.9.1) 仍然近似适用。再将上式对 $\mathcal{H}$ 求偏导数，可得**磁化率**

$$\chi = \left(\frac{\partial \mathcal{M}}{\partial \mathcal{H}}\right)_T = \frac{\mu_0}{a + 3b\mathcal{M}^2} = \begin{cases} \dfrac{\mu_0}{a_0}t^{-1}, & t>0 \\[3mm] \dfrac{\mu_0}{2a_0}(-t)^{-1}, & t<0 \end{cases} \tag{3.9.12}$$

最后一步已将式中的 $\mathcal{M}^2$ 近似地用自发磁化强度值式 (3.9.7) 代入。与式 (3.8.2′) 比较可知，**$\gamma=\gamma'=1$**。

最后，在 $T=T_c$ 即 $a=0$ 时，式 (3.9.11) 给出

$$\mathcal{H} \propto \mathcal{M}^3 \tag{3.9.13}$$

与式 (3.8.3′) 比较可知，**临界指数 $\delta=3$**。

综上所述，对于单轴铁磁体，朗道理论导出了描述其临界行为的式 (3.8.1′) 至式 (3.8.4′)，所得临界指数为

$$\alpha = \alpha' = 0, \qquad \beta = \frac{1}{2}, \qquad \gamma = \gamma' = 1, \qquad \delta = 3 \tag{3.9.14}$$

我们看到朗道自由能的展开表达式 (3.9.1) 含有与具体系统特性有关的参量 $a$、$b$。但式 (3.9.14) 给出的临界指数却与参量 $a$、$b$ 无关，**显示物质系统在临界点邻域行为的普适性**。不过目前的实验结果和理论分析表明，**临界指数与物质系统的空间维数 $d$ 有关，而朗道理论的临界指数与空间维数无关，普适性显得过高了**。另外，与 §3.8 中介绍的实验值比较，式 (3.9.14) 的临界指数数值上也存在差异。产生这些差异的根本原因是，**作为热力学理论，朗道理论将序参量看作是在整个系统中均匀而没有涨落的**。实际上，在连续相变临界点的邻域，序参量往往有强烈的涨落。我们将在 §10.2 至 §10.4 中以铁磁体为例，介绍序参量的涨落和涨落的关联，引入新的与涨落有关的临界指数，并在此基础上导出临界指数间存在的关系——**标度关系**。尽管存在不足，朗道理论仍然是非常重要和成功的理论，在定性分析物质的相结构时有很大价值；由于数学处理简单，面对新的问题，人们往往用朗道理论进行初步的分析。值得强调，在相对涨落较小的情形，例如低温零磁场下金属正常态与超导态的转变，以及液晶中的某种转变，朗道理论（或推广的朗道理论）得到的结果与实验结果是相符的。

> **脚注**（教材 p.86）：Chaikin P M, Lubensky T C. *Principles of Condensed Matter Physics*[M]. Cambridge: Cambridge University Press, 1997; Chap 4.

朗道理论提出了一些重要概念如**序参量、对称破缺、普适性**等，这些概念在临界现象的现代理论中仍然起着重要作用。

---

## 本章小结

### 一、公式速查

**三种平衡判据**

| 判据 | 适用条件 | 稳定平衡的充要条件 | 平衡条件 | 稳定性条件 |
| --- | --- | --- | --- | --- |
| 熵判据 | 孤立系（体积、内能不变） | $\Delta S<0$ (3.1.1) | $\delta S=0$ | $\delta^2 S<0$ |
| 自由能判据 | 等温等体 | $\Delta F>0$ (3.1.3) | $\delta F=0$ | $\delta^2 F>0$ |
| 吉布斯函数判据 | 等温等压 | $\Delta G>0$ (3.1.5) | $\delta G=0$ | $\delta^2 G>0$ |

**均匀系的平衡与稳定性**

| 内容 | 公式 |
| --- | --- |
| 平衡条件 | $T=T_0$，$p=p_0$ (3.1.9) |
| $\delta^2 S$ 的平方和形式 | $\delta^2 S = -\dfrac{C_V}{T^2}(\delta T)^2 + \dfrac{1}{T}\left(\dfrac{\partial p}{\partial V}\right)_T(\delta V)^2 < 0$ (3.1.13) |
| 稳定性条件 | $C_V>0$，$\left(\dfrac{\partial p}{\partial V}\right)_T<0$ (3.1.14) |

**开系、化学势与巨热力势**

| 内容 | 公式 |
| --- | --- |
| 开系的吉布斯函数 | $\mathrm{d}G = -S\mathrm{d}T + V\mathrm{d}p + \mu\mathrm{d}n$ (3.2.2) |
| 化学势 | $\mu = \left(\dfrac{\partial G}{\partial n}\right)_{T,p} = G_m$ (3.2.3)(3.2.5) |
| $G(T,p,n)$ 的偏导数 | $S=-\left(\dfrac{\partial G}{\partial T}\right)_{p,n}$，$V=\left(\dfrac{\partial G}{\partial p}\right)_{T,n}$，$\mu=\left(\dfrac{\partial G}{\partial n}\right)_{T,p}$ (3.2.6) |
| 开系基本方程 | $\mathrm{d}U = T\mathrm{d}S - p\mathrm{d}V + \mu\mathrm{d}n$ (3.2.7)；$\mathrm{d}H$ (3.2.8)；$\mathrm{d}F$ (3.2.9) |
| 巨热力势 | $J = F - \mu n$ (3.2.10)；$\mathrm{d}J = -S\mathrm{d}T - p\mathrm{d}V - n\mathrm{d}\mu$ (3.2.11)；$J = F - G = -pV$ (3.2.13) |

**单元复相系的平衡与相平衡曲线**

| 内容 | 公式 |
| --- | --- |
| 两相平衡条件 | $T^\alpha=T^\beta$（热）、$p^\alpha=p^\beta$（力学）、$\mu^\alpha=\mu^\beta$（相变）(3.3.6) |
| 两相平衡曲线 | $\mu^\alpha(T,p)=\mu^\beta(T,p)$ (3.4.1) |
| 三相平衡 | $T^\alpha=T^\beta=T^\gamma$，$p^\alpha=p^\beta=p^\gamma$，$\mu^\alpha=\mu^\beta=\mu^\gamma$ (3.4.2) |
| 克拉珀龙方程 | $\dfrac{\mathrm{d}p}{\mathrm{d}T} = \dfrac{L}{T\left(V_m^\beta - V_m^\alpha\right)}$ (3.4.6)，$L=T\left(S_m^\beta-S_m^\alpha\right)$ (3.4.5) |
| 蒸气压方程 | $\dfrac{1}{p}\dfrac{\mathrm{d}p}{\mathrm{d}T} = \dfrac{L}{RT^2}$ (3.4.7)；$\ln p = -\dfrac{L}{RT} + A$ (3.4.8) |

**临界点、范德瓦耳斯方程与对应态定律**

| 内容 | 公式 |
| --- | --- |
| 杠杆规则 | $v = xv_l + (1-x)v_g$ (3.5.1) |
| 范德瓦耳斯方程 | $\left(p+\dfrac{a}{V_m^2}\right)(V_m-b) = RT$ (3.5.2) |
| 化学势差 | $\mu-\mu_0 = \displaystyle\int_{p_0}^{p} V_m\,\mathrm{d}p$ (3.5.4) |
| 麦克斯韦等面积法则 | 面积$(BND)$ = 面积$(DJA)$ (3.5.5)，即 $\displaystyle\int_{BNDJA}V_m\mathrm{d}p = 0$ |
| 临界点方程 | $\left(\dfrac{\partial p}{\partial V_m}\right)_T = 0$，$\left(\dfrac{\partial^2 p}{\partial V_m^2}\right)_T = 0$ (3.5.6) |
| 临界参量 | $T_c = \dfrac{8a}{27Rb}$，$p_c = \dfrac{a}{27b^2}$，$V_{mc} = 3b$ |
| 临界系数 | $\dfrac{RT_c}{p_cV_{mc}} = \dfrac{8}{3} = 2.667$ (3.5.7) |
| 对比方程 | $\left(p^*+\dfrac{3}{v^{*2}}\right)\left(v^*-\dfrac{1}{3}\right) = \dfrac{8}{3}t^*$ (3.5.8) |
| 临界态稳定性 | $\left(\dfrac{\partial p}{\partial V_m}\right)_T = 0$，$\left(\dfrac{\partial^2 p}{\partial V_m^2}\right)_T = 0$，$\left(\dfrac{\partial^3 p}{\partial V_m^3}\right)_T < 0$ (3.5.9) |

**液滴与曲面分界面（*§3.6）**

| 内容 | 公式 |
| --- | --- |
| 力学平衡条件 | $p^\alpha = p^\beta + \dfrac{2\sigma}{r}$ (3.6.6)；相变平衡条件 $\mu^\alpha=\mu^\beta$ (3.6.7) |
| 曲面平衡蒸气压 | $\left(p'-p+\dfrac{2\sigma}{r}\right)v^\alpha = RT\ln\dfrac{p'}{p}$ (3.6.13)；$\ln\dfrac{p'}{p} = \dfrac{2\sigma v^\alpha}{RTr}$ (3.6.14) |
| 临界半径 | $r_c = \dfrac{2\sigma v^\alpha}{RT\ln\frac{p'}{p}}$ (3.6.15) |
| 气泡（$r\to-r$） | $p^\beta = p^\alpha + \dfrac{2\sigma}{r}$ (3.6.16)；$\ln\dfrac{p}{p'} = \dfrac{2\sigma v^\alpha}{RTr}$ (3.6.17) |

**相变的分类（§3.7）**

| 内容 | 公式 |
| --- | --- |
| 一级相变 | $\mu^{(1)}=\mu^{(2)}$，$\dfrac{\partial \mu^{(1)}}{\partial T}\neq\dfrac{\partial \mu^{(2)}}{\partial T}$，$\dfrac{\partial \mu^{(1)}}{\partial p}\neq\dfrac{\partial \mu^{(2)}}{\partial p}$ (3.7.1) |
| 二级相变的特征量 | $c_p=-T\dfrac{\partial^2\mu}{\partial T^2}$，$\alpha=\dfrac{1}{v}\dfrac{\partial^2\mu}{\partial T\partial p}$，$\kappa_T=-\dfrac{1}{v}\dfrac{\partial^2\mu}{\partial p^2}$ (3.7.2) |
| 埃伦菲斯特方程 | $\dfrac{\mathrm{d}p}{\mathrm{d}T} = \dfrac{\alpha^{(2)}-\alpha^{(1)}}{\kappa_T^{(2)}-\kappa_T^{(1)}} = \dfrac{c_p^{(2)}-c_p^{(1)}}{Tv\left(\alpha^{(2)}-\alpha^{(1)}\right)}$ (3.7.3) |

**临界指数（\*§3.8）**

| 系统  | 规律                                                                                                                                                                                                                     | 实验值                                                                                                                                                             |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 液气  | $\rho_l-\rho_g \propto (-t)^\beta$ (3.8.1)；$\kappa_T \propto t^{-\gamma},\ (-t)^{-\gamma'}$ (3.8.2)；$p-p_c \propto \pm\lvert\rho-\rho_c\rvert^\delta$ (3.8.3)；$c_V \propto t^{-\alpha},\ (-t)^{-\alpha'}$ (3.8.4)      | $\beta\approx0.32\text{–}0.35$，$\gamma\approx1.2\text{–}1.3$，$\gamma'\approx1.1\text{–}1.2$，$\delta\approx4.6\text{–}5.0$，$\alpha=\alpha'\approx0.1\text{–}0.2$ |
| 铁磁  | $\mathcal{M} \propto (-t)^\beta$ (3.8.1′)；$\chi \propto t^{-\gamma},\ (-t)^{-\gamma'}$ (3.8.2′)；$\mathcal{M} \propto \mathcal{H}^{1/\delta}$ (3.8.3′)；$c_{\mathcal{H}} \propto t^{-\alpha},\ (-t)^{-\alpha'}$ (3.8.4′) | $\beta\approx0.30\text{–}0.36$，$\gamma\approx1.2\text{–}1.4$，$\gamma'\approx1.0\text{–}1.2$，$\delta\approx4.2\text{–}4.8$，$\alpha=\alpha'\approx0.0\text{–}0.2$ |

**朗道连续相变理论（§3.9）**

| 内容 | 公式 |
| --- | --- |
| 朗道自由能 | $F = F_0(T) + \dfrac{1}{2}a(T)\mathcal{M}^2 + \dfrac{1}{4}b(T)\mathcal{M}^4 + \cdots$ (3.9.1) |
| 极值条件 | $\dfrac{\partial F}{\partial \mathcal{M}} = \mathcal{M}\left(a+b\mathcal{M}^2\right)=0$ (3.9.2)；$\dfrac{\partial^2 F}{\partial \mathcal{M}^2} = a+3b\mathcal{M}^2>0$ (3.9.3) |
| 序参量 | $\mathcal{M}=0$ 或 $\pm\sqrt{-\dfrac{a}{b}}$ (3.9.4)；$a=a_0t$ (3.9.5)；$b>0$ 为常量 (3.9.6)；$\mathcal{M}\propto(-t)^{1/2}$ (3.9.7) |
| 两个分支的自由能 | $F=F_0$（$T>T_c$）(3.9.8a)；$F=F_0-\dfrac{a^2}{4b}$（$T<T_c$）(3.9.8b) |
| 比热容突变 | $c(t\to-0)-c(t\to+0) = \dfrac{a_0^2}{2bT_c}$ (3.9.9) |
| 状态方程与磁化率 | $\mu_0\mathcal{H} = a\mathcal{M}+b\mathcal{M}^3$ (3.9.11)；$\chi = \dfrac{\mu_0}{a+3b\mathcal{M}^2}$ (3.9.12)；$\mathcal{H}\propto\mathcal{M}^3$ (3.9.13) |
| 临界指数 | $\alpha=\alpha'=0$，$\beta=\dfrac{1}{2}$，$\gamma=\gamma'=1$，$\delta=3$ (3.9.14) |

### 二、概念脉络

```
§1.16 熵增加原理 / §1.18 F、G 的极小性质
        │
        ▼
§3.1 热动平衡判据（熵判据 / 自由能判据 / 吉布斯函数判据）
        │  δ(…) = 0        → 平衡条件
        │  δ²(…) > 0 或 < 0 → 稳定性条件
        ▼
   均匀系：T = T0，p = p0；C_V > 0，(∂p/∂V)_T < 0
        │
        ▼
§3.2 开系基本方程（引入 μ）：dU = T dS − p dV + μ dn
        │  μ = G_m；巨热力势 J = F − μn = −pV
        ▼
§3.3 单元系复相平衡条件
     T^α = T^β（热平衡）、p^α = p^β（力学平衡）、μ^α = μ^β（相变平衡）
        │
        ▼
§3.4 单元复相系的平衡性质
     相图：汽化线 / 熔化线 / 升华线 / 三相点 / 临界点
     两相平衡曲线方程 μ^α(T,p) = μ^β(T,p)
        ├──► 克拉珀龙方程 dp/dT = L / [T(V_m^β − V_m^α)]
        └──► 蒸气压方程 ln p = −L/(RT) + A
        │
        ▼
§3.5 临界点和气液两相的转变
     p-V 等温线、杠杆规则、麦克斯韦等面积法则
     临界点 (∂p/∂V_m)_T = (∂²p/∂V_m²)_T = 0，(∂³p/∂V_m³)_T < 0
     范德瓦耳斯对比方程 → 对应态定律
        │
        ├──► *§3.6 液滴的形成：p^α = p^β + 2σ/r
        │        → ln(p′/p) = 2σv^α/(RTr) → 临界半径 r_c
        │        过饱和蒸气（缺凝结核）、过热液体（缺汽化核）
        │
        ├──► §3.7 相变的分类
        │        一级相变（有 L、Δv，可有亚稳态）
        │        连续相变（无 L、Δv；c_p、κ_T、χ 跃变或发散）
        │        → 埃伦菲斯特方程
        │
        └──► *§3.8 临界现象和临界指数 α、α′、β、γ、γ′、δ
                 （液气系统 ↔ 铁磁系统 一一对应，具有普适性）
                    │
                    ▼
             §3.9 朗道连续相变理论
                   序参量、对称破缺、朗道自由能
                   α = α′ = 0，β = 1/2，γ = γ′ = 1，δ = 3（平均场）
                   → §10.2–§10.4 序参量涨落与标度关系
```

### 三、易错点

1. **三个判据的适用条件不能混**：熵判据用于孤立系，自由能判据用于等温等体，吉布斯函数判据用于等温等压。判据的形式（$\Delta<0$ 还是 $\Delta>0$）由"极大"还是"极小"决定，不要死记。
2. **平衡条件由 $\delta(\cdots)=0$ 给出，稳定性条件由 $\delta^2(\cdots)$ 的符号给出**，两者必须同时满足才是稳定平衡；$\delta^2=0$ 属中性平衡，亚稳平衡只是局域（而非全局）的极值。
3. $\delta^2\widetilde{S}\approx\delta^2 S$ 用的是 **$n_0\gg n$**：二阶偏导数与物质的量成反比，所以介质的二阶项可以忽略。
4. **相变平衡条件写成 $\mu^\alpha/T^\alpha=\mu^\beta/T^\beta$**（由熵判据得到），只有当两相温度相等时才退化为 $\mu^\alpha=\mu^\beta$；力学平衡条件同理。
5. **化学势对单元系等于摩尔吉布斯函数** $G_m$；因为 $G$ 是广延量，$G=nG_m$，这是 $\mu=G_m$ 的根据。
6. **克拉珀龙方程的次序**：$L$ 是 1 mol 物质由 $\alpha$ 相变到 $\beta$ 相吸收的潜热，分母必须是 $V_m^\beta-V_m^\alpha$，二者次序一旦颠倒就会差一个负号。
7. **平衡曲线上 $T$、$p$ 只有一个独立**，这是复相系自由度减少的体现（对比 §4.3 吉布斯相律）；两相平衡属于 §3.1 所说的**中性平衡**。
8. **杠杆规则 $v=xv_l+(1-x)v_g$** 中 $x$、$1-x$ 分别是液相、气相的"比例"，谁乘谁不要弄反。
9. **麦克斯韦等面积法则的实质是 $\mu_A=\mu_B$**，即 $\displaystyle\int_{BNDJA}V_m\,\mathrm{d}p=0$（两块面积相等、符号相反）。
10. 范德瓦耳斯等温线上 $\left(\partial p/\partial V_m\right)_T>0$ 的一段（$JDN$）**不满足稳定性条件、物理上不可实现**；$BN$ 与 $AJ$ 段分别是过饱和蒸气与过热液体的亚稳态。
11. **临界点上稳定性条件失效**：$\left(\partial p/\partial V_m\right)_T=0$ 使 (3.1.14) 的第二式不成立，必须用高级微分（教材用到 $\delta^4S_m$）才能判定，这正是 §3.5.5 的来由。
12. **临界系数是量纲一的量**；范德瓦耳斯理论给出 $8/3$，实测值（3.27~4.37）都比它大，说明范德瓦耳斯方程只在定性上正确。
13. **曲面分界面上两相的压强不相等**：$p^\alpha=p^\beta+2\sigma/r$，因此 $\mu^\alpha=\mu^\beta$ 中的两个压强取值不同——这是液滴问题最关键、最容易被忽略的一点。
14. 液滴公式中**把 $r$ 换成 $-r$** 即得气泡的公式；$r_c$ 是临界半径，$r>r_c$ 液滴凝结长大，$r<r_c$ 液滴汽化消失。
15. **一级相变与连续相变的判据是"化学势的哪一级偏导数出现突变"**，而不是"有没有潜热"（后者只是结果）。连续相变不存在两相共存，也不存在亚稳态，相变点是热力学函数的奇异性点。
16. **朗道理论是平均场理论**：它忽略序参量的涨落，所得临界指数（$\beta=1/2$，$\gamma=1$，$\delta=3$，$\alpha=0$）既与实验值有差别，又与空间维数无关，普适性"过头"了；改进的方向是引入序参量涨落与标度关系。

### 四、名词中英对照

| 中文 | 英文 | 中文 | 英文 |
| --- | --- | --- | --- |
| 虚变动 | virtual variation | 热动平衡判据 | criterion of thermodynamic equilibrium |
| 熵判据 | entropy criterion | 自由能判据 | free energy criterion |
| 吉布斯函数判据 | Gibbs function criterion | 平衡条件 | equilibrium condition |
| 稳定性条件 | stability condition | 中性平衡 | neutral equilibrium |
| 亚稳平衡 | metastable equilibrium | 单元系 | one-component system |
| 组元 | component | 复相系 | multiphase system |
| 开系 | open system | 化学势 | chemical potential |
| 摩尔吉布斯函数 | molar Gibbs function | 巨热力势 | grand potential |
| 热平衡条件 | thermal equilibrium condition | 力学平衡条件 | mechanical equilibrium condition |
| 相变平衡条件 | phase equilibrium condition | 相图 | phase diagram |
| 汽化线 | vaporization curve | 熔化线 | melting curve |
| 升华线 | sublimation curve | 三相点 | triple point |
| 临界点 | critical point | 临界温度 / 临界压强 | critical temperature / pressure |
| 平衡相变 | equilibrium phase transition | 克拉珀龙方程 | Clapeyron equation |
| 相变潜热 | latent heat | 饱和蒸气 | saturated vapor |
| 蒸气压方程 | vapor-pressure equation | 等温线 | isotherm |
| 杠杆规则 | lever rule | 麦克斯韦等面积法则 | Maxwell's equal-area rule |
| 过饱和蒸气 | supersaturated vapor | 过热液体 | superheated liquid |
| 对比温度 / 压强 / 体积 | reduced temperature / pressure / volume | 对应态定律 | law of corresponding states |
| 临界系数 | critical coefficient | 表面张力系数 | surface tension coefficient |
| 液滴 | droplet | 凝结核 | condensation nucleus |
| 临界半径 | critical radius | 气泡 | bubble |
| 一级相变 | first-order phase transition | 二级相变 | second-order phase transition |
| 连续相变 | continuous phase transition | 埃伦菲斯特方程 | Ehrenfest equation |
| 序参量 | order parameter | 对称破缺 | symmetry breaking |
| 临界现象 | critical phenomena | 临界指数 | critical exponent |
| 约化温度 | reduced temperature | 自发磁化强度 | spontaneous magnetization |
| 居里温度 | Curie temperature | 磁化率 | magnetic susceptibility |
| 朗道自由能 | Landau free energy | 标度关系 | scaling relation |
| 普适性 | universality | 平均场理论 | mean-field theory |
