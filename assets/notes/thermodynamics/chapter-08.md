---
tags:
  - 热统
  - 统计物理
  - 第八章
created: 2026-10-22
---

# 第八章　玻色统计和费米统计

> **教材**：汪志诚《热力学·统计物理》（第六版）
> **范围**：§8.1–§8.5（教材 pp. 186–201；其中 *§8.2 为带星号的选读内容）

## 本章地图

| 节 | 标题 | 一句话 |
| --- | --- | --- |
| §8.1 | 热力学量的统计表达式 | 以巨配分函数 $\ln\Xi$ 为特性函数，统一导出玻色、费米气体的 $\bar N,U,Y,S,J$ |
| *§8.2 | 弱简并理想玻色气体和费米气体 | 对量子分布作弱简并展开，得到量子统计关联对粒子数密度和内能的首个修正 |
| §8.3 | 玻色-爱因斯坦凝聚 | 温度降到 $T_c$ 后，宏观数目的玻色子占据最低能级；自由气体满足 $n\lambda^3\geqslant2.612$ |
| §8.4 | 光子气体 | 光子数不守恒使化学势为零；玻色统计给出普朗克公式、斯特藩-玻耳兹曼定律和辐射压强 |
| §8.5 | 金属中的自由电子气体 | 强简并费米气体在 $T=0$ 仍有费米能和简并压，低温电子热容与 $T$ 成正比 |

**主线**：巨配分函数 $\Xi$ → 玻色与费米分布的统计热力学 → 弱简并展开 → 玻色凝聚 → 光子气体 → 强简并电子气体。量子统计关联在弱简并区只是小修正，在强简并区则改变基态、压强、熵和热容。

---

## §8.1　热力学量的统计表达式

### 8.1.1 非简并与简并气体

第七章的经典极限条件可写成

$$
e^\alpha=\frac{V}{N}\left(\frac{2\pi mkT}{h^2}\right)^{3/2}\gg1,
$$

或

$$
n\lambda^3=\frac{N}{V}\left(\frac{h^2}{2\pi mkT}\right)^{3/2}\ll1.
$$

其中 $\lambda=h/\sqrt{2\pi mkT}$ 是热德布罗意波长。满足这些条件的气体称为**非简并气体**，无论粒子自旋满足玻色还是费米统计，都可用玻耳兹曼统计处理。

不满足上述条件的气体称为**简并气体**，必须分别使用玻色或费米分布。粒子全同性造成的量子统计关联会改变宏观性质；本节先用平均分布处理开放系统，第九章将用巨正则系综重新建立其分布。

在开放系统中，粒子可从粒子源进出。达到平衡时，系统与源具有相同的温度 $T$ 和化学势 $\mu$，因而具有相同的 $\alpha,\beta$，且

$$
\bar N=\sum_l\frac{\omega_l}{e^{\alpha+\beta\varepsilon_l}\pm1},
\qquad
U=\sum_l\frac{\omega_l\varepsilon_l}{e^{\alpha+\beta\varepsilon_l}\pm1}.
$$

本书约定：式中上面的“$+$”属于费米气体，下面的“$-$”属于玻色气体。平均分布的表达式与第六章最概然分布完全相同；区别在于这里把源确定的 $\alpha,\beta$ 和外参量 $y$ 看作已知量，热力学量表现为 $\alpha,\beta,y$ 的函数。

### 8.1.2 玻色气体的巨配分函数

对玻色系统，平均总粒子数为

$$
\bar N=\sum_l a_l
=\sum_l\frac{\omega_l}{e^{\alpha+\beta\varepsilon_l}-1}. \tag{8.1.1}
$$

定义**巨配分函数** $\Xi$ 为

$$
\Xi=\prod_l\Xi_l
=\prod_l\left(1-e^{-\alpha-\beta\varepsilon_l}\right)^{-\omega_l}. \tag{8.1.2}
$$

取对数：

$$
\ln\Xi
=-\sum_l\omega_l\ln\left(1-e^{-\alpha-\beta\varepsilon_l}\right). \tag{8.1.3}
$$

对式 (8.1.3) 关于 $\alpha$ 求导：

$$
\begin{aligned}
-\frac{\partial}{\partial\alpha}\ln\Xi
&=-\sum_l\omega_l
\frac{e^{-\alpha-\beta\varepsilon_l}}
{1-e^{-\alpha-\beta\varepsilon_l}}\\
&=\sum_l\frac{\omega_l}{e^{\alpha+\beta\varepsilon_l}-1}
=\bar N.
\end{aligned}
$$

故

$$
\bar N=-\frac{\partial}{\partial\alpha}\ln\Xi. \tag{8.1.4}
$$

系统的平均总能量为

$$
U=\sum_l\varepsilon_l a_l
=\sum_l\frac{\varepsilon_l\omega_l}{e^{\alpha+\beta\varepsilon_l}-1}. \tag{8.1.5}
$$

因为能级 $\varepsilon_l$ 与 $\beta$ 无关，

$$
-\frac{\partial}{\partial\beta}\ln\Xi
=\sum_l\frac{\varepsilon_l\omega_l}
{e^{\alpha+\beta\varepsilon_l}-1}
=U,
$$

即

$$
U=-\frac{\partial}{\partial\beta}\ln\Xi. \tag{8.1.6}
$$

### 8.1.3 广义力、熵和巨热力势

若能级依赖外参量 $y$，外界对位于第 $l$ 个能级的单个粒子所做的广义力为 $-\partial\varepsilon_l/\partial y$，因此系统所受平均广义力为

$$
\begin{aligned}
Y
&=-\sum_l a_l\frac{\partial\varepsilon_l}{\partial y}
=-\sum_l\frac{\omega_l}{e^{\alpha+\beta\varepsilon_l}-1}
\frac{\partial\varepsilon_l}{\partial y}\\
&=-\frac1\beta\frac{\partial}{\partial y}
\ln\Xi.
\end{aligned}
$$

于是

$$
Y=-\frac1\beta\frac{\partial}{\partial y}\ln\Xi. \tag{8.1.7}
$$

对于体积功 $Y=-p$、$y=V$，有

$$
p=\frac1\beta\frac{\partial}{\partial V}\ln\Xi. \tag{8.1.8}
$$

由式 (8.1.4)–(8.1.7) 可得

$$
\begin{aligned}
\beta\left(\mathrm dU-Y\,\mathrm dy+\frac{\alpha}{\beta}\mathrm d\bar N\right)
={}&-\beta\,\mathrm d\left(\frac{\partial\ln\Xi}{\partial\beta}\right)
+\frac{\partial\ln\Xi}{\partial y}\mathrm dy
-\alpha\,\mathrm d\left(\frac{\partial\ln\Xi}{\partial\alpha}\right).
\end{aligned}
$$

又因为 $\ln\Xi=\ln\Xi(\alpha,\beta,y)$，

$$
\mathrm d\ln\Xi
=\frac{\partial\ln\Xi}{\partial\alpha}\mathrm d\alpha
+\frac{\partial\ln\Xi}{\partial\beta}\mathrm d\beta
+\frac{\partial\ln\Xi}{\partial y}\mathrm dy.
$$

因此前式成为

$$
\beta\left(\mathrm dU-Y\,\mathrm dy+\frac{\alpha}{\beta}\mathrm d\bar N\right)
=\mathrm d\left(
\ln\Xi-\alpha\frac{\partial\ln\Xi}{\partial\alpha}
-\beta\frac{\partial\ln\Xi}{\partial\beta}
\right).
$$

开放系统的热力学基本关系为

$$
\frac1T\left(\mathrm dU-Y\,\mathrm dy-\mu\,\mathrm d\bar N\right)=\mathrm dS.
$$

比较两式的积分因子与化学势项，得

$$
\beta=\frac1{kT},
\qquad
\alpha=-\frac{\mu}{kT}. \tag{8.1.9}
$$

于是

$$
\mathrm dS
=k\,\mathrm d\left(
\ln\Xi-\alpha\frac{\partial\ln\Xi}{\partial\alpha}
-\beta\frac{\partial\ln\Xi}{\partial\beta}
\right).
$$

积分后

$$
\begin{aligned}
S
&=k\left(
\ln\Xi-\alpha\frac{\partial\ln\Xi}{\partial\alpha}
-\beta\frac{\partial\ln\Xi}{\partial\beta}
\right)\\
&=k\left(\ln\Xi+\alpha\bar N+\beta U\right).
\end{aligned}
\tag{8.1.10}
$$

将式 (8.1.3) 代入，并与玻耳兹曼关系比较，可得

$$
S=k\ln\Omega. \tag{8.1.11}
$$

### 8.1.4 费米系统

费米系统的巨配分函数只需改为

$$
\Xi=\prod_l\Xi_l
=\prod_l\left(1+e^{-\alpha-\beta\varepsilon_l}\right)^{\omega_l}, \tag{8.1.12}
$$

其对数

$$
\ln\Xi
=\sum_l\omega_l\ln\left(1+e^{-\alpha-\beta\varepsilon_l}\right). \tag{8.1.13}
$$

前面导出的式 (8.1.4)、(8.1.6)–(8.1.8)、(8.1.10) 对费米系统同样适用。因此，只要知道能级、简并度和 $\ln\Xi$，便能通过求导与代数运算得到全部平衡热力学函数。

以 $T,V,\mu$ 为自然变量的特性函数是**巨热力势**

$$
J=U-TS-\bar N\mu.
$$

与式 (8.1.10) 比较得

$$
J=-kT\ln\Xi. \tag{8.1.14}
$$

实际常给定 $N,T,V$。此时由式 (8.1.4) 令 $\bar N=N$，得到 $N=\bar N(\alpha,T,V)$，由此隐式确定 $\alpha$；再代回热力学量。不同简并程度要求不同的近似。

> 教材脚注：那里指出得到 $N,T,V$ 表象需要分级近似，并引用王竹溪《统计物理学导论》第二版第 64、65 节。

---

## *§8.2　弱简并理想玻色气体和费米气体

### 8.2.1 能态密度与基本积分

本节讨论 $e^{-\alpha}$ 或 $n\lambda^3$ 很小但不为零的弱简并气体。忽略内部自由度，只保留平动。单粒子能量为

$$
\varepsilon=\frac1{2m}\left(p_x^2+p_y^2+p_z^2\right). \tag{8.2.1}
$$

体积 $V$ 内，能量在 $\varepsilon$ 到 $\varepsilon+\mathrm d\varepsilon$ 间的微观状态数为

$$
D(\varepsilon)\mathrm d\varepsilon
=\frac{2\pi V}{h^3}(2m)^{3/2}\varepsilon^{1/2}\mathrm d\varepsilon, \tag{8.2.2}
$$

其中已暂取自旋简并度 $g=1$。若粒子自旋简并度为 $g$，右端再乘 $g$。总粒子数和内能分别为

$$
N=g\frac{2\pi V}{h^3}(2m)^{3/2}
\int_0^\infty\frac{\varepsilon^{1/2}\mathrm d\varepsilon}
{e^{\alpha+\beta\varepsilon}\pm1}, \tag{8.2.3}
$$

$$
U=g\frac{2\pi V}{h^3}(2m)^{3/2}
\int_0^\infty\frac{\varepsilon^{3/2}\mathrm d\varepsilon}
{e^{\alpha+\beta\varepsilon}\pm1}. \tag{8.2.4}
$$

引入 $x=\beta\varepsilon$，利用 $\beta=1/(kT)$：

$$
N=g\frac{2\pi V}{h^3}(2mkT)^{3/2}
\int_0^\infty\frac{x^{1/2}\mathrm dx}{e^{\alpha+x}\pm1}, \tag{8.2.3'}
$$

$$
U=g\frac{2\pi V}{h^3}(2mkT)^{3/2}kT
\int_0^\infty\frac{x^{3/2}\mathrm dx}{e^{\alpha+x}\pm1}. \tag{8.2.4'}
$$

### 8.2.2 弱简并展开

当 $e^{-\alpha}\ll1$ 时，被积函数中的分布可展开：

$$
\begin{aligned}
\frac1{e^{\alpha+x}\pm1}
&=\frac1{e^{\alpha+x}(1\pm e^{-\alpha-x})}\\
&=e^{-\alpha-x}\left(1\mp e^{-\alpha-x}+\cdots\right).
\end{aligned}
$$

只取前两项：

$$
\frac1{e^{\alpha+x}\pm1}
=e^{-\alpha-x}\left(1\mp e^{-\alpha-x}\right). \tag{8.2.5}
$$

第一项对应玻耳兹曼分布；第二项是粒子全同性造成的首个量子统计修正。

利用

$$
\int_0^\infty x^{1/2}e^{-x}\mathrm dx=\frac{\sqrt\pi}{2},
\qquad
\int_0^\infty x^{3/2}e^{-x}\mathrm dx=\frac{3\sqrt\pi}{4},
$$

以及 $e^{-2x}$ 带来的因子 $2^{-3/2}$、$2^{-5/2}$，有

$$
N=g\left(\frac{2\pi mkT}{h^2}\right)^{3/2}V e^{-\alpha}
\left(1\mp\frac1{2^{3/2}}e^{-\alpha}\right), \tag{8.2.6}
$$

$$
U=\frac32g\left(\frac{2\pi mkT}{h^2}\right)^{3/2}VkTe^{-\alpha}
\left(1\mp\frac1{2^{5/2}}e^{-\alpha}\right). \tag{8.2.7}
$$

将两式相除，并在分母中保留到 $e^{-\alpha}$ 的一阶：

$$
\begin{aligned}
U
&=\frac32NkT
\frac{1\mp2^{-5/2}e^{-\alpha}}
{1\mp2^{-3/2}e^{-\alpha}}\\
&\approx\frac32NkT
\left[1\pm\left(2^{-3/2}-2^{-5/2}\right)e^{-\alpha}\right]\\
&=\frac32NkT\left(1\pm\frac1{4\sqrt2}e^{-\alpha}\right).
\end{aligned}
$$

第二项只需用零级近似。令式 (8.2.6) 的修正为零，得

$$
e^{-\alpha}
=\frac{N}{gV}\left(\frac{h^2}{2\pi mkT}\right)^{3/2}.
$$

代回上式：

$$
U=\frac32NkT\left[
1\pm\frac1{4\sqrt2}\frac1g\frac NV
\left(\frac{h^2}{2\pi mkT}\right)^{3/2}
\right], \tag{8.2.8}
$$

或写成

$$
U=\frac32NkT\left(1\pm\frac1{4\sqrt2g}n\lambda^3\right). \tag{8.2.8'}
$$

**符号约定**：上面的“$+$”对应费米气体，“$-$”对应玻色气体。因此费米气体的附加内能为正，表现为等效排斥；玻色气体的附加内能为负，表现为等效吸引。

---

## §8.3　玻色-爱因斯坦凝聚

### 8.3.1 化学势的上限

考虑由 $N$ 个自旋为零的全同近独立玻色子组成的理想气体。由式 (8.1.9)，玻色分布可写为

$$
a_l=\frac{\omega_l}{e^{(\varepsilon_l-\mu)/(kT)}-1}. \tag{8.3.1}
$$

占据数非负要求对所有 $\varepsilon_l$ 均为正，因此

$$
\varepsilon_l-\mu>0.
$$

设最低能级为 $\varepsilon_0$，则

$$
\varepsilon_0>\mu. \tag{8.3.2}
$$

若取能量零点 $\varepsilon_0=0$，则

$$
\mu<0. \tag{8.3.3}
$$

化学势由

$$
\frac1V\sum_l\frac{\omega_l}
{e^{(\varepsilon_l-\mu)/(kT)}-1}
=\frac NV=n \tag{8.3.4}
$$

确定。由于 $\varepsilon_l,\omega_l$ 与温度无关，给定数密度 $n$ 时，温度越低，$\mu$ 越大且越接近 $0$。

### 8.3.2 临界温度

在热力学极限且能级间距远小于 $kT$ 时，可将求和换成积分：

$$
\frac{2\pi}{h^3}(2m)^{3/2}
\int_0^\infty\frac{\varepsilon^{1/2}\mathrm d\varepsilon}
{e^{(\varepsilon-\mu)/(kT)}-1}
=n. \tag{8.3.5}
$$

当温度降到某一临界温度 $T_c$ 时，$\mu\to-0$。此时

$$
\frac{2\pi}{h^3}(2m)^{3/2}
\int_0^\infty\frac{\varepsilon^{1/2}\mathrm d\varepsilon}
{e^{\varepsilon/(kT_c)}-1}
=n. \tag{8.3.6}
$$

令 $x=\varepsilon/(kT_c)$：

$$
\frac{2\pi}{h^3}(2mkT_c)^{3/2}
\int_0^\infty\frac{x^{1/2}\mathrm dx}{e^x-1}
=n. \tag{8.3.7}
$$

利用

$$
\int_0^\infty\frac{x^{1/2}\mathrm dx}{e^x-1}
=\frac{\sqrt\pi}{2}\times2.612,
$$

得到

$$
T_c=\frac{2\pi}{(2.612)^{2/3}}
\frac{h^2}{mk}n^{2/3}. \tag{8.3.8}
$$

### 8.3.3 基态宏观占据

若 $T<T_c$，仍有 $\mu<0$，但 $\mu$ 只能无限接近 $0$。若继续使用式 (8.3.5)，激发态积分给出的粒子数将小于 $n$，与固定总粒子数矛盾。

原因是态密度含因子 $\varepsilon^{1/2}$，把求和换成积分时漏掉了 $\varepsilon=0$ 的基态。$T>T_c$ 时基态占据数很小；$T<T_c$ 时基态可容纳宏观数目的粒子，不能再把基态并入积分。正确关系为

$$
n_0(T)+\frac{2\pi}{h^3}(2m)^{3/2}
\int_0^\infty\frac{\varepsilon^{1/2}\mathrm d\varepsilon}
{e^{\varepsilon/(kT)}-1}
=n, \tag{8.3.9}
$$

其中 $n_0(T)=a_0/V$ 是基态粒子数密度；第二项是 $\varepsilon>0$ 的激发态粒子数密度 $n_{\varepsilon>0}$，其中已取 $\mu\to-0$。

令 $x=\varepsilon/(kT)$，并利用式 (8.3.7) 的同一积分：

$$
\begin{aligned}
n_{\varepsilon>0}
&=\frac{2\pi}{h^3}(2m)^{3/2}
\int_0^\infty\frac{\varepsilon^{1/2}\mathrm d\varepsilon}
{e^{\varepsilon/(kT)}-1}\\
&=\frac{2\pi}{h^3}(2mkT)^{3/2}
\int_0^\infty\frac{x^{1/2}\mathrm dx}{e^x-1}\\
&=n\left(\frac T{T_c}\right)^{3/2}.
\end{aligned}
\tag{8.3.10}
$$

代入式 (8.3.9)：

$$
n_0(T)=n\left[1-\left(\frac T{T_c}\right)^{3/2}\right]. \tag{8.3.11}
$$

温度低于 $T_c$ 时，宏观数目的粒子占据 $\varepsilon=0$ 的基态，这称为**玻色-爱因斯坦凝聚**，简称玻色凝聚；$T_c$ 称为**凝聚温度**。凝聚体的能量、动量和熵均为零，对压强也没有贡献。

> 教材图的含义：图 8.3.1 给出 $n_0/n$ 随 $T/T_c$ 的变化，在 $T=T_c$ 处降为零；图 8.3.2 给出 $C_V/(Nk)$ 在 $T_c$ 处出现导数不连续的尖峰。

### 8.3.4 凝聚体的内能与热容

$T<T_c$ 时，只有激发态粒子对内能有贡献：

$$
\begin{aligned}
U
&=\frac{2\pi V}{h^3}(2m)^{3/2}
\int_0^\infty\frac{\varepsilon^{3/2}\mathrm d\varepsilon}
{e^{\varepsilon/(kT)}-1}\\
&=\frac{2\pi V}{h^3}(2m)^{3/2}(kT)^{5/2}
\int_0^\infty\frac{x^{3/2}\mathrm dx}{e^x-1}.
\end{aligned}
$$

利用积分并代入式 (8.3.8)：

$$
U=0.770NkT\left(\frac T{T_c}\right)^{3/2}. \tag{8.3.12}
$$

定容热容

$$
\begin{aligned}
C_V
&=\left(\frac{\partial U}{\partial T}\right)_V\\
&=\frac52\frac UT\\
&=1.925Nk\left(\frac T{T_c}\right)^{3/2}.
\end{aligned}
\tag{8.3.13}
$$

$T<T_c$ 时 $C_V\propto T^{3/2}$；在 $T=T_c$ 处达到 $1.925Nk$；高温时应趋于经典值 $3Nk/2$。$C_V$ 在 $T_c$ 连续，但其对温度的导数发生跳变。

式 (8.3.8) 可改写为

$$
n\left(\frac{h}{\sqrt{2\pi mkT_c}}\right)^3
=n\lambda^3=2.612. \tag{8.3.14}
$$

因而理想玻色气体发生凝聚的条件为

$$
n\lambda^3\geqslant2.612. \tag{8.3.15}
$$

降低温度和增大数密度都有利于凝聚。

### 8.3.5 实验联系与磁光阱

液 $^{4}\mathrm{He}$ 在 $T_\lambda=2.17\ \mathrm K$ 发生 $\lambda$ 相变：高温相为 He I，低温相为 He II 并具有超流性。把 $m=6.65\times10^{-27}\ \mathrm{kg}$、$V_m=27.6\times10^{-6}\ \mathrm{m^3\cdot mol^{-1}}$ 代入式 (8.3.8)，得 $T_c=3.13\ \mathrm K$，与 $T_\lambda$ 接近。伦敦在 1938 年提出 $\lambda$ 相变可能同玻色凝聚有关；但液氦原子间有强相互作用，并非理想玻色气体。

对磁光阱中的稀薄原子气体，原子在三维各向异性谐振子势中运动，能级为

$$
\begin{aligned}
\varepsilon_{n_x,n_y,n_z}
={}&\hbar\omega_x\left(n_x+\frac12\right)
+\hbar\omega_y\left(n_y+\frac12\right)
+\hbar\omega_z\left(n_z+\frac12\right),\\
&n_x,n_y,n_z=0,1,2,\cdots.
\end{aligned}
\tag{8.3.16}
$$

$T<T_c$ 时，凝聚在基态

$$
\varepsilon_0=\frac{\hbar}{2}(\omega_x+\omega_y+\omega_z) \tag{8.3.17}
$$

上。临界温度由

$$
kT_c=\hbar\bar\omega\left(\frac N{1.202}\right)^{1/3},
\qquad
\bar\omega=(\omega_x\omega_y\omega_z)^{1/3} \tag{8.3.18}
$$

确定，凝聚比例

$$
\frac{N_0}{N}=1-\left(\frac T{T_c}\right)^3. \tag{8.3.19}
$$

式 (8.3.18) 和式 (8.3.19) 的推导在教材中留作习题。自由二维理想玻色气体不会凝聚，但处在磁光阱中的二维原子可凝聚。

> 教材脚注（实验史料）：1995 年 Anderson、Davis、Bradley 等人在 $^{87}\mathrm{Rb}$、$^{23}\mathrm{Na}$ 中实现玻色凝聚；随后 $^{7}\mathrm{Li}$ 亦实现凝聚。教材给出原子密度、数目和温度量级，并列出物质波干涉及理论综述文献。稀薄原子仍有相互作用，三种原子的相互作用符号并不相同。

---

## §8.4　光子气体

### 8.4.1 光子数不守恒与分布

空腔平衡辐射可视为光子气体。光子能量和动量满足

$$
\left\{
\begin{aligned}
p&=\hbar k,\\
\varepsilon&=\hbar\omega,
\end{aligned}
\right. \tag{8.4.1}
$$

利用 $\omega=ck$，得

$$
\varepsilon=cp. \tag{8.4.2}
$$

> 教材脚注：也可从相对论质能关系 $\varepsilon^2=c^2p^2+m^2c^4$ 令光子静质量 $m=0$ 得到。

光子是玻色子。腔壁不断发射和吸收光子，平衡时只有能量守恒，没有固定粒子数，因此只引入一个拉格朗日乘子 $\beta$。对应 $\alpha=0$，也就是

$$
\mu=0.
$$

光子的统计分布为

$$
a_l=\frac{\omega_l}{e^{\beta\varepsilon_l}-1}. \tag{8.4.3}
$$

光子自旋量子数为 $1$，自旋沿动量方向的两个投影 $\pm\hbar$ 对应两个独立偏振。体积 $V$ 中，动量在 $p$ 到 $p+\mathrm dp$ 内的量子态数为

$$
\frac{8\pi V}{h^3}p^2\mathrm dp. \tag{8.4.4}
$$

由 $p=\hbar\omega/c$、$\mathrm dp=\hbar\,\mathrm d\omega/c$，频率在 $\omega$ 到 $\omega+\mathrm d\omega$ 内的态数为

$$
\frac{V}{\pi^2c^3}\omega^2\mathrm d\omega. \tag{8.4.5}
$$

因此平均光子数为

$$
\frac{V}{\pi^2c^3}
\frac{\omega^2\mathrm d\omega}{e^{\hbar\omega/(kT)}-1}. \tag{8.4.6}
$$

### 8.4.2 普朗克公式及极限

每个光子能量为 $\hbar\omega$，所以辐射场的内能分布为

$$
U(\omega,T)\mathrm d\omega
=\frac{V}{\pi^2c^3}
\frac{\hbar\omega^3}{e^{\hbar\omega/(kT)}-1}\mathrm d\omega. \tag{8.4.7}
$$

这就是**普朗克公式**。它同时消除了经典能量均分导致的低频不准确和高频发散。

低频极限 $\hbar\omega/(kT)\ll1$ 时，

$$
e^{\hbar\omega/(kT)}-1
\approx\frac{\hbar\omega}{kT},
$$

故

$$
U(\omega,T)\mathrm d\omega
=\frac{V}{\pi^2c^3}\omega^2kT\,\mathrm d\omega. \tag{8.4.8}
$$

这是**瑞利-金斯公式**。

高频极限 $\hbar\omega/(kT)\gg1$ 时，分母中的 $-1$ 可略，

$$
U(\omega,T)\mathrm d\omega
=\frac{V}{\pi^2c^3}\hbar\omega^3
e^{-\hbar\omega/(kT)}\mathrm d\omega. \tag{8.4.9}
$$

这是**维恩公式**。指数因子使高频光子极少出现。

从波动图像看，一个确定圆频率、波矢和偏振的平面波是一个振动自由度，其量子化能量为

$$
\varepsilon_n=\hbar\omega\left(n+\frac12\right),
\qquad n=0,1,2,\cdots.
$$

玻色分布给出平均量子数

$$
\bar n=\frac1{e^{\hbar\omega/(kT)}-1}.
$$

$\bar n$ 既是平均光子数，也是振动自由度的平均量子数。低频模式 $\hbar\omega\ll kT$ 近似恢复经典值 $kT$；高频模式被冻结在 $n=0$ 的基态。

### 8.4.3 总能量、维恩位移和热力学函数

由式 (8.4.7)：

$$
U=\frac{V}{\pi^2c^3}
\int_0^\infty\frac{\hbar\omega^3}{e^{\hbar\omega/(kT)}-1}\mathrm d\omega.
$$

令 $x=\hbar\omega/(kT)$，则

$$
U=\frac{V\hbar}{\pi^2c^3}
\left(\frac{kT}{\hbar}\right)^4
\int_0^\infty\frac{x^3\mathrm dx}{e^x-1}.
$$

利用

$$
\int_0^\infty\frac{x^3\mathrm dx}{e^x-1}
=\frac{\pi^4}{15},
$$

得

$$
U=\frac{\pi^2k^4}{15c^3\hbar^3}VT^4. \tag{8.4.10}
$$

这说明平衡辐射的内能密度与 $T^4$ 成正比。热力学只能把比例常量留给实验；统计物理给出了该常量。

谱能量密度极大处满足

$$
\frac{\mathrm d}{\mathrm dx}\left(\frac{x^3}{e^x-1}\right)=0.
$$

求导后

$$
3-3e^{-x}=x.
$$

数值解为

$$
x=\frac{\hbar\omega_m}{kT}\approx2.822. \tag{8.4.11}
$$

因此极大频率 $\omega_m\propto T$，即**维恩位移定律**。

### 8.4.4 巨配分函数与熵

光子的巨配分函数对数为

$$
\begin{aligned}
\ln\Xi
&=-\sum_l\omega_l
\ln\left(1-e^{-\alpha-\beta\varepsilon_l}\right)\\
&=-\frac{V}{\pi^2c^3}
\int_0^\infty\omega^2
\ln\left(1-e^{-\beta\hbar\omega}\right)\mathrm d\omega.
\end{aligned}
$$

令 $x=\beta\hbar\omega$：

$$
\ln\Xi
=-\frac{V}{\pi^2c^3}\frac1{(\beta\hbar)^3}
\int_0^\infty x^2\ln(1-e^{-x})\mathrm dx. \tag{8.4.12}
$$

分部积分：

$$
\begin{aligned}
\int_0^\infty x^2\ln(1-e^{-x})\mathrm dx
&=\left[\frac{x^3}{3}\ln(1-e^{-x})\right]_0^\infty
-\frac13\int_0^\infty\frac{x^3\mathrm dx}{e^x-1}\\
&=0-\frac13\frac{\pi^4}{15}
=-\frac{\pi^4}{45}.
\end{aligned}
$$

故

$$
\ln\Xi
=\frac{\pi^2V}{45c^3}\frac1{(\beta\hbar)^3}. \tag{8.4.13}
$$

内能

$$
U=-\frac{\partial}{\partial\beta}\ln\Xi
=\frac{\pi^2k^4}{15c^3\hbar^3}VT^4, \tag{8.4.14}
$$

与式 (8.4.10) 一致。压强

$$
p=\frac1\beta\frac{\partial}{\partial V}\ln\Xi
=\frac{\pi^2k^4}{45c^3\hbar^3}T^4. \tag{8.4.15}
$$

比较得

$$
p=\frac13\frac UV. \tag{8.4.16}
$$

熵为

$$
\begin{aligned}
S
&=k\left(\ln\Xi-\beta\frac{\partial\ln\Xi}{\partial\beta}\right)
=k(\ln\Xi+\beta U)\\
&=\frac{4\pi^2k^4}{45c^3\hbar^3}T^3V.
\end{aligned}
\tag{8.4.17}
$$

$T\to0$ 时 $S\to0$，满足热力学第三定律。辐射通量密度

$$
J_u=\frac c4\frac UV
=\frac{\pi^2k^4}{60c^2\hbar^3}T^4. \tag{8.4.18}
$$

---

## §8.5　金属中的自由电子气体

### 8.5.1 强简并判据

金属中的价电子可近似为在金属内部自由运动的全同费米子。若按经典能量均分，每个电子应贡献 $3k/2$ 的热容，但实验上电子热容通常远小于离子振动热容。索末菲在 1928 年用费米分布解释了这一点。

以铜为例，电子数密度约为

$$
\begin{aligned}
n
&=\frac{8.9\times10^3\ \mathrm{kg\cdot m^{-3}}}
{63\ \mathrm{g\cdot mol^{-1}}}N_A\\
&=8.5\times10^{28}\ \mathrm{m^{-3}}.
\end{aligned}
$$

取 $m=9.1\times10^{-31}\ \mathrm{kg}$，在 $T=300\ \mathrm K$ 时

$$
n\lambda^3
=\frac NV\left(\frac{h^2}{2\pi mkT}\right)^{3/2}
=\frac{3.53\times10^7\ \mathrm{K^{3/2}}}{T^{3/2}}
\approx6800\gg1.
$$

因此自由电子形成强简并费米气体。

费米分布中，单个量子态的平均占据数为

$$
f=\frac1{e^{(\varepsilon-\mu)/(kT)}+1}. \tag{8.5.1}
$$

考虑到自旋投影有两个值，能量在 $\varepsilon$ 到 $\varepsilon+\mathrm d\varepsilon$ 内的量子态数为

$$
D(\varepsilon)\mathrm d\varepsilon
=\frac{4\pi V}{h^3}(2m)^{3/2}\varepsilon^{1/2}\mathrm d\varepsilon.
$$

所以该能壳内平均电子数为

$$
\frac{4\pi V}{h^3}(2m)^{3/2}
\frac{\varepsilon^{1/2}\mathrm d\varepsilon}
{e^{(\varepsilon-\mu)/(kT)}+1}. \tag{8.5.2}
$$

给定 $N,T,V$，化学势由

$$
\frac{4\pi V}{h^3}(2m)^{3/2}
\int_0^\infty\frac{\varepsilon^{1/2}\mathrm d\varepsilon}
{e^{(\varepsilon-\mu)/(kT)}+1}
=N \tag{8.5.3}
$$

确定。

### 8.5.2 绝对零度的分布

$T\to0$ 时：

$$
f=
\begin{cases}
1,&\varepsilon<\mu(0),\\
0,&\varepsilon>\mu(0).
\end{cases} \tag{8.5.4}
$$

泡利不相容原理允许每个量子态最多容纳一个电子。绝对零度下电子从最低能级依次填充到 $\mu(0)$，故 $\mu(0)$ 是 $0\ \mathrm K$ 时电子的最大能量：

$$
\frac{4\pi V}{h^3}(2m)^{3/2}
\int_0^{\mu(0)}\varepsilon^{1/2}\mathrm d\varepsilon=N. \tag{8.5.5}
$$

完成积分：

$$
\frac{4\pi V}{h^3}(2m)^{3/2}
\frac23\mu(0)^{3/2}=N,
$$

所以

$$
\mu(0)
=\frac{\hbar^2}{2m}\left(3\pi^2\frac NV\right)^{2/3}. \tag{8.5.6}
$$

$\mu(0)$ 称为**费米能**。令 $\mu(0)=p_F^2/(2m)$，则

$$
p_F=(3\pi^2n)^{1/3}\hbar,
$$

$p_F$ 和 $v_F=p_F/m$ 分别称为**费米动量**和**费米速率**。对铜，$\mu(0)=1.12\times10^{-18}\ \mathrm J=7.0\ \mathrm{eV}$。定义费米温度

$$
kT_F=\mu(0),
$$

得 $T_F\approx8.1\times10^4\ \mathrm K\gg300\ \mathrm K$。

$0\ \mathrm K$ 时内能为

$$
\begin{aligned}
U(0)
&=\frac{4\pi V}{h^3}(2m)^{3/2}
\int_0^{\mu(0)}\varepsilon^{3/2}\mathrm d\varepsilon\\
&=\frac35N\mu(0).
\end{aligned}
\tag{8.5.7}
$$

因此每个电子的平均能量为 $3\mu(0)/5$。相应的简并压为

$$
p(0)=\frac23\frac{U(0)}V
=\frac25n\mu(0). \tag{8.5.8}
$$

铜中该压强约为 $3.8\times10^{10}\ \mathrm{Pa}$，由电子与离子的静电吸引所平衡。

### 8.5.3 有限温度下的费米面

由式 (8.5.1)：

$$
\begin{cases}
f>\dfrac12,&\varepsilon<\mu,\\
f=\dfrac12,&\varepsilon=\mu,\\
f<\dfrac12,&\varepsilon>\mu.
\end{cases} \tag{8.5.9}
$$

只有能量在 $\mu$ 附近、宽度约为 $kT$ 的壳层内的占据情况与 $0\ \mathrm K$ 明显不同。这个边界称为**费米面**。$kT\ll\mu(0)$ 时，$T\ll T_F$，仍是强简并费米气体。

对热容作粗略估计。对热容有贡献的有效电子数约为

$$
N_{\text{有效}}\approx\frac{kT}{\mu}N.
$$

若每个有效电子贡献约 $3k/2$，则

$$
C_V\approx\frac32Nk\frac{kT}{\mu}
=\frac32Nk\frac T{T_F}. \tag{8.5.10}
$$

在室温下 $T/T_F\approx1/270$，所以电子热容远小于经典值。低温下离子振动热容以 $T^3$ 衰减，电子热容以 $T$ 衰减，后者最终占优。

### 8.5.4 索末菲展开与低温热容

由式 (8.5.3) 与能量积分，

$$
N=\frac{4\pi V}{h^3}(2m)^{3/2}
\int_0^\infty\frac{\varepsilon^{1/2}\mathrm d\varepsilon}
{e^{(\varepsilon-\mu)/(kT)}+1}, \tag{8.5.11}
$$

$$
U=\frac{4\pi V}{h^3}(2m)^{3/2}
\int_0^\infty\frac{\varepsilon^{3/2}\mathrm d\varepsilon}
{e^{(\varepsilon-\mu)/(kT)}+1}. \tag{8.5.12}
$$

两式都有形式

$$
I=\int_0^\infty\frac{\eta(\varepsilon)\mathrm d\varepsilon}
{e^{(\varepsilon-\mu)/(kT)}+1}, \tag{8.5.13}
$$

其中 $\eta(\varepsilon)=C\varepsilon^{1/2}$ 或 $C\varepsilon^{3/2}$，且

$$
C=\frac{4\pi V}{h^3}(2m)^{3/2}.
$$

令 $\varepsilon-\mu=kTx$：

$$
\begin{aligned}
I
&=kT\int_{-\mu/(kT)}^\infty
\frac{\eta(\mu+kTx)}{e^x+1}\mathrm dx\\
&=kT\int_0^{\mu/(kT)}
\frac{\eta(\mu-kTx)}{e^{-x}+1}\mathrm dx
+kT\int_0^\infty
\frac{\eta(\mu+kTx)}{e^x+1}\mathrm dx.
\end{aligned}
$$

对第一项使用

$$
\frac1{e^{-x}+1}=1-\frac1{e^x+1},
$$

得到

$$
I=\int_0^\mu\eta(\varepsilon)\mathrm d\varepsilon
+kT\int_0^\infty
\frac{\eta(\mu+kTx)-\eta(\mu-kTx)}{e^x+1}\mathrm dx.
$$

上部积分上限已取为 $\infty$，因为 $\mu/(kT)\gg1$ 且费米因子使主要贡献来自 $x\sim O(1)$。将分子作泰勒展开：

$$
\eta(\mu+kTx)-\eta(\mu-kTx)
=2\eta'(\mu)kTx+\cdots.
$$

因此

$$
\begin{aligned}
I
&=\int_0^\mu\eta(\varepsilon)\mathrm d\varepsilon
+2(kT)^2\eta'(\mu)
\int_0^\infty\frac{x\mathrm dx}{e^x+1}+\cdots\\
&=\int_0^\mu\eta(\varepsilon)\mathrm d\varepsilon
+\frac{\pi^2}{6}(kT)^2\eta'(\mu)+\cdots.
\end{aligned}
\tag{8.5.14}
$$

其中用到

$$
\int_0^\infty\frac{x\mathrm dx}{e^x+1}=\frac{\pi^2}{12}.
$$

代入 $\eta=C\varepsilon^{1/2}$：

$$
N=\frac C3\mu^{3/2}
+\frac{\pi^2}{6}(kT)^2\frac C2\mu^{-1/2}
=\frac23C\mu^{3/2}
\left[1+\frac{\pi^2}{8}\left(\frac{kT}{\mu}\right)^2\right]. \tag{8.5.15}
$$

代入 $\eta=C\varepsilon^{3/2}$：

$$
U=\frac C4\mu^{5/2}
+\frac{\pi^2}{6}(kT)^2\frac{3C}{2}\mu^{1/2}
=\frac25C\mu^{5/2}
\left[1+\frac{5\pi^2}{8}\left(\frac{kT}{\mu}\right)^2\right]. \tag{8.5.16}
$$

由式 (8.5.15)：

$$
\mu=\left(\frac{3N}{2C}\right)^{2/3}
\left[1+\frac{\pi^2}{8}\left(\frac{kT}{\mu}\right)^2\right]^{-2/3}.
$$

修正项中可作 $\mu\approx\mu(0)$ 的零级替换，并展开到 $(kT/\mu(0))^2$：

$$
\begin{aligned}
\mu
&\approx\mu(0)
\left\{1+\frac{\pi^2}{8}
\left[\frac{kT}{\mu(0)}\right]^2\right\}^{-2/3}\\
&\approx\mu(0)
\left\{1-\frac{\pi^2}{12}
\left[\frac{kT}{\mu(0)}\right]^2\right\}.
\end{aligned}
\tag{8.5.17}
$$

代回式 (8.5.16)：

$$
\begin{aligned}
U
&=\frac25C\mu(0)^{5/2}
\left\{1-\frac{\pi^2}{12}
\left[\frac{kT}{\mu(0)}\right]^2\right\}^{5/2}
\left\{1+\frac{5\pi^2}{8}
\left[\frac{kT}{\mu(0)}\right]^2\right\}\\
&=\frac35N\mu(0)
\left\{1+\frac{5\pi^2}{12}
\left[\frac{kT}{\mu(0)}\right]^2\right\}.
\end{aligned}
\tag{8.5.18}
$$

因此

$$
\begin{aligned}
C_V
&=\left(\frac{\partial U}{\partial T}\right)_V\\
&=Nk\frac{\pi^2}{2}\frac{kT}{\mu(0)}
=\gamma_0T.
\end{aligned}
\tag{8.5.19}
$$

低温下金属的定容热容为电子与离子振动两部分之和：

$$
C_V=\gamma T+AT^3. \tag{8.5.20}
$$

作图时以 $C_V/T$ 对 $T^2$ 为坐标，应得到直线，纵轴截距为 $\gamma$。对铜，实验截距约为 $(0.688\pm0.002)\ \mathrm{mJ\cdot mol^{-1}\cdot K^{-2}}$，由式 (8.5.19) 算出的 $\gamma_0=0.50\ \mathrm{mJ\cdot mol^{-1}\cdot K^{-2}}$，存在一定差异。

---

## 本章小结

### 一、公式速查

| 内容 | 公式 |
| --- | --- |
| 非简并判据 | $e^\alpha\gg1$；$n\lambda^3\ll1$ |
| 巨配分函数（玻色 / 费米） | (8.1.2)(8.1.12) |
| 粒子数、内能、广义力 | $\bar N=-\partial_\alpha\ln\Xi$ (8.1.4)；$U=-\partial_\beta\ln\Xi$ (8.1.6)；$Y=-\beta^{-1}\partial_y\ln\Xi$ (8.1.7) |
| 压强、熵、巨热力势 | $p=\beta^{-1}\partial_V\ln\Xi$ (8.1.8)；$S=k(\ln\Xi+\alpha\bar N+\beta U)$ (8.1.10)；$J=-kT\ln\Xi$ (8.1.14) |
| 弱简并内能 | $U=\frac32NkT[1\pm n\lambda^3/(4\sqrt2g)]$ (8.2.8′) |
| 玻色化学势 | $\mu<0$；临界时 $\mu\to-0$ (8.3.3) |
| 凝聚温度与判据 | $T_c=\frac{2\pi}{(2.612)^{2/3}}\frac{h^2}{mk}n^{2/3}$ (8.3.8)；$n\lambda^3\geqslant2.612$ (8.3.15) |
| 凝聚比例、内能、热容 | (8.3.11)(8.3.12)(8.3.13) |
| 谐振子势中的凝聚 | (8.3.18)(8.3.19) |
| 光子分布 | $\mu=0$；$a_l=\omega_l/(e^{\beta\varepsilon_l}-1)$ (8.4.3) |
| 普朗克公式 | (8.4.7) |
| 维恩位移 | $\hbar\omega_m/(kT)\approx2.822$ (8.4.11) |
| 辐射热力学 | $U=\frac{\pi^2k^4}{15c^3\hbar^3}VT^4$ (8.4.14)；$p=U/(3V)$ (8.4.16)；$S=\frac{4\pi^2k^4}{45c^3\hbar^3}T^3V$ (8.4.17) |
| 费米分布 | $f=1/[e^{(\varepsilon-\mu)/(kT)}+1]$ (8.5.1) |
| 费米能、费米动量 | (8.5.6)；$p_F=(3\pi^2n)^{1/3}\hbar$ |
| 零温能量与压强 | $U(0)=\frac35N\mu(0)$ (8.5.7)；$p(0)=\frac25n\mu(0)$ (8.5.8) |
| 有限温度化学势 | (8.5.17) |
| 电子热容与金属低温热容 | (8.5.19)(8.5.20) |

### 二、概念脉络

```text
量子全同粒子
    │
    ├──► 玻色气体：nλ³ ≪ 1 时近似经典
    │       │
    │       ├──► 弱简并：U 修正为负
    │       │        （等效吸引）
    │       └──► 强简并：nλ³ ≥ 2.612
    │                │
    │                ├──► μ → 0⁻，基态宏观占据
    │                ├──► BEC：n₀/n = 1-(T/Tc)^(3/2)
    │                └──► U ∝ T^(5/2), C_V ∝ T^(3/2)
    │
    └──► 费米气体：每个量子态占据 0 或 1
            │
            ├──► 弱简并：U 修正为正
            │        （等效排斥）
            └──► 强简并：T ≪ T_F
                     │
                     ├──► 费米面，μ(T) 缓变
                     ├──► 零温仍有 U 和简并压
                     └──► 低温激发仅 kT 壳层
                          C_V ∝ T

巨配分函数 lnΞ
    ├──► N̄, U, Y, p, S, J 的统一统计表达
    ├──► 光子：μ=0，光子数不守恒
    │       ├──► Planck 公式
    │       ├──► U/V ∝ T⁴
    │       ├──► p=U/(3V)
    │       └──► Wien 位移
    └──► 自由电子：强简并 Fermi 气体
            └──► 金属低温 C_V = γT + AT³
```

### 三、易错点

1. **$\alpha$ 与 $\mu$ 的符号关系是 $\alpha=-\mu/(kT)$**：玻色分布的分母应写成 $e^{(\varepsilon-\mu)/(kT)}-1$，费米分布则分母为 $+1$。
2. **巨配分函数不是玻耳兹曼配分函数**：$\ln\Xi$ 的自然变量是 $T,V,\mu$；$J=-kT\ln\Xi$ 才是对应的巨热力势。
3. **弱简并展开中的正负号**：费米气体附加内能为正，玻色气体附加内能为负。不能把“$+$”简单理解成总取增强。
4. **玻色化学势必须有上限**：$\mu<\varepsilon_0$；取 $\varepsilon_0=0$ 后 $\mu<0$。$\mu\to-0$ 表示从负侧趋近，不是成为正值。
5. **$T_c$ 不是发生凝聚的唯一温度旋钮**：凝聚判据为 $n\lambda^3\geqslant2.612$，升温与降低密度都有利于离开凝聚区。
6. **把求和改成积分时必须单独处理基态**：$T<T_c$ 时基态占据宏观量级，态密度积分漏掉的正是 $n_0$；这正是式 (8.3.9) 的物理来源。
7. **凝聚体熵为零但不等于整个系统熵为零**：基态宏观占据只是粒子数现象；激发态仍可携带熵和热容。
8. **玻色-爱因斯坦凝聚不等于“粒子都静止不动”的经典液滴**：讨论的是理想量子气体在动量空间的基态占据；液氦超流还包含强相互作用。
9. **光子化学势为零源于光子数不守恒**：不要因为 $\mu=0$ 就把光子当成“没有能量”；能量分布仍满足普朗克公式。
10. **普朗克公式与态密度、玻色因子缺一不可**：$U(\omega)\mathrm d\omega$ 同时包含 $\omega^2$ 模式密度与 $1/(e^{\beta\hbar\omega}-1)$。
11. **瑞利-金斯与维恩公式只是极限**：低频极限丢失量子能隙，高频极限丢失受激玻色增强；完整公式是普朗克式。
12. **辐射压强不是 $p=U/V$**：光子气体由式 (8.4.16) 给出 $p=U/(3V)$，这个 1/3 来自三维无质量色散。
13. **费米能不是热能**：铜的 $T_F\sim10^5\ \mathrm K$，而室温 $T\ll T_F$；因此金属电子仍是强简并的。
14. **电子热容不是 $3Nk/2$**：只有费米面附近约 $kT/\mu$ 比例的电子可被热激发，故 $C_V\propto T$。
15. **费米气体的零温压强不为零**：泡利不相容原理迫使电子填充到费米能，产生简并压；不能把它误认为热压强。

### 四、名词中英对照

| 中文 | 英文 | 中文 | 英文 |
| --- | --- | --- | --- |
| 玻色统计 | Bose statistics | 费米统计 | Fermi statistics |
| 玻色-爱因斯坦分布 | Bose-Einstein distribution | 费米-狄拉克分布 | Fermi-Dirac distribution |
| 巨配分函数 | grand partition function | 巨热力势 | grand thermodynamic potential |
| 弱简并气体 | weakly degenerate gas | 强简并气体 | strongly degenerate gas |
| 玻色-爱因斯坦凝聚 | Bose-Einstein condensation | 凝聚温度 | condensation temperature |
| 凝聚体 | condensate | 基态占据 | ground-state occupation |
| 光子气体 | photon gas | 普朗克公式 | Planck formula |
| 瑞利-金斯公式 | Rayleigh-Jeans formula | 维恩公式 | Wien formula |
| 维恩位移定律 | Wien displacement law | 辐射通量密度 | radiative flux density |
| 自由电子气体 | free-electron gas | 费米面 | Fermi surface |
| 费米能 | Fermi energy | 费米动量 / 费米速率 | Fermi momentum / Fermi velocity |
| 费米温度 | Fermi temperature | 简并压 | degeneracy pressure |
| 索末菲展开 | Sommerfeld expansion | 电子热容系数 | electronic heat-capacity coefficient |
