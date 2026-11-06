---
tags:
  - 热统
  - 热力学
  - 第二章
created: 2026-09-24
---

# 第二章　均匀物质的热力学性质

> **教材**：汪志诚《热力学·统计物理》（第六版）
> **范围**：§2.1–§2.8（教材 pp. 42–61；其中 §2.8 为带星号的选读内容）

## 本章地图

| 节 | 标题 | 一句话 |
| --- | --- | --- |
| §2.1 | 内能、焓、自由能和吉布斯函数的全微分 | 四个全微分式与自然变量，导出麦克斯韦关系 |
| §2.2 | 麦克斯韦关系的简单应用 | 把不可直接测量的量用物态方程与热容表示 |
| §2.3 | 气体的节流过程和绝热膨胀过程 | 焦汤系数 $\mu$；两种获得低温的常用手段 |
| §2.4 | 基本热力学函数的确定 | 由物态方程 + 热容积分求 $U$、$S$、$H$、$G$ |
| §2.5 | 特性函数 | 马休定理：一个特性函数确定全部热力学函数 |
| §2.6 | 热辐射的热力学理论 | $u=aT^4$；斯特藩-玻耳兹曼定律；基尔霍夫定律 |
| §2.7 | 磁介质的热力学 | 代换 $p\to-\mu_0\mathcal{H}$、$V\to\mathfrak{m}$；绝热去磁制冷 |
| *§2.8 | 获得低温的方法（选读） | 节流、绝热膨胀、磁冷却、稀释制冷、激光制冷 |

**主线**：§2.1 导出全部微分关系 → §2.2 用它把不可测量化为可测量 → §2.3 应用（节流、绝热膨胀）→ §2.4–§2.5 由少量实验数据确定全部热力学函数 → §2.6–§2.7 把这套方法推广到辐射场与磁介质。

---

## §2.1　内能、焓、自由能和吉布斯函数的全微分

### 2.1.1 内能的全微分

第一章已引入物态方程、内能和熵，并导出热力学基本方程

$$\mathrm{d}U = T\,\mathrm{d}S - p\,\mathrm{d}V \tag{2.1.1}$$

它给出相邻两平衡态的内能、熵、体积之间的关系，可看作 $U$ 作为 $S$、$V$ 的函数时的全微分。$U$ 作为 $S$、$V$ 的函数，其全微分为

$$\mathrm{d}U = \left(\frac{\partial U}{\partial S}\right)_V \mathrm{d}S + \left(\frac{\partial U}{\partial V}\right)_S \mathrm{d}V$$

与 (2.1.1) 比较，可知

$$\left(\frac{\partial U}{\partial S}\right)_V = T, \qquad \left(\frac{\partial U}{\partial V}\right)_S = -p \tag{2.1.2}$$

考虑到求偏导数的次序可以交换，即 $\dfrac{\partial^2 U}{\partial V\partial S} = \dfrac{\partial^2 U}{\partial S\partial V}$，可得

$$\left(\frac{\partial T}{\partial V}\right)_S = -\left(\frac{\partial p}{\partial S}\right)_V \tag{2.1.3}$$

### 2.1.2 焓的全微分

焓的定义是 $H = U + pV$。对 $H$ 求微分，并将 (2.1.1) 代入，得

$$\mathrm{d}H = T\,\mathrm{d}S + V\,\mathrm{d}p \tag{2.1.4}$$

$H$ 作为 $S$、$p$ 的函数，其全微分为 $\mathrm{d}H = \left(\dfrac{\partial H}{\partial S}\right)_p \mathrm{d}S + \left(\dfrac{\partial H}{\partial p}\right)_S \mathrm{d}p$，与 (2.1.4) 比较得

$$\left(\frac{\partial H}{\partial S}\right)_p = T, \qquad \left(\frac{\partial H}{\partial p}\right)_S = V \tag{2.1.5}$$

再由二阶偏导数可交换，易得

$$\left(\frac{\partial T}{\partial p}\right)_S = \left(\frac{\partial V}{\partial S}\right)_p \tag{2.1.6}$$

### 2.1.3 自由能的全微分

自由能的定义是 $F = U - TS$。对 $F$ 求微分，并将 (2.1.1) 代入，得

$$\mathrm{d}F = -S\,\mathrm{d}T - p\,\mathrm{d}V \tag{2.1.7}$$

类似可得

$$\left(\frac{\partial F}{\partial T}\right)_V = -S, \qquad \left(\frac{\partial F}{\partial V}\right)_T = -p \tag{2.1.8}$$

以及

$$\left(\frac{\partial S}{\partial V}\right)_T = \left(\frac{\partial p}{\partial T}\right)_V \tag{2.1.9}$$

### 2.1.4 吉布斯函数的全微分

吉布斯函数的定义是 $G = U - TS + pV$。对 $G$ 求微分，并将 (2.1.1) 代入，得

$$\mathrm{d}G = -S\,\mathrm{d}T + V\,\mathrm{d}p \tag{2.1.10}$$

类似可得

$$\left(\frac{\partial G}{\partial T}\right)_p = -S, \qquad \left(\frac{\partial G}{\partial p}\right)_T = V \tag{2.1.11}$$

以及

$$\left(\frac{\partial S}{\partial p}\right)_T = -\left(\frac{\partial V}{\partial T}\right)_p \tag{2.1.12}$$

### 2.1.5 两组公式的性质

| 热力学函数 | 自然变量 | 全微分 | 偏导数关系 |
| --- | --- | --- | --- |
| 内能 $U$ | $S$、$V$ | $\mathrm{d}U = T\mathrm{d}S - p\mathrm{d}V$ (2.1.1) | (2.1.2) |
| 焓 $H = U + pV$ | $S$、$p$ | $\mathrm{d}H = T\mathrm{d}S + V\mathrm{d}p$ (2.1.4) | (2.1.5) |
| 自由能 $F = U - TS$ | $T$、$V$ | $\mathrm{d}F = -S\mathrm{d}T - p\mathrm{d}V$ (2.1.7) | (2.1.8) |
| 吉布斯函数 $G = U - TS + pV$ | $T$、$p$ | $\mathrm{d}G = -S\mathrm{d}T + V\mathrm{d}p$ (2.1.10) | (2.1.11) |

$U(S,V)$、$H(S,p)$、$F(T,V)$、$G(T,p)$ 是 §2.5 要讲的**特性函数**的例子，其自变量称为各特性函数的**自然变量**。

- 式 (2.1.2)、(2.1.5)、(2.1.8)、(2.1.11) 四式将 $S$、$T$、$p$、$V$ 四个变量用热力学函数 $U$、$H$、$F$、$G$ 的偏导数表达出来；
- 式 (2.1.3)、(2.1.6)、(2.1.9)、(2.1.12) 则给出 $S$、$T$、$p$、$V$ 这四个变量的偏导数之间的关系。

> **上述两组关系都来源于热力学基本方程 (2.1.1)，是根据热力学基本方程的全微分性质得到的。** 本章将利用它们通过数学推演得出简单系统平衡性质的关系，并导出简单系统热力学函数的一般表达式。所得的热力学关系非常普遍，适用于处于平衡态的任何简单系统；对于其他热力学均匀系统，通过直接代换也可得到类似关系（参阅 §2.7）。

---

## §2.2　麦克斯韦关系的简单应用

### 2.2.1 麦克斯韦关系

上节导出了 $S$、$T$、$p$、$V$ 四个变量的偏导数之间的关系：

$$\left(\frac{\partial T}{\partial V}\right)_S = -\left(\frac{\partial p}{\partial S}\right)_V \tag{2.2.1}$$

$$\left(\frac{\partial T}{\partial p}\right)_S = \left(\frac{\partial V}{\partial S}\right)_p \tag{2.2.2}$$

$$\left(\frac{\partial S}{\partial V}\right)_T = \left(\frac{\partial p}{\partial T}\right)_V \tag{2.2.3}$$

$$\left(\frac{\partial S}{\partial p}\right)_T = -\left(\frac{\partial V}{\partial T}\right)_p \tag{2.2.4}$$

以上偏导数关系是由麦克斯韦（Maxwell）最先导出来的，称为**麦克斯韦关系**。利用麦克斯韦关系，**可以把一些不能直接通过实验测量的物理量以物态方程（或 $\alpha$、$\kappa_T$）和热容等可以直接通过实验测量的物理量表达出来**。

### 2.2.2 内能随体积的变化

选 $T$、$V$ 为状态参量，$U$ 的全微分为 $\mathrm{d}U = \left(\dfrac{\partial U}{\partial T}\right)_V \mathrm{d}T + \left(\dfrac{\partial U}{\partial V}\right)_T \mathrm{d}V$；又 $\mathrm{d}U = T\mathrm{d}S - p\mathrm{d}V$，而 $T$、$V$ 为自变量时 $\mathrm{d}S = \left(\dfrac{\partial S}{\partial T}\right)_V \mathrm{d}T + \left(\dfrac{\partial S}{\partial V}\right)_T \mathrm{d}V$，代入可得

$$\mathrm{d}U = T\left(\frac{\partial S}{\partial T}\right)_V \mathrm{d}T + \left[T\left(\frac{\partial S}{\partial V}\right)_T - p\right]\mathrm{d}V$$

比较得

$$C_V = \left(\frac{\partial U}{\partial T}\right)_V = T\left(\frac{\partial S}{\partial T}\right)_V \tag{2.2.5}$$

及

$$\left(\frac{\partial U}{\partial V}\right)_T = T\left(\frac{\partial S}{\partial V}\right)_T - p \tag{2.2.6}$$

（(2.2.5) 给出定容热容的另一表达式。）将 (2.2.3) 代入 (2.2.6) 得

$$\left(\frac{\partial U}{\partial V}\right)_T = T\left(\frac{\partial p}{\partial T}\right)_V - p \tag{2.2.7}$$

**式 (2.2.7) 给出在温度保持不变时内能随体积的变化率与物态方程的关系。**

例如：

- **理想气体** $pV_m = RT$，由 (2.2.7) 得 $\left(\dfrac{\partial U_m}{\partial V_m}\right)_T = 0$——这正是焦耳定律的结果。
- **范德瓦耳斯气体** $\left(p + \dfrac{a}{V_m^2}\right)(V_m - b) = RT$，由 (2.2.7) 得

$$\left(\frac{\partial U_m}{\partial V_m}\right)_T = \frac{RT}{V_m - b} - p = \frac{a}{V_m^2}$$

这是在温度保持不变时范德瓦耳斯气体的内能随体积的变化率。

### 2.2.3 焓随压强的变化

选 $T$、$p$ 为独立变量，$\mathrm{d}H = \left(\dfrac{\partial H}{\partial T}\right)_p \mathrm{d}T + \left(\dfrac{\partial H}{\partial p}\right)_T \mathrm{d}p$；又 $\mathrm{d}H = T\mathrm{d}S + V\mathrm{d}p$，以 $T$、$p$ 为自变量时 $\mathrm{d}S = \left(\dfrac{\partial S}{\partial T}\right)_p \mathrm{d}T + \left(\dfrac{\partial S}{\partial p}\right)_T \mathrm{d}p$，代入可得

$$\mathrm{d}H = T\left(\frac{\partial S}{\partial T}\right)_p \mathrm{d}T + \left[T\left(\frac{\partial S}{\partial p}\right)_T + V\right]\mathrm{d}p$$

比较得

$$C_p = \left(\frac{\partial H}{\partial T}\right)_p = T\left(\frac{\partial S}{\partial T}\right)_p \tag{2.2.8}$$

$$\left(\frac{\partial H}{\partial p}\right)_T = T\left(\frac{\partial S}{\partial p}\right)_T + V \tag{2.2.9}$$

（(2.2.8) 给出定压热容的另一表达式。）将 (2.2.4) 代入 (2.2.9)：

$$\left(\frac{\partial H}{\partial p}\right)_T = V - T\left(\frac{\partial V}{\partial T}\right)_p \tag{2.2.10}$$

**式 (2.2.10) 给出在温度保持不变时焓随压强的变化率与物态方程的关系。**

### 2.2.4 任意简单系统的 $C_p - C_V$

由 (2.2.5) 和 (2.2.8)：

$$C_p - C_V = T\left(\frac{\partial S}{\partial T}\right)_p - T\left(\frac{\partial S}{\partial T}\right)_V$$

而由函数关系 $S(T,p) = S[T, V(T,p)]$ 可得

$$\left(\frac{\partial S}{\partial T}\right)_p = \left(\frac{\partial S}{\partial T}\right)_V + \left(\frac{\partial S}{\partial V}\right)_T \left(\frac{\partial V}{\partial T}\right)_p$$

因此

$$C_p - C_V = T\left(\frac{\partial S}{\partial V}\right)_T \left(\frac{\partial V}{\partial T}\right)_p$$

利用麦克斯韦关系 (2.2.3)，可写成

$$C_p - C_V = T\left(\frac{\partial p}{\partial T}\right)_V \left(\frac{\partial V}{\partial T}\right)_p \tag{2.2.11}$$

**上式给出两热容之差与物态方程的关系。** 例如对理想气体可得 $C_p - C_V = nR$，即 (1.7.8)。

> **注意**：§1.7 中根据热力学第一定律只求得了理想气体的 $C_p$ 与 $C_V$ 之差，而 (2.2.11) 适用于**任意的简单系统**。

也可以将 (2.2.11) 表示为

$$C_p - C_V = \frac{VT\alpha^2}{\kappa_T} \tag{2.2.12}$$

式中 $\alpha$ 和 $\kappa_T$ 分别为体胀系数和等温压缩系数。上式右方不可能取负值，因此**恒有 $C_p - C_V \geqslant 0$**。水的密度在 4 °C 时具有极大值，此时 $\alpha = 0$，因此在 4 °C 时水的 $C_p = C_V$。

实验中难以测量固体和液体的定容热容，但可以根据 (2.2.12) 由定压热容及 $\alpha$、$\kappa_T$ 计算出来。例如，在 0 °C 和 1 atm 下，测得水银的下列数据：$C_p = 28.0\ \mathrm{J\cdot K^{-1}}$，$V = 1.47\times 10^{-5}\ \mathrm{m^3}$，$\alpha = 181\times 10^{-6}\ \mathrm{K^{-1}}$，$\kappa_T = 3.94\times 10^{-6}\ \mathrm{atm^{-1}} = 3.89\times 10^{-11}\ \mathrm{Pa^{-1}}$，$T = 273\ \mathrm{K}$，代入 (2.2.12) 即求得 $C_V = 24.6\ \mathrm{J\cdot K^{-1}}$。

### 2.2.5 雅可比行列式的应用

在热力学中往往要进行导数变换的运算，雅可比（Jacobi）行列式是进行导数变换运算的一个有用的工具（雅可比行列式的性质可参阅附录 A）。举两个例子。

**例 2.2.1**　求证：绝热压缩系数 $\kappa_S$ 与等温压缩系数 $\kappa_T$ 之比等于定容热容与定压热容之比。

**证明**　$\kappa_S$ 和 $\kappa_T$ 的定义分别是

$$\kappa_S = -\frac{1}{V}\left(\frac{\partial V}{\partial p}\right)_S, \qquad \kappa_T = -\frac{1}{V}\left(\frac{\partial V}{\partial p}\right)_T \tag{2.2.13}$$

因此

$$\frac{\kappa_S}{\kappa_T} = \frac{-\frac{1}{V}\left(\frac{\partial V}{\partial p}\right)_S}{-\frac{1}{V}\left(\frac{\partial V}{\partial p}\right)_T} = \frac{\frac{\partial (V,S)}{\partial (p,S)}}{\frac{\partial (V,T)}{\partial (p,T)}} = \frac{\frac{\partial (V,S)}{\partial (V,T)}}{\frac{\partial (p,S)}{\partial (p,T)}} = \frac{\left(\frac{\partial S}{\partial T}\right)_V}{\left(\frac{\partial S}{\partial T}\right)_p} = \frac{C_V}{C_p} \tag{2.2.14}$$

**例 2.2.2**　求证：

$$C_p - C_V = -T\frac{\left(\frac{\partial p}{\partial T}\right)_V^2}{\left(\frac{\partial p}{\partial V}\right)_T} \tag{2.2.15}$$

**证明**

$$C_p = T\left(\frac{\partial S}{\partial T}\right)_p = T\frac{\partial (S,p)}{\partial (T,p)} = T\frac{\frac{\partial (S,p)}{\partial (T,V)}}{\frac{\partial (T,p)}{\partial (T,V)}}$$

$$= T\frac{\left(\frac{\partial S}{\partial T}\right)_V\left(\frac{\partial p}{\partial V}\right)_T - \left(\frac{\partial S}{\partial V}\right)_T\left(\frac{\partial p}{\partial T}\right)_V}{\left(\frac{\partial p}{\partial V}\right)_T} = C_V - T\frac{\left(\frac{\partial p}{\partial T}\right)_V^2}{\left(\frac{\partial p}{\partial V}\right)_T}$$

**补充（麦克斯韦关系的雅可比行列式形式）**

> **性质说明**：教材 §2.2.5 只给出上面两个应用例题，并没有把麦克斯韦关系本身写成雅可比行列式。下面的统一形式由附录 A 的性质直接推出，属笔记补充内容，故不占教材式号。

所引附录 A 的性质：

$$\frac{\partial (u,v)}{\partial (x,y)} = \begin{vmatrix} \dfrac{\partial u}{\partial x} & \dfrac{\partial u}{\partial y} \\[2mm] \dfrac{\partial v}{\partial x} & \dfrac{\partial v}{\partial y} \end{vmatrix} = \frac{\partial u}{\partial x}\frac{\partial v}{\partial y} - \frac{\partial u}{\partial y}\frac{\partial v}{\partial x} \tag{A.10}$$

（式 (A.10) 即雅可比行列式的定义：设 $u$、$v$ 是独立变量 $x$、$y$ 的函数 $u=u(x,y)$、$v=v(x,y)$；下面出现的 $\dfrac{\partial (T,S)}{\partial (p,V)}$ 一类记号即按该定义理解。）

$$\left(\frac{\partial u}{\partial x}\right)_y = \frac{\partial (u,y)}{\partial (x,y)} \tag{A.11}$$

$$\frac{\partial (u,v)}{\partial (x,y)} = -\frac{\partial (v,u)}{\partial (x,y)} \tag{A.12}$$

$$\frac{\partial (u,v)}{\partial (x,y)} = \frac{\partial (u,v)}{\partial (r,s)}\frac{\partial (r,s)}{\partial (x,y)} \tag{A.13}$$

由 (A.11)，四个麦克斯韦关系 (2.2.1)–(2.2.4) 可以用**一个**雅可比行列式等式统一表示：

$$\frac{\partial (T,S)}{\partial (p,V)} = 1$$

**等价说法**：对任意一对独立变量 $(x,y)$，

$$\frac{\partial (T,S)}{\partial (x,y)} = \frac{\partial (p,V)}{\partial (x,y)}$$

（由 (A.13) 立得：$\dfrac{\partial (T,S)}{\partial (x,y)} = \dfrac{\partial (T,S)}{\partial (p,V)}\dfrac{\partial (p,V)}{\partial (x,y)}$。）也就是说，**$(T,S)$ 与 $(p,V)$ 在雅可比行列式中的地位完全对等，可以整对地互相替换**——这正是 $\mathrm{d}U = T\,\mathrm{d}S - p\,\mathrm{d}V$ 中 $S$ 与 $V$、$T$ 与 $-p$ 互为共轭变量的数学表述。

**等价性的验证**：分别取 $(x,y)=(S,V)$、$(S,p)$、$(T,V)$、$(T,p)$，按定义 (A.10) 展开，即逐一还原出四个麦克斯韦关系：

| $(x,y)$ | $\dfrac{\partial (T,S)}{\partial (x,y)}$ | $\dfrac{\partial (p,V)}{\partial (x,y)}$ | 还原出的关系 |
| --- | --- | --- | --- |
| $(S,V)$ | $-\left(\dfrac{\partial T}{\partial V}\right)_S$ | $\left(\dfrac{\partial p}{\partial S}\right)_V$ | (2.2.1) |
| $(S,p)$ | $-\left(\dfrac{\partial T}{\partial p}\right)_S$ | $-\left(\dfrac{\partial V}{\partial S}\right)_p$ | (2.2.2) |
| $(T,V)$ | $\left(\dfrac{\partial S}{\partial V}\right)_T$ | $\left(\dfrac{\partial p}{\partial T}\right)_V$ | (2.2.3) |
| $(T,p)$ | $\left(\dfrac{\partial S}{\partial p}\right)_T$ | $-\left(\dfrac{\partial V}{\partial T}\right)_p$ | (2.2.4) |

以第一行为例。取 $(x,y)=(S,V)$，由定义 (A.10)

$$\frac{\partial (T,S)}{\partial (S,V)} = \left(\frac{\partial T}{\partial S}\right)_V\left(\frac{\partial S}{\partial V}\right)_S - \left(\frac{\partial T}{\partial V}\right)_S\left(\frac{\partial S}{\partial S}\right)_V = -\left(\frac{\partial T}{\partial V}\right)_S$$

$$\frac{\partial (p,V)}{\partial (S,V)} = \left(\frac{\partial p}{\partial S}\right)_V\left(\frac{\partial V}{\partial V}\right)_S - \left(\frac{\partial p}{\partial V}\right)_S\left(\frac{\partial V}{\partial S}\right)_V = \left(\frac{\partial p}{\partial S}\right)_V$$

两者相等即给出 (2.2.1)；其余三行同法可得。

**用法示例**：把待求偏导数补成配对的雅可比行列式后，用 $\partial (T,S)=\partial (p,V)$ 整对替换，可省去逐项导数变换。例如重新导出 (2.2.3)：

$$\left(\frac{\partial S}{\partial V}\right)_T = \frac{\partial (S,T)}{\partial (V,T)} = -\frac{\partial (T,S)}{\partial (V,T)} = -\frac{\partial (p,V)}{\partial (V,T)} = \left(\frac{\partial p}{\partial T}\right)_V$$

其中第二步用 (A.12)，第三步对 $(x,y)=(V,T)$ 用上面的一般等式，最后一步用了 $\dfrac{\partial (p,V)}{\partial (V,T)} = -\left(\dfrac{\partial p}{\partial T}\right)_V$。

---

## §2.3　气体的节流过程和绝热膨胀过程

热力学中往往用偏导数描述一个物理效应。例如，可逆绝热过程中熵保持不变，可逆绝热过程中温度随压强的变化率用 $\left(\dfrac{\partial T}{\partial p}\right)_S$ 描述；绝热自由膨胀过程中内能保持不变，绝热自由膨胀过程中温度随体积的变化率用 $\left(\dfrac{\partial T}{\partial V}\right)_U$ 描述。

为了求出某一效应的变化率，可以将描述该效应的偏导数用 $C_p$、$\alpha$、$\kappa_T$ 表示出来，或者求出描述该效应的偏导数与描述另一效应的偏导数之间的关系。作为例子，本节讨论气体的**节流过程**和**绝热膨胀过程**——这两种过程都是获得低温的常用方法。

### 2.3.1 节流过程与焦汤系数

如图 2.3.1 所示，管子用不导热的材料包着，管子中间有一个多孔塞或节流阀。多孔塞两边各维持着较高的压强 $p_1$ 和较低的压强 $p_2$，于是气体从高压的一边经多孔塞不断地流到低压的一边，并达到定常状态。这个过程就叫做**节流过程**。

测量气体在多孔塞两边的温度表明，节流过程前后气体的温度发生了变化，该效应称为**焦耳-汤姆孙效应**（简称焦汤效应），是焦耳和汤姆孙（W. Thomson，后来被封为开尔文男爵）在 1852 年用多孔塞实验研究气体内能时发现的。

**理论分析**：设过程中有一定数量的气体通过了多孔塞。通过多孔塞前，其压强为 $p_1$、体积为 $V_1$、内能为 $U_1$；通过多孔塞后，压强为 $p_2$、体积为 $V_2$、内能为 $U_2$。在过程中，外界对这部分气体所做的功是 $p_1V_1 - p_2V_2$。因为过程是绝热的，根据热力学第一定律有

$$U_2 - U_1 = p_1 V_1 - p_2 V_2$$

即 $U_2 + p_2V_2 = U_1 + p_1V_1$，或

$$H_1 = H_2 \tag{2.3.1}$$

**在节流过程前后，气体的焓相等。**

定义

$$\mu = \left(\frac{\partial T}{\partial p}\right)_H \tag{2.3.2}$$

表示在焓不变的条件下气体温度随压强的变化率，称为**焦耳-汤姆孙系数**，简称**焦汤系数**。

取 $T$、$p$ 为状态参量，$H = H(T,p)$，其偏导数间存在下述关系：

$$\left(\frac{\partial T}{\partial p}\right)_H \left(\frac{\partial p}{\partial H}\right)_T \left(\frac{\partial H}{\partial T}\right)_p = -1$$

或

$$\left(\frac{\partial T}{\partial p}\right)_H = -\frac{\left(\frac{\partial H}{\partial p}\right)_T}{\left(\frac{\partial H}{\partial T}\right)_p} \tag{2.3.3}$$

将 (2.2.8) 和 (2.2.10) 代入，得

$$\mu = \left(\frac{\partial T}{\partial p}\right)_H = \frac{1}{C_p}\left[T\left(\frac{\partial V}{\partial T}\right)_p - V\right] \tag{2.3.4}$$

或

$$\mu = \frac{V}{C_p}(T\alpha - 1) \tag{2.3.5}$$

**式 (2.3.4) 和 (2.3.5) 给出了焦汤系数与物态方程和热容的关系。**

- **理想气体**：$\alpha = \dfrac{1}{T}$，所以 $\mu = 0$。即理想气体在节流过程前后温度不变。
- **实际气体**：若 $\alpha T > 1$，有 $\mu > 0$；若 $\alpha T < 1$，有 $\mu < 0$。

一般说来 $\alpha$ 是 $T$、$p$ 的函数，所以 $\alpha = \dfrac{1}{T}$ 对应于 $T$-$p$ 图上的一条曲线，称为**反转曲线**。曲线给出使 $\mu = 0$ 的温度（称为**反转温度**）与压强的关系。图 2.3.2 中的虚线是实测的 N₂ 的反转曲线：在曲线的一侧 $\mu > 0$，气体经节流过程后降温，称为**制冷区**；在曲线的另一侧 $\mu < 0$，气体在节流过程后升温，称为**制温区**。利用节流过程的降温效应可以使气体降温而液化。

### 2.3.2 由昂内斯方程分析焦汤效应

只保留到第二位力系数，昂内斯方程近似为

$$p = \frac{nRT}{V}\left[1 + \frac{n}{V}B(T)\right] \tag{2.3.6}$$

注意修正项 $\dfrac{n}{V}B(T)$ 远小于 1，可以把零级近似 $\dfrac{n}{V} = \dfrac{p}{RT}$ 代入该修正项，而将上式表示为

$$p = \frac{nRT}{V}\left(1 + \frac{p}{RT}B\right)$$

或

$$V = n\left(\frac{RT}{p} + B\right)$$

将根据上式算得的 $\alpha$ 代入 (2.3.5)，得

$$\mu = \frac{n}{C_p}\left(T\frac{\mathrm{d}B}{\mathrm{d}T} - B\right) \tag{2.3.7}$$

由图 1.3.1 可以看出，$T\dfrac{\mathrm{d}B}{\mathrm{d}T}$ 是正的。**在足够低的温度下**，分子间吸引力的影响显著，使 $B$ 取负值，因此 (2.3.7) 给出的 $\mu > 0$。**温度足够高时**，斥力的影响显著，使 $B$ 取正值，有可能使 $\mu < 0$。**反转温度的存在是分子间吸引力和斥力的影响相互竞争的表现。**

> 图 2.3.2 中的实线是等焓线。不过在节流过程中气体的状态并非沿等焓线变化——**节流过程是一个不可逆过程**，我们只是将节流过程前后达到定常状态的气体的初态和终态近似看作平衡态，利用其焓相等即 $H(T_1,p_1) = H(T_2,p_2)$ 的条件，确定压强降低 $p_1 - p_2$ 后所能获得的温度降低。

### 2.3.3 绝热膨胀

如果把过程近似地看作准静态，则在准静态绝热过程中气体的熵保持不变。由

$$\mathrm{d}S = \left(\frac{\partial S}{\partial T}\right)_p \mathrm{d}T + \left(\frac{\partial S}{\partial p}\right)_T \mathrm{d}p = 0$$

可得

$$\left(\frac{\partial T}{\partial p}\right)_S = -\frac{\left(\frac{\partial S}{\partial p}\right)_T}{\left(\frac{\partial S}{\partial T}\right)_p} = \frac{T}{C_p}\left(\frac{\partial V}{\partial T}\right)_p = \frac{VT\alpha}{C_p} \tag{2.3.8}$$

式 (2.3.8) 给出在准静态绝热过程中气体的温度随压强的变化率。**上式右方是恒正的**，所以随着体积膨胀、压强降低，气体的温度必然下降。

从能量转化的角度看：气体在绝热膨胀过程中减少其内能而对外做功；膨胀后气体分子间的平均距离增大，吸力的影响减弱而使分子间的相互作用能有所增加；内能减少，相互作用能又增加，分子的平均动能必然减少，因而气体的温度下降。

气体的绝热膨胀过程也被用来使气体降温并液化。§2.8 中将进一步讨论利用节流过程和绝热膨胀过程获得低温的问题。

---

## §2.4　基本热力学函数的确定

前面所引进的热力学函数中，最基本的是**物态方程、内能和熵**，其他热力学函数均可由这三个基本函数导出。本节导出简单系统的基本热力学函数的一般表达式，即这三个函数与状态参量的函数关系。

### 2.4.1 以 $T$、$V$ 为状态参量

物态方程为

$$p = p(T, V) \tag{2.4.1}$$

（前面已经说过，热力学中物态方程要由实验测定。）

根据 (2.2.5) 和 (2.2.7)，内能的全微分为

$$\mathrm{d}U = C_V\,\mathrm{d}T + \left[T\left(\frac{\partial p}{\partial T}\right)_V - p\right]\mathrm{d}V \tag{2.4.2}$$

沿一条任意的积分路线求积分，可得

$$U = \int\left\{C_V\,\mathrm{d}T + \left[T\left(\frac{\partial p}{\partial T}\right)_V - p\right]\mathrm{d}V\right\} + U_0 \tag{2.4.3}$$

（内能的积分表达式。）根据 (2.2.5) 和 (2.2.3)，熵的全微分为

$$\mathrm{d}S = \frac{C_V}{T}\mathrm{d}T + \left(\frac{\partial p}{\partial T}\right)_V \mathrm{d}V \tag{2.4.4}$$

求线积分得

$$S = \int\left[\frac{C_V}{T}\mathrm{d}T + \left(\frac{\partial p}{\partial T}\right)_V \mathrm{d}V\right] + S_0 \tag{2.4.5}$$

（熵的积分表达式。）

由 (2.4.3) 和 (2.4.5) 可知，**只要测得物质的 $C_V$ 和物态方程，即可求得物质的内能函数和熵函数**。还可以证明，只要测得在某一体积下的定容热容 $C_V^0$，则任意体积下的定容热容都可根据物态方程求出来（习题 2.9）。因此，只需知道物态方程和某一体积下的定容热容，就可以求得内能和熵。

### 2.4.2 以 $T$、$p$ 为状态参量

物态方程是

$$V = V(T, p) \tag{2.4.6}$$

关于内能函数，在选 $T$、$p$ 为独立变量时，先求焓较为方便。由 (2.2.8) 和 (2.2.10) 得焓的全微分

$$\mathrm{d}H = C_p\,\mathrm{d}T + \left[V - T\left(\frac{\partial V}{\partial T}\right)_p\right]\mathrm{d}p \tag{2.4.7}$$

求线积分，得

$$H = \int\left\{C_p\,\mathrm{d}T + \left[V - T\left(\frac{\partial V}{\partial T}\right)_p\right]\mathrm{d}p\right\} + H_0 \tag{2.4.8}$$

（焓的积分表达式。）再由 $U = H - pV$ 即可求得内能。

关于熵函数，由 (2.2.8) 和 (2.2.4) 得熵的全微分

$$\mathrm{d}S = \frac{C_p}{T}\mathrm{d}T - \left(\frac{\partial V}{\partial T}\right)_p \mathrm{d}p \tag{2.4.9}$$

求线积分得

$$S = \int\left[\frac{C_p}{T}\mathrm{d}T - \left(\frac{\partial V}{\partial T}\right)_p \mathrm{d}p\right] + S_0 \tag{2.4.10}$$

（熵的积分表达式。）

由 (2.4.8) 和 (2.4.10) 可知，**只要测得物质的定压热容 $C_p$ 和物态方程，即可求得物质的内能函数（先求焓）和熵函数**。还可以证明，只要测得某一压强下的定压热容 $C_p^0$，则任意压强下的 $C_p$ 都可根据物态方程求出来（习题 2.9）。因此，只需知道物态方程和某一压强下的定压热容，就可以求得内能和熵。

对于固体和液体，定容热容在实验上难以直接测定，选 $T$、$p$ 为自变量比较方便。

### 2.4.3 例 2.4.1：理想气体的焓、熵和吉布斯函数

以 $T$、$p$ 为状态参量，1 mol 理想气体的物态方程为 $pV_m = RT$。由物态方程得

$$\left(\frac{\partial V_m}{\partial T}\right)_p = \frac{R}{p}, \qquad V_m - T\left(\frac{\partial V_m}{\partial T}\right)_p = 0$$

代入 (2.4.8)，得理想气体的摩尔焓为

$$H_m = \int C_{p,m}\,\mathrm{d}T + H_{m0} \tag{2.4.11}$$

如果摩尔定压热容 $C_{p,m}$ 可以看作常量，则有 $H_m = C_{p,m}T + H_{m0}$ 　(2.4.11')

将 $\left(\dfrac{\partial V_m}{\partial T}\right)_p = \dfrac{R}{p}$ 代入 (2.4.10)，得理想气体的摩尔熵为

$$S_m = \int\frac{C_{p,m}}{T}\mathrm{d}T - R\ln p + S_{m0} \tag{2.4.12}$$

如果 $C_{p,m}$ 可以看作常量，则有 $S_m = C_{p,m}\ln T - R\ln p + S_{m0}$ 　(2.4.12')

式 (2.4.11) 和式 (2.4.12) 就是式 (1.7.7) 和式 (1.15.6)。

根据 (1.18.7)，摩尔吉布斯函数为 $G_m = H_m - TS_m$。将 (2.4.11) 和 (2.4.12) 代入，得

$$G_m = \int C_{p,m}\mathrm{d}T - T\int\frac{C_{p,m}}{T}\mathrm{d}T + RT\ln p + H_{m0} - TS_{m0} \tag{2.4.13}$$

如果热容可以看作常量，则有

$$G_m = C_{p,m}T - C_{p,m}T\ln T + RT\ln p + H_{m0} - TS_{m0} \tag{2.4.13'}$$

式 (2.4.13) 可以表示为另一形式。利用分部积分公式 $\displaystyle\int x\,\mathrm{d}y = xy - \int y\,\mathrm{d}x$，令其中的 $x = \dfrac{1}{T}$、$y = \displaystyle\int C_{p,m}\mathrm{d}T$，(2.4.13) 可表示为

$$G_m = -T\int\frac{\mathrm{d}T}{T^2}\int C_{p,m}\mathrm{d}T + RT\ln p + H_{m0} - TS_{m0} \tag{2.4.14}$$

通常将 $G_m$ 写成

$$G_m = RT(\varphi + \ln p) \tag{2.4.15}$$

其中 $\varphi$ 是温度的函数：

$$\varphi = \frac{H_{m0}}{RT} - \int\frac{\mathrm{d}T}{RT^2}\int C_{p,m}\mathrm{d}T - \frac{S_{m0}}{R} \tag{2.4.16}$$

如果热容可以看作常量，则有

$$\varphi = \frac{H_{m0}}{RT} - \frac{C_{p,m}\ln T}{R} + \frac{C_{p,m} - S_{m0}}{R} \tag{2.4.16'}$$

以后我们要用到上述理想气体热力学函数的表达式，特别是式 (2.4.15)。

### 2.4.4 例 2.4.2：范德瓦耳斯气体的内能和熵

1 mol 范德瓦耳斯气体的物态方程为 $\left(p + \dfrac{a}{V_m^2}\right)(V_m - b) = RT$。由范德瓦耳斯方程可得

$$\left(\frac{\partial p}{\partial T}\right)_V = \frac{R}{V_m - b}, \qquad T\left(\frac{\partial p}{\partial T}\right)_V - p = \frac{a}{V_m^2}$$

代入 (2.4.3) 和 (2.4.5)，并注意范德瓦耳斯气体的 $C_V$ 只是 $T$ 的函数、与体积无关（习题 2.10），可分别得

$$U_m = \int C_{V,m}\mathrm{d}T - \frac{a}{V_m} + U_{m0} \tag{2.4.17}$$

$$S_m = \int\frac{C_{V,m}}{T}\mathrm{d}T + R\ln(V_m - b) + S_{m0} \tag{2.4.18}$$

### 2.4.5 例 2.4.3：简单固体的内能和熵

简单固体的物态方程为

$$V(T,p) = V_0(T_0,0)\left[1 + \alpha(T - T_0) - \kappa_T p\right]$$

引入符号 $V_1 = V_0 - \alpha V_0 T_0$，可将物态方程表示为

$$V = V_1 + V_0(\alpha T - \kappa_T p)$$

由此可得

$$\left(\frac{\partial p}{\partial T}\right)_V = \frac{\alpha}{\kappa_T}, \qquad T\left(\frac{\partial p}{\partial T}\right)_V - p = \frac{V - V_1}{\kappa_T V_0}$$

代入 (2.4.3) 和 (2.4.5)。注意在物态方程 (1.3.14) 中，$p$ 是 $T$ 的线性函数，因此简单固体的定容热容 $C_V$ 与体积无关、只是 $T$ 的函数（习题 2.9）。由此可得

$$U = \int C_V\,\mathrm{d}T + \frac{1}{2}\frac{(V - V_1)^2}{\kappa_T V_0} + U_0 \tag{2.4.19}$$

$$S = \int\frac{C_V}{T}\mathrm{d}T + \frac{\alpha}{\kappa_T}V + S_0 \tag{2.4.20}$$

---

## §2.5　特性函数

马休（Massieu）在 1869 年证明：**如果适当选择独立变量（称为自然变量），只要知道一个热力学函数，就可以通过求偏导数而求得均匀系统的全部热力学函数，从而把均匀系统的平衡性质完全确定。** 这个已知的热力学函数称为**特性函数**，表明它是表征均匀系统的特性的。

在 §2.1 中说过：内能 $U$ 作为 $S$、$V$ 的函数，焓 $H$ 作为 $S$、$p$ 的函数，自由能 $F$ 作为 $T$、$V$ 的函数，吉布斯函数 $G$ 作为 $T$、$p$ 的函数，都是特性函数。还有其他热力学函数也可以作为特性函数，本书不作讨论。

在应用上，最重要的特性函数是**自由能和吉布斯函数**。

### 2.5.1 自由能作为特性函数

式 (2.1.7) 给出了自由能的全微分表达式

$$\mathrm{d}F = -S\,\mathrm{d}T - p\,\mathrm{d}V \tag{2.5.1}$$

因此

$$S = -\frac{\partial F}{\partial T}, \qquad p = -\frac{\partial F}{\partial V} \tag{2.5.2}$$

如果已知 $F(T,V)$：

- 求 $F$ 对 $T$ 的偏导数即得出熵 $S(T,V)$；
- 求 $F$ 对 $V$ 的偏导数即得出压强 $p(T,V)$，这就是物态方程。

根据自由能的定义 $F = U - TS$，有

$$U = F + TS = F - T\frac{\partial F}{\partial T} \tag{2.5.3}$$

上式给出内能 $U(T,V)$。这样，三个基本的热力学函数便都可由 $F(T,V)$ 求出来了。**式 (2.5.3) 称为吉布斯-亥姆霍兹（Gibbs-Helmholtz）方程。**

在第七章玻耳兹曼统计和第九章正则系综理论中，就是用统计物理学方法求出自由能作为 $N$、$T$、$V$ 的函数，再进而求其他热力学函数的。

### 2.5.2 吉布斯函数作为特性函数

根据式 (2.1.10)，吉布斯函数的全微分为

$$\mathrm{d}G = -S\,\mathrm{d}T + V\,\mathrm{d}p \tag{2.5.4}$$

因此

$$S = -\frac{\partial G}{\partial T}, \qquad V = \frac{\partial G}{\partial p} \tag{2.5.5}$$

如果已知 $G(T,p)$，求 $G$ 对 $T$ 的偏导数即得出 $S(T,p)$；求 $G$ 对 $p$ 的偏导数即得出 $V(T,p)$，这就是物态方程。

由吉布斯函数的定义，有

$$U = G + TS - pV = G - T\frac{\partial G}{\partial T} - p\frac{\partial G}{\partial p} \tag{2.5.6}$$

上式给出 $U(T,p)$。这样，三个基本的热力学函数便都可由 $G(T,p)$ 求出来了。由焓的定义 $H = U + pV$，得

$$H = G - T\frac{\partial G}{\partial T} \tag{2.5.7}$$

**式 (2.5.7) 是吉布斯-亥姆霍兹方程的另一表达形式。**

### 2.5.3 例 2.5.1：表面系统的热力学函数

**解**　表面指液体和其他相的交界面。它实际上是很薄的一层，在垂直于分界面的方向上，表面的性质有急剧的变化。现在把它理想化为一个几何面，并把分界面两侧的两相都看成均匀的。

将表面看作一个热力学系统，描述表面系统的状态参量是表面张力系数 $\sigma$ 和表面积 $A$（相当于流体的 $p$ 和 $V$），表面系统的物态方程是 $\sigma$、$A$ 和 $T$ 的关系：

$$f(\sigma, A, T) = 0 \tag{2.5.8}$$

实验指出，**表面张力系数只是 $T$ 的函数，与表面积 $A$ 无关**。所以物态方程简化为 $\sigma = \sigma(T)$。

根据式 (1.4.4)，当表面的面积有 $\mathrm{d}A$ 的改变时，外界所做的功为 $\mathrm{d}W = \sigma\,\mathrm{d}A$。因此表面系统自由能的全微分为

$$\mathrm{d}F = -S\,\mathrm{d}T + \sigma\,\mathrm{d}A \tag{2.5.9}$$

由此得

$$S = -\frac{\partial F}{\partial T}, \qquad \sigma = \frac{\partial F}{\partial A} \tag{2.5.10}$$

将 (2.5.10) 的第二式积分，注意 $\sigma$ 与 $A$ 无关，即得

$$F = \sigma A \tag{2.5.11}$$

当 $A\to 0$ 时，表面系统就不再存在，其自由能也应为零。所以式 (2.5.11) 不含积分常量。**式 (2.5.11) 指出，$\sigma$ 是单位面积的自由能。**

将式 (2.5.11) 代入式 (2.5.10) 的第一式，得表面系统的熵为

$$S = -A\frac{\mathrm{d}\sigma}{\mathrm{d}T} \tag{2.5.12}$$

由 $U = F + TS$，得表面系统的内能为

$$U = A\left(\sigma - T\frac{\mathrm{d}\sigma}{\mathrm{d}T}\right) \tag{2.5.13}$$

因此，只要测得表面张力随温度的变化关系 $\sigma = \sigma(T)$，就可以求得表面系统的热力学函数。

---

## §2.6　热辐射的热力学理论

受热的固体会辐射电磁波，称为**热辐射**。一般情形下，热辐射的强度、强度按频率的分布与辐射体的温度、性质都有关。如果辐射体对电磁波的吸收和辐射达到平衡，热辐射的特性将只取决于温度，与辐射体的其他特性无关，称为**平衡辐射**。

考虑一个封闭的空窖，窖壁保持一定的温度 $T$。窖壁将不断向空窖发射并吸收电磁波，窖内辐射场与窖壁达到平衡后，二者具有共同的温度，显然空窖内的辐射就是平衡辐射，也称为**黑体辐射**（原因将在后面提到）。

平衡辐射包含各种频率、沿各个方向传播的电磁波，这些电磁波的振幅和相位是无规的。由热力学的一般理论可以证明，**窖内平衡辐射是空间均匀和各向同性的，它的内能密度和内能密度按频率的分布只取决于温度，与空窖的其他特性无关**。下面仅就最后一点进行论证。

> **论证**：设想有两个空窖，温度相同，但形状、体积和窖壁材料不同。开一个小窗把两个空窖连通起来，在窗上放上滤光片，滤光片只允许圆频率在 $\omega$ 到 $\omega+\mathrm{d}\omega$ 范围的电磁波通过，如图 2.6.1 所示。如果辐射场在 $\omega$ 到 $\omega+\mathrm{d}\omega$ 范围的内能密度在两窖不等，能量将通过小窗从内能密度较高的空窖辐射到内能密度较低的空窖，使前者温度降低、后者温度升高。这样就在温度相同的两个空窖间自发地产生温度差，可以利用该温度差来获得有用的功，这违背热力学第二定律，显然是不可能的。所以空窖辐射的内能密度和内能密度按频率的分布只可能是温度的函数。

### 2.6.1 内能密度与温度的关系

现在根据热力学理论导出窖内平衡辐射的热力学函数。这里要用到电磁理论关于辐射压强 $p$ 与辐射能量密度 $u$ 之间的关系：

$$p = \frac{1}{3}u \tag{2.6.1}$$

该关系于 1901 年得到列别杰夫（Lebedev）的实验证明，根据统计物理理论可以导出上式 [参阅式 (8.4.16)]。

将窖内平衡辐射看作热力学系统，选温度 $T$ 和体积 $V$ 为状态参量。由于空窖辐射是均匀的，其内能密度只是温度 $T$ 的函数，空窖辐射的内能 $U(T,V)$ 可以表示为

$$U(T,V) = u(T)V \tag{2.6.2}$$

利用热力学公式 (2.2.7) 和式 (2.6.1) 可得

$$u = \frac{T}{3}\frac{\mathrm{d}u}{\mathrm{d}T} - \frac{u}{3}$$

即

$$T\frac{\mathrm{d}u}{\mathrm{d}T} = 4u$$

积分得

$$u = aT^4 \tag{2.6.3}$$

其中 $a$ 是积分常量。**式 (2.6.3) 指出，空窖辐射的内能密度与热力学温度 $T$ 的四次方成正比。**

### 2.6.2 熵、可逆绝热关系与吉布斯函数

现在求空窖辐射的熵。将式 (2.6.3) 的 $u$ 和式 (2.6.1) 的 $p$ 代入热力学基本方程 $\mathrm{d}S = \dfrac{\mathrm{d}U + p\,\mathrm{d}V}{T}$，有

$$\mathrm{d}S = \frac{1}{T}\mathrm{d}(aT^4V) + \frac{1}{3}aT^3\mathrm{d}V = 4aT^2V\,\mathrm{d}T + \frac{4}{3}aT^3\mathrm{d}V = \frac{4}{3}a\,\mathrm{d}(VT^3)$$

积分得

$$S = \frac{4}{3}aT^3V \tag{2.6.4}$$

式 (2.6.4) 中没有积分常量，因为 $VT^3 = 0$ 时就不存在辐射场了。

在可逆绝热过程中辐射场的熵不变，这时有

$$T^3V = \text{常量} \tag{2.6.5}$$

将式 (2.6.1)、式 (2.6.3) 和式 (2.6.4) 代入吉布斯函数 $G = U - TS + pV$，可得**平衡辐射的吉布斯函数为零**，即

$$G = 0 \tag{2.6.6}$$

根据统计物理学可以导出平衡辐射的热力学函数（§8.4）。我们将看到，**式 (2.6.6) 是平衡辐射光子数不守恒的结果**。

### 2.6.3 辐射通量密度与斯特藩-玻耳兹曼定律

如果在窖壁开一小孔，电磁辐射将从小孔射出。假设小孔足够小，使窖内辐射场的平衡状态不受显著破坏。以 $J_u$ 表示单位时间内通过小孔单位面积向一侧辐射的辐射能量，称为**辐射通量密度**。辐射通量密度 $J_u$ 与辐射内能密度 $u$ 之间存在以下关系：

$$J_u = \frac{1}{4}cu \tag{2.6.7}$$

**证明**：考虑在单位时间内通过面积元 $\mathrm{d}A$ 向一侧辐射的能量。如果投射到 $\mathrm{d}A$ 上的是一束传播方向与 $\mathrm{d}A$ 的法线方向平行的平面电磁波，则单位时间内通过 $\mathrm{d}A$ 向一侧辐射的能量为 $cu\,\mathrm{d}A$。各向同性的辐射场包含各种传播方向，因此传播方向在 $\mathrm{d}\Omega$ 立体角内的辐射内能密度为 $\dfrac{u\,\mathrm{d}\Omega}{4\pi}$。单位时间内，传播方向在 $\mathrm{d}\Omega$ 立体角内通过 $\mathrm{d}A$ 向一侧辐射的能量为 $\dfrac{cu\,\mathrm{d}\Omega}{4\pi}\cos\theta\,\mathrm{d}A$，其中 $\theta$ 是传播方向与 $\mathrm{d}A$ 法线方向的夹角，如图 2.6.2 所示。对所有传播方向求积分，就可以得到单位时间内通过 $\mathrm{d}A$ 向一侧辐射的总辐射能量：

$$J_u\,\mathrm{d}A = \frac{cu\,\mathrm{d}A}{4\pi}\int \cos\theta\,\mathrm{d}\Omega = \frac{cu\,\mathrm{d}A}{4\pi}\int_0^{\pi/2}\sin\theta\cos\theta\,\mathrm{d}\theta\int_0^{2\pi}\mathrm{d}\varphi = \frac{1}{4}cu\,\mathrm{d}A$$

这就证明了式 (2.6.7)。将式 (2.6.3) 代入式 (2.6.7)，得

$$J_u = \frac{1}{4}caT^4 = \sigma T^4 \tag{2.6.8}$$

式 (2.6.8) 称为**斯特藩-玻耳兹曼定律**。这是 1879 年斯特藩（Stefan）在观测中发现、1884 年玻耳兹曼（Boltzmann）用热力学理论导出的。$\sigma$ 称为**斯特藩常量**，它的数值为

$$\sigma = 5.670\,374\,419\cdots\times 10^{-8}\ \mathrm{W\cdot m^{-2}\cdot K^{-4}}$$

在热力学中，$\sigma$ 的数值要由实验测定；统计物理学可以将 $\sigma$ 用其他基本物理常量表示出来 [式 (8.4.18)]。

### 2.6.4 基尔霍夫定律与绝对黑体

空窖内的辐射场与窖壁达到平衡后，其内能密度按频率的分布 $u(\omega)$ 只是温度的函数，与窖壁物质的特性无关。**这一事实说明，物质对各种频率电磁波的发射特性和吸收特性必然有某种联系**：例如，某一频率范围发射较强，吸收也必然较强。

将物体置于辐射场中，根据式 (2.6.7)，单位时间内投射到物体的单位面积上、圆频率在 $\mathrm{d}\omega$ 范围的辐射能量为 $\dfrac{c}{4}u(\omega)\mathrm{d}\omega$。以 $\alpha_\omega$ 表示其中被物体吸收的百分比，称为物体对频率在 $\omega$ 附近的辐射能量的**吸收因数**。这就是说，单位时间内被物体的单位面积所吸收、频率在 $\mathrm{d}\omega$ 范围的辐射能量为 $\dfrac{c}{4}\alpha_\omega u(\omega)\mathrm{d}\omega$，其余被物体反射。

用 $e_\omega\,\mathrm{d}\omega$ 表示单位时间内从物体的单位面积发射、频率在 $\mathrm{d}\omega$ 范围的辐射能量，$e_\omega$ 称为物体对频率在 $\omega$ 附近的电磁波的**面辐射强度**。注意 $\alpha_\omega$ 和 $e_\omega$ 表征物体的固有属性，与辐射场是否与物体达到平衡无关。

如果吸收与发射达到平衡，必有 $e_\omega\,\mathrm{d}\omega = \dfrac{c}{4}\alpha_\omega u(\omega,T)\mathrm{d}\omega$，或

$$\frac{e_\omega}{\alpha_\omega} = \frac{c}{4}u(\omega, T) \tag{2.6.9}$$

式中 $u(\omega,T)$ 是平衡辐射在 $\omega$ 处的能量密度。式 (2.6.9) 称为**基尔霍夫（Kirchhoff）定律**，它指出：**物体在任何频率处的面辐射强度与吸收因数之比对所有物体都相同，是频率和温度的普遍函数。**

吸收因数等于 1 的物体称为**绝对黑体**，它把投射到其表面的任何频率的电磁波完全吸收。绝对黑体是最好的吸收体，由式 (2.6.9) 知它也是最好的辐射体。对于绝对黑体，式 (2.6.9) 简化为

$$e_\omega = \frac{c}{4}u(\omega, T) \tag{2.6.10}$$

因此，**黑体的面辐射强度与平衡辐射的辐射通量密度完全相同。由于这个原因，平衡辐射也称为黑体辐射。**

我们日常接触到的物体都不是绝对黑体。开有小孔的空窖非常接近于绝对黑体，如图 2.6.3 所示：通过小孔射入空窖的电磁波，要经过窖壁多次反射才有可能从小孔射出，每经一次反射，窖壁都要吸收一部分电磁波，多次反射后从小孔射出的电磁波就极其微弱了。也就是说，**从小孔射入空窖的任何频率的电磁波实际上都将被空窖完全吸收**。这从另一角度说明空窖辐射是黑体辐射。

---

## §2.7　磁介质的热力学

### 2.7.1 基本方程与代换规则

在 §1.4 中我们求得了磁介质中磁场强度和磁化强度发生改变时外界所做的功

$$\mathrm{d}W = V\,\mathrm{d}\!\left(\frac{1}{2}\mu_0\mathcal{H}^2\right) + \mu_0 V\mathcal{H}\,\mathrm{d}\mathcal{M} \tag{2.7.1}$$

式中右方第一项是激发磁场所做的功，第二项是使介质磁化所做的功。**当热力学系统只包括介质而不包括磁场时，功的表达式只取式 (2.7.1) 右方的第二项**。这一项也可表示为

$$\mathrm{d}W = \mu_0\mathcal{H}\,\mathrm{d}\mathfrak{m} \tag{2.7.2}$$

其中 $\mathfrak{m} = \mathcal{M}V$ 是介质的总磁矩，我们假设介质是均匀磁化的。

如果忽略磁介质的体积变化，根据式 (1.14.7)，磁介质的热力学基本方程为

$$\mathrm{d}U = T\,\mathrm{d}S + \mu_0\mathcal{H}\,\mathrm{d}\mathfrak{m} \tag{2.7.3}$$

式 (2.7.3) 可由式 (1.14.6) 通过代换

$$p \to -\mu_0\mathcal{H}, \qquad V \to \mathfrak{m} \tag{2.7.4}$$

而得到。从这里可以看出，**由热力学基本方程 (1.14.6) 出发、通过数学推演而得到的关于简单系统的一般热力学关系，经代换式 (2.7.4) 后同样适用于磁介质**。例如，可以类似地定义磁介质的焓、自由能和吉布斯函数。吉布斯函数 $G$ 为

$$G = U - TS - \mu_0\mathcal{H}\mathfrak{m} \tag{2.7.5}$$

对式 (2.7.5) 求微分，将式 (2.7.3) 代入，可得 $G$ 的全微分：

$$\mathrm{d}G = -S\,\mathrm{d}T - \mu_0\mathfrak{m}\,\mathrm{d}\mathcal{H} \tag{2.7.6}$$

由全微分条件可得

$$\left(\frac{\partial S}{\partial \mathcal{H}}\right)_T = \mu_0\left(\frac{\partial \mathfrak{m}}{\partial T}\right)_{\mathcal{H}} \tag{2.7.7}$$

式 (2.7.7) 是磁介质的一个麦克斯韦关系。

### 2.7.2 绝热去磁制冷

由于存在函数关系 $S = S(T,\mathcal{H})$，故有

$$\left(\frac{\partial S}{\partial \mathcal{H}}\right)_T \left(\frac{\partial \mathcal{H}}{\partial T}\right)_S \left(\frac{\partial T}{\partial S}\right)_{\mathcal{H}} = -1$$

或

$$\left(\frac{\partial T}{\partial \mathcal{H}}\right)_S = -\left(\frac{\partial S}{\partial \mathcal{H}}\right)_T \left(\frac{\partial T}{\partial S}\right)_{\mathcal{H}} \tag{2.7.8}$$

在磁场不变时，磁介质的热容 $C_{\mathcal{H}}$ 为

$$C_{\mathcal{H}} = T\left(\frac{\partial S}{\partial T}\right)_{\mathcal{H}} \tag{2.7.9}$$

（上式与 (2.2.8) 相当。）将式 (2.7.7) 和式 (2.7.9) 代入式 (2.7.8)，得

$$\left(\frac{\partial T}{\partial \mathcal{H}}\right)_S = -\frac{\mu_0 T}{C_{\mathcal{H}}}\left(\frac{\partial \mathfrak{m}}{\partial T}\right)_{\mathcal{H}} \tag{2.7.10}$$

假设磁介质遵从居里定律

$$\mathfrak{m} = \frac{CV}{T}\mathcal{H} \tag{2.7.11}$$

代入式 (2.7.10)，得

$$\left(\frac{\partial T}{\partial \mathcal{H}}\right)_S = \frac{CV}{C_{\mathcal{H}}T}\mu_0\mathcal{H} \tag{2.7.12}$$

式 (2.7.12) 的右方是正数，这说明**在绝热条件下减小磁场时，磁介质的温度将降低**。这个效应称为**绝热去磁制冷效应**，是获得 1 K 以下低温的有效方法。在 §2.8 中我们将进一步讨论利用绝热去磁制冷效应获得低温的方法。

### 2.7.3 计及体积变化的推广

如果考虑磁介质体积的变化，根据式 (1.14.7)，热力学基本方程应为

$$\mathrm{d}U = T\,\mathrm{d}S - p\,\mathrm{d}V + \mu_0\mathcal{H}\,\mathrm{d}\mathfrak{m} \tag{2.7.13}$$

吉布斯函数为

$$G = U - TS + pV - \mu_0\mathcal{H}\mathfrak{m} \tag{2.7.14}$$

其全微分为

$$\mathrm{d}G = -S\,\mathrm{d}T + V\,\mathrm{d}p - \mu_0\mathfrak{m}\,\mathrm{d}\mathcal{H} \tag{2.7.15}$$

由全微分条件可得

$$\left(\frac{\partial V}{\partial \mathcal{H}}\right)_{T,p} = -\mu_0\left(\frac{\partial \mathfrak{m}}{\partial p}\right)_{T,\mathcal{H}} \tag{2.7.16}$$

这也是磁介质的一个麦克斯韦关系。式 (2.7.16) 左方的偏导数给出在温度和压强保持不变时体积随磁场的变化率，它描述**磁致伸缩效应**；右方给出在温度和磁场保持不变时介质磁矩随压强的变化率，它描述**压磁效应**。**式 (2.7.16) 给出磁致伸缩效应与压磁效应的关系。**

### 2.7.4 磁力功的另一表达式

除了式 (2.7.1) 和式 (2.7.2) 两个表达式以外，磁化功往往还使用另一个表达式。设空间中存在给定的不均匀磁场，例如由永久磁铁产生的磁场。将样品从无穷远处移入磁场内，例如从 $x = -\infty$ 处沿 $x$ 轴移到 $x = a$ 处，介质将被磁化。当样品在 $x$ 处时，所受磁场的力为

$$\mu_0\mathfrak{m}(x)\frac{\mathrm{d}\mathcal{H}(x)}{\mathrm{d}x}$$

> 脚注：磁矩 $\mathfrak{m}$ 在外磁场 $\mathcal{H}$ 中所受的力为 $\mu_0(\mathfrak{m}\cdot\nabla)\mathcal{H}$，可参阅 Jackson J D. Classical Electrodynamics, Second Edition, §5.7。为简单起见，我们考虑一维的情形。

移动样品时，外界必须克服此力而做功：

$$W = -\mu_0\int_{-\infty}^{a}\mathfrak{m}(x)\frac{\mathrm{d}\mathcal{H}(x)}{\mathrm{d}x}\mathrm{d}x = -\mu_0\int_0^{\mathcal{H}(a)}\mathfrak{m}\,\mathrm{d}\mathcal{H} \tag{2.7.17}$$

其中 $\mathcal{H}(a)$ 是在 $x = a$ 处的磁场强度，在 $x = -\infty$ 处磁场为零。上式可通过分部积分表示为

$$W = -\mu_0\mathfrak{m}(a)\mathcal{H}(a) + \int_0^{\mathfrak{m}(a)}\mu_0\mathcal{H}\,\mathrm{d}\mathfrak{m} \tag{2.7.18}$$

右方第一项是磁矩 $\mathfrak{m}(a)$ 在磁场 $\mathcal{H}(a)$ 中的势能；由式 (2.7.2) 可知，式 (2.7.18) 右方第二项是将介质磁化所做的功。与式 (2.7.17) 相应的微功是

$$\mathrm{d}W = -\mu_0\mathfrak{m}\,\mathrm{d}\mathcal{H} \tag{2.7.19}$$

它不但包含外场改变 $\mathrm{d}\mathcal{H}$ 时为使样品磁矩发生改变所做的功，而且包括样品在外磁场中势能的改变。

由于两态内能之差是通过绝热过程的功定义的，使用式 (2.7.2) 和式 (2.7.19) 两种不同的功的表达式时，内能的含义显然也将不同。以 $U$ 和 $U^*$ 表示相应的内能，由式 (2.7.18) 知

$$U^* = -\mu_0\mathfrak{m}\mathcal{H} + U \tag{2.7.20}$$

即 **$U^*$ 包含样品在外磁场中的势能**。其他热力学函数也有相应的变化。

---

## *§2.8　获得低温的方法（选读）

低温技术在现代科学技术中有重要的应用。本节对获得低温的方法作简略的介绍。

### 2.8.1 气体液化

将气体液化可以获得低至 1 K 的低温。目前常用节流过程或者节流过程与绝热膨胀相结合的方法来液化气体。

在 §2.3 中讲过，令气体在制冷区节流膨胀可使气体降温。用节流过程制冷有两个优点：

1. 装置没有移动的部分（低温下移动部分的润滑是技术上十分困难的问题）；
2. 在一定的压强降落 $\Delta$ 下，温度越低所获得的温度降落越大 [参阅式 (2.3.7) 和图 1.3.1]。焦汤效应的典型值是 $10^{-1}\sim 1\ \mathrm{K\cdot atm^{-1}}$。

为了使气体的温度降至临界温度以下而液化，可以令节流过程重复进行，并通过**逆流热交换器**使经节流膨胀降温后的气体对后来进入的气体进行预冷，从而把各次节流膨胀所获得的冷却效应积累起来。人们用这种方法成功地实现了沸点最低的氢（杜瓦，Dewar，1898 年）和氦（昂内斯，Onnes，1908 年）的液化。

但是用节流过程降温，**气体的初温必须低于反转温度**。表 2.8.1 列出几种气体的反转温度。早期的工作用液氮将氢预冷，液氢容易爆炸是一个严重缺点。

| 气体 | 最高反转温度/K | 临界温度/K | 1 atm 下的沸点/K |
| --- | --- | --- | --- |
| $\mathrm{CO_2}$ | 1 275 | 304 | — |
| $\mathrm{N_2}$ | 607 | 126 | 77.3 |
| $\mathrm{H_2}$ | 204 | 33.1 | 20.4 |
| He | 43 | 5.2 | 4.2 |

### 2.8.2 绝热膨胀

气体经绝热膨胀后温度总是降低的，因此**用绝热膨胀过程降温不必经过预冷**。缺点是膨胀机有移动的部分，而且温度越低降温效应越小 [式 (2.3.8)]。

卡皮查（Kapitza，1934 年）将绝热膨胀过程和节流过程结合使用：先用绝热膨胀过程使氦降温到反转温度以下，再用节流过程将氦液化。

1 atm 下氦的沸点是 4.2 K。用抽气机将氦的蒸气抽走，使液氦迅速蒸发或低压沸腾可进一步降温。不过氦的饱和蒸气压随温度降低而迅速减小，降温效应随之下降。用这种方法一般可获得 1 K 的低温。

### 2.8.3 磁冷却法

产生 1 K 以下低温的一个有效方法是**磁冷却法**，这是德拜（Debye）在 1926 年提出来的。

在 §2.7 中讲过，在绝热过程中顺磁性固体的温度随磁场的减小而下降。图 2.8.1 是实验装置的示意图。图 2.8.2 画出了不同磁场下顺磁体的熵 $S$ 随温度 $T$ 变化的曲线 [式 (7.8.7)]。

实验中将顺磁体放在装有低压氦气的容器内，通过低压氦气与液氦的接触而保持在 1 K 左右的低温 $T_i$，加上磁场 $\mathcal{H}_i$（量级为 $10^8\ \mathrm{A\cdot m^{-1}}$）使顺磁体磁化。磁化过程释放出的热量（习题 2.20）由液氦吸收，从而保证磁化过程是**等温**的。在这个等温磁化过程中，顺磁体的状态由 $a$ 点变到 $b$ 点。顺磁体磁化后，抽出低压氦气而使顺磁体**绝热**，然后准静态地使磁场减小为 $\mathcal{H}_f$（$\mathcal{H}_f$ 一般为零）。在这个**绝热去磁**过程中，顺磁体的熵保持不变，其状态由 $b$ 点变到 $c$ 点，温度降为 $T_f$。

利用固体中顺磁离子的绝热去磁效应可以产生 1 K 以下至 mK（$10^{-3}\ \mathrm{K}$）量级的低温。例如从 0.5 K 出发，使硝酸铈镁 $2\mathrm{Ce(NO_3)_3\cdot 3Mg(NO_3)_2\cdot 24H_2O}$ 绝热去磁可降温到 2 mK。**当温度降到 mK 量级时，顺磁离子磁矩间的相互作用便不能忽略**——磁矩间的相互作用相当于产生一个等效的磁场（大小为 $10^3\sim 10^4\ \mathrm{A\cdot m^{-1}}$），使磁矩的分布有序化，该方法便不再有效。

### 2.8.4 核绝热去磁

核磁矩（来自核子磁矩）的大小约为原子磁矩（来自电子磁矩）的 1/2 000。因此，核磁矩间的相互作用较顺磁离子间的相互作用要弱得多，利用核绝热去磁可以获得更低的温度。

- 戈特（Gorter）在 1934 年提出这种方法。
- 1956 年库尔蒂（Kurti）和西蒙（Simon）等成功地将铜的核自旋温度降到 $\mu\mathrm{K}$（$10^{-6}\ \mathrm{K}$）量级。
- 20 世纪 80 年代中期，芬兰赫尔辛基技术大学小组将铜的核自旋温度降到 25 nK（$1\ \mathrm{nK} = 10^{-9}\ \mathrm{K}$），其中导电电子和晶格的温度为 $50\ \mu\mathrm{K}$。核自旋系统的弛豫时间为 ms 量级，在 $50\ \mu\mathrm{K}$ 的低温下自旋-晶格的弛豫时间约为 3 h。在 20 世纪 80 年代末，该小组更将银的核自旋温度降至 2 nK。

### 2.8.5 稀释制冷与激光制冷

绝热去磁制冷过程是单一循环，不能连续工作。1951 年，伦敦（London）提出 $^3\mathrm{He}$-$^4\mathrm{He}$ **稀释制冷**的制冷新方法。在 20 世纪 70 年代末，已实现温度为 2 mK 的 $^3\mathrm{He}$-$^4\mathrm{He}$ 稀释制冷的连续运转，代替绝热去磁在许多研究领域得到应用。我们将在 §4.4 中介绍 $^3\mathrm{He}$-$^4\mathrm{He}$ 的相图，在 §9.8 中介绍 $^3\mathrm{He}$-$^4\mathrm{He}$ 稀释制冷的原理。

20 世纪 80 年代发展了一种新的制冷方法——**激光制冷**，应用该方法在 20 世纪 90 年代中期获得了低至 170 nK 的低温。我们将在 §10.7 中对激光制冷作初步的介绍。

---

## 本章小结

### 一、公式速查

**四个全微分与自然变量**

| 热力学函数 | 自然变量 | 全微分 | 偏导数关系 |
| --- | --- | --- | --- |
| $U$ | $S$、$V$ | $\mathrm{d}U = T\mathrm{d}S - p\mathrm{d}V$ (2.1.1) | $\left(\dfrac{\partial U}{\partial S}\right)_V = T$，$\left(\dfrac{\partial U}{\partial V}\right)_S = -p$ (2.1.2) |
| $H = U + pV$ | $S$、$p$ | $\mathrm{d}H = T\mathrm{d}S + V\mathrm{d}p$ (2.1.4) | $\left(\dfrac{\partial H}{\partial S}\right)_p = T$，$\left(\dfrac{\partial H}{\partial p}\right)_S = V$ (2.1.5) |
| $F = U - TS$ | $T$、$V$ | $\mathrm{d}F = -S\mathrm{d}T - p\mathrm{d}V$ (2.1.7) | $\left(\dfrac{\partial F}{\partial T}\right)_V = -S$，$\left(\dfrac{\partial F}{\partial V}\right)_T = -p$ (2.1.8) |
| $G = U - TS + pV$ | $T$、$p$ | $\mathrm{d}G = -S\mathrm{d}T + V\mathrm{d}p$ (2.1.10) | $\left(\dfrac{\partial G}{\partial T}\right)_p = -S$，$\left(\dfrac{\partial G}{\partial p}\right)_T = V$ (2.1.11) |

**麦克斯韦关系**

| 来源 | 关系 |
| --- | --- |
| $U(S,V)$ | $\left(\dfrac{\partial T}{\partial V}\right)_S = -\left(\dfrac{\partial p}{\partial S}\right)_V$ (2.2.1) |
| $H(S,p)$ | $\left(\dfrac{\partial T}{\partial p}\right)_S = \left(\dfrac{\partial V}{\partial S}\right)_p$ (2.2.2) |
| $F(T,V)$ | $\left(\dfrac{\partial S}{\partial V}\right)_T = \left(\dfrac{\partial p}{\partial T}\right)_V$ (2.2.3) |
| $G(T,p)$ | $\left(\dfrac{\partial S}{\partial p}\right)_T = -\left(\dfrac{\partial V}{\partial T}\right)_p$ (2.2.4) |

**把不可测量化为可测量**

| 内容 | 公式 |
| --- | --- |
| 热容的另一表达式 | $C_V = T\left(\dfrac{\partial S}{\partial T}\right)_V$ (2.2.5)；$C_p = T\left(\dfrac{\partial S}{\partial T}\right)_p$ (2.2.8) |
| 内能随体积的变化 | $\left(\dfrac{\partial U}{\partial V}\right)_T = T\left(\dfrac{\partial p}{\partial T}\right)_V - p$ (2.2.7) |
| 焓随压强的变化 | $\left(\dfrac{\partial H}{\partial p}\right)_T = V - T\left(\dfrac{\partial V}{\partial T}\right)_p$ (2.2.10) |
| 两热容之差 | $C_p - C_V = T\left(\dfrac{\partial p}{\partial T}\right)_V\left(\dfrac{\partial V}{\partial T}\right)_p$ (2.2.11) $= \dfrac{VT\alpha^2}{\kappa_T}$ (2.2.12) $= -T\dfrac{(\partial p/\partial T)_V^2}{(\partial p/\partial V)_T}$ (2.2.15) |
| 绝热与等温压缩系数之比 | $\dfrac{\kappa_S}{\kappa_T} = \dfrac{C_V}{C_p}$ (2.2.14) |

**节流与绝热膨胀**

| 内容 | 公式 |
| --- | --- |
| 节流前后焓相等 | $H_1 = H_2$ (2.3.1) |
| 焦汤系数 | $\mu = \left(\dfrac{\partial T}{\partial p}\right)_H = \dfrac{1}{C_p}\left[T\left(\dfrac{\partial V}{\partial T}\right)_p - V\right]$ (2.3.4) $= \dfrac{V}{C_p}(T\alpha - 1)$ (2.3.5) |
| 昂内斯方程下 | $\mu = \dfrac{n}{C_p}\left(T\dfrac{\mathrm{d}B}{\mathrm{d}T} - B\right)$ (2.3.7) |
| 绝热膨胀温变率 | $\left(\dfrac{\partial T}{\partial p}\right)_S = \dfrac{VT\alpha}{C_p}$ (2.3.8) |

**基本热力学函数的确定**

| 自变量 | 公式 |
| --- | --- |
| $p = p(T,V)$ | $U = \displaystyle\int\left\{C_V\mathrm{d}T + \left[T\left(\dfrac{\partial p}{\partial T}\right)_V - p\right]\mathrm{d}V\right\} + U_0$ (2.4.3)；$S = \displaystyle\int\left[\dfrac{C_V}{T}\mathrm{d}T + \left(\dfrac{\partial p}{\partial T}\right)_V\mathrm{d}V\right] + S_0$ (2.4.5) |
| $V = V(T,p)$ | $H = \displaystyle\int\left\{C_p\mathrm{d}T + \left[V - T\left(\dfrac{\partial V}{\partial T}\right)_p\right]\mathrm{d}p\right\} + H_0$ (2.4.8)；$S = \displaystyle\int\left[\dfrac{C_p}{T}\mathrm{d}T - \left(\dfrac{\partial V}{\partial T}\right)_p\mathrm{d}p\right] + S_0$ (2.4.10) |
| 理想气体 | $H_m = \displaystyle\int C_{p,m}\mathrm{d}T + H_{m0}$ (2.4.11)；$S_m = \displaystyle\int\dfrac{C_{p,m}}{T}\mathrm{d}T - R\ln p + S_{m0}$ (2.4.12)；$G_m = RT(\varphi + \ln p)$ (2.4.15) |
| 范德瓦耳斯气体 | $U_m = \displaystyle\int C_{V,m}\mathrm{d}T - \dfrac{a}{V_m} + U_{m0}$ (2.4.17)；$S_m = \displaystyle\int\dfrac{C_{V,m}}{T}\mathrm{d}T + R\ln(V_m - b) + S_{m0}$ (2.4.18) |
| 简单固体 | $U = \displaystyle\int C_V\mathrm{d}T + \dfrac{1}{2}\dfrac{(V - V_1)^2}{\kappa_T V_0} + U_0$ (2.4.19)；$S = \displaystyle\int\dfrac{C_V}{T}\mathrm{d}T + \dfrac{\alpha}{\kappa_T}V + S_0$ (2.4.20) |

**特性函数**

| 内容 | 公式 |
| --- | --- |
| $F(T,V)$ | $S = -\dfrac{\partial F}{\partial T}$，$p = -\dfrac{\partial F}{\partial V}$ (2.5.2)；$U = F - T\dfrac{\partial F}{\partial T}$ (2.5.3，吉布斯-亥姆霍兹方程) |
| $G(T,p)$ | $S = -\dfrac{\partial G}{\partial T}$，$V = \dfrac{\partial G}{\partial p}$ (2.5.5)；$U = G - T\dfrac{\partial G}{\partial T} - p\dfrac{\partial G}{\partial p}$ (2.5.6)；$H = G - T\dfrac{\partial G}{\partial T}$ (2.5.7) |
| 表面系统 | $F = \sigma A$ (2.5.11)；$S = -A\dfrac{\mathrm{d}\sigma}{\mathrm{d}T}$ (2.5.12)；$U = A\left(\sigma - T\dfrac{\mathrm{d}\sigma}{\mathrm{d}T}\right)$ (2.5.13) |

**热辐射**

| 内容 | 公式 |
| --- | --- |
| 辐射压强 | $p = \dfrac{1}{3}u$ (2.6.1) |
| 内能密度 | $U = u(T)V$ (2.6.2)；$u = aT^4$ (2.6.3) |
| 熵与可逆绝热 | $S = \dfrac{4}{3}aT^3V$ (2.6.4)；$T^3V = \text{常量}$ (2.6.5) |
| 吉布斯函数 | $G = 0$ (2.6.6) |
| 辐射通量密度 | $J_u = \dfrac{1}{4}cu$ (2.6.7) $= \sigma T^4$ (2.6.8)，$\sigma = 5.670\,374\,419\cdots\times 10^{-8}\ \mathrm{W\cdot m^{-2}\cdot K^{-4}}$ |
| 基尔霍夫定律 | $\dfrac{e_\omega}{\alpha_\omega} = \dfrac{c}{4}u(\omega,T)$ (2.6.9)；黑体 $e_\omega = \dfrac{c}{4}u(\omega,T)$ (2.6.10) |

**磁介质**

| 内容 | 公式 |
| --- | --- |
| 代换规则 | $p \to -\mu_0\mathcal{H}$，$V \to \mathfrak{m}$ (2.7.4) |
| 基本方程与吉布斯函数 | $\mathrm{d}U = T\mathrm{d}S + \mu_0\mathcal{H}\mathrm{d}\mathfrak{m}$ (2.7.3)；$G = U - TS - \mu_0\mathcal{H}\mathfrak{m}$ (2.7.5)；$\mathrm{d}G = -S\mathrm{d}T - \mu_0\mathfrak{m}\mathrm{d}\mathcal{H}$ (2.7.6) |
| 麦克斯韦关系 | $\left(\dfrac{\partial S}{\partial \mathcal{H}}\right)_T = \mu_0\left(\dfrac{\partial \mathfrak{m}}{\partial T}\right)_{\mathcal{H}}$ (2.7.7)；$\left(\dfrac{\partial V}{\partial \mathcal{H}}\right)_{T,p} = -\mu_0\left(\dfrac{\partial \mathfrak{m}}{\partial p}\right)_{T,\mathcal{H}}$ (2.7.16) |
| 绝热去磁 | $\left(\dfrac{\partial T}{\partial \mathcal{H}}\right)_S = -\dfrac{\mu_0 T}{C_{\mathcal{H}}}\left(\dfrac{\partial \mathfrak{m}}{\partial T}\right)_{\mathcal{H}}$ (2.7.10)；居里定律下 $= \dfrac{CV}{C_{\mathcal{H}}}\mu_0\mathcal{H}$ (2.7.12) |
| 另一种磁功与内能 | $\mathrm{d}W = -\mu_0\mathfrak{m}\,\mathrm{d}\mathcal{H}$ (2.7.19)；$U^* = U - \mu_0\mathfrak{m}\mathcal{H}$ (2.7.20) |

### 二、概念脉络

```
第一章的热力学基本方程 dU = T dS - p dV
        │
        ▼
§2.1 四个全微分与自然变量 U(S,V)、H(S,p)、F(T,V)、G(T,p)
        │
        ├──► 偏导数关系 (2.1.2)(2.1.5)(2.1.8)(2.1.11)
        └──► 麦克斯韦关系 (2.1.3)(2.1.6)(2.1.9)(2.1.12)
        │
        ▼
§2.2 把不可测量化为可测量
     (∂U/∂V)_T、(∂H/∂p)_T、C_p - C_V、κ_S/κ_T
        │
        ▼
§2.3 应用：节流过程（H 不变，μ = (V/C_p)(Tα-1)）
            绝热膨胀（S 不变，(∂T/∂p)_S = VTα/C_p）
        │
        ▼
§2.4 基本热力学函数的确定：物态方程 + 热容 → U、S、H、G
        │
        ▼
§2.5 特性函数（马休定理）：一个 F(T,V) 或 G(T,p) 决定全部热力学函数
        │
        ├──► §2.6 热辐射：u = aT^4，J_u = σT^4，基尔霍夫定律
        ├──► §2.7 磁介质：代换 p → -μ0H、V → m，绝热去磁制冷
        └──► *§2.8 获得低温的方法（节流 / 绝热膨胀 / 磁冷却 / 稀释制冷 / 激光制冷）
```

### 三、易错点

1. **全微分条件是对"自然变量"而言的**。只有写成 $U(S,V)$、$H(S,p)$、$F(T,V)$、$G(T,p)$，交叉偏导数相等才直接给出麦克斯韦关系；自变量选错会导出错误关系。
2. **四个麦克斯韦关系的记忆**：分母是 $S$、$V$ 的来自 $U$ 和 $H$，分母是 $T$、$p$ 的来自 $F$ 和 $G$；两个式子带负号 (2.2.1)(2.2.4)，两个不带 (2.2.2)(2.2.3)。
3. **$C_p - C_V \geqslant 0$ 恒成立**，因为 (2.2.12) 右方 $VT\alpha^2/\kappa_T$ 不可能为负；只有当 $\alpha = 0$（如水在 4 °C）时两者相等。
4. **理想气体的焦汤系数为零**（$\alpha = 1/T$），节流前后温度不变；只有实际气体才有制冷区与制温区，二者以**反转曲线**（$\alpha T = 1$）分界。
5. **节流过程不是等焓过程**，只是初、终态焓相等；它本身是**不可逆过程**，$p$-$V$ 图上的等焓线只是初、终态的连线。
6. **节流降温要求初温低于反转温度**，绝热膨胀降温则不必预冷——这是两者最重要的工程差别。
7. **固体和液体的 $C_V$ 难以直接测量**，所以选 $T$、$p$ 为自变量、用 $C_p$ 更实用；$C_V$ 可由 (2.2.12) 反算。
8. **特性函数的"特性"在于：一个函数 + 自然变量 → 全部热力学函数**。由 $F$ 求对 $V$ 的偏导即得物态方程，求对 $T$ 的偏导即得熵，再由 $F = U - TS$ 求内能。
9. **理想气体的 $U$、$H$ 只依赖 $T$**，故例 2.4.1 中 $V_m - T(\partial V_m/\partial T)_p = 0$；范德瓦耳斯气体则不然，内能中含 $-a/V_m$ 项。
10. **热辐射的 $G = 0$** 是光子数不守恒的结果，不是巧合。
11. **$\sigma$ 是"一词两用"**：§2.5 中是表面张力系数，§2.6 中是斯特藩常量，看上下文区分。
12. **磁介质的代换是 $p \to -\mu_0\mathcal{H}$、$V \to \mathfrak{m}$**，负号不能漏，否则吉布斯函数与麦克斯韦关系的符号会错。
13. **两种磁功表达式对应两种内能**：$\mathrm{d}W = \mu_0\mathcal{H}\mathrm{d}\mathfrak{m}$ 对应不含外场势能的 $U$；$\mathrm{d}W = -\mu_0\mathfrak{m}\mathrm{d}\mathcal{H}$ 对应包含外场势能的 $U^* = U - \mu_0\mathfrak{m}\mathcal{H}$ (2.7.20)，混用会得到不同的内能表达式。

### 四、名词中英对照

| 中文 | 英文 | 中文 | 英文 |
| --- | --- | --- | --- |
| 自然变量 | natural variable | 特性函数 | characteristic function |
| 全微分 | total differential | 麦克斯韦关系 | Maxwell relations |
| 雅可比行列式 | Jacobian determinant | 体胀系数 | expansion coefficient ($\alpha$) |
| 等温压缩系数 | isothermal compressibility ($\kappa_T$) | 绝热压缩系数 | adiabatic compressibility ($\kappa_S$) |
| 焦耳-汤姆孙系数（焦汤系数） | Joule-Thomson coefficient ($\mu$) | 节流过程 | throttling process |
| 反转曲线 / 反转温度 | inversion curve / inversion temperature | 制冷区 / 制温区 | cooling / heating region |
| 绝热膨胀 | adiabatic expansion | 逆流热交换器 | counter-current heat exchanger |
| 表面张力系数 | surface tension ($\sigma$) | 表面系统 | surface system |
| 吉布斯-亥姆霍兹方程 | Gibbs-Helmholtz equation | 马休定理 | Massieu's theorem |
| 热辐射 / 平衡辐射 | thermal radiation / equilibrium radiation | 黑体辐射 | black-body radiation |
| 辐射能量密度 | radiation energy density ($u$) | 辐射通量密度 | radiation flux density ($J_u$) |
| 斯特藩-玻耳兹曼定律 | Stefan-Boltzmann law | 斯特藩常量 | Stefan constant |
| 吸收因数 | absorptivity ($\alpha_\omega$) | 面辐射强度 | emissive power ($e_\omega$) |
| 基尔霍夫定律 | Kirchhoff's law | 绝对黑体 | absolute black body |
| 磁介质 | magnetic medium | 磁化强度 | magnetization ($\mathcal{M}$) |
| 总磁矩 | total magnetic moment ($\mathfrak{m}$) | 居里定律 | Curie's law |
| 绝热去磁制冷 | adiabatic demagnetization | 磁致伸缩效应 / 压磁效应 | magnetostriction / piezomagnetism |
| 核绝热去磁 | nuclear adiabatic demagnetization | 稀释制冷 | dilution refrigeration |
| 磁冷却法 | magnetic cooling | 激光制冷 | laser cooling |

