---
tags:
  - 热统
  - 统计物理
  - 第七章
created: 2026-10-22
---

# 第七章　玻耳兹曼统计

> **教材**：汪志诚《热力学·统计物理》（第六版）
> **范围**：§7.1–§7.9（教材 pp. 157–181；其中 *§7.9 为带星号的选读内容）

## 本章地图

| 节 | 标题 | 一句话 |
| --- | --- | --- |
| §7.1 | 热力学量的统计表达式 | 由粒子配分函数 $Z_1$ 求出内能、广义力、熵和自由能；$\alpha,\beta$ 与温度和化学势相联系 |
| §7.2 | 理想气体的物态方程 | 对自由粒子的配分函数作相空间积分，导出 $pV=NkT$，并给出经典极限 $n\lambda^3\ll1$ |
| §7.3 | 麦克斯韦速度分布律 | 由玻耳兹曼分布导出速度、速率分布，以及最概然速率、平均速率、方均根速率和碰壁数 |
| §7.4 | 能量均分定理 | 经典系统中每个独立平方项平均贡献 $kT/2$；解释理想气体与固体热容，但在低温及黑体辐射处失效 |
| §7.5 | 理想气体的内能和热容 | 用配分函数分解计算平动、振动、转动贡献，解决经典理论无法解释的低温“冻结”问题 |
| §7.6 | 理想气体的熵 | 经典熵不满足广延性；计入 $1/N!$ 后得到绝对熵、萨克尔-特特罗德公式及负化学势 |
| §7.7 | 固体热容的爱因斯坦理论 | 把固体看成频率相同的 $3N$ 个量子振子，得到高温杜隆-珀蒂极限和低温指数衰减 |
| §7.8 | 顺磁性固体 | 二能级磁矩模型给出磁化强度 $n\mu\tanh(\mu B/kT)$，并解释极化饱和、居里定律和熵 |
| *§7.9 | 负温度状态 | 能级有上界且与正温系统弱耦合的系统中，熵随内能下降时可出现负绝对温度 |

**主线**：由玻耳兹曼分布引入 $Z_1$ → 用配分函数统一计算热力学量 → 应用于理想气体、速度分布、热容和顺磁体 → 在量子能级离散性不可忽略时，以振动、转动和磁矩的量子模型修正经典结论 → *讨论有能级上界系统的负温度。

---

## §7.1　热力学量的统计表达式

### 7.1.1 粒子配分函数

本章讨论两类遵从玻耳兹曼分布的系统：

- 可分辨的定域粒子系统；
- 满足经典极限条件的玻色或费米气体。

系统的内能是粒子无规则运动总能量的统计平均值：

$$
U=\sum_l a_l\varepsilon_l
=\sum_l\varepsilon_l\omega_l e^{-\alpha-\beta\varepsilon_l} \tag{7.1.1}
$$

引入

$$
Z_1=\sum_l\omega_l e^{-\beta\varepsilon_l} \tag{7.1.2}
$$

称为**粒子配分函数**，简称配分函数。由式 (6.6.9) 得

$$
N=e^{-\alpha}\sum_l\omega_l e^{-\beta\varepsilon_l}
=e^{-\alpha}Z_1 \tag{7.1.3}
$$

因此

$$
e^{-\alpha}=\frac{N}{Z_1}.
$$

此式给出了参数 $\alpha$、粒子数 $N$ 与配分函数的关系，可用来消去 $U$ 中的 $e^{-\alpha}$。

### 7.1.2 内能与广义力的统计表达式

由式 (7.1.1)：

$$
\begin{aligned}
U
&=e^{-\alpha}\sum_l\varepsilon_l
\omega_l e^{-\beta\varepsilon_l}\\
&=e^{-\alpha}\left(-\frac{\partial}{\partial\beta}\right)
\sum_l\omega_l e^{-\beta\varepsilon_l}\\
&=\frac{N}{Z_1}\left(-\frac{\partial}{\partial\beta}\right)Z_1
=-N\frac{\partial}{\partial\beta}\ln Z_1.
\end{aligned}
\tag{7.1.4}
$$

在无穷小过程中，热力学第一定律为

$$
\mathrm dU=\mathrm dW+\mathrm dQ, \tag{7.1.5}
$$

其中 $\mathrm dW$ 是外界对系统做的功，$\mathrm dQ$ 是系统从外界吸收的热量。

若过程准静态，外界对系统做功可写成

$$
\mathrm dW=Y\,\mathrm dy,
$$

其中 $y$ 是外参量，$Y$ 是与 $y$ 共轭的广义力。例如体积改变 $\mathrm dV$ 时，

$$
\mathrm dW=-p\,\mathrm dV.
$$

粒子的能量通常是外参量的函数。外参量改变时，外界作用于处于能级 $\varepsilon_l$ 的单个粒子的广义力为

$$
-\frac{\partial\varepsilon_l}{\partial y}.
$$

系统所受总广义力为

$$
\begin{aligned}
Y
&=-\sum_l\frac{\partial\varepsilon_l}{\partial y}a_l
=-\sum_l\frac{\partial\varepsilon_l}{\partial y}
\omega_l e^{-\alpha-\beta\varepsilon_l}\\
&=e^{-\alpha}\left(-\frac{1}{\beta}\frac{\partial}{\partial y}\right)
\sum_l\omega_l e^{-\beta\varepsilon_l}\\
&=\frac{N}{Z_1}
\left(-\frac{1}{\beta}\frac{\partial}{\partial y}\right)Z_1
=-\frac{N}{\beta}\frac{\partial}{\partial y}\ln Z_1.
\end{aligned}
\tag{7.1.6}
$$

对体积功，$y=V$，故

$$
p=\frac{N}{\beta}\frac{\partial}{\partial V}\ln Z_1. \tag{7.1.7}
$$

### 7.1.3 功、热与熵

外参量改变 $\mathrm dy$ 时，外界对系统做功

$$
Y\,\mathrm dy
=\mathrm dy\sum_l\frac{\partial\varepsilon_l}{\partial y}a_l
=\sum_l a_l\,\mathrm d\varepsilon_l. \tag{7.1.8}
$$

对

$$
U=\sum_l\varepsilon_l a_l
$$

求全微分：

$$
\mathrm dU
=\sum_l a_l\,\mathrm d\varepsilon_l
+\sum_l\varepsilon_l\,\mathrm da_l. \tag{7.1.9}
$$

第一项是分布不变、仅因能级随外参量改变而产生的内能变化，恰好等于外界做功；第二项是能级不变、因粒子重新分布产生的内能变化，等于系统吸热：

$$
\mathrm dQ
=\mathrm dU-Y\,\mathrm dy
=\sum_l\varepsilon_l\,\mathrm da_l.
$$

热力学第二定律给出

$$
\frac{1}{T}\mathrm dQ
=\frac{1}{T}\left(\mathrm dU-Y\,\mathrm dy\right)
=\mathrm dS. \tag{7.1.10}
$$

利用式 (7.1.4)、式 (7.1.6)：

$$
\mathrm dQ
=\mathrm dU-Y\,\mathrm dy
=-N\,\mathrm d\left(\frac{\partial\ln Z_1}{\partial\beta}\right)
+\frac{N}{\beta}\frac{\partial\ln Z_1}{\partial y}\mathrm dy.
$$

乘以 $\beta$，并把 $Z_1$ 的全微分

$$
\mathrm d\ln Z_1
=\frac{\partial\ln Z_1}{\partial\beta}\mathrm d\beta
+\frac{\partial\ln Z_1}{\partial y}\mathrm dy
$$

代入，可得

$$
\beta(\mathrm dU-Y\,\mathrm dy)
=N\,\mathrm d\left(
\ln Z_1-\beta\frac{\partial\ln Z_1}{\partial\beta}
\right). \tag{7.1.11}
$$

因此 $\beta$ 也是 $\mathrm dQ$ 的积分因子。与 $1/T$ 比较，令

$$
\beta=\frac{1}{kT}. \tag{7.1.12}
$$

对两个互为热平衡的系统，总能量守恒使两系统具有共同的 $\beta$，这与它们具有共同温度一致。常数 $k$ 称为**玻耳兹曼常量**：

$$
k=1.381\times10^{-23}\ \mathrm{J\cdot K^{-1}}.
$$

比较式 (7.1.10) 与式 (7.1.11)，积分得

$$
S=Nk\left(
\ln Z_1-\beta\frac{\partial\ln Z_1}{\partial\beta}
\right). \tag{7.1.13}
$$

积分常量已取为零。

### 7.1.4 熵的统计意义与系统类型

由式 (7.1.3)，

$$
\ln Z_1=\ln N+\alpha.
$$

代入式 (7.1.13)：

$$
\begin{aligned}
S
&=k\left(N\ln N+\alpha N+\beta U\right)\\
&=k\left[N\ln N+\sum_l(\alpha+\beta\varepsilon_l)a_l\right].
\end{aligned}
$$

由玻耳兹曼分布，

$$
\alpha+\beta\varepsilon_l
=\ln\frac{\omega_l}{a_l}.
$$

因此

$$
S=k\left(
N\ln N+\sum_l a_l\ln\omega_l
-\sum_l a_l\ln a_l
\right). \tag{7.1.14}
$$

与式 (6.6.4) 比较：

$$
S=k\ln\Omega. \tag{7.1.15}
$$

这就是**玻耳兹曼关系**：宏观状态的熵等于玻耳兹曼常量乘该系统微观状态数的对数。

> **教材脚注**（p.159）：式 (7.1.15) 适用于近独立玻耳兹曼系统；玻色、费米系统见式 (8.1.11)，有相互作用系统见式 (9.3.10)。

式 (7.1.15) 中的 $\Omega$ 是 $\Omega_{\mathrm{M.B.}}$，故式 (7.1.13) 和式 (7.1.15) 直接适用于粒子可分辨的定域系统。对满足经典极限条件的玻色或费米气体，内能与广义力的表达式 (7.1.4)、式 (7.1.7) 仍然适用，但微观状态数须取 $\Omega_{\mathrm{M.B.}}/N!$。于是熵改为

$$
S=Nk\left(
\ln Z_1-\beta\frac{\partial\ln Z_1}{\partial\beta}
\right)-k\ln N! \tag{7.1.13'}
$$

或

$$
S=k\ln\frac{\Omega_{\mathrm{M.B.}}}{N!}. \tag{7.1.15'}
$$

§7.6 将证明：只有后一种形式才满足熵的广延性。

### 7.1.5 自由能与经典配分函数

自由能

$$
F=U-TS.
$$

代入式 (7.1.4) 与式 (7.1.13)：

$$
\begin{aligned}
F
&=-N\frac{\partial\ln Z_1}{\partial\beta}
-NkT\left(
\ln Z_1-\beta\frac{\partial\ln Z_1}{\partial\beta}
\right)\\
&=-NkT\ln Z_1.
\end{aligned}
\tag{7.1.16}
$$

对满足经典极限条件的玻色、费米气体，

$$
F=-NkT\ln Z_1+kT\ln N!. \tag{7.1.16'}
$$

因此，求得配分函数 $Z_1$ 后，可由式 (7.1.4)、式 (7.1.7)、式 (7.1.13) 或式 (7.1.13′) 得到系统的内能、物态方程和熵，从而确定全部平衡性质。

经典统计中，

$$
Z_1=\sum_l e^{-\beta\varepsilon_l}
\frac{\Delta\omega_l}{h_0^r}, \tag{7.1.17}
$$

连续化后

$$
Z_1=\int e^{-\beta\varepsilon}\frac{\mathrm d\omega}{h_0^r}
=\int\cdots\int e^{-\beta\varepsilon(p,q)}
\frac{\mathrm dq_1\cdots\mathrm dq_r
\mathrm dp_1\cdots\mathrm dp_r}{h_0^r}. \tag{7.1.18}
$$

经典分布

$$
a_l=e^{-\alpha-\beta\varepsilon_l}
\frac{\Delta\omega_l}{h_0^r} \tag{7.1.19}
$$

可写为

$$
a_l=\frac{N}{Z_1}e^{-\beta\varepsilon_l}
\frac{\Delta\omega_l}{h_0^r}. \tag{7.1.20}
$$

$h_0$ 在式 (7.1.19) 与 $Z_1$ 的相消中消失，故内能和物态方程与 $h_0$ 的选择无关；熵则可能相差一个常量。绝对熵需要量子理论。

---

## §7.2　理想气体的物态方程

### 7.2.1 配分函数

单原子理想气体可近似看作无内部结构的质点。无外场时，分子能量为

$$
\varepsilon=\frac{1}{2m}\left(p_x^2+p_y^2+p_z^2\right). \tag{7.2.1}
$$

在

$$
\mathrm dx\,\mathrm dy\,\mathrm dz\,
\mathrm dp_x\,\mathrm dp_y\,\mathrm dp_z
$$

范围内的微观状态数为

$$
\frac{\mathrm dx\,\mathrm dy\,\mathrm dz\,
\mathrm dp_x\,\mathrm dp_y\,\mathrm dp_z}{h^3}. \tag{7.2.2}
$$

代入式 (7.1.2)：

$$
Z_1=\frac{1}{h^3}
\int\cdots\int e^{-\frac{\beta}{2m}
(p_x^2+p_y^2+p_z^2)}
\mathrm dx\,\mathrm dy\,\mathrm dz\,
\mathrm dp_x\,\mathrm dp_y\,\mathrm dp_z. \tag{7.2.3}
$$

积分分解为位置积分和三个动量积分：

$$
\begin{aligned}
Z_1
&=\frac{V}{h^3}
\left(\int_{-\infty}^{+\infty}
e^{-\frac{\beta}{2m}p_x^2}\mathrm dp_x\right)
\left(\int_{-\infty}^{+\infty}
e^{-\frac{\beta}{2m}p_y^2}\mathrm dp_y\right)\\
&\quad\times
\left(\int_{-\infty}^{+\infty}
e^{-\frac{\beta}{2m}p_z^2}\mathrm dp_z\right).
\end{aligned}
$$

利用高斯积分可得

$$
Z_1=V\left(\frac{2\pi m}{h^2\beta}\right)^{3/2}. \tag{7.2.4}
$$

### 7.2.2 物态方程与经典极限

由式 (7.1.7)：

$$
p=\frac{N}{\beta}\frac{\partial}{\partial V}\ln Z_1
=\frac{NkT}{V}. \tag{7.2.5}
$$

即

$$
pV=NkT.
$$

这就是理想气体物态方程。对双原子或多原子分子，加上转动、振动等内部能量后，配分函数对 $V$ 的依赖关系仍只来自平动因子，因此物态方程仍是式 (7.2.5)。

由式 (7.1.3)，

$$
e^\alpha=\frac{Z_1}{N}.
$$

将式 (7.2.4) 代入，经典极限条件 $e^\alpha\gg1$ 化为

$$
e^\alpha
=\frac{V}{N}\left(\frac{2\pi mkT}{h^2}\right)^{3/2}\gg1. \tag{7.2.6}
$$

因此，下列情形更易满足经典极限：

- 粒子数密度 $N/V$ 小，即气体稀薄；
- 温度高；
- 分子质量大。

| 气体 | $1\ \mathrm{atm}$ 下的沸点/K | $e^\alpha$ |
| --- | ---: | ---: |
| He | $4.2$ | $7.5$ |
| H₂ | $20.3$ | $1.4\times10^2$ |
| Ne | $27.2$ | $9.3\times10^3$ |
| Ar | $87.4$ | $4.7\times10^5$ |

除极低温的 He 外，普通气体通常满足经典极限。低温下 He 的量子统计偏离会被较强的原子间相互作用掩盖。

引入分子热运动平均动量对应的德布罗意波长

$$
\lambda=\frac{h}{p}
=\frac{h}{\sqrt{2\pi mkT}},
\qquad n=\frac{N}{V}.
$$

式 (7.2.6) 可写成

$$
n\lambda^3\ll1. \tag{7.2.7}
$$

这表示分子德布罗意波的平均热波长远小于分子平均间距，或每个热波长立方体中的平均粒子数远小于 $1$。

---

## §7.3　麦克斯韦速度分布律

### 7.3.1 从玻耳兹曼分布到速度分布

玻耳兹曼分布为

$$
a_l=\omega_l e^{-\alpha-\beta\varepsilon_l}. \tag{7.3.1}
$$

对自由分子的质心平动，

$$
\varepsilon=\frac{1}{2m}\left(p_x^2+p_y^2+p_z^2\right).
$$

在

$$
\mathrm dp_x\,\mathrm dp_y\,\mathrm dp_z
$$

范围内的状态数为

$$
\frac{V}{h^3}\mathrm dp_x\,\mathrm dp_y\,\mathrm dp_z.
$$

故动量在此范围内的分子数为

$$
\frac{V}{h^3}e^{-\alpha-\frac{1}{2mkT}
(p_x^2+p_y^2+p_z^2)}
\mathrm dp_x\,\mathrm dp_y\,\mathrm dp_z. \tag{7.3.2}
$$

由总分子数为 $N$ 定 $\alpha$：

$$
\frac{V}{h^3}
\iiint e^{-\alpha-\frac{1}{2mkT}
(p_x^2+p_y^2+p_z^2)}
\mathrm dp_x\,\mathrm dp_y\,\mathrm dp_z=N. \tag{7.3.3}
$$

积分得

$$
e^{-\alpha}
=\frac{N}{V}\left(\frac{h^2}{2\pi mkT}\right)^{3/2}. \tag{7.3.4}
$$

因此动量在给定范围内的分子数为

$$
N\left(\frac{1}{2\pi mkT}\right)^{3/2}
e^{-\frac{1}{2mkT}(p_x^2+p_y^2+p_z^2)}
\mathrm dp_x\,\mathrm dp_y\,\mathrm dp_z. \tag{7.3.5}
$$

令

$$
p_x=mv_x,\qquad p_y=mv_y,\qquad p_z=mv_z,
$$

则速度在

$$
\mathrm dv_x\,\mathrm dv_y\,\mathrm dv_z
$$

范围内的分子数为

$$
N\left(\frac{m}{2\pi kT}\right)^{3/2}
e^{-\frac{m}{2kT}(v_x^2+v_y^2+v_z^2)}
\mathrm dv_x\,\mathrm dv_y\,\mathrm dv_z. \tag{7.3.6}
$$

以 $n=N/V$ 表示数密度，单位体积内速度在给定范围内的分子数为

$$
f(v_x,v_y,v_z)\,\mathrm dv_x\,\mathrm dv_y\,\mathrm dv_z
=n\left(\frac{m}{2\pi kT}\right)^{3/2}
e^{-\frac{m}{2kT}(v_x^2+v_y^2+v_z^2)}
\mathrm dv_x\,\mathrm dv_y\,\mathrm dv_z. \tag{7.3.7}
$$

归一化条件为

$$
\int_{-\infty}^{+\infty}
\int_{-\infty}^{+\infty}
\int_{-\infty}^{+\infty}
f(v_x,v_y,v_z)
\mathrm dv_x\,\mathrm dv_y\,\mathrm dv_z=n. \tag{7.3.8}
$$

### 7.3.2 速率分布

在速度空间采用球坐标 $(v,\theta,\varphi)$，体积元为

$$
v^2\sin\theta\,\mathrm dv\,\mathrm d\theta\,\mathrm d\varphi.
$$

对方向积分后，单位体积内速率在 $v\sim v+\mathrm dv$ 的分子数为

$$
4\pi n\left(\frac{m}{2\pi kT}\right)^{3/2}
e^{-\frac{m}{2kT}v^2}v^2\,\mathrm dv. \tag{7.3.9}
$$

归一化：

$$
4\pi n\left(\frac{m}{2\pi kT}\right)^{3/2}
\int_0^\infty e^{-\frac{m}{2kT}v^2}v^2\,\mathrm dv=n. \tag{7.3.10}
$$

### 7.3.3 三种特征速率

最概然速率 $v_p$ 使速率分布函数取极大值：

$$
\frac{\mathrm d}{\mathrm dv}
\left(e^{-\frac{m}{2kT}v^2}v^2\right)=0. \tag{7.3.11}
$$

由此

$$
v_p=\sqrt{\frac{2kT}{m}}. \tag{7.3.12}
$$

平均速率：

$$
\begin{aligned}
\bar v
&=4\pi\left(\frac{m}{2\pi kT}\right)^{3/2}
\int_0^\infty
ve^{-\frac{m}{2kT}v^2}v^2\,\mathrm dv\\
&=\sqrt{\frac{8kT}{\pi m}}.
\end{aligned}
\tag{7.3.13}
$$

方均根速率：

$$
\overline{v^2}
=4\pi\left(\frac{m}{2\pi kT}\right)^{3/2}
\int_0^\infty
v^2e^{-\frac{m}{2kT}v^2}v^2\,\mathrm dv
=\frac{3kT}{m},
$$

故

$$
v_{\mathrm s}=\sqrt{\overline{v^2}}
=\sqrt{\frac{3kT}{m}}. \tag{7.3.14}
$$

三者之比为

$$
v_{\mathrm s}:\bar v:v_p
=\sqrt{\frac{3}{2}}:\frac{2}{\sqrt\pi}:1
=1.225:1.128:1.
$$

若 $M=N_{\mathrm A}m$ 是摩尔质量，则 $k/m=R/M$，方均根速率也可写成

$$
v_{\mathrm s}=\sqrt{\frac{3RT}{M}}. \tag{7.3.15}
$$

### 7.3.4 碰壁数与泻流

取器壁面元 $\mathrm dA$，外法线沿 $x$ 轴。在 $\mathrm dt$ 内，能撞到面元且速度为 $(v_x,v_y,v_z)$ 的分子位于底面积 $\mathrm dA$、高 $v_x\mathrm dt$ 的柱体内。碰壁数的微分为

$$
\mathrm d\Gamma=fv_x\,\mathrm dv_x\,\mathrm dv_y\,\mathrm dv_z. \tag{7.3.16}
$$

单位时间单位面积的总碰壁数为

$$
\Gamma
=\int_{-\infty}^{+\infty}\mathrm dv_y
\int_{-\infty}^{+\infty}\mathrm dv_z
\int_0^{+\infty}v_xf\,\mathrm dv_x.
$$

代入麦克斯韦分布：

$$
\Gamma
=n\left(\frac{m}{2\pi kT}\right)^{1/2}
\int_0^\infty
v_xe^{-\frac{m}{2kT}v_x^2}\mathrm dv_x
=n\sqrt{\frac{kT}{2\pi m}}. \tag{7.3.17}
$$

因 $\bar v=\sqrt{8kT/(\pi m)}$，也可写成

$$
\Gamma=\frac14n\bar v. \tag{7.3.18}
$$

若器壁有小孔且孔足够小，容器内平衡分布几乎不受影响，则逸出分子数等于撞到孔上的分子数。这个过程称为**泻流**。

---

## §7.4　能量均分定理

### 7.4.1 定理的表述与证明

> **能量均分定理**：对于处在温度为 $T$ 的平衡状态的经典系统，粒子能量中每一个独立的平方项的平均值等于 $\frac12 kT$。

粒子能量可写成动能与势能之和。动能为动量平方项之和：

$$
\varepsilon_p=\frac12\sum_{i=1}^r a_ip_i^2, \tag{7.4.1}
$$

其中 $a_i>0$，可为广义坐标的函数，但与相应动量无关。

考察 $\frac12a_1p_1^2$ 的平均值：

$$
\begin{aligned}
\overline{\frac12a_1p_1^2}
&=\frac1N\int
\frac12a_1p_1^2e^{-\alpha-\beta\varepsilon}
\frac{\mathrm dq_1\cdots\mathrm dq_r
\mathrm dp_1\cdots\mathrm dp_r}{h_0^r}\\
&=\frac1{Z_1}\int
\frac12a_1p_1^2e^{-\beta\varepsilon}
\frac{\mathrm dq_1\cdots\mathrm dq_r
\mathrm dp_1\cdots\mathrm dp_r}{h_0^r}.
\end{aligned}
$$

对 $p_1$ 分部积分：

$$
\begin{aligned}
\int_{-\infty}^{+\infty}
\frac12a_1p_1^2e^{-\frac{\beta}{2}a_1p_1^2}\mathrm dp_1
&=\left.
-\frac{p_1}{2\beta}
e^{-\frac{\beta}{2}a_1p_1^2}
\right|_{-\infty}^{+\infty}\\
&\quad+\frac{1}{2\beta}
\int_{-\infty}^{+\infty}
e^{-\frac{\beta}{2}a_1p_1^2}\mathrm dp_1.
\end{aligned}
$$

边界项为零，故

$$
\overline{\frac12a_1p_1^2}
=\frac{1}{2\beta}
\frac1{Z_1}\int e^{-\beta\varepsilon}
\frac{\mathrm dq_1\cdots\mathrm dq_r
\mathrm dp_1\cdots\mathrm dp_r}{h_0^r}
=\frac{1}{2\beta}
=\frac12kT. \tag{7.4.2}
$$

若势能中有平方项

$$
\varepsilon_q=\frac12\sum_{i=1}^{r'}b_iq_i^2
+\varepsilon'_q(q_{r'+1},\cdots,q_r), \tag{7.4.3}
$$

其中 $b_i>0$ 且不依赖于相应的 $q_i$，对 $q_i$ 作同样的分部积分，可得

$$
\overline{\frac12b_iq_i^2}=\frac12kT. \tag{7.4.4}
$$

### 7.4.2 单原子与双原子理想气体

单原子分子的能量

$$
\varepsilon=\frac{1}{2m}
(p_x^2+p_y^2+p_z^2) \tag{7.4.5}
$$

有三个平方项，故

$$
\bar\varepsilon=\frac32kT.
$$

于是

$$
U=\frac32NkT,\qquad
C_V=\frac32Nk,\qquad
C_p=\frac52Nk,
$$

$$
\gamma=\frac{C_p}{C_V}=\frac53=1.667. \tag{7.4.6}
$$

| 气体 | 温度/K | $\gamma$ |
| --- | ---: | ---: |
| He | 291 | 1.660 |
| He | 93 | 1.673 |
| Ne | 292 | 1.642 |
| Ar | 288 | 1.650 |
| Ar | 93 | 1.690 |
| Kr | 292 | 1.689 |
| Xe | 292 | 1.666 |
| Na 蒸气 | 750–926 | 1.680 |
| K 蒸气 | 660–1000 | 1.640 |
| Hg 蒸气 | 548–629 | 1.666 |

双原子分子能量可写成

$$
\begin{aligned}
\varepsilon
&=\frac{1}{2m}(p_x^2+p_y^2+p_z^2)
+\frac{1}{2I}
\left(p_\theta^2+\frac{p_\varphi^2}{\sin^2\theta}\right)\\
&\quad+\frac{1}{2m_\mu}p_r^2+u(r).
\end{aligned}
\tag{7.4.7}
$$

第一项为质心平动，第二项为绕质心转动，第三项为两原子相对运动的动能，$u(r)$ 为相互作用势。

若不考虑相对运动，双原子分子有 $5$ 个平方项，故

$$
\bar\varepsilon=\frac52kT,
\qquad
U=\frac52NkT,
\qquad
C_V=\frac52Nk,
\qquad
C_p=\frac72Nk,
$$

$$
\gamma=\frac{C_p}{C_V}=\frac75=1.40. \tag{7.4.8}
$$

多数双原子气体在常温下与 $1.40$ 相符；低温氢气等偏离则来自量子转动、振动自由度冻结。

### 7.4.3 固体热容与黑体辐射的困难

假设每个原子在平衡位置附近作三个独立的一维简谐振动。一个自由度上的能量为

$$
\varepsilon=\frac1{2m}p^2+\frac12m\omega^2q^2, \tag{7.4.9}
$$

有两个平方项。每个原子有 $3$ 个自由度，故平均能量为

$$
\bar\varepsilon=3kT,
\qquad
U=3NkT,
\qquad
C_V=3Nk. \tag{7.4.10}
$$

这就是杜隆-珀蒂定律，在室温和高温下与实验符合；低温下实验热容趋于零，经典理论失效。

对空腔平衡辐射，把电磁场分解为单色平面波：

$$
\boldsymbol E=\boldsymbol E_0e^{i(\boldsymbol k\cdot\boldsymbol r-\omega t)}, \tag{7.4.11}
$$

其中

$$
\begin{cases}
k_x=\dfrac{2\pi}{L}n_x,\\
k_y=\dfrac{2\pi}{L}n_y,\\
k_z=\dfrac{2\pi}{L}n_z,
\end{cases}
\qquad n_x,n_y,n_z=0,\pm1,\pm2,\cdots \tag{7.4.12}
$$

每个波矢有两个独立偏振方向。把式 (7.4.11) 代入

$$
\nabla^2\boldsymbol E-\frac{1}{c^2}
\frac{\partial^2\boldsymbol E}{\partial t^2}=0
$$

得

$$
\omega=ck. \tag{7.4.13}
$$

在 $\mathrm dk_x\mathrm dk_y\mathrm dk_z$ 内，计入两个偏振方向的振动模式数为

$$
\frac{V\,\mathrm dk_x\mathrm dk_y\mathrm dk_z}{4\pi^3}
\times2.
$$

换成圆频率后：

$$
D(\omega)\,\mathrm d\omega
=\frac{V}{\pi^2c^3}\omega^2\,\mathrm d\omega. \tag{7.4.14}
$$

由能量均分，每个振动模式平均能量为 $kT$，故

$$
U_\omega\,\mathrm d\omega
=\frac{V}{\pi^2c^3}\omega^2kT\,\mathrm d\omega. \tag{7.4.15}
$$

总能量

$$
U=\int_0^\infty U_\omega\,\mathrm d\omega
=\frac{VkT}{\pi^2c^3}
\int_0^\infty\omega^2\,\mathrm d\omega
\to\infty. \tag{7.4.16}
$$

这就是经典理论的**紫外灾难**。它表明高频自由度不能继续按 $kT$ 均分，必须用量子理论处理。

---

## §7.5　理想气体的内能和热容

### 7.5.1 配分函数的分解

忽略原子内电子运动，双原子分子的能量可分成平动、振动和转动：

$$
\varepsilon=\varepsilon^{\mathrm t}
+\varepsilon^{\mathrm v}
+\varepsilon^{\mathrm r}. \tag{7.5.1}
$$

相应简并度乘积为

$$
\omega_l=\omega_l^{\mathrm t}
\omega_l^{\mathrm v}
\omega_l^{\mathrm r}.
$$

配分函数因此分解：

$$
\begin{aligned}
Z_1
&=\sum_{t,v,r}
\omega_t^{\mathrm t}\omega_v^{\mathrm v}\omega_r^{\mathrm r}
e^{-\beta(\varepsilon_t^{\mathrm t}
+\varepsilon_v^{\mathrm v}
+\varepsilon_r^{\mathrm r})}\\
&=\left(\sum_t\omega_t^{\mathrm t}
e^{-\beta\varepsilon_t^{\mathrm t}}\right)
\left(\sum_v\omega_v^{\mathrm v}
e^{-\beta\varepsilon_v^{\mathrm v}}\right)
\left(\sum_r\omega_r^{\mathrm r}
e^{-\beta\varepsilon_r^{\mathrm r}}\right)\\
&=Z_1^{\mathrm t}Z_1^{\mathrm v}Z_1^{\mathrm r}.
\end{aligned}
\tag{7.5.2}
$$

内能：

$$
U=-N\frac{\partial}{\partial\beta}\ln Z_1
=U^{\mathrm t}+U^{\mathrm v}+U^{\mathrm r}. \tag{7.5.3}
$$

热容：

$$
C_V=C_V^{\mathrm t}+C_V^{\mathrm v}+C_V^{\mathrm r}. \tag{7.5.4}
$$

### 7.5.2 平动

平动配分函数由式 (7.2.4) 给出：

$$
Z_1^{\mathrm t}
=V\left(\frac{2\pi m}{h^2\beta}\right)^{3/2}.
$$

因此

$$
U^{\mathrm t}
=-N\frac{\partial}{\partial\beta}
\ln Z_1^{\mathrm t}
=\frac{3N}{2\beta}
=\frac32NkT,
\qquad
C_V^{\mathrm t}=\frac32Nk. \tag{7.5.5}
$$

与能量均分定理一致。

### 7.5.3 振动

双原子分子相对振动近似为线性谐振子：

$$
\varepsilon_n=\left(n+\frac12\right)\hbar\omega,
\qquad n=0,1,2,\cdots.
$$

振动配分函数

$$
Z_1^{\mathrm v}
=\sum_{n=0}^{\infty}
e^{-\beta\hbar\omega(n+1/2)}. \tag{7.5.6}
$$

利用

$$
1+x+x^2+\cdots=\frac1{1-x},\qquad |x|<1,
$$

令 $x=e^{-\beta\hbar\omega}$：

$$
Z_1^{\mathrm v}
=\frac{e^{-\beta\hbar\omega/2}}
{1-e^{-\beta\hbar\omega}}. \tag{7.5.7}
$$

内能：

$$
\begin{aligned}
U^{\mathrm v}
&=-N\frac{\partial}{\partial\beta}
\ln Z_1^{\mathrm v}\\
&=\frac{N\hbar\omega}{2}
+\frac{N\hbar\omega}{e^{\beta\hbar\omega}-1}.
\end{aligned}
\tag{7.5.8}
$$

第一项是零点能，第二项是温度 $T$ 下的热激发能。

热容：

$$
C_V^{\mathrm v}
=\left(\frac{\partial U^{\mathrm v}}{\partial T}\right)_V
=Nk\left(\frac{\hbar\omega}{kT}\right)^2
\frac{e^{\hbar\omega/kT}}
{\left(e^{\hbar\omega/kT}-1\right)^2}. \tag{7.5.9}
$$

定义振动特征温度

$$
k\theta_{\mathrm v}=\hbar\omega, \tag{7.5.10}
$$

则

$$
U^{\mathrm v}
=\frac{Nk\theta_{\mathrm v}}2
+\frac{Nk\theta_{\mathrm v}}
{e^{\theta_{\mathrm v}/T}-1},
$$

$$
C_V^{\mathrm v}
=Nk\left(\frac{\theta_{\mathrm v}}T\right)^2
\frac{e^{\theta_{\mathrm v}/T}}
{\left(e^{\theta_{\mathrm v}/T}-1\right)^2}.
$$

| 分子 | $\theta_{\mathrm v}/(10^3\ \mathrm K)$ | 分子 | $\theta_{\mathrm v}/(10^3\ \mathrm K)$ |
| --- | ---: | --- | ---: |
| H₂ | 6.10 | CO | 3.07 |
| N₂ | 3.34 | NO | 2.69 |
| O₂ | 2.23 | HCl | 4.14 |

常温 $T\ll\theta_{\mathrm v}$ 时，

$$
U^{\mathrm v}
\approx\frac{Nk\theta_{\mathrm v}}2
+Nk\theta_{\mathrm v}e^{-\theta_{\mathrm v}/T}, \tag{7.5.8'}
$$

$$
C_V^{\mathrm v}
\approx Nk\left(\frac{\theta_{\mathrm v}}T\right)^2
e^{-\theta_{\mathrm v}/T}. \tag{7.5.9'}
$$

振动自由度被冻结在基态，对常温热容几乎无贡献。

### 7.5.4 异核双原子分子的转动

转动能级

$$
\varepsilon_l^{\mathrm r}
=\frac{l(l+1)\hbar^2}{2I},
\qquad l=0,1,2,\cdots \tag{7.5.11}
$$

简并度 $2l+1$。转动配分函数

$$
Z_1^{\mathrm r}
=\sum_{l=0}^{\infty}(2l+1)
e^{-\frac{l(l+1)\hbar^2}{2IkT}}. \tag{7.5.12}
$$

定义转动特征温度

$$
\frac{\hbar^2}{2I}=k\theta_{\mathrm r}, \tag{7.5.13}
$$

则

$$
Z_1^{\mathrm r}
=\sum_{l=0}^{\infty}(2l+1)
e^{-\frac{\theta_{\mathrm r}}T l(l+1)}. \tag{7.5.12'}
$$

常温通常 $\theta_{\mathrm r}/T\ll1$，求和可由积分代替。令

$$
x=\frac{\theta_{\mathrm r}}T l(l+1),
\qquad
\mathrm dx=(2l+1)\frac{\theta_{\mathrm r}}T\mathrm dl,
$$

得

$$
Z_1^{\mathrm r}
\approx\frac{T}{\theta_{\mathrm r}}
\int_0^\infty e^{-x}\mathrm dx
=\frac{T}{\theta_{\mathrm r}}
=\frac{2I}{\beta\hbar^2}. \tag{7.5.14}
$$

因此

$$
U^{\mathrm r}
=-N\frac{\partial}{\partial\beta}
\ln Z_1^{\mathrm r}=NkT,
\qquad
C_V^{\mathrm r}=Nk. \tag{7.5.15}
$$

| 分子 | $\theta_{\mathrm r}/\mathrm K$ | 分子 | $\theta_{\mathrm r}/\mathrm K$ |
| --- | ---: | --- | ---: |
| H₂ | 85.4 | CO | 2.77 |
| N₂ | 2.86 | NO | 2.42 |
| O₂ | 2.70 | HCl | 15.1 |

### 7.5.5 同核双原子分子与氢的正、仲态

同核双原子分子必须考虑核自旋与全同性。以 H₂ 为例：

| 核自旋状态 | 允许的转动量子数 | 名称 | 正常比例 |
| --- | --- | --- | ---: |
| 平行 | $l=1,3,5,\cdots$（奇数） | 正氢 | $3/4$ |
| 反平行 | $l=0,2,4,\cdots$（偶数） | 仲氢 | $1/4$ |

正氢和仲氢相互转变概率很小，通常可看作非平衡混合物。定义

$$
\begin{cases}
Z_{1o}^{\mathrm r}
=\displaystyle\sum_{l=1,3,\cdots}(2l+1)
e^{-l(l+1)\theta_{\mathrm r}/T},\\[8pt]
Z_{1p}^{\mathrm r}
=\displaystyle\sum_{l=0,2,\cdots}(2l+1)
e^{-l(l+1)\theta_{\mathrm r}/T}.
\end{cases} \tag{7.5.16}
$$

氢的转动内能为

$$
U^{\mathrm r}
=-\frac34N\frac{\partial}{\partial\beta}
\ln Z_{1o}^{\mathrm r}
-\frac14N\frac{\partial}{\partial\beta}
\ln Z_{1p}^{\mathrm r}. \tag{7.5.17}
$$

高温 $T\gg\theta_{\mathrm r}$ 时，奇、偶求和各约占全体的一半，可近似用积分，恢复

$$
C_V^{\mathrm r}=Nk. \tag{7.5.18}
$$

低温时要直接求和。H₂ 的转动热容从低温零逐渐增长到 $Nk$，与实验曲线符合。

### 7.5.6 电子对热容的贡献

原子基态无自旋、轨道角动量且无精细结构时，电子激发能通常为 eV 量级，相应特征温度为 $10^4\sim10^5\ \mathrm K$。常温下电子冻结在基态，对热容无贡献。

若基态发生自旋-轨道分裂，且特征温度与实验温区相当，则电子自由度可对热容有贡献。例如 O 原子的特征温度约 $230\ \mathrm K$ 和 $320\ \mathrm K$，NO 分子存在约 $178\ \mathrm K$ 的能级分裂。

### 7.5.7 经典配分函数的对照

异核双原子分子在经典近似下能量为

$$
\begin{aligned}
\varepsilon
&=\frac{1}{2m}(p_x^2+p_y^2+p_z^2)
+\frac1{2I}\left(p_\theta^2+\frac{p_\varphi^2}{\sin^2\theta}\right)\\
&\quad+\frac1{2m_\mu}
(p_r^2+m_\mu^2\omega^2r^2).
\end{aligned}
\tag{7.5.19}
$$

总配分函数分解为

$$
Z_1=Z_1^{\mathrm t}Z_1^{\mathrm v}Z_1^{\mathrm r},
$$

其中

$$
\begin{cases}
Z_1^{\mathrm t}
=\displaystyle\int
e^{-\frac{\beta}{2m}(p_x^2+p_y^2+p_z^2)}
\frac{\mathrm dx\,\mathrm dy\,\mathrm dz\,\mathrm dp_x\,\mathrm dp_y\,\mathrm dp_z}
{h_0^3},\\[8pt]
Z_1^{\mathrm v}
=\displaystyle\int
e^{-\frac{\beta}{2m_\mu}(p_r^2+m_\mu^2\omega^2r^2)}
\frac{\mathrm dp_r\,\mathrm dr}{h_0},\\[8pt]
Z_1^{\mathrm r}
=\displaystyle\int
e^{-\frac\beta2\left(p_\theta^2+p_\varphi^2/\sin^2\theta\right)}
\frac{\mathrm dp_\theta\,\mathrm dp_\varphi
\mathrm d\theta\,\mathrm d\varphi}{h_0^2}.
\end{cases} \tag{7.5.20}
$$

分别积分：

$$
Z_1^{\mathrm t}
=V\left(\frac{2\pi m}{h_0^2\beta}\right)^{3/2}, \tag{7.5.21}
$$

$$
Z_1^{\mathrm v}
=\left(\frac{2\pi m_\mu}{h_0\beta}\right)^{1/2}
\left(\frac{2\pi m_\mu}{h_0\beta m_\mu^2\omega^2}\right)^{1/2}
=\frac{2\pi}{h_0\beta\omega}, \tag{7.5.22}
$$

$$
\begin{aligned}
Z_1^{\mathrm r}
&=\frac1{h_0^2}
\int_0^{2\pi}\mathrm d\varphi
\int_0^\pi\mathrm d\theta
\int_{-\infty}^{+\infty}
e^{-\frac{\beta}{2I}p_\theta^2}\mathrm dp_\theta\\
&\quad\times
\int_{-\infty}^{+\infty}
e^{-\frac{\beta}{2I\sin^2\theta}p_\varphi^2}
\mathrm dp_\varphi\\
&=\frac1{h_0^2}
\int_0^{2\pi}\mathrm d\varphi
\int_0^\pi\mathrm d\theta
\left(\frac{2\pi I}{\beta}\right)^{1/2}
\left(\frac{2\pi I\sin^2\theta}{\beta}\right)^{1/2}\\
&=\frac{2\pi I}{h_0^2\beta}
\int_0^{2\pi}\mathrm d\varphi
\int_0^\pi\sin\theta\,\mathrm d\theta\\
&=\frac{8\pi^2I}{h_0^2\beta}.
\end{aligned}
\tag{7.5.23}
$$

因此

$$
C_V^{\mathrm t}=\frac32Nk,
\qquad
C_V^{\mathrm v}=Nk,
\qquad
C_V^{\mathrm r}=Nk, \tag{7.5.24}
$$

与能量均分定理的结果一致，且与 $h_0$ 无关。

---

## §7.6　理想气体的熵

### 7.6.1 经典熵与广延性矛盾

由经典式 (7.1.13) 与平动配分函数 (7.2.4)，单原子理想气体的熵为

$$
S=\frac32Nk\ln T+Nk\ln V
+\frac32Nk\left[
1+\ln\left(\frac{2\pi mk}{h_0^2}\right)
\right]. \tag{7.6.1}
$$

其中含有与 $h_0$ 有关的常数。更严重的是，当 $N,V$ 同时加倍时，式 (7.6.1) 中的 $N\ln V$ 变为原来的 $2$ 倍再加交叉项，不满足广延性。

历史上，吉布斯要求把式 (7.6.1) 减去 $k\ln N!$。量子统计建立后，这一修正得到解释：理想气体满足经典极限时，

$$
\Omega_{\mathrm{B.E.}}
\approx\Omega_{\mathrm{F.D.}}
\approx\frac{\Omega_{\mathrm{M.B.}}}{N!},
$$

故应使用

$$
S=k\ln\frac{\Omega_{\mathrm{M.B.}}}{N!}.
$$

### 7.6.2 量子修正后的绝对熵

由单原子气体

$$
Z_1=V\left(\frac{2\pi mkT}{h^2}\right)^{3/2},
\qquad
\ln N!\approx N(\ln N-1),
$$

代入式 (7.1.13′)，得

$$
S=\frac32Nk\ln T
+Nk\ln\frac{V}{N}
+\frac32Nk\left[
\frac53+\ln\left(\frac{2\pi mk}{h^2}\right)
\right]. \tag{7.6.2}
$$

它满足熵的广延性，且给出绝对熵（不含任意相加常量）。

### 7.6.3 萨克尔-特特罗德公式与化学势

把与凝聚相平衡的饱和蒸气近似为理想气体，用物态方程消去 $V/N$，并对每分子使用 $k$，可把式 (7.6.2) 写成

$$
\ln p=\frac52\ln T+\frac52
+\ln\left[
k^{5/2}\left(\frac{2\pi m}{h^2}\right)^{3/2}
\right]
-\frac{S_{\mathrm{vap}}}{Nk}. \tag{7.6.3}
$$

其中 $S_{\mathrm{vap}}$ 是饱和蒸气的熵。凝聚相与蒸气相的熵差满足

$$
S_{\mathrm{vap}}-S_{\mathrm{con}}=\frac{L}{T}, \tag{7.6.4}
$$

$L$ 是相变潜热。低温下 $S_{\mathrm{con}}\ll L/T$，可忽略，于是

$$
\ln p
=-\frac{L}{RT}
+\frac52\ln T+\frac52
+\ln\left[
k^{5/2}\left(\frac{2\pi m}{h^2}\right)^{3/2}
\right]. \tag{7.6.5}
$$

这就是**萨克尔-特特罗德公式**，其与实验蒸气压曲线相符，是式 (7.6.2) 的有力支持。

定义单个分子的化学势

$$
\mu=\left(\frac{\partial F}{\partial N}\right)_{T,V}. \tag{7.6.6}
$$

由式 (7.1.16′)：

$$
\mu=-kT\ln\frac{Z_1}{N}. \tag{7.6.7}
$$

代入式 (7.2.4)：

$$
\mu=kT\ln\left[
\frac NV\left(\frac{h^2}{2\pi mkT}\right)^{3/2}
\right]. \tag{7.6.8}
$$

由经典极限式 (7.2.6)，方括号中的量远小于 $1$，故理想气体的化学势为负。

---

## §7.7　固体热容的爱因斯坦理论

固体中的 $N$ 个原子可看成 $3N$ 个定域量子振子。爱因斯坦假设这 $3N$ 个振子的频率都相同，记为 $\omega$。振子能级为

$$
\varepsilon_n=\hbar\omega\left(n+\frac12\right),
\qquad n=0,1,2,\cdots. \tag{7.7.1}
$$

振子在平衡位置附近定域，可以区分，故遵从玻耳兹曼分布。一个振子的配分函数为

$$
Z_1=\sum_{n=0}^{\infty}
e^{-\beta\hbar\omega(n+1/2)}
=\frac{e^{-\beta\hbar\omega/2}}
{1-e^{-\beta\hbar\omega}}. \tag{7.7.2}
$$

由式 (7.1.4)，$3N$ 个振子的总内能为

$$
U=-3N\frac{\partial}{\partial\beta}\ln Z_1
=3N\frac{\hbar\omega}{2}
+\frac{3N\hbar\omega}
{e^{\beta\hbar\omega}-1}. \tag{7.7.3}
$$

定容热容

$$
C_V
=\left(\frac{\partial U}{\partial T}\right)_V
=3Nk\left(\frac{\hbar\omega}{kT}\right)^2
\frac{e^{\hbar\omega/kT}}
{\left(e^{\hbar\omega/kT}-1\right)^2}. \tag{7.7.4}
$$

定义爱因斯坦特征温度

$$
k\theta_{\mathrm E}=\hbar\omega. \tag{7.7.5}
$$

于是

$$
C_V
=3Nk\left(\frac{\theta_{\mathrm E}}T\right)^2
\frac{e^{\theta_{\mathrm E}/T}}
{\left(e^{\theta_{\mathrm E}/T}-1\right)^2}. \tag{7.7.6}
$$

高温 $T\gg\theta_{\mathrm E}$：

$$
e^{\theta_{\mathrm E}/T}-1\approx\frac{\theta_{\mathrm E}}T,
$$

故

$$
C_V\approx3Nk. \tag{7.7.7}
$$

恢复杜隆-珀蒂定律。

低温 $T\ll\theta_{\mathrm E}$：

$$
e^{\theta_{\mathrm E}/T}-1\approx e^{\theta_{\mathrm E}/T},
$$

故

$$
C_V
\approx3Nk\left(\frac{\theta_{\mathrm E}}T\right)^2
e^{-\theta_{\mathrm E}/T}. \tag{7.7.8}
$$

$T\to0$ 时 $C_V\to0$，定性解释了热容随温度降低而减小。

但实验低温热容按 $T^3$ 趋于零，而式 (7.7.8) 为指数衰减，量上不符。原因是爱因斯坦理论把 $3N$ 个振子都设成同一频率，过于简化；§9.7 的德拜模型将进一步改进。

---

## §7.8　顺磁性固体

### 7.8.1 二能级磁矩模型

磁性离子定域在晶体格点上，彼此距离足够远。只讨论总角动量量子数 $J=1/2$ 的简单情形。磁矩 $\mu$ 在外磁场中的两个能量为

$$
-\mu B\quad\text{和}\quad+\mu B.
$$

配分函数

$$
Z_1=e^{\beta\mu B}+e^{-\beta\mu B}
=2\cosh\left(\frac{\mu B}{kT}\right). \tag{7.8.1}
$$

体积内单位体积的磁矩数记为 $n$。由式 (7.1.6) 的磁学对应关系，磁化强度为

$$
\mathcal M
=\frac n\beta\frac{\partial}{\partial B}\ln Z_1. \tag{7.8.2}
$$

代入式 (7.8.1)：

$$
\begin{aligned}
\mathcal M
&=n\mu
\frac{e^{\beta\mu B}-e^{-\beta\mu B}}
{e^{\beta\mu B}+e^{-\beta\mu B}}
\\
&=n\mu\tanh\left(\frac{\mu B}{kT}\right).
\end{aligned}
\tag{7.8.3}
$$

### 7.8.2 弱场与强场极限

弱场或高温，$\mu B/kT\ll1$：

$$
\tanh\left(\frac{\mu B}{kT}\right)
\approx\frac{\mu B}{kT}.
$$

因此

$$
\mathcal M
=\frac{n\mu^2}{kT}B
=\chi\mathcal H, \tag{7.8.4}
$$

其中 $B=\mu_0\mathcal H$，磁化率为

$$
\chi=\frac{n\mu^2\mu_0}{kT}.
$$

这就是**居里定律**。

强场或低温，$\mu B/kT\gg1$：

$$
\tanh\left(\frac{\mu B}{kT}\right)\to1,
$$

故

$$
\mathcal M=n\mu. \tag{7.8.5}
$$

几乎所有磁矩都沿外场方向，磁化达到饱和。

### 7.8.3 内能与熵

顺磁固体单位体积的内能为

$$
u=-n\frac{\partial}{\partial\beta}\ln Z_1
=-n\mu B\tanh\left(\frac{\mu B}{kT}\right)
=-\mathcal M B. \tag{7.8.6}
$$

这是磁矩在外场中的势能。

由式 (7.1.13)，单位体积的熵为

$$
s=nk\left[
\ln2+\ln\cosh\left(\frac{\mu B}{kT}\right)
-\frac{\mu B}{kT}
\tanh\left(\frac{\mu B}{kT}\right)
\right]. \tag{7.8.7}
$$

弱场或高温 $x=\mu B/kT\ll1$：

$$
\ln\cosh x\approx\frac12x^2,
\qquad
\tanh x\approx x.
$$

代入得

$$
s\approx nk\ln2
=k\ln2^n. \tag{7.8.8}
$$

这对应每个磁矩有向上、向下两个等概率状态，系统微观状态数为 $2^n$。

强场或低温 $x\gg1$：

$$
\cosh x\approx\frac12e^x,
\qquad
\tanh x\approx1,
$$

故

$$
s\approx0. \tag{7.8.9}
$$

所有磁矩都沿外场方向，微观状态数为 $1$。

> **教材脚注**（p.179）：磁介质中磁矩的势能为 $-\boldsymbol\mu\cdot\boldsymbol B$，其中 $\boldsymbol B$ 是磁矩所在位置的局域场。顺磁体其他磁矩的贡献很小，故可近似取外部磁感应强度。

对核自旋系统，结论形式相同，但核磁矩约为电子磁矩的 $m_e/m_p\sim1/2000$。要达到相同的磁化比例，核自旋系统需要更低的温度。

---

## *§7.9　负温度状态

### 7.9.1 温度与熵的热力学关系

根据热力学基本方程，当外参量 $y$ 不变时，

$$
\frac1T=\left(\frac{\partial S}{\partial U}\right)_y. \tag{7.9.1}
$$

多数系统能级无上界，微观状态数随内能单调增加，故 $S$ 随 $U$ 单调增加，$T>0$。但若系统的能级有上界，熵可能在某段内能范围内随 $U$ 增加而减小，于是

$$
\left(\frac{\partial S}{\partial U}\right)_y<0,
$$

对应**负绝对温度**。

### 7.9.2 核自旋二能级模型

取核自旋量子数 $1/2$。在外磁场中，磁矩的两个能级为

$$
\pm\varepsilon,
\qquad
\varepsilon=\frac{e\hbar B}{2m_p}.
$$

设 $N_+$、$N_-$ 分别为处于 $+\varepsilon$、$-\varepsilon$ 的磁矩数，则

$$
N_+ + N_-=N, \tag{7.9.2}
$$

$$
E=(N_+-N_-)\varepsilon. \tag{7.9.3}
$$

解出

$$
N_+=\frac N2\left(1+\frac{E}{N\varepsilon}\right),
\qquad
N_-=\frac N2\left(1-\frac{E}{N\varepsilon}\right). \tag{7.9.4}
$$

微观状态数

$$
\Omega=\frac{N!}{N_+!N_-!},
$$

故

$$
S=k\ln\Omega
=k\ln\frac{N!}{N_+!N_-!}. \tag{7.9.5}
$$

利用 $\ln m!\approx m(\ln m-1)$，代入式 (7.9.4)：

$$
\begin{aligned}
S
&=Nk\left[
\ln2
-\frac12\left(1+\frac{E}{N\varepsilon}\right)
\ln\left(1+\frac{E}{N\varepsilon}\right)\right.\\
&\qquad\qquad\left.
-\frac12\left(1-\frac{E}{N\varepsilon}\right)
\ln\left(1-\frac{E}{N\varepsilon}\right)
\right].
\end{aligned}
\tag{7.9.6}
$$

### 7.9.3 负温度的出现

对式 (7.9.6) 求导：

$$
\frac1T
=\left(\frac{\partial S}{\partial E}\right)_B
=\frac{k}{2\varepsilon}
\ln\frac{N\varepsilon-E}{N\varepsilon+E}. \tag{7.9.7}
$$

$S(E)$ 关于 $E=0$ 对称。

- $E<0$：$\partial S/\partial E>0$，系统处于正温状态；
- $E=0$：$\partial S/\partial E=0$，$T=\pm\infty$；
- $E>0$：$\partial S/\partial E<0$，系统处于负温状态。

| 状态 | 磁矩取向 | 内能 $E$ | 熵 | 温度 |
| --- | --- | ---: | ---: | --- |
| 基态 | 全部沿场 | $-N\varepsilon$ | $0$ | $+0\ \mathrm K$ |
| 等概率 | 沿场、逆场各半 | $0$ | $Nk\ln2$ | $\pm\infty$ |
| 激发态 | 全部逆场 | $+N\varepsilon$ | $0$ | $-0\ \mathrm K$ |

负温度系统的能量高于正温度系统。它与正温系统接触时，热量从负温系统流向正温系统。按热力学温标的“热”到“冷”顺序：

$$
+0\ \mathrm K,\cdots,+273\ \mathrm K,\cdots,
+\infty\ \mathrm K,\cdots,-273\ \mathrm K,\cdots,-0\ \mathrm K.
$$

因此负温度比任何正温度都“热”，但负温度不是“比绝对零度还低”。

### 7.9.4 出现负温度的条件

1. **粒子能级必须有上界。** 平动、转动、振动自由度能级无上界，不会出现负温；核自旋二能级系统满足条件。
2. **系统与正温系统弱耦合，或内部平衡远快于与外界交换。** 核自旋与晶格的自旋-晶格弛豫时间要长，才能暂时维持负温。

一个系统不能经准静态过程由正温连续变为负温；负温态通常通过快速反转强磁场、使核自旋来不及跟随而制备。

---

## 本章小结

### 一、公式速查

**配分函数与热力学量**

| 内容 | 公式 |
| --- | --- |
| 粒子配分函数 | $Z_1=\sum_l\omega_l e^{-\beta\varepsilon_l}$ (7.1.2) |
| 粒子数与 $\alpha$ | $N=e^{-\alpha}Z_1$ (7.1.3) |
| 内能 | $U=-N\partial_\beta\ln Z_1$ (7.1.4) |
| 广义力 | $Y=-(N/\beta)\partial_y\ln Z_1$ (7.1.6)；$p=(N/\beta)\partial_V\ln Z_1$ (7.1.7) |
| 温度参数 | $\beta=1/(kT)$ (7.1.12) |
| 熵 | $S=Nk(\ln Z_1-\beta\partial_\beta\ln Z_1)$ (7.1.13)；气体须再减 $k\ln N!$ (7.1.13′) |
| 玻耳兹曼关系 | $S=k\ln\Omega$ (7.1.15)；气体为 $S=k\ln(\Omega_{\mathrm{M.B.}}/N!)$ (7.1.15′) |
| 自由能 | $F=-NkT\ln Z_1$ (7.1.16)；气体再加 $kT\ln N!$ (7.1.16′) |
| 经典配分函数 | 式 (7.1.17)–(7.1.20) |

**理想气体**

| 内容 | 公式 |
| --- | --- |
| 单原子分子配分函数 | $Z_1=V(2\pi m/(h^2\beta))^{3/2}$ (7.2.4) |
| 物态方程 | $p=NkT/V$ (7.2.5) |
| 经典极限 | $e^\alpha=(V/N)(2\pi mkT/h^2)^{3/2}\gg1$ (7.2.6)；$n\lambda^3\ll1$ (7.2.7) |
| 速度分布 | 式 (7.3.6)；数密度形式 (7.3.7) |
| 速率分布 | $4\pi n(m/(2\pi kT))^{3/2}e^{-mv^2/(2kT)}v^2\mathrm dv$ (7.3.9) |
| 特征速率 | $v_p=\sqrt{2kT/m}$ (7.3.12)；$\bar v=\sqrt{8kT/(\pi m)}$ (7.3.13)；$v_s=\sqrt{3kT/m}$ (7.3.14) |
| 碰壁数 | $\Gamma=n\sqrt{kT/(2\pi m)}=\frac14n\bar v$ (7.3.17)(7.3.18) |

**能量均分与理想气体热容**

| 内容 | 公式或结果 |
| --- | --- |
| 能量均分 | 每个独立平方项平均值 $kT/2$ (7.4.2)(7.4.4) |
| 单原子气体 | $U=3NkT/2,\ C_V=3Nk/2,\ C_p=5Nk/2,\ \gamma=5/3$ (7.4.5)(7.4.6) |
| 双原子气体（忽略相对运动） | $U=5NkT/2,\ C_V=5Nk/2,\ C_p=7Nk/2,\ \gamma=7/5$ (7.4.7)(7.4.8) |
| 经典固体 | $U=3NkT,\ C_V=3Nk$ (7.4.10) |
| 黑体辐射 | $D(\omega)=V\omega^2/(\pi^2c^3)$ (7.4.14)；瑞利-金斯公式 (7.4.15)；总能量发散 (7.4.16) |

**振动、转动与熵**

| 内容 | 公式 |
| --- | --- |
| 配分函数分解 | $Z_1=Z_1^{\mathrm t}Z_1^{\mathrm v}Z_1^{\mathrm r}$ (7.5.2) |
| 振动配分函数 | $Z_1^{\mathrm v}=e^{-\beta\hbar\omega/2}/(1-e^{-\beta\hbar\omega})$ (7.5.7) |
| 振动内能 / 热容 | (7.5.8)(7.5.9)；特征温度 $k\theta_{\mathrm v}=\hbar\omega$ (7.5.10) |
| 转动配分函数 | 式 (7.5.12)；高温 $Z_1^{\mathrm r}\approx T/\theta_{\mathrm r}=2I/(\beta\hbar^2)$ (7.5.14) |
| 转动内能 / 热容 | $U^{\mathrm r}=NkT,\ C_V^{\mathrm r}=Nk$ (7.5.15) |
| 氢的正、仲态 | 式 (7.5.16)(7.5.17) |
| 经典理想气体熵 | 式 (7.6.1)（不广延）；量子修正后 (7.6.2) |
| 萨克尔-特特罗德公式 | (7.6.5) |
| 化学势 | $\mu=-kT\ln(Z_1/N)$ (7.6.7)；理想气体式 (7.6.8) |

**爱因斯坦固体与顺磁固体**

| 内容 | 公式 |
| --- | --- |
| 爱因斯坦配分函数 | $Z_1=e^{-\beta\hbar\omega/2}/(1-e^{-\beta\hbar\omega})$ (7.7.2) |
| 内能 | $U=3N\hbar\omega/2+3N\hbar\omega/(e^{\beta\hbar\omega}-1)$ (7.7.3) |
| 热容 | (7.7.4)；特征温度 (7.7.5)；无量纲式 (7.7.6) |
| 高温 / 低温极限 | $3Nk$ (7.7.7)；指数衰减 (7.7.8) |
| 顺磁配分函数 | $Z_1=e^{\beta\mu B}+e^{-\beta\mu B}$ (7.8.1) |
| 磁化强度 | $\mathcal M=n\mu\tanh(\mu B/kT)$ (7.8.3) |
| 弱场 / 强场 | 居里定律 (7.8.4)；饱和磁化 (7.8.5) |
| 能量与熵 | (7.8.6)(7.8.7)；弱场 $s=nk\ln2$ (7.8.8)；强场 $s=0$ (7.8.9) |

**负温度（*§7.9）**

| 内容 | 公式 |
| --- | --- |
| 温度定义 | $1/T=(\partial S/\partial U)_y$ (7.9.1) |
| 二能级约束 | (7.9.2)(7.9.3)；占据数 (7.9.4) |
| 熵 | (7.9.5)(7.9.6) |
| 温度 | $1/T=(k/2\varepsilon)\ln[(N\varepsilon-E)/(N\varepsilon+E)]$ (7.9.7) |

### 二、概念脉络

```text
玻耳兹曼分布 a_l = ω_l exp(-α-βε)
        │
        ▼
引入粒子配分函数 Z₁ = Σ ω_l exp(-βε_l)
        │
        ├──► N = exp(-α) Z₁
        ├──► U = -N ∂lnZ₁/∂β
        ├──► Y = -(N/β) ∂lnZ₁/∂y
        ├──► S = Nk(lnZ₁ - β∂lnZ₁/∂β)
        └──► F = -NkT lnZ₁
        │
        ├──► 理想气体：Z₁ ∝ V → pV = NkT
        │     经典极限 nλ³ ≪ 1
        │     └──► Maxwell velocity distribution
        │           ├── v_p < v̄ < v_s
        │           └── wall flux Γ = n v̄/4
        │
        ├──► 经典能量均分：每个平方项 kT/2
        │     ├── 单原子：3/2 Nk
        │     ├── 双原子：5/2 Nk
        │     ├── 固体：3Nk
        │     └── 黑体：ω² 模式无穷 → 紫外灾难
        │
        └──► 量子配分函数处理离散能级
              ├── 振动：θ_v ≫ T 时冻结
              ├── 转动：θ_r ≪ T 时恢复 Nk
              ├── 爱因斯坦固体：低温指数衰减
              └── 二能级顺磁体：tanh(μB/kT)
        │
        ▼
熵与微观状态数
    S = k lnΩ_MB      （定域系统）
    S = k ln(Ω_MB/N!) （经典极限气体）
        │
        ├──► 经典熵不广延；1/N! 修正后广延
        ├──► 绝对熵、Sackur-Tetrode 公式
        └──► 理想气体化学势为负
        │
        ▼
能级有上界 + 弱耦合 → 负温度
    S(E) 在 E>0 后随 E 减小
    -0 K 比 +∞ K 更热
```

### 三、易错点

1. **$Z_1$ 是单粒子配分函数**：$Z_1=\sum_l\omega_l e^{-\beta\varepsilon_l}$，不是整个 $N$ 粒子系统的配分函数。定域系统的总配分函数为 $Z_1^N$，而经典极限气体还要除以 $N!$。
2. **$\beta=1/(kT)$ 由热力学积分因子确定**：不能只凭量纲认定；式 (7.1.11) 表明 $\beta$ 是 $đQ$ 的积分因子，与 $1/T$ 只差普适常数。
3. **广义力公式带负号**：$Y=-(N/\beta)\partial_y\ln Z_1$；对体积功 $Y=-p$，因此 $p=(N/\beta)\partial_V\ln Z_1$。符号错误会把压强变成负值。
4. **定域系统与经典极限气体的熵不同**：两者内能、物态方程相同，但熵分别用 $S=Nk(\ln Z_1-\beta\partial_\beta\ln Z_1)$ 与再减 $k\ln N!$ 的式子。不要漏掉不可分辨性修正。
5. **经典熵不是绝对熵**：$h_0$ 的选择只影响熵的相加常量，但不影响内能、广义力和物态方程；只有量子配分函数才给出绝对熵。
6. **理想气体经典极限不是只由密度决定**：判据是 $n\lambda^3\ll1$，同时依赖数密度、温度和质量。稀薄、高温和较大质量都有利于经典极限；轻粒子在低温下更容易显示量子统计效应。
7. **麦克斯韦速度分布与速率分布不同**：前者是三维速度矢量分布，含 $e^{-mv^2/(2kT)}$；后者还要乘速度空间球壳因子 $4\pi v^2\mathrm dv$。
8. **三种特征速率关系固定**：$v_p< \bar v<v_s$，比为 $\sqrt{3/2}:2/\sqrt\pi:1$。不要漏掉 $\pi$ 或把 $v_s$ 写成 $\sqrt{2kT/m}$。
9. **碰壁数只对 $v_x>0$ 积分**：这是 $\Gamma=\frac14n\bar v$ 中因子 $1/4$ 的来源；对全部速度方向积分会得到零净粒子流。
10. **能量均分只适用于经典平方项**：量子能级间距远大于 $kT$ 时该自由度会“冻结”；黑体高频模式正因如此不能继续贡献 $kT$。
11. **$C_V$ 与 $C_p$ 的差别不能忽略**：理想气体 $C_p=C_V+Nk$，固体实验常给 $C_p$，与理论 $C_V$ 比较前需作热力学修正。
12. **振动零点能不进入热容**：式 (7.5.8) 的第一项 $N\hbar\omega/2$ 是零点能，与温度无关，求导后不出现在 $C_V$ 中。
13. **振动“冻结”是能级间距效应**：常温下 $\theta_{\mathrm v}=10^3\ \mathrm K\gg T$，振子主要处于基态；这不是振子没有零点能。
14. **转动特征温度很小，常温通常可经典处理**：但 H₂ 的 $\theta_{\mathrm r}=85.4\ \mathrm K$，低温必须区分正氢、仲氢并直接求和。
15. **同核双原子分子不能忽略核自旋**：正氢与仲氢的转动量子数选择不同，低温热容不能用简单的 $Nk$。
16. **经典熵的 $N\ln V$ 不广延**：必须减 $k\ln N!$ 或改用 $Nk\ln(V/N)$；这正是全同粒子不可分辨性在热力学中的体现。
17. **爱因斯坦模型只定性解释低温趋势**：它给出指数衰减，而实验固体热容按 $T^3$ 趋于零；需要德拜的频谱模型。
18. **负温度不是“低温”**：负温系统的能量高于正温系统，接触时向正温系统放热；$-0\ \mathrm K$ 才接近负温侧的最高有序极限。
19. **负温度要求能级有上界**：无限能级系统不会出现负温度；二能级核自旋系统的熵在 $E>0$ 时随能量增加而减小，才有负温。

### 四、名词中英对照

| 中文 | 英文 | 中文 | 英文 |
| --- | --- | --- | --- |
| 玻耳兹曼统计 | Boltzmann statistics | 粒子配分函数 | single-particle partition function |
| 配分函数 | partition function | 玻耳兹曼常量 | Boltzmann constant |
| 广义力 | generalized force | 外参量 | external parameter |
| 热力学第一定律 | first law of thermodynamics | 广延性 | extensivity |
| 玻耳兹曼关系 | Boltzmann relation | 自由能 | Helmholtz free energy |
| 经典极限 | classical limit | 德布罗意热波长 | thermal de Broglie wavelength |
| 麦克斯韦速度分布 | Maxwell velocity distribution | 速率分布 | speed distribution |
| 最概然速率 | most probable speed | 平均速率 | mean speed |
| 方均根速率 | root-mean-square speed | 碰壁数 | wall collision rate |
| 泻流 | effusion | 能量均分定理 | equipartition theorem |
| 杜隆-珀蒂定律 | Dulong-Petit law | 紫外灾难 | ultraviolet catastrophe |
| 瑞利-金斯公式 | Rayleigh-Jeans formula | 平动配分函数 | translational partition function |
| 振动配分函数 | vibrational partition function | 转动配分函数 | rotational partition function |
| 振动特征温度 | vibrational characteristic temperature | 转动特征温度 | rotational characteristic temperature |
| 正氢 / 仲氢 | orthohydrogen / parahydrogen | 零点能 | zero-point energy |
| 萨克尔-特特罗德公式 | Sackur-Tetrode formula | 化学势 | chemical potential |
| 爱因斯坦固体模型 | Einstein model of a solid | 爱因斯坦特征温度 | Einstein temperature |
| 顺磁性 | paramagnetism | 顺磁固体 | paramagnetic solid |
| 磁化强度 | magnetization | 磁化率 | magnetic susceptibility |
| 居里定律 | Curie's law | 饱和磁化 | saturation magnetization |
| 负温度 | negative temperature | 自旋-晶格弛豫 | spin-lattice relaxation |
