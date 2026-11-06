---
tags:
  - 热统
  - 统计物理
  - 第十章
created: 2026-10-22
---

# 第十章　涨落理论

> **教材**：汪志诚《热力学·统计物理》（第六版）
> **范围**：§10.1–§10.7（教材 pp. 248–269；其中 *§10.2–*§10.4 为带星号的选读内容）

## 本章地图

| 节 | 标题 | 一句话 |
| --- | --- | --- |
| §10.1 | 涨落的准热力学理论 | 由复合孤立系统的熵偏差得到正则、巨正则系统热力学量涨落的高斯分布 |
| *§10.2 | 临界点邻域序参量的涨落 | 金兹堡-朗道模型把序参量写成空间场，由高斯近似求长波涨落与热容发散 |
| *§10.3 | 序参量涨落的空间关联 | 关联函数及其傅里叶分量描述临界涨落；关联长度 $\xi$ 在临界点发散 |
| *§10.4 | 临界指数的标度关系　普适性 | 标度假设把热力学奇异性归结为关联长度，导出临界指数之间的标度律和普适类 |
| §10.5 | 布朗运动理论 | 朗之万方程给出均方位移 $\overline{x^2}=2kTt/\alpha$；扩散图像导出爱因斯坦关系 |
| §10.6 | 布朗颗粒动量的扩散和时间关联 | 涨落力的 $\delta$ 关联导致动量扩散；动量关联时间由阻尼系数决定 |
| §10.7 | 布朗运动简例 | 热噪声、光学黏胶中的多普勒冷却和磁光陷阱都可化为朗之万方程 |

**主线**：§10.1 用宏观热力学计算平衡涨落 → *§10.2–§10.4 把涨落推广到临界点的空间场与标度分析 → §10.5 以朗之万方程研究随机力驱动下的位移扩散 → §10.6 建立涨落力、动量扩散与耗散的联系 → §10.7 给出电子线路和激光冷却中的实例。

---

## §10.1　涨落的准热力学理论

### 10.1.1 涨落的统计定义

统计物理中的宏观量是相应微观量在给定宏观条件下的平均值：

$$
\overline B=\sum_s B_s\rho_s.
$$

其中 $\rho_s$ 是系统处于微观状态 $s$ 的概率，$B_s$ 是微观量 $B$ 在状态 $s$ 的取值。偏差

$$
\delta B_s=B_s-\overline B
$$

的平均值为零。用偏差平方的平均值

$$
\overline{(B_s-\overline B)^2}
=\sum_s\rho_s(B_s-\overline B)^2
=\overline{B^2}-(\overline B)^2
$$

表示 $B$ 对平均值的涨落。

本章讨论**准热力学涨落理论**：直接求给定宏观条件下热力学量取各种涨落值的概率分布。温度、熵等量不是独立微观量，其涨落应理解为：先取能量的涨落值 $E$ 和体积的涨落值 $V$，再用函数关系 $S=S(E,V)$ 得到相应的熵值，并与 $S(\overline E,\overline V)$ 比较。

### 10.1.2 正则系统的能量和体积涨落

考虑系统与热源接触并达到平衡，二者共同构成孤立复合系统。系统和源的能量、体积变化满足

$$
\Delta E+\Delta E_{\mathrm r}=0,
\qquad
\Delta V+\Delta V_{\mathrm r}=0. \tag{10.1.1}
$$

用 $\overline\Omega^{(0)}$ 表示系统处于平均能量 $\overline E$、平均体积 $\overline V$ 时复合系统的微观状态数。复合系统平均熵为

$$
\overline S^{(0)}
=k\ln\overline\Omega^{(0)}. \tag{10.1.2}
$$

当系统能量和体积偏离平均值 $\Delta E,\Delta V$ 时，复合系统微观状态数和熵分别为 $\Omega^{(0)}$、$S^{(0)}$：

$$
S^{(0)}=k\ln\Omega^{(0)}. \tag{10.1.3}
$$

孤立复合系统的各可达微观状态等概率，因此该涨落出现的概率正比于相应微观状态数：

$$
W\propto\Omega^{(0)}
=\overline\Omega^{(0)}
\exp\left[
\frac{S^{(0)}-\overline S^{(0)}}k
\right],
$$

即

$$
W\propto e^{\Delta S^{(0)}/k}. \tag{10.1.4}
$$

熵具有广延性：

$$
\Delta S^{(0)}
=\Delta S+\Delta S_{\mathrm r}, \tag{10.1.5}
$$

其中 $\Delta S$、$\Delta S_{\mathrm r}$ 分别是系统和源的熵偏差。

源很大，其温度 $T$ 和压强 $p$ 在涨落中可视为常数。对源应用热力学基本方程，

$$
\Delta S_{\mathrm r}
=\frac{\Delta E_{\mathrm r}+p\Delta V_{\mathrm r}}T.
$$

利用式 (10.1.1)：

$$
\Delta S_{\mathrm r}
=-\frac{\Delta E+p\Delta V}{T}. \tag{10.1.6}
$$

代回式 (10.1.4)：

$$
W\propto
\exp\left[
-\frac{\Delta E-T\Delta S+p\Delta V}{kT}
\right]. \tag{10.1.7}
$$

将 $E=E(S,V)$ 在 $(\overline S,\overline V)$ 附近作二阶泰勒展开：

$$
\begin{aligned}
E
&=\overline E
+\left(\frac{\partial E}{\partial S}\right)_0\Delta S
+\left(\frac{\partial E}{\partial V}\right)_0\Delta V\\
&\quad+\frac12\left[
\left(\frac{\partial^2E}{\partial S^2}\right)_0(\Delta S)^2
+2\left(\frac{\partial^2E}{\partial S\partial V}\right)_0\Delta S\Delta V
+\left(\frac{\partial^2E}{\partial V^2}\right)_0(\Delta V)^2
\right].
\end{aligned} \tag{10.1.8}
$$

下标 $0$ 表示在平均值处取值。由于

$$
\left(\frac{\partial E}{\partial S}\right)_0=T,
\qquad
\left(\frac{\partial E}{\partial V}\right)_0=-p,
$$

一阶项互相抵消。整理二阶项：

$$
\begin{aligned}
\Delta E-T\Delta S+p\Delta V
&=\frac12\left[
\Delta S\,\delta T-\Delta V\,\delta p
\right],
\end{aligned} \tag{10.1.9}
$$

其中 $\delta T,\delta p$ 是相对于平均温度和平均压强的偏差。于是

$$
W\propto
\exp\left[
-\frac{\Delta S\Delta T-\Delta V\Delta p}{2kT}
\right]. \tag{10.1.10}
$$

四个偏差量中只有两个独立。若取 $\Delta T,\Delta V$ 为自变量，

$$
\Delta S
=\left(\frac{\partial S}{\partial T}\right)_V\Delta T
+\left(\frac{\partial S}{\partial V}\right)_T\Delta V
=\frac{C_V}{T}\Delta T
+\left(\frac{\partial p}{\partial T}\right)_V\Delta V,
$$

$$
\Delta p
=\left(\frac{\partial p}{\partial T}\right)_V\Delta T
+\left(\frac{\partial p}{\partial V}\right)_T\Delta V.
$$

注意 $\Delta S\Delta T$、$\Delta V\Delta p$ 中交叉项在最后相减时抵消，得到

$$
W\propto
\exp\left[
-\frac{C_V}{2kT^2}(\Delta T)^2
+\frac1{2kT}
\left(\frac{\partial p}{\partial V}\right)_T(\Delta V)^2
\right]. \tag{10.1.11}
$$

> **脚注**（教材 p.250）：若在式 (10.1.4) 中直接用
> $\Delta S^{(0)}=\frac12\delta^2S$，也可直接得到式 (10.1.11)。

因此 $\Delta T,\Delta V$ 是两个独立高斯变量，均值为零，且

$$
\overline{\Delta T\,\Delta V}
=\overline{\Delta T}\;\overline{\Delta V}=0,
\qquad
\overline{(\Delta T)^2}
=\frac{kT^2}{C_V},
\qquad
\overline{(\Delta V)^2}
=-kT\left(\frac{\partial V}{\partial p}\right)_T. \tag{10.1.12}
$$

因为方差恒正，所以平衡稳定性要求

$$
C_V>0,
\qquad
\left(\frac{\partial p}{\partial V}\right)_T<0.
$$

广延量的均方涨落正比于粒子数 $N$，强度量的均方涨落反比于 $N$；但它们各自的**相对涨落**都正比于 $1/\sqrt N$，所以通常极小。临界点附近是例外。

### 10.1.3 磁介质中的涨落

对磁介质，用磁矩 $\mathscr m$ 代替体积，用外场 $\mathscr H$ 代替负压强。类似推导给出

$$
W\propto
\exp\left[
-\frac{C_{\mathscr m}}{2kT^2}(\Delta T)^2
-\frac{\mu_0}{2kT}
\left(\frac{\partial\mathscr H}{\partial\mathscr m}\right)_T
(\Delta\mathscr m)^2
\right]. \tag{10.1.13}
$$

因此

$$
\overline{\Delta T\,\Delta\mathscr m}=0,
\qquad
\overline{(\Delta T)^2}
=\frac{kT^2}{C_{\mathscr m}},
\qquad
\overline{(\Delta\mathscr m)^2}
=\frac{kT}{\mu_0}
\left(\frac{\partial\mathscr m}{\partial\mathscr H}\right)_T. \tag{10.1.14}
$$

### 10.1.4 巨正则系统的涨落

考虑系统、热源和粒子源共同构成孤立复合系统。开系涨落概率为

$$
W\propto
\exp\left[
-\frac{\Delta S\Delta T-\Delta V\Delta p+\Delta N\Delta\mu}{2kT}
\right]. \tag{10.1.15}
$$

固定 $V$，取 $T,N$ 为自变量：

$$
\Delta S
=\left(\frac{\partial S}{\partial T}\right)_{V,N}\Delta T
+\left(\frac{\partial S}{\partial N}\right)_{V,T}\Delta N
=\frac{C_V}{T}\Delta T
-\left(\frac{\partial\mu}{\partial T}\right)_{V,N}\Delta N,
$$

$$
\Delta\mu
=\left(\frac{\partial\mu}{\partial T}\right)_{V,N}\Delta T
+\left(\frac{\partial\mu}{\partial N}\right)_{V,T}\Delta N.
$$

交叉项抵消：

$$
W\propto
\exp\left[
-\frac{C_V}{2kT^2}(\Delta T)^2
+\frac1{2kT}
\left(\frac{\partial\mu}{\partial N}\right)_{V,T}
(\Delta N)^2
\right]. \tag{10.1.16}
$$

因此

$$
\overline{(\Delta T)^2}
=\frac{kT^2}{C_V},
\qquad
\overline{\Delta N\Delta T}=0,
\qquad
\overline{(\Delta N)^2}
=kT\left(\frac{\partial N}{\partial\mu}\right)_{V,T}. \tag{10.1.17}
$$

第二式与巨正则分布直接计算的粒子数涨落一致。

能量偏差为

$$
\Delta E
=\left(\frac{\partial E}{\partial T}\right)_{V,N}\Delta T
+\left(\frac{\partial E}{\partial N}\right)_{V,T}\Delta N.
$$

取平方平均，交叉项为零：

$$
\begin{aligned}
\overline{(\Delta E)^2}
&=\left(\frac{\partial E}{\partial T}\right)_{V,N}^2
\overline{(\Delta T)^2}
+\left(\frac{\partial E}{\partial N}\right)_{V,T}^2
\overline{(\Delta N)^2}\\
&=kT^2C_V
+\left(\frac{\partial E}{\partial N}\right)_{V,T}^2
\overline{(\Delta N)^2}.
\end{aligned}
\tag{10.1.18}
$$

第一项来自热源交换能量；第二项来自粒子源交换粒子。粒子数涨落很大时第二项不可忽略。

---

## *§10.2　临界点邻域序参量的涨落　补充材料

### 10.2.1 金兹堡-朗道模型

朗道理论把序参量看成全空间均匀量，忽略了临界点邻域的强烈空间涨落。金兹堡-朗道模型引入**局域序参量** $\mathscr M(\boldsymbol r)$，其定义是 $\boldsymbol r$ 附近体积元内原子自旋磁矩之和除以体积。体积元须远大于原子间距，同时足够小以保留空间变化。

> 教材从 §10.2 到 §10.4 将玻耳兹曼常量记为 $k_{\mathrm B}$。

以 $k_{\mathrm B}T$ 为能量单位的等效哈密顿量为

$$
\frac{\mathscr H}{k_{\mathrm B}T}
=\int\mathrm d\boldsymbol r\left[
a_0+a_2\mathscr M^2+a_4\mathscr M^4
+c(\nabla\mathscr M)^2
-\mu_0 h(\boldsymbol r)\cdot\mathscr M(\boldsymbol r)
\right]. \tag{10.2.1}
$$

其中 $a_0,a_2,a_4,c$ 是温度的函数，

$$
h(\boldsymbol r)=\frac{\mathscr H(\boldsymbol r)}{k_{\mathrm B}T}.
$$

- $a_0+a_2\mathscr M^2+a_4\mathscr M^4$ 是序参量的局域势能密度；
- $-\mu_0h\cdot\mathscr M$ 是序参量在外场中的势能；
- $c(\nabla\mathscr M)^2$ 表示空间不均匀和近邻交换作用的能量代价。

某一序参量场 $\mathscr M(\boldsymbol r)$ 的概率为

$$
\rho
=\frac1Z
\exp\left[
-\frac{\mathscr H[\mathscr M(\boldsymbol r)]}{k_{\mathrm B}T}
\right], \tag{10.2.2}
$$

配分函数是对全部场构型作泛函积分：

$$
Z
=\int\exp\left[
-\frac{\mathscr H[\mathscr M(\boldsymbol r)]}{k_{\mathrm B}T}
\right]\mathcal D\mathscr M. \tag{10.2.3}
$$

无外场且序参量均匀时，$\nabla\mathscr M=0$。最概然均匀态由

$$
\left\{
\begin{aligned}
&\frac{\partial}{\partial\mathscr M}
\left(a_0+a_2\mathscr M^2+a_4\mathscr M^4\right)=0,\\[2mm]
&\frac{\partial^2}{\partial\mathscr M^2}
\left(a_0+a_2\mathscr M^2+a_4\mathscr M^4\right)>0
\end{aligned}
\right. \tag{10.2.4}
$$

确定，与朗道连续相变理论的稳定性条件相同。

### 10.2.2 序参量的傅里叶分解

令

$$
\mathscr M(\boldsymbol r)
=\frac1{L^{3/2}}\sum_{\boldsymbol k}
e^{i\boldsymbol k\cdot\boldsymbol r}
\mathscr M_{\boldsymbol k}, \tag{10.2.5}
$$

其中

$$
k_{\alpha}=\frac{2\pi}{L}n_{\alpha},
\qquad
\alpha=x,y,z,
\qquad
n_{\alpha}=0,\pm1,\pm2,\cdots.
$$

局域序参量不分辨原子间距尺度，故只保留

$$
k\lesssim\Lambda\sim\frac1{d_0}
$$

的模式。逆变换为

$$
\mathscr M_{\boldsymbol k}
=\frac1{L^{3/2}}
\int\mathrm d\boldsymbol r\,
\mathscr M(\boldsymbol r)e^{-i\boldsymbol k\cdot\boldsymbol r}. \tag{10.2.6}
$$

实场满足

$$
\mathscr M_{\boldsymbol k}^{*}
=\mathscr M_{-\boldsymbol k}. \tag{10.2.7}
$$

零波矢分量给出系统平均磁矩密度：

$$
\mathscr M_{\boldsymbol 0}
=\frac1{L^{3/2}}
\int\mathscr M(\boldsymbol r)\,\mathrm d\boldsymbol r
=L^{3/2}\overline{\mathscr M}. \tag{10.2.8}
$$

### 10.2.3 傅里叶空间中的哈密顿量

将式 (10.2.5) 代入式 (10.2.1)，无外场时

$$
\begin{aligned}
\frac{\mathscr H}{k_{\mathrm B}T}
&=a_0L^3
+\sum_{\boldsymbol k}
\mathscr M_{\boldsymbol k}\cdot\mathscr M_{-\boldsymbol k}
(a_2+ck^2)\\
&\quad+\frac1{L^3}\sum_{\boldsymbol k,\boldsymbol k_1,\boldsymbol k_2}
a_4
(\mathscr M_{\boldsymbol k}\cdot\mathscr M_{\boldsymbol k_1})
(\mathscr M_{\boldsymbol k_2}
\cdot\mathscr M_{-\boldsymbol k-\boldsymbol k_1-\boldsymbol k_2}).
\end{aligned} \tag{10.2.9}
$$

四阶项耦合不同傅里叶模式。精确求配分函数很困难；**高斯近似**忽略四阶项，将各模式看成独立高斯变量。

对全部 $\mathscr M_{\boldsymbol k}$ 作泛函积分：

$$
Z
=\prod_{|\boldsymbol k|<\Lambda}
\int\exp\left[
-\frac{\mathscr H}{k_{\mathrm B}T}
\right]
\mathrm d(\operatorname{Re}\mathscr M_{\boldsymbol k})
\mathrm d(\operatorname{Im}\mathscr M_{\boldsymbol k}). \tag{10.2.10}
$$

利用实场条件，只需对半空间积分：

$$
Z
=\prod_{\substack{|\boldsymbol k|<\Lambda\\k_z\geqslant0}}
\int\exp\left[
-\frac{\mathscr H}{k_{\mathrm B}T}
\right]
\mathrm d(\operatorname{Re}\mathscr M_{\boldsymbol k})
\mathrm d(\operatorname{Im}\mathscr M_{\boldsymbol k}). \tag{10.2.11}
$$

若对完整球体积分，则最后取平方根。

### 10.2.4 顺磁相中的序参量涨落

在 $T>T_c$ 的顺磁相，$\mathscr M_{\boldsymbol 0}=0$。高斯近似下

$$
\frac{\mathscr H}{k_{\mathrm B}T}
=\frac{\mathscr H(\overline{\mathscr M})}{k_{\mathrm B}T}
+\sum_{|\boldsymbol k|<\Lambda}
|\mathscr M_{\boldsymbol k}|^2(a_2+ck^2), \tag{10.2.12}
$$

且

$$
\frac{\mathscr H(\overline{\mathscr M})}{k_{\mathrm B}T}
=a_0L^3.
$$

各模式平均值为零：

$$
\overline{\mathscr M_{\boldsymbol k}}=0. \tag{10.2.13}
$$

高斯积分给出

$$
\overline{|\mathscr M_{\boldsymbol k}|^2}
=\frac1{2(a_2+ck^2)}. \tag{10.2.14}
$$

由朗道理论

$$
a_2=a_2'(T-T_c),
$$

因此

$$
\overline{|\mathscr M_{\boldsymbol k}|^2}
=\frac1{2\left[a_2'(T-T_c)+ck^2\right]}. \tag{10.2.14'}
$$

当 $T\to T_c$ 且 $k\to0$ 时，长波序参量涨落异常增强。

完成高斯积分：

$$
Z
=e^{-a_0L^3}
\prod_{\substack{|\boldsymbol k|<\Lambda\\k_z\geqslant0}}
\frac{\pi}{2(a_2+ck^2)}. \tag{10.2.15}
$$

自由能：

$$
F
=-k_{\mathrm B}T\ln Z
=k_{\mathrm B}Ta_0L^3
-\frac{k_{\mathrm B}T}{2}
\sum_{|\boldsymbol k|<\Lambda}
\ln\frac{\pi}{2(a_2+ck^2)}. \tag{10.2.16}
$$

热容

$$
C=-T\frac{\partial^2F}{\partial T^2}.
$$

把求和换成积分并只保留最奇异项：

$$
C
=\frac{k_{\mathrm B}(a_2'T)^2L^3}{4\pi^2}
\int_0^\Lambda
\frac{k^2\,\mathrm dk}{(a_2+ck^2)^2}
+\text{奇异性较低的项}. \tag{10.2.17}
$$

三维 $k$ 空间求和与积分的关系为

$$
\sum_{|\boldsymbol k|<\Lambda}
\longrightarrow
\frac{L^3}{8\pi^3}
\int_0^\Lambda4\pi k^2\,\mathrm dk. \tag{10.2.18}
$$

### 10.2.5 关联长度与热容发散

令

$$
k=\frac{k'}\xi, \tag{10.2.19}
$$

其中

$$
\xi
=\left[
\frac{c}{a_2'(T-T_c)}
\right]^{1/2}. \tag{10.2.20}
$$

则

$$
C
\equiv\frac{k_{\mathrm B}(a_2'T/c)^2L^3}{4\pi^2}\xi
\int_0^{\Lambda\xi}
\frac{k'^2\,\mathrm dk'}{(1+k'^2)^2}
+\text{奇异性较低的项}
\equiv C_0\xi
+\text{奇异性较低的项}. \tag{10.2.21}
$$

当 $T\to T_c$ 时 $\xi\to\infty$，因此热容发散。物理原因是：模式 $\boldsymbol k$ 的涨落正比于 $(a_2+ck^2)^{-1}$，长波模式在临界点附近异常增强。

在 $d$ 维空间，

$$
\sum_{|\boldsymbol k|<\Lambda}
\longrightarrow
\frac{L^d}{(2\pi)^d}
\int_0^\Lambda k^{d-1}\mathrm dk
\int\mathrm d\Omega_d. \tag{10.2.22}
$$

由式 (10.2.21) 的量纲和积分行为：

$$
C\sim\xi^{4-d}. \tag{10.2.23}
$$

顺磁相 $T>T_c$：

$$
C\sim(T-T_c)^{-\alpha},
\qquad
\alpha=2-\frac d2. \tag{10.2.23'}
$$

铁磁相 $T<T_c$ 得到相同幂次，比例系数不同：

$$
C\sim\xi^{4-d}, \tag{10.2.24}
$$

$$
C\sim(T_c-T)^{-\alpha'},
\qquad
\alpha'=2-\frac d2. \tag{10.2.24'}
$$

因此 $d<4$ 时高斯近似已经显示热容在临界点发散；$d\geqslant4$ 时高斯近似给出有限或较弱奇异。

---

## *§10.3　序参量涨落的空间关联　补充材料

### 10.3.1 关联函数

局域序参量的偏差为

$$
\mathscr M(\boldsymbol r)
-\overline{\mathscr M(\boldsymbol r)}.
$$

其平均值为零：

$$
\overline{\mathscr M(\boldsymbol r)-\overline{\mathscr M(\boldsymbol r)}}=0. \tag{10.3.1}
$$

定义关联函数

$$
C(\boldsymbol r,\boldsymbol r')
=\overline{
\left[\mathscr M(\boldsymbol r)
-\overline{\mathscr M(\boldsymbol r)}\right]
\cdot
\left[\mathscr M(\boldsymbol r')
-\overline{\mathscr M(\boldsymbol r')}\right]
}. \tag{10.3.2}
$$

同一点：

$$
C(\boldsymbol r,\boldsymbol r)
=\overline{
\left[\mathscr M(\boldsymbol r)
-\overline{\mathscr M(\boldsymbol r)}\right]^2
}. \tag{10.3.3}
$$

若不同地点的涨落独立，

$$
C(\boldsymbol r,\boldsymbol r')=0. \tag{10.3.4}
$$

对均匀各向同性系统，$C$ 只依赖

$$
R=|\boldsymbol R|=|\boldsymbol r-\boldsymbol r'|.
$$

### 10.3.2 从傅里叶分量求关联函数

把序参量偏差展开为

$$
\mathscr M(\boldsymbol r)-\overline{\mathscr M}
=\frac1{L^{3/2}}
\sum_{\substack{\boldsymbol k\\k\neq0}}
\mathscr M_{\boldsymbol k}e^{i\boldsymbol k\cdot\boldsymbol r}. \tag{10.3.5}
$$

逆变换：

$$
\mathscr M_{\boldsymbol k}
=\frac1{L^{3/2}}
\int\left[
\mathscr M(\boldsymbol r)-\overline{\mathscr M}
\right]e^{-i\boldsymbol k\cdot\boldsymbol r}
\mathrm d\boldsymbol r. \tag{10.3.6}
$$

实场给出 $\mathscr M_{-\boldsymbol k}=\mathscr M_{\boldsymbol k}^{*}$，所以

$$
|\mathscr M_{\boldsymbol k}|^2
=\frac1{L^3}
\int\mathrm d\boldsymbol r\,\mathrm d\boldsymbol r'
\left[\mathscr M(\boldsymbol r)-\overline{\mathscr M}\right]
\left[\mathscr M(\boldsymbol r')-\overline{\mathscr M}\right]
e^{-i\boldsymbol k\cdot(\boldsymbol r-\boldsymbol r')}. \tag{10.3.7}
$$

取系综平均：

$$
\begin{aligned}
\overline{|\mathscr M_{\boldsymbol k}|^2}
&=\frac1{L^3}
\int\mathrm d\boldsymbol r\,\mathrm d\boldsymbol r'\,
C(|\boldsymbol r-\boldsymbol r'|)
e^{-i\boldsymbol k\cdot(\boldsymbol r-\boldsymbol r')}\\
&=\frac1{L^3}
\int\mathrm d\boldsymbol r\int\mathrm d\boldsymbol R\,
C(R)e^{-i\boldsymbol k\cdot\boldsymbol R}
=L^{3/2}C(\boldsymbol k).
\end{aligned}
\tag{10.3.8}
$$

其中

$$
C(\boldsymbol k)
=\frac1{L^{3/2}}
\int C(\boldsymbol R)e^{-i\boldsymbol k\cdot\boldsymbol R}
\mathrm d\boldsymbol R. \tag{10.3.9}
$$

逆变换：

$$
C(\boldsymbol R)
=\frac1{L^{3/2}}
\sum_{\boldsymbol k}
C(\boldsymbol k)e^{i\boldsymbol k\cdot\boldsymbol R}. \tag{10.3.10}
$$

代入式 (10.3.8)：

$$
C(\boldsymbol R)
=\frac1{L^3}
\sum_{\boldsymbol k}
\overline{|\mathscr M_{\boldsymbol k}|^2}
e^{i\boldsymbol k\cdot\boldsymbol R}. \tag{10.3.11}
$$

### 10.3.3 高斯近似下的关联函数

将式 (10.2.14) 代入式 (10.3.11)：

$$
C(\boldsymbol R)
=\frac1{L^3}\sum_{\boldsymbol k}
\frac1{2(a_2+ck^2)}e^{i\boldsymbol k\cdot\boldsymbol R}. \tag{10.3.12}
$$

利用

$$
\frac1{L^3}\sum_{\boldsymbol k}
\longrightarrow
\frac1{(2\pi)^3}\int\mathrm d\boldsymbol k,
$$

得到

$$
C(\boldsymbol R)
=\frac1{(2\pi)^3}
\int\frac{e^{i\boldsymbol k\cdot\boldsymbol R}}
{2(a_2+ck^2)}
\mathrm d\boldsymbol k. \tag{10.3.13}
$$

三维积分给出

$$
C(\boldsymbol R)
=\frac1{4\pi c}\frac{e^{-R/\xi}}R, \tag{10.3.14}
$$

其中

$$
\xi
=\left(\frac c{a_2}\right)^{1/2}
=\left[
\frac{c}{a_2'(T-T_c)}
\right]^{1/2}. \tag{10.3.15}
$$

当 $R\ll\xi$ 时涨落显著关联；当 $R\gg\xi$ 时关联按指数衰减。$\xi$ 称为**关联长度**。在临界点 $T\to T_c$ 时 $\xi\to\infty$，系统出现长程关联：

$$
C(R)\sim\frac1R. \tag{10.3.16}
$$

### 10.3.4 关联长度和关联函数的临界指数

实验上关联长度随约化温度

$$
t=\frac{T-T_c}{T_c}
$$

发散：

$$
\xi\sim|t|^{-\nu}. \tag{10.3.17}
$$

三维实验值 $\nu\approx0.6\sim0.7$。

在 $d$ 维空间，关联函数的标度形式为

$$
C(R)\sim\frac1{R^{d-2+\eta}}e^{-R/\xi}, \tag{10.3.18}
$$

其中 $\eta$ 是另一个临界指数，典型值约 $0.05\sim0.1$。临界点处

$$
C(R)\sim\frac1{R^{d-2+\eta}}. \tag{10.3.19}
$$

至此已引入六个临界指数

$$
\alpha,\ \beta,\ \gamma,\ \delta,\ \nu,\ \eta.
$$

下一节讨论它们是否独立。

---

## *§10.4　临界指数的标度关系　普适性　补充材料

### 10.4.1 标度假设

临界点邻域只有一个发散的宏观特征长度——关联长度 $\xi$。当 $\xi$ 远大于晶格常数等微观长度时，微观细节被抹去，热力学奇异部分应由 $\xi$ 唯一决定。这就是**标度假设**。

以下用自然单位，把能量的量纲用 $k_{\mathrm B}T$ 消去。以 $L$ 表示长度量纲。$d$ 维系统的自由能密度

$$
f=\frac{F}{k_{\mathrm B}TL^d}
$$

的量纲为

$$
[f]=L^{-d}. \tag{10.4.1}
$$

由式 (10.3.18)：

$$
[C(R)]=L^{2-d-\eta}. \tag{10.4.2}
$$

由式 (10.3.2)，序参量满足

$$
[\mathscr M]=[C(R)]^{1/2}
=L^{(2-d-\eta)/2}. \tag{10.4.3}
$$

### 10.4.2 磁化率的关联函数表示

磁化率可写成关联函数的积分：

$$
\chi
=\frac{\mu_0}{k_{\mathrm B}T}
\int C(R)\,\mathrm d^dR. \tag{10.4.4}
$$

根据式 (10.1.14)：

$$
\overline{(\Delta\mathscr m)^2}
=\frac{k_{\mathrm B}T}{\mu_0}
\left(\frac{\partial\mathscr m}{\partial\mathscr H}\right)_T. \tag{10.4.5}
$$

左边可写成

$$
\begin{aligned}
\overline{(\Delta\mathscr m)^2}
&=\int
\left[\mathscr M(\boldsymbol r)-\overline{\mathscr M}\right]
\cdot
\left[\mathscr M(\boldsymbol r')-\overline{\mathscr M}\right]
\mathrm d^d r\,\mathrm d^d r'\\
&=\int C(|\boldsymbol r-\boldsymbol r'|)
\mathrm d^d r\,\mathrm d^d r'\\
&=\int\mathrm d^d r\int C(R)\,\mathrm d^dR
=L^d\int C(R)\,\mathrm d^dR.
\end{aligned}
$$

右边为

$$
\frac{k_{\mathrm B}T}{\mu_0}
\left(\frac{\partial\mathscr m}{\partial\mathscr H}\right)_T
=\frac{k_{\mathrm B}T}{\mu_0}L^d
\left(\frac{\partial\overline{\mathscr M}}
{\partial\mathscr H}\right)_T
=\frac{k_{\mathrm B}T}{\mu_0}L^d\chi.
$$

两式相等即得式 (10.4.4)。因此

$$
[k_{\mathrm B}T\chi]=L^{2-\eta}. \tag{10.4.6}
$$

### 10.4.3 各热力学量的长度量纲

磁介质的自由能微分为

$$
\mathrm dF=-S\mathrm dT+\mu_0\mathscr H\,\mathrm d\mathscr m.
$$

因此

$$
\mathscr H
=\frac1{\mu_0}\frac{\partial F}{\partial\mathscr m}
=\frac1{\mu_0L^d}
\frac{\partial F}{\partial\overline{\mathscr M}}.
$$

于是

$$
\left[\frac{\mathscr H}{k_{\mathrm B}T}\right]
=\frac{L^{-d}}{L^{(2-d-\eta)/2}}
=L^{-(2+d-\eta)/2}. \tag{10.4.7}
$$

临界点附近 $\xi$ 是唯一进入量纲分析的长度，自由能密度奇异部分必须满足

$$
f\sim\xi^{-d}\sim|t|^{\nu d}. \tag{10.4.8}
$$

类似地，由式 (10.4.3)、(10.4.6)、(10.4.7)：

$$
\mathscr M\sim\xi^{-(2-d-\eta)/2}
\sim|t|^{-\nu(2-d-\eta)/2}, \tag{10.4.9}
$$

$$
k_{\mathrm B}T\chi
\sim\xi^{2-\eta}
\sim|t|^{-\nu(2-\eta)}, \tag{10.4.10}
$$

$$
\frac{\mathscr H}{k_{\mathrm B}T}
\sim\xi^{-(2+d-\eta)/2}
\sim|t|^{\nu(2+d-\eta)/2}. \tag{10.4.11}
$$

### 10.4.4 标度律

将热容实验规律

$$
C\sim|t|^{-\alpha}
$$

与式 (10.4.8) 比较：

$$
-\alpha=\nu d-2. \tag{10.4.12}
$$

将

$$
\chi\sim|t|^{-\gamma}
$$

与式 (10.4.10) 比较：

$$
\gamma=\nu(2-\eta). \tag{10.4.13}
$$

将

$$
\mathscr M\sim t^\beta
$$

与式 (10.4.9) 比较：

$$
\beta=-\frac{\nu(2-d-\eta)}2. \tag{10.4.14}
$$

由 $\mathscr H\sim\mathscr M^\delta$ 及 $\mathscr M\sim t^\beta$，再与式 (10.4.11) 比较：

$$
\beta\delta=\frac{\nu(2+d-\eta)}2. \tag{10.4.15}
$$

整理后得到四个标度律：

$$
\gamma=\nu(2-\eta). \tag{10.4.16}
$$

$$
\alpha+2\beta+\gamma=2. \tag{10.4.17}
$$

$$
\gamma=\beta(\delta-1). \tag{10.4.18}
$$

$$
\nu d=2-\alpha. \tag{10.4.19}
$$

| 名称 | 标度关系 |
| --- | --- |
| 费希尔（Fisher）标度律 | (10.4.16) |
| 拉什布鲁克（Rushbrooke）标度律 | (10.4.17) |
| 维多姆（Widom）标度律 | (10.4.18) |
| 约瑟夫森（Josephson）标度律 | (10.4.19) |

因此六个临界指数只有两个独立。实验值、平均场、二维伊辛模型和三维伊辛模型的结果满足这些关系。

### 10.4.5 普适性与序参量维数

表 10.4.1　临界指数

| 指数 | 实验值 | 平均场 | 二维伊辛 | 三维伊辛 |
| --- | --- | --- | --- | --- |
| $\alpha$ | $0\sim0.2$ | $0$ | $0$ | $0.12$ |
| $\beta$ | $0.3\sim0.4$ | $1/2$ | $1/8$ | $0.31$ |
| $\gamma$ | $1.2\sim1.4$ | $1$ | $1.75$ | $1.25$ |
| $\delta$ | $4\sim5$ | $3$ | $15$ | $5$ |
| $\nu$ | $0.6\sim0.7$ | $1/2$ | $1$ | $0.64$ |
| $\eta$ | $0\sim0.1$ | $0$ | $0.25$ | $0.05$ |

不同物质系统的临界指数相近，说明临界行为由少数全局性质决定。**普适性假设**：临界行为只由空间维数 $d$ 和序参量维数 $n$ 决定；$(d,n)$ 相同的系统属于同一普适类。

| 系统 | 序参量 | $n$ |
| --- | --- | --- |
| 单轴各向异性铁磁体 | 标量磁化强度 | $1$ |
| 平面各向异性铁磁体 | 平面内二维磁矩 | $2$ |
| 三维各向同性磁体 | 三维磁矩 | $3$ |
| 液-气临界点 | 密度差 $\rho_l-\rho_g$ | $1$ |
| 超导、超流 | 复宏观波函数 $\psi=\psi_0e^{i\varphi}$ | $2$ |

液-气系统的对应关系需注意：顺磁相对应液气不分的状态，铁磁相的两个磁化方向对应液相和气相。超导序参量的模 $\psi_0$ 取非零值、相位 $\varphi$ 取定，破坏了原来的旋转对称性。

---

## §10.5　布朗运动理论

### 10.5.1 朗之万方程

布朗颗粒受到介质分子碰撞。在水平一维投影中，设颗粒质量为 $m$，位置为 $x(t)$，介质分子的净作用力为 $f(t)$，其他外力为 $\mathscr F(t)$。牛顿方程为

$$
m\frac{\mathrm d^2x}{\mathrm d t^2}
=f(t)+\mathscr F(t). \tag{10.5.1}
$$

将 $f(t)$ 分解为平均黏性阻力和随机涨落力。当颗粒速度不大时，黏性阻力为 $-\alpha v$。若布朗颗粒近似为半径 $a$ 的小球，介质黏度为 $\eta$，斯托克斯公式为

$$
\alpha=6\pi a\eta. \tag{10.5.2}
$$

剩余涨落力记作 $F(t)$，满足

$$
\overline{F(t)}=0.
$$

因此运动方程写成

$$
m\frac{\mathrm d^2x}{\mathrm d t^2}
=-\alpha\frac{\mathrm dx}{\mathrm dt}
+F(t)+\mathscr F(t). \tag{10.5.3}
$$

这就是**朗之万方程**。无其他外力时

$$
m\frac{\mathrm d^2x}{\mathrm d t^2}
=-\alpha\frac{\mathrm dx}{\mathrm dt}+F(t). \tag{10.5.4}
$$

### 10.5.2 均方位移

用 $x$ 乘式 (10.5.4)，并利用

$$
x\ddot x
=\frac12\frac{\mathrm d^2}{\mathrm d t^2}x^2-\dot x^2,
$$

得到

$$
\frac12\frac{\mathrm d^2}{\mathrm d t^2}
\overline{mx^2}
-m\overline{\dot x^2}
=-\frac\alpha2\frac{\mathrm d}{\mathrm d t}
\overline{x^2}
+\overline{xF(t)}. \tag{10.5.5}
$$

随机涨落力与位置无关，且均值为零，所以

$$
\overline{xF(t)}
=\overline x\,\overline{F(t)}=0.
$$

达到热平衡后，能量均分给出

$$
\frac12m\overline{\dot x^2}
=\frac12kT.
$$

于是

$$
\frac{\mathrm d^2}{\mathrm d t^2}\overline{x^2}
+\frac\alpha m
\frac{\mathrm d}{\mathrm d t}\overline{x^2}
-\frac{2kT}{m}=0. \tag{10.5.6}
$$

其通解为

$$
\overline{x^2}
=\frac{2kT}{\alpha}t
+C_1e^{-\alpha t/m}+C_2. \tag{10.5.7}
$$

颗粒为半径 $a$、密度 $\rho$ 的小球时，

$$
m=\frac{4\pi a^3\rho}{3},
\qquad
\frac\alpha m=\frac{9\eta}{2a^2\rho}.
$$

典型水中布朗颗粒的 $\alpha/m\sim10^7\ \mathrm{s^{-1}}$，瞬态项在 $t>10^{-6}\ \mathrm s$ 后即可忽略。若 $t=0$ 时全部粒子位于 $x=0$，则 $C_2=0$，长时行为为

$$
\overline{x^2}
=\frac{2kT}{\alpha}t. \tag{10.5.8}
$$

均方位移正比于 $t$，这是随机过程的典型结果。若粒子以恒定速率作弹道运动，则会有 $\overline{x^2}\propto t^2$。

### 10.5.3 扩散图像与爱因斯坦关系

布朗颗粒的密度为 $n(x,t)$，流量为 $J(x,t)$。菲克定律：

$$
J=-D\nabla n. \tag{10.5.9}
$$

连续性方程：

$$
\frac{\partial n}{\partial t}
+\nabla\cdot J=0. \tag{10.5.10}
$$

联立得到扩散方程：

$$
\frac{\partial n}{\partial t}
=D\nabla^2n. \tag{10.5.11}
$$

初始条件

$$
n(x,0)=N\delta(x). \tag{10.5.12}
$$

一维解为高斯分布：

$$
n(x,t)
=\frac{N}{2\sqrt{\pi Dt}}
e^{-x^2/(4Dt)}. \tag{10.5.13}
$$

由此

$$
\overline{x^2}
=\frac1N\int_{-\infty}^{+\infty}
x^2n(x,t)\,\mathrm dx
=2Dt. \tag{10.5.14}
$$

与式 (10.5.8) 比较：

$$
D=\frac{kT}{\alpha}. \tag{10.5.15}
$$

这就是**爱因斯坦关系**：扩散系数与温度成正比，与黏性阻力系数成反比。

---

## §10.6　布朗颗粒动量的扩散和时间关联

### 10.6.1 涨落力的时间关联函数

令

$$
\gamma=\frac\alpha m.
$$

朗之万方程写成动量形式：

$$
\frac{\mathrm dp}{\mathrm d t}
=-\gamma p+F(t). \tag{10.6.1}
$$

定义涨落力的时间关联函数：

$$
\overline{F(t)F(t+\tau)}
=\frac1N\sum_{i=1}^N
F_i(t)F_i(t+\tau). \tag{10.6.2}
$$

若 $\tau$ 大于涨落力的关联时间 $\tau_c$，两个时刻的力互不相关，关联函数趋于零。$\tau_c$ 是极短的微观时间，因此可近似为

$$
\overline{F(t)F(t+\tau)}
=2D_p\delta(\tau). \tag{10.6.3}
$$

这里 $D_p$ 是动量扩散系数。长时间平均值

$$
\langle F(t)F(t+\tau)\rangle
=\lim_{T_0\to\infty}
\frac1{T_0}\int_0^{T_0}
F(t)F(t+\tau)\,\mathrm dt \tag{10.6.4}
$$

等于系综平均值：

$$
\langle F(t)F(t+\tau)\rangle
=\overline{F(t)F(t+\tau)}. \tag{10.6.5}
$$

### 10.6.2 平均动量与动量涨落

式 (10.6.1) 两边乘 $e^{\gamma t}$：

$$
\frac{\mathrm d}{\mathrm d t}
\left(pe^{\gamma t}\right)
=e^{\gamma t}F(t).
$$

从 $0$ 到 $t$ 积分：

$$
p(t)=p(0)e^{-\gamma t}
+e^{-\gamma t}\int_0^tF(\xi)e^{\gamma\xi}\,\mathrm d\xi. \tag{10.6.6}
$$

取系综平均，利用 $\overline{F(\xi)}=0$：

$$
\overline{p(t)}
=\overline{p(0)}e^{-\gamma t}. \tag{10.6.7}
$$

平均动量以时间 $1/\gamma$ 指数衰减。因多次碰撞才显著改变平均动量，通常

$$
\frac1\gamma\gg\tau_c.
$$

动量偏差：

$$
\overline{(\Delta p)^2}
=\overline{\left[p(t)-\overline{p(t)}\right]^2}. \tag{10.6.8}
$$

由式 (10.6.6)、(10.6.7)：

$$
p(t)-\overline{p(t)}
=\int_0^t\mathrm d\xi\,
F(\xi)e^{-\gamma(t-\xi)}. \tag{10.6.9}
$$

因此

$$
\overline{(\Delta p)^2}
=\int_0^t\mathrm d\xi
\int_0^t\mathrm d\xi'\,
\overline{F(\xi)F(\xi')}
e^{-\gamma(t-\xi)}
e^{-\gamma(t-\xi')}. \tag{10.6.10}
$$

代入式 (10.6.3)：

$$
\begin{aligned}
\overline{(\Delta p)^2}
&=2D_p\int_0^t
e^{-2\gamma(t-\xi)}\,\mathrm d\xi\\
&=\frac{D_p}{\gamma}
\left(1-e^{-2\gamma t}\right).
\end{aligned}
\tag{10.6.11}
$$

当 $\tau_c\ll t\ll1/\gamma$ 时，

$$
\overline{(\Delta p)^2}
\approx2D_pt. \tag{10.6.12}
$$

这正是动量空间的扩散。

当 $t\gg1/\gamma$ 时，初始动量影响消失，$\overline{p(t)}=0$，故

$$
\overline{p^2}
=\frac{D_p}{\gamma}. \tag{10.6.13}
$$

由能量均分定理：

$$
\frac{\overline{p^2}}{2m}
=\frac12kT. \tag{10.6.14}
$$

因此

$$
D_p=m\gamma kT=\alpha kT. \tag{10.6.15}
$$

这联系了涨落力的强度、动量扩散系数和黏性耗散，是涨落-耗散定理在布朗运动中的特例。

### 10.6.3 动量时间关联与均方位移

当 $t,t'>1/\gamma$ 时，初始项可忽略。由式 (10.6.6)，

$$
\begin{aligned}
\overline{p(t)p(t')}
&=e^{-\gamma(t+t')}
\int_0^t\mathrm d\xi
\int_0^{t'}\mathrm d\xi'\,
\overline{F(\xi)F(\xi')}
e^{\gamma(\xi+\xi')}\\
&=2D_p
\int_0^t\mathrm d\xi
\int_0^{t'}\mathrm d\xi'\,
\delta(\xi-\xi')
e^{-\gamma(t-\xi)}
e^{-\gamma(t'-\xi')}.
\end{aligned}
\tag{10.6.16}
$$

若 $t>t'$，先对 $\xi$ 积分：

$$
\int_0^t\mathrm d\xi\,
e^{-\gamma(t-\xi)}
\delta(\xi-\xi')
=e^{-\gamma(t-\xi')}.
$$

再对 $\xi'$ 积分：

$$
\overline{p(t)p(t')}
=\frac{D_p}{\gamma}
\left[
e^{-\gamma(t-t')}
-e^{-\gamma(t+t')}
\right].
$$

若 $t<t'$，结果把 $t,t'$ 对调。合并为

$$
\overline{p(t)p(t')}
=\frac{D_p}{\gamma}
\left[
e^{-\gamma|t-t'|}
-e^{-\gamma(t+t')}
\right].
$$

当 $t,t'\gg1/\gamma$：

$$
\overline{p(t)p(t')}
=\frac{D_p}{\gamma}
e^{-\gamma|t-t'|}
=mkTe^{-\gamma|t-t'|}. \tag{10.6.17}
$$

因此动量关联时间为 $1/\gamma$，关联强度为 $D_p/\gamma=mkT$。

位移为

$$
x(t)=\frac1m\int_0^tp(\xi)\,\mathrm d\xi.
$$

利用式 (10.6.17)：

$$
\begin{aligned}
\overline{x^2(t)}
&=\frac1{m^2}\int_0^t\mathrm d\xi
\int_0^t\mathrm d\xi'\,
\overline{p(\xi)p(\xi')}\\
&=\frac{kT}{m}\int_0^t\mathrm d\xi
\int_0^t\mathrm d\xi'\,
e^{-\gamma|\xi-\xi'|}.
\end{aligned}
$$

把积分区域分成 $\xi>\xi'$ 与 $\xi'<\xi$ 两块，或直接利用被积函数关于对角线对称：

$$
\overline{x^2(t)}
=\frac{2kT}{m^2}\int_0^t\mathrm d\xi
\int_0^\xi\mathrm d\xi'\,
e^{-\gamma(\xi-\xi')}.
$$

完成积分：

$$
\overline{x^2(t)}
=2D_pt
\equiv\frac{2kT}{m\gamma}t
=\frac{2kT}{\alpha}t. \tag{10.6.18}
$$

与式 (10.5.8) 完全一致。

---

## §10.7　布朗运动简例

### 10.7.1 热噪声

电路电感为 $L$，外加电动势为 $\mathscr E$，电流为 $i$。导体中电子散射造成等效电压，可分成平均阻尼 $-Ri$ 和涨落电压 $V(t)$：

$$
L\frac{\mathrm di}{\mathrm dt}
=\mathscr E-Ri+V(t). \tag{10.7.1}
$$

把朗之万方程写成速度形式：

$$
m\frac{\mathrm dv}{\mathrm dt}
=\mathscr F-\alpha v+F(t). \tag{10.7.2}
$$

两式在对应关系

$$
L\leftrightarrow m,
\qquad
i\leftrightarrow v,
\qquad
\mathscr E\leftrightarrow\mathscr F,
\qquad
R\leftrightarrow\alpha,
\qquad
V\leftrightarrow F \tag{10.7.3}
$$

下完全等价。

涨落电压的相关时间约等于电子自由碰撞时间，室温下约 $10^{-14}\ \mathrm s$，远小于无线电周期，故可写成

$$
\overline{V(t)V(t+\tau)}
=2RkT\delta(\tau). \tag{10.7.4}
$$

傅里叶变换：

$$
V(t)=\int_{-\infty}^{+\infty}
V(\omega)e^{i\omega t}\,\mathrm d\omega,
$$

$$
V(\omega)=\frac1{2\pi}
\int_{-\infty}^{+\infty}
V(t)e^{-i\omega t}\,\mathrm dt. \tag{10.7.5}
$$

由 $V^*(t)=V(t)$ 得 $V^*(\omega)=V(-\omega)$。计算频谱关联并用式 (10.7.4)：

$$
\overline{V(\omega)V^*(\omega')}
=\frac{kTR}{\pi}
\delta(\omega-\omega'). \tag{10.7.6}
$$

这就是**奈奎斯特定理**。不同频率分量统计独立，每个频率分量的均方涨落与频率无关，表现为白噪声。

### 10.7.2 光学黏胶与多普勒冷却

二能级原子的共振跃迁频率为

$$
\omega_A=\frac{\varepsilon_2-\varepsilon_1}{\hbar}.
$$

频率略低于 $\omega_A$ 的激光从相反方向照射原子。由于多普勒效应，逆着原子运动方向传播的光在原子参考系中频率升高，更接近共振，吸收概率更大。多次吸收-自发辐射循环后，原子受到与速度反向的平均力：

$$
\overline{f_z}=-\alpha v_z.
$$

自发辐射反冲力的平均值为零。于是

$$
m\frac{\mathrm dv_z}{\mathrm dt}
=-\alpha v_z+F_z. \tag{10.7.7}
$$

若三个方向都加入反向传播的激光束，

$$
m\frac{\mathrm dv_i}{\mathrm dt}
=-\alpha v_i+F_i,
\qquad i=x,y,z. \tag{10.7.8}
$$

这组方程就是三维布朗运动的朗之万方程。激光场对原子形成由光子构成的黏性介质，称为**光学黏胶**。

由式 (10.6.15)，光学黏胶中原子达到的平衡温度为

$$
kT=\frac{D_p}{\alpha}. \tag{10.7.9}
$$

估算最低温度得到多普勒极限：

$$
kT_D=\frac{\hbar\Gamma}{2}. \tag{10.7.10}
$$

其中 $\Gamma$ 是激发态自发光速率。卤金属原子的 $T_D$ 约为 $100\ \mu\mathrm K$，如 Na 约 $240\ \mu\mathrm K$、Cs 约 $125\ \mu\mathrm K$。

### 10.7.3 磁光陷阱

光学黏胶能降低速度，却不能长期把原子限制在空间小区域。磁光陷阱加入随坐标线性变化的磁场

$$
B_z=\lambda z.
$$

塞曼效应对激发态 $j=1$ 的三个磁子能级产生能量移动

$$
\Delta\varepsilon_{2,m_j}
=\frac{e\hbar}{2m}B_zm_j.
$$

在 $z>0$ 区域，$m_j=-1$ 的子能级最低；在 $z<0$ 区域，$m_j=+1$ 的子能级最低。两束反向传播、圆偏振相反的激光与不同塞曼子能级共振，使原子总被推回 $z=0$。在小位移区域磁力可写成 $-Kz$，再加多普勒黏性阻力，一维平均方程为

$$
m\frac{\mathrm dv_z}{\mathrm dt}
=-Kz-\alpha v_z+F_z.
$$

三维推广为

$$
m\frac{\mathrm dv_i}{\mathrm dt}
=-K_i-\alpha v_i+F_i(t),
\qquad i=x,y,z. \tag{10.7.11}
$$

其中恢复力把原子束缚在陷阱中心，黏性力降低其平均动能，随机反冲力造成布朗涨落。
这里 $K_i$ 表示第 $i$ 个方向的位置恢复力分量，并不是统一的三维弹簧常量。

---

## 本章小结

### 一、公式速查

**准热力学涨落**

| 内容 | 公式 |
| --- | --- |
| 正则系统涨落概率 | $W\propto e^{-(\Delta E-T\Delta S+p\Delta V)/(kT)}$ (10.1.7)；$W\propto e^{-(\Delta S\Delta T-\Delta V\Delta p)/(2kT)}$ (10.1.10) |
| 温度、体积涨落 | $\overline{(\Delta T)^2}=kT^2/C_V$；$\overline{(\Delta V)^2}=-kT(\partial V/\partial p)_T$ (10.1.12) |
| 巨正则涨落 | $\overline{(\Delta N)^2}=kT(\partial N/\partial\mu)_{V,T}$ (10.1.17) |
| 巨正则能量涨落 | (10.1.18) |

**临界涨落与标度**

| 内容 | 公式 |
| --- | --- |
| 金兹堡-朗道哈密顿量 | (10.2.1) |
| 高斯模式涨落 | $\overline{\|\mathscr M_k\|^2}=[2(a_2+ck^2)]^{-1}$ (10.2.14) |
| 关联长度 | $\xi=[c/a_2'(T-T_c)]^{1/2}$ (10.2.20)(10.3.15) |
| 三维关联函数 | $C(R)=(4\pi cR)^{-1}e^{-R/\xi}$ (10.3.14) |
| 临界关联函数 | $C(R)\sim R^{-(d-2+\eta)}$ (10.3.19) |
| 标度律 | $\gamma=\nu(2-\eta)$；$\alpha+2\beta+\gamma=2$；$\gamma=\beta(\delta-1)$；$\nu d=2-\alpha$ (10.4.16)–(10.4.19) |

**布朗运动**

| 内容 | 公式 |
| --- | --- |
| 朗之万方程 | $m\ddot x=-\alpha\dot x+F(t)$ (10.5.4) |
| 斯托克斯阻力系数 | $\alpha=6\pi a\eta$ (10.5.2) |
| 均方位移 | $\overline{x^2}=2kTt/\alpha=2Dt$ (10.5.8)(10.5.14) |
| 爱因斯坦关系 | $D=kT/\alpha$ (10.5.15) |
| 涨落力关联 | $\overline{F(t)F(t+\tau)}=2D_p\delta(\tau)$ (10.6.3) |
| 动量扩散 | $\overline{(\Delta p)^2}\approx2D_pt$ (10.6.12) |
| 涨落-耗散关系 | $D_p=\alpha kT$ (10.6.15) |
| 动量时间关联 | $\overline{p(t)p(t')}=mkTe^{-\gamma\|t-t'\|}$ (10.6.17) |

### 二、概念脉络

```text
§10.1 平衡态涨落
    复合孤立系统 + 等概率原理
        ├── 正则系综：ΔE, ΔV 的高斯分布
        ├── 巨正则系综：加入 ΔN
        └── 涨落方差与 Cv、压缩率、磁化率相联系
                    │
                    ▼
§10.2 临界点涨落
    均匀朗道理论 → 局域序参量 M(r)
        ├── 金兹堡-朗道泛函
        ├── 傅里叶模式 Mk
        ├── 高斯近似：忽略四阶模式耦合
        └── 长波模式异常增强 → 热容发散
                    │
                    ▼
§10.3 空间关联
    C(r,r') = <δM(r)·δM(r')>
        ├── 傅里叶关联 C(k)
        ├── ξ ~ |t|^(-ν)
        ├── C(R) ~ R^(-d+2-η)
        └── ξ → ∞：临界长程关联
                    │
                    ▼
§10.4 标度与普适性
    ξ 是唯一宏观长度
        ├── 量纲分析 → 指数关系
        ├── 四个标度律
        └── 普适类由 (d,n) 决定

§10.5–§10.7 布朗运动
    朗之万方程 m dv/dt = -αv + F(t)
        ├── 平均粘性阻力 → 耗散
        ├── 随机力 F(t) → 扩散
        ├── 位移均方：2Dt
        └── 爱因斯坦关系：D = kT/α
                    │
                    ▼
    涨落力时间关联
        ├── F(t)F(t+τ) = 2Dp δ(τ)
        ├── 动量扩散 Dp = αkT
        ├── 动量关联时间 1/γ
        └── 位移涨落重新得到 2Dt
                    │
                    ▼
    物理实例
        ├── 热噪声 ↔ LR 电路
        ├── 光学黏胶 ↔ 多普勒冷却
        └── 磁光陷阱 ↔ 恢复力 + 粘性 + 随机力
```

### 三、易错点

1. **涨落不是偏离平均值的简单差**：统计物理中的涨落用均方偏差或概率分布描述；线性偏差的平均值恒为零。
2. **准热力学理论仍假设局域平衡**：它把热力学基本方程用于发生小涨落的状态，因此适合小偏差的二阶展开。
3. **不是所有偏差都独立**：式 (10.1.10) 中有四个偏差量，但只有两个独立。
4. **$C_V>0$ 和 $(\partial p/\partial V)_T<0$ 是涨落概率可归一化的稳定性条件**：方差不能为负。
5. **临界点附近相对涨落不一定小**：通常的 $1/\sqrt N$ 估计在大数极限下成立，但临界点关联长度发散，长波涨落异常增强。
6. **金兹堡-朗道模型不是严格的微观理论**：它是序参量场的唯象泛函模型。
7. **高斯近似忽略的是模式耦合**：四阶项包含四个 $\mathscr M_k$，忽略它就是忽略不同傅里叶模式之间的相互作用。
8. **傅里叶分量不是独立物理粒子**：$\mathscr M_k$ 是序参量空间分布的展开系数；其独立性只在高斯近似下成立。
9. **$k=0$ 是均匀模式**：有限 $k$ 描述空间非均匀涨落；临界长波增强来自 $k\to0$。
10. **关联长度和关联函数不是同一个概念**：$\xi$ 是关联随距离衰减的特征长度；$C(R)$ 是实际的关联强度随距离的函数。
11. **标度律不是说所有系统的临界指数完全相同**：普适性只在同一 $(d,n)$ 普适类内成立。
12. **液-气与铁磁对应不能把液相直接对应顺磁相**：液相和气相对应铁磁的两个磁化方向；不分相的气液状态对应顺磁相。
13. **布朗运动的均方位移是 $2Dt$，不是 $(kT/m)t^2$**：后者是弹道运动；随机碰撞导致 $t$ 的一次幂。
14. **$D=kT/\alpha$ 联系涨落和耗散**：扩散系数小不代表随机力弱；它还取决于阻尼。
15. **朗之万方程不是一件确定力学方程**：$F(t)$ 是随机函数，必须讨论其统计平均值和关联函数。
16. **动量扩散系数 $D_p$ 与空间扩散系数 $D$ 不同**：量纲和物理意义不同，不能混用。
17. **多普勒冷却的最低温度不是绝对零度**：自发辐射反冲给出多普勒极限 $kT_D=\hbar\Gamma/2$。
18. **光学黏胶不是普通流体**：它是激光场形成的有效黏性介质，原子仍受随机反冲力驱动。
19. **磁光陷阱需要速度阻尼和位置恢复力同时存在**：只有阻尼会降低速度，不能长期把原子束缚在空间区域。

### 四、名词中英对照

| 中文 | 英文 | 中文 | 英文 |
| --- | --- | --- | --- |
| 准热力学涨落 | quasi-thermodynamic fluctuation | 广延量 | extensive quantity |
| 强度量 | intensive quantity | 正则系综 | canonical ensemble |
| 巨正则系综 | grand canonical ensemble | 涨落概率 | fluctuation probability |
| 高斯分布 | Gaussian distribution | 金兹堡-朗道模型 | Ginzburg-Landau model |
| 序参量 | order parameter | 局域序参量 | local order parameter |
| 泛函积分 | functional integral | 傅里叶分量 | Fourier component |
| 高斯近似 | Gaussian approximation | 关联函数 | correlation function |
| 空间关联 | spatial correlation | 关联长度 | correlation length |
| 临界指数 | critical exponent | 标度假设 | scaling hypothesis |
| 标度律 | scaling law | 普适性 | universality |
| 普适类 | universality class | 布朗运动 | Brownian motion |
| 朗之万方程 | Langevin equation | 黏性阻力 | viscous drag |
| 涨落力 | random force | 斯托克斯公式 | Stokes formula |
| 均方位移 | mean-square displacement | 扩散系数 | diffusion coefficient |
| 菲克定律 | Fick law | 爱因斯坦关系 | Einstein relation |
| 涨落-耗散定理 | fluctuation-dissipation theorem | 动量扩散 | momentum diffusion |
| 时间关联函数 | time-correlation function | 白噪声 | white noise |
| 热噪声 | thermal noise | 奈奎斯特定理 | Nyquist theorem |
| 光学黏胶 | optical molasses | 多普勒冷却 | Doppler cooling |
| 多普勒极限 | Doppler limit | 磁光陷阱 | magneto-optical trap |
