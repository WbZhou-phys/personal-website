---
tags:
  - 热统
  - 统计物理
  - 第九章
created: 2026-10-22
---

# 第九章　系综理论

> **教材**：汪志诚《热力学·统计物理》（第六版）
> **范围**：§9.1–§9.12（教材 pp. 205–245；其中 *§9.3、*§9.8、*§9.9、*§9.12 为带星号或补充材料）

## 本章地图

| 节 | 标题 | 一句话 |
| --- | --- | --- |
| §9.1 | 相空间　刘维尔定理 | 相空间代表点密度沿哈密顿轨道运动守恒，代表点流满足连续性方程 |
| §9.2 | 微正则系综 | 孤立系统在给定能量壳层内各微观状态等概率，宏观量取系综平均 |
| *§9.3 | 微正则系综理论的热力学公式　补充材料 | 由 $\Omega(N,E,V)$ 确定 $\beta,p/T,\mu/T$，得到微正则理论的热力学公式；补充材料推出理想气体 |
| §9.4 | 正则系综 | 恒温闭系与热源弱耦合，微观状态概率为 $e^{-\beta E_s}/Z$ |
| §9.5 | 正则系综理论的热力学公式 | 由配分函数 $Z(N,V,T)$ 求 $U,Y,p,S,F$，并导出能量涨落 |
| §9.6 | 实际气体的物态方程 | 用位形积分和迈耶集团展开得到第二位力系数与范德瓦耳斯方程 |
| §9.7 | 固体的热容 | 将原子振动分解为简正模式，德拜频谱给出高温 $3Nk$ 与低温 $T^3$ 热容 |
| *§9.8 | 液 $^4\mathrm{He}$ 的性质和朗道超流理论　补充材料 | 二流体模型、声子与旋子元激发解释超流、第二声和低温热容 |
| *§9.9 | 伊辛模型的平均场理论　补充材料 | 用平均场把相互作用自旋化为独立自旋，得到自发磁化、临界温度和临界指数 |
| §9.10 | 巨正则系综 | 恒温、恒化学势开系的状态概率为 $e^{-\alpha N-\beta E_s}/\Xi$ |
| §9.11 | 巨正则系综理论的热力学公式 | 由 $\Xi(\mu,V,T)$ 求平均粒子数、内能、压强、熵和粒子数涨落 |
| *§9.12 | 巨正则系综理论的简单应用　补充材料 | 处理吸附，并严格重导玻色、费米分布及分析其涨落 |

**主线**：相空间与刘维尔定理 → 孤立系用微正则系综 → 恒温闭系用正则系综 → 恒温开放系用巨正则系综。三种系综给出相同的宏观热力学信息；系综理论可处理粒子间相互作用，而不仅是近独立粒子系统。

---

## §9.1　相空间　刘维尔定理

### 9.1.1 系统微观状态的相空间描述

设系统有 $N$ 个全同粒子，每个粒子有 $r$ 个自由度，则系统总自由度为

$$
f=Nr.
$$

若系统含多种粒子，第 $i$ 种粒子的粒子数和自由度为 $N_i,r_i$，则

$$
f=\sum_i N_ir_i.
$$

系统的微观状态由 $f$ 个广义坐标 $q_1,\cdots,q_f$ 和 $f$ 个共轭动量 $p_1,\cdots,p_f$ 确定。以这 $2f$ 个变量构成 **$\Gamma$ 空间**，系统的一个微观状态对应其中一点，即代表点。

保守系统的哈密顿量 $H(q,p)$ 不显含时间，运动方程为

$$
\dot q_i=\frac{\partial H}{\partial p_i},
\qquad
\dot p_i=-\frac{\partial H}{\partial q_i},
\qquad i=1,2,\cdots,f. \tag{9.1.1}
$$

孤立系统能量守恒：

$$
H(q_1,\cdots,q_f;p_1,\cdots,p_f)=E. \tag{9.1.2}
$$

式 (9.1.2) 是 $\Gamma$ 空间中的**能量曲面**，代表点始终位于该曲面上。由一阶常微分方程组 (9.1.1) 的初值唯一性，轨道于同一点只有一条，轨道之间互不相交。

### 9.1.2 系综与分布函数

设想大量结构相同、宏观条件相同而初态不同的系统，各自独立运动。它们的代表点在 $\Gamma$ 空间形成统计系综。体积元

$$
\mathrm d\Omega=\mathrm dq_1\cdots\mathrm dq_f\,
\mathrm dp_1\cdots\mathrm dp_f
$$

内的代表点数，在 $t$ 时刻为

$$
\rho(q_1,\cdots,q_f;p_1,\cdots,p_f;t)\mathrm d\Omega, \tag{9.1.3}
$$

其中 $\rho$ 为**分布函数**。若系综总系统数为 $N$，则

$$
\int\rho(q,p;t)\mathrm d\Omega=N, \tag{9.1.4}
$$

其中 $N$ 不随时间改变。

代表点从 $(q_i,p_i)$ 经 $\mathrm dt$ 到

$$
\left(q_i+\dot q_i\mathrm dt,\ p_i+\dot p_i\mathrm dt\right).
$$

沿轨道跟踪代表点时，密度随时间的全导数为

$$
\frac{\mathrm d\rho}{\mathrm dt}
=\frac{\partial\rho}{\partial t}
+\sum_i\left(
\frac{\partial\rho}{\partial q_i}\dot q_i
+\frac{\partial\rho}{\partial p_i}\dot p_i
\right). \tag{9.1.5}
$$

刘维尔定理为

$$
\frac{\mathrm d\rho}{\mathrm dt}=0. \tag{9.1.6}
$$

### 9.1.3 连续性方程与刘维尔定理的证明

固定体积元 $\mathrm d\Omega$ 内代表点数的增加率为

$$
\frac{\partial\rho}{\partial t}\mathrm dt\,\mathrm d\Omega. \tag{9.1.7}
$$

在 $q_i$ 与 $q_i+\mathrm dq_i$ 一对面上，粒子流分别为 $\rho\dot q_i\,\mathrm dt\,\mathrm dA$ 和

$$
\left[
\rho\dot q_i+\frac{\partial(\rho\dot q_i)}{\partial q_i}\mathrm dq_i
\right]\mathrm dt\,\mathrm dA.
$$

净流入为

$$
-\frac{\partial(\rho\dot q_i)}{\partial q_i}\mathrm dt\,\mathrm d\Omega.
$$

动量方向同样给出

$$
-\frac{\partial(\rho\dot p_i)}{\partial p_i}\mathrm dt\,\mathrm d\Omega.
$$

全部 $2f$ 个面求和并与式 (9.1.7) 比较：

$$
\frac{\partial\rho}{\partial t}
+\sum_i\left[
\frac{\partial(\rho\dot q_i)}{\partial q_i}
+\frac{\partial(\rho\dot p_i)}{\partial p_i}
\right]=0. \tag{9.1.8}
$$

由哈密顿方程，

$$
\frac{\partial\dot q_i}{\partial q_i}
+\frac{\partial\dot p_i}{\partial p_i}
=\frac{\partial^2H}{\partial q_i\partial p_i}
-\frac{\partial^2H}{\partial p_i\partial q_i}
=0.
$$

故

$$
\frac{\partial\rho}{\partial t}
+\sum_i\left(
\frac{\partial\rho}{\partial q_i}\dot q_i
+\frac{\partial\rho}{\partial p_i}\dot p_i
\right)=0. \tag{9.1.9}
$$

代入式 (9.1.5) 即证式 (9.1.6)。刘维尔定理也可写成

$$
\frac{\partial\rho}{\partial t}
=-\sum_i\left(
\frac{\partial\rho}{\partial q_i}\frac{\partial H}{\partial p_i}
-\frac{\partial\rho}{\partial p_i}\frac{\partial H}{\partial q_i}
\right). \tag{9.1.10}
$$

在 $t\to-t$ 下该方程不变，因此纯力学演化是可逆的。刘维尔定理只来自哈密顿力学，尚未引入统计假设。

---

## §9.2　微正则系综

### 9.2.1 系综平均

宏观量是相应微观量在给定宏观条件下全部可能微观状态上的平均。经典情形中，微观状态落在 $\mathrm d\Omega$ 内的概率为

$$
\rho(q,p,t)\mathrm d\Omega, \tag{9.2.1}
$$

归一化条件为

$$
\int\rho(q,p,t)\mathrm d\Omega=1. \tag{9.2.2}
$$

微观量 $B(q,p)$ 的系综平均为

$$
\overline{B(t)}
=\int B(q,p)\rho(q,p,t)\mathrm d\Omega. \tag{9.2.3}
$$

量子情形中，系统处在状态 $s$ 的概率为 $\rho_s(t)$：

$$
\sum_s\rho_s(t)=1, \tag{9.2.4}
$$

$$
\overline{B(t)}=\sum_s\rho_s(t)B_s. \tag{9.2.5}
$$

### 9.2.2 微正则分布

平衡孤立系的概率分布不显含时间。由式 (9.1.10)，

$$
\sum_i\left(
\frac{\partial\rho}{\partial q_i}\frac{\partial H}{\partial p_i}
-\frac{\partial\rho}{\partial p_i}\frac{\partial H}{\partial q_i}
\right)=0. \tag{9.2.6}
$$

若孤立系统的能量在 $E$ 到 $E+\Delta E$ 之间，$\rho$ 只通过能量依赖 $q,p$，则式 (9.2.6) 自动满足。结合等概率原理，微正则分布为

$$
\begin{cases}
\rho(q,p)=\text{常量},&E\leqslant H(q,p)\leqslant E+\Delta E,\\[2mm]
\rho(q,p)=0,&H(q,p)<E,\quad E+\Delta E<H(q,p).
\end{cases} \tag{9.2.7}
$$

即在能量壳层内的相体积上均匀分布。量子微正则分布为

$$
\rho_s=\frac1\Omega, \tag{9.2.8}
$$

其中 $\Omega$ 是能量壳层内可能的微观状态数。

### 9.2.3 微观状态数与等概率假设

对含 $N$ 个粒子、每个粒子自由度 $r$ 的全同粒子系统：

$$
\Omega=\frac1{N!h^{Nr}}
\int_{E\leqslant H(q,p)\leqslant E+\Delta E}\mathrm d\Omega. \tag{9.2.9}
$$

相空间微小格胞体积取 $h^{Nr}$；$N!$ 去除全同粒子交换重复计数。若含多类粒子，

$$
\Omega=
\frac1{\prod_iN_i!\,h^{N_ir_i}}
\int_{E\leqslant H(q,p)\leqslant E+\Delta E}\mathrm d\Omega. \tag{9.2.10}
$$

刘维尔定理本身不能从力学唯一推出微正则分布。等概率原理是平衡统计物理的基本假设，其正确性由推论与实验相符来确认。

若宏观量的相对涨落满足

$$
\frac{\overline{B^2}-(\overline B)^2}{(\overline B)^2}\ll1, \tag{9.2.11}
$$

则最概然值与系综平均值趋于相同。对宏观系统，相对涨落通常为 $N^{-1/2}$ 量级。

长时间平均可写为

$$
\langle B\rangle
=\lim_{T\to\infty}\frac1T
\int_0^TB[q(t),p(t)]\mathrm dt. \tag{9.2.12}
$$

时间平均与系综平均相等的**各态遍历假说**并非对理想保守系统的严格定理。现代处理通常把系统与外界不可避免的弱随机作用纳入考虑。

---

## *§9.3　微正则系综理论的热力学公式　补充材料

### 9.3.1 两个弱耦合系统的热平衡

孤立复合系统 $A^{(0)}=A_1+A_2$。若两子系统微观状态数分别为 $\Omega_1,\Omega_2$，则

$$
\Omega^{(0)}(E_1,E_2)
=\Omega_1(E_1)\Omega_2(E_2). \tag{9.3.1}
$$

两系统只交换能量，保持粒子数和体积不变：

$$
E_1+E_2=E^{(0)}. \tag{9.3.2}
$$

因此

$$
\Omega^{(0)}
=\Omega_1(E_1)\Omega_2(E^{(0)}-E_1). \tag{9.3.3}
$$

最概然能量分配 $\overline E_1,\overline E_2$ 使 $\Omega^{(0)}$ 极大：

$$
\frac{\partial\Omega^{(0)}}{\partial E_1}=0.
$$

利用式 (9.3.3) 和 $\partial E_2/\partial E_1=-1$，除以 $\Omega_1\Omega_2$：

$$
\left(\frac{\partial\ln\Omega_1}{\partial E_1}\right)_{N_1,V_1}
=
\left(\frac{\partial\ln\Omega_2}{\partial E_2}\right)_{N_2,V_2}. \tag{9.3.4}
$$

定义

$$
\beta=\left(\frac{\partial\ln\Omega(N,E,V)}{\partial E}\right)_{N,V}, \tag{9.3.5}
$$

则热平衡条件为

$$
\beta_1=\beta_2. \tag{9.3.6}
$$

热力学给出

$$
\left(\frac{\partial S_1}{\partial U_1}\right)_{N_1,V_1}
=
\left(\frac{\partial S_2}{\partial U_2}\right)_{N_2,V_2}, \tag{9.3.7}
$$

$$
\left(\frac{\partial S}{\partial U}\right)_{N,V}=\frac1T. \tag{9.3.8}
$$

比较可知 $\beta\propto1/T$。令比例常数为 $1/k$：

$$
\beta=\frac1{kT}. \tag{9.3.9}
$$

再比较式 (9.3.4) 和式 (9.3.7)，得到

$$
S=k\ln\Omega. \tag{9.3.10}
$$

此式也称为玻耳兹曼关系。由 $\Omega=\prod_i\Omega_i$ 与 $S=\sum_iS_i$ 可看出，它对由弱耦合部分组成的非平衡系统也成立。

### 9.3.2 体积与粒子交换平衡

若两系统还能交换体积和粒子，则

$$
\left(\frac{\partial\ln\Omega_1}{\partial E_1}\right)_{N_1,V_1}
=
\left(\frac{\partial\ln\Omega_2}{\partial E_2}\right)_{N_2,V_2}, \tag{9.3.11}
$$

$$
\left(\frac{\partial\ln\Omega_1}{\partial V_1}\right)_{N_1,E_1}
=
\left(\frac{\partial\ln\Omega_2}{\partial V_2}\right)_{N_2,E_2}, \tag{9.3.12}
$$

$$
\left(\frac{\partial\ln\Omega_1}{\partial N_1}\right)_{E_1,V_1}
=
\left(\frac{\partial\ln\Omega_2}{\partial N_2}\right)_{E_2,V_2}. \tag{9.3.13}
$$

定义

$$
\gamma=\left(\frac{\partial\ln\Omega}{\partial V}\right)_{N,E}, \tag{9.3.14}
$$

$$
\alpha=\left(\frac{\partial\ln\Omega}{\partial N}\right)_{E,V}. \tag{9.3.15}
$$

则一般平衡条件为

$$
\beta_1=\beta_2,\qquad
\gamma_1=\gamma_2,\qquad
\alpha_1=\alpha_2. \tag{9.3.16}
$$

由式 (9.3.10)：

$$
\mathrm d\ln\Omega
=\beta\,\mathrm dE+\gamma\,\mathrm dV+\alpha\,\mathrm dN.
$$

与开系热力学基本方程

$$
\mathrm dS=\frac{\mathrm dU}{T}+\frac{p}{T}\mathrm dV-\frac{\mu}{T}\mathrm dN
$$

比较，得

$$
\gamma=\frac{p}{kT},\qquad
\alpha=-\frac{\mu}{kT}. \tag{9.3.17}
$$

因此式 (9.3.16) 等价于

$$
T_1=T_2,\qquad p_1=p_2,\qquad \mu_1=\mu_2. \tag{9.3.18}
$$

### 9.3.3 微正则系综的热力学程序

经典理想气体满足 $\Omega(N,E,V)\propto V^N$，故

$$
\frac p{kT}
=\frac{\partial\ln\Omega}{\partial V}
=\frac{\partial}{\partial V}\ln V^N
=\frac NV. \tag{9.3.19}
$$

与 $pV=nRT$ 比较，$k=nR/N=R/N_A$，即玻耳兹曼常量。

微正则系综的计算步骤：

1. 求微观状态数 $\Omega(N,E,V)$；
2. 由

$$
S(N,E,V)=k\ln\Omega(N,E,V);
$$

3. 反解出 $E=E(S,V,N)$；
4. 利用

$$
\mathrm dE=T\,\mathrm dS-p\,\mathrm dV, \tag{9.3.20}
$$

得到

$$
T=\left(\frac{\partial E}{\partial S}\right)_{V,N},
\qquad
p=-\left(\frac{\partial E}{\partial V}\right)_{S,N}. \tag{9.3.21}
$$

### 9.3.4 补充材料：微正则理论求理想气体

单原子理想气体的哈密顿量为

$$
H=\sum_{i=1}^{3N}\frac{p_i^2}{2m}.
$$

由式 (9.2.9)，先计算能量不超过 $E$ 的累积态数

$$
\Sigma(E)=\frac1{N!h^{3N}}
\int_{H(q,p)\leqslant E}
\mathrm dq_1\cdots\mathrm dq_{3N}\,
\mathrm dp_1\cdots\mathrm dp_{3N}.
$$

坐标积分给出 $V^N$。令 $p_i=\sqrt{2mE}\,x_i$：

$$
\Sigma(E)=K\frac{V^N}{N!h^{3N}}(2mE)^{3N/2}, \tag{9.3.22}
$$

其中

$$
K=\int_{\sum_i x_i^2\leqslant1}
\mathrm dx_1\cdots\mathrm dx_{3N}
$$

是 $3N$ 维单位球体积。高斯积分给出

$$
K=\frac{\pi^{3N/2}}{\Gamma(3N/2+1)}
=\frac{\pi^{3N/2}}{(3N/2)!}. \tag{9.3.23}
$$

因此

$$
\Sigma(E)=
\left(\frac V{h^3}\right)^N
\frac{(2\pi mE)^{3N/2}}{N!(3N/2)!}. \tag{9.3.24}
$$

能量在 $E$ 到 $E+\Delta E$ 间的态数为

$$
\Omega(E)
=\frac{\partial\Sigma}{\partial E}\Delta E
=\frac{3N}{2}\frac{\Delta E}{E}\Sigma(E). \tag{9.3.25}
$$

所以熵为

$$
S=k\ln\Omega
=Nk\ln\left[
\frac V{h^3N}
\left(\frac{4\pi mE}{3N}\right)^{3/2}
\right]+\frac52Nk
+k\left[\ln\left(\frac{3N}{2}\right)+\ln\left(\frac{\Delta E}{E}\right)\right].
$$

热力学极限中最后一项相对前项趋于零。利用斯特林公式，得

$$
S=Nk\ln\left[
\frac V{h^3N}
\left(\frac{4\pi mE}{3N}\right)^{3/2}
\right]+\frac52Nk. \tag{9.3.26}
$$

这是广延量，且与能量壳层宽度 $\Delta E$ 无关。反解：

$$
E(S,V,N)
=\frac{3h^2N^{5/3}}{4\pi mV^{2/3}}
\exp\left(\frac{2S}{3Nk}-\frac53\right). \tag{9.3.27}
$$

由式 (9.3.21)：

$$
T=\left(\frac{\partial E}{\partial S}\right)_{V,N}
=\frac2{3Nk}E,
$$

故

$$
E=\frac32NkT,
$$

$$
p=-\left(\frac{\partial E}{\partial V}\right)_{S,N}
=\frac23\frac EV,
$$

两式联立给出

$$
pV=NkT.
$$

代回式 (9.3.26)，得到熟知的萨克尔-特特罗德型熵。

$K$ 与能量无关的证明也可用配分积分：一方面

$$
\int e^{-\beta E}\frac{\mathrm d\Omega}{N!h^{3N}}
=\frac{V^N}{N!h^{3N}}
\left(\frac{2\pi m}{\beta}\right)^{3N/2}.
$$

另一方面用态数积分

$$
\int e^{-\beta E}\frac{\mathrm d\Omega}{N!h^{3N}}
=\int_0^\infty e^{-\beta E}
\frac{\mathrm d\Sigma}{\mathrm dE}\mathrm dE
=K\frac{V^N}{N!h^{3N}}(2m)^{3N/2}
\left(\frac{3N}{2}\right)!\beta^{-3N/2}.
$$

比较两式即得式 (9.3.23)。

---

## §9.4　正则系综

### 9.4.1 系统与热源

考虑具有确定 $N,V,T$ 的闭系。系统与大热源热接触，只交换能量。复合系统孤立，满足

$$
E+E_r=E^{(0)}, \tag{9.4.1}
$$

且 $E\ll E^{(0)}$。

系统处于能量 $E_s$ 的状态 $s$ 时，热源可处在任意能量为 $E^{(0)}-E_s$ 的微观状态。概率正比于热源微观状态数：

$$
\rho_s\propto\Omega_r(E^{(0)}-E_s). \tag{9.4.2}
$$

将热源熵展开：

$$
\begin{aligned}
\ln\Omega_r(E^{(0)}-E_s)
&=\ln\Omega_r(E^{(0)})
+\left(\frac{\partial\ln\Omega_r}{\partial E_r}\right)_{E_r=E^{(0)}}(-E_s)\\
&=\ln\Omega_r(E^{(0)})-\beta E_s.
\end{aligned}
\tag{9.4.3}
$$

其中

$$
\beta=\left(\frac{\partial\ln\Omega_r}{\partial E_r}\right)_{E_r=E^{(0)}}
=\frac1{kT}.
$$

第一项对系统是常量，故 $\rho_s\propto e^{-\beta E_s}$。归一化：

$$
\rho_s=\frac1Z e^{-\beta E_s}. \tag{9.4.4}
$$

### 9.4.2 配分函数

$$
Z=\sum_s e^{-\beta E_s}. \tag{9.4.5}
$$

若能级 $E_l$ 的简并度为 $\Omega_l$，系统处于该能级的概率和配分函数分别为

$$
\rho_l=\frac1Z\Omega_l e^{-\beta E_l}, \tag{9.4.6}
$$

$$
Z=\sum_l\Omega_l e^{-\beta E_l}. \tag{9.4.7}
$$

经典正则分布与配分函数为

$$
\rho(q,p)\mathrm d\Omega
=\frac1{N!h^{Nr}}
\frac{e^{-\beta E(q,p)}}Z\mathrm d\Omega, \tag{9.4.8}
$$

$$
Z=\frac1{N!h^{Nr}}
\int e^{-\beta E(q,p)}\mathrm d\Omega. \tag{9.4.9}
$$

因子 $1/N!$ 对应全同粒子不可分辨性。

---

## §9.5　正则系综理论的热力学公式

### 9.5.1 内能、广义力、压强与熵

内能为

$$
\begin{aligned}
U=\overline E
&=\frac1Z\sum_sE_se^{-\beta E_s}\\
&=\frac1Z\left(-\frac{\partial}{\partial\beta}\right)
\sum_se^{-\beta E_s}\\
&=-\frac{\partial}{\partial\beta}\ln Z.
\end{aligned}
\tag{9.5.1}
$$

广义力为

$$
\begin{aligned}
Y
&=\frac1Z\sum_s\frac{\partial E_s}{\partial y}e^{-\beta E_s}\\
&=\frac1Z\left(-\frac1\beta\frac{\partial}{\partial y}\right)
\sum_se^{-\beta E_s}\\
&=-\frac1\beta\frac{\partial}{\partial y}\ln Z.
\end{aligned}
\tag{9.5.2}
$$

对体积功 $Y=-p$：

$$
p=\frac1\beta\frac{\partial}{\partial V}\ln Z. \tag{9.5.3}
$$

由式 (9.5.1) 和 (9.5.2)：

$$
\beta(\mathrm dU-Y\,\mathrm dy)
=-\beta\,\mathrm d\left(\frac{\partial\ln Z}{\partial\beta}\right)
+\frac{\partial\ln Z}{\partial y}\mathrm dy.
$$

又

$$
\mathrm d\ln Z
=\frac{\partial\ln Z}{\partial\beta}\mathrm d\beta
+\frac{\partial\ln Z}{\partial y}\mathrm dy,
$$

故

$$
\beta(\mathrm dU-Y\,\mathrm dy)
=\mathrm d\left(
\ln Z-\beta\frac{\partial\ln Z}{\partial\beta}
\right).
$$

与 $T^{-1}(\mathrm dU-Y\,\mathrm dy)=\mathrm dS$ 比较，得到

$$
\beta=\frac1{kT},
$$

$$
S=k\left(\ln Z-\beta\frac{\partial\ln Z}{\partial\beta}\right). \tag{9.5.4}
$$

亥姆霍兹自由能为

$$
F=U-TS=-kT\ln Z. \tag{9.5.5}
$$

所以 $Z(N,V,T)$ 通过特性函数 $F(N,V,T)$ 给出全部平衡性质。

### 9.5.2 能量涨落

定义能量涨落：

$$
\begin{aligned}
\overline{(E-\overline E)^2}
&=\sum_s\rho_s(E_s-\overline E)^2\\
&=\overline{E^2}-(\overline E)^2.
\end{aligned}
$$

由正则分布

$$
\frac{\partial\overline E}{\partial\beta}
=-\left[\overline{E^2}-(\overline E)^2\right],
$$

所以

$$
\overline{(E-\overline E)^2}
=-\frac{\partial\overline E}{\partial\beta}
=kT^2\frac{\partial\overline E}{\partial T}
=kT^2C_V. \tag{9.5.6}
$$

相对涨落为

$$
\frac{\overline{(E-\overline E)^2}}{(\overline E)^2}
=\frac{kT^2C_V}{(\overline E)^2}. \tag{9.5.7}
$$

$C_V$ 与 $\overline E$ 都是广延量，故相对涨落 $\propto1/N$。热力学极限下正则系综与微正则系综给出相同宏观结果。

---

## §9.6　实际气体的物态方程

### 9.6.1 位形积分

单原子经典气体的能量为

$$
E=\sum_{i=1}^N\frac{p_i^2}{2m}
+\sum_{i<j}\phi(r_{ij}), \tag{9.6.1}
$$

其中 $\phi$ 是二体相互作用势。正则配分函数为

$$
Z=\frac1{N!h^{3N}}
\int\cdots\int
e^{-\beta E}\mathrm dq_1\cdots\mathrm dq_{3N}
\mathrm dp_1\cdots\mathrm dp_{3N}. \tag{9.6.2}
$$

动量积分与坐标积分分离：

$$
\prod_{i=1}^{3N}\int_{-\infty}^{\infty}
e^{-\beta p_i^2/(2m)}\mathrm dp_i
=\left(\frac{2\pi m}{\beta}\right)^{3N/2}.
$$

因此

$$
Z=\frac1{N!}\left(\frac{2\pi m}{\beta h^2}\right)^{3N/2}Q, \tag{9.6.3}
$$

其中

$$
Q=\int\cdots\int
e^{-\beta\sum_{i<j}\phi(r_{ij})}
\mathrm dr_1\cdots\mathrm dr_N \tag{9.6.4}
$$

称为**位形积分**。

### 9.6.2 迈耶集团展开

定义迈耶函数

$$
f_{ij}=e^{-\beta\phi(r_{ij})}-1. \tag{9.6.5}
$$

短程势 $\phi$ 只在分子力程内非零，所以 $f_{ij}$ 也只在很小范围内非零。由 $1+f_{ij}=e^{-\beta\phi(r_{ij})}$，

$$
\begin{aligned}
Q
&=\int\cdots\int\prod_{i<j}(1+f_{ij})
\mathrm dr_1\cdots\mathrm dr_N\\
&=\int\cdots\int\left(
1+\sum_{i<j}f_{ij}
+\sum_{i<j}\sum_{i'<j'}f_{ij}f_{i'j'}+\cdots
\right)\mathrm dr_1\cdots\mathrm dr_N.
\end{aligned}
\tag{9.6.6}
$$

第一项给出理想气体 $V^N$。第二项描述一对分子同时处于力程内的二体集团；更高项描述三体及多体集团。只保留二体项，即第二位力系数近似：

$$
Q=\int\cdots\int\left(1+\sum_{i<j}f_{ij}\right)
\mathrm dr_1\cdots\mathrm dr_N. \tag{9.6.7}
$$

共有 $N(N-1)/2$ 个相同积分。对 $i,j$ 之外 $N-2$ 个粒子积分得 $V^{N-2}$，且

$$
\int\cdots\int f_{12}\mathrm dr_1\cdots\mathrm dr_N
=V^{N-2}\iint f_{12}\mathrm dr_1\mathrm dr_2.
$$

忽略边界效应，先对 $\mathrm dr_2$ 积分，再对 $\mathrm dr_1$ 积分，得到

$$
\iint f_{12}\mathrm dr_1\mathrm dr_2
=V\int f_{12}\mathrm dr.
$$

所以

$$
Q=V^N+\frac{N^2}{2}V^{N-1}\int f_{12}\mathrm dr
=V^N\left(1+\frac{N^2}{2V}\int f_{12}\mathrm dr\right).
$$

取对数：

$$
\ln Q=N\ln V
+\ln\left(1+\frac{N^2}{2V}\int f_{12}\mathrm dr\right). \tag{9.6.8}
$$

把对数作级数展开并保留第一位力修正：

$$
\ln Q=N\ln V+\frac{N^2}{2V}\int f_{12}\mathrm dr. \tag{9.6.9}
$$

### 9.6.3 第二位力系数

由式 (9.5.3)，式 (9.6.3) 中只与 $V$ 有关的部分是 $Q$，故

$$
p=\frac1\beta\frac{\partial\ln Z}{\partial V}
=\frac1\beta\frac{\partial\ln Q}{\partial V}.
$$

代入式 (9.6.9)：

$$
p=\frac1\beta\frac NV\left(
1-\frac N{2V}\int f_{12}\mathrm dr
\right).
$$

即

$$
pV=NkT\left(1+\frac{nB}{V}\right), \tag{9.6.10}
$$

其中

$$
B=-\frac{N_A}{2}\int f_{12}\mathrm dr. \tag{9.6.11}
$$

$B$ 称为**第二位力系数**。

### 9.6.4 伦纳德-琼斯势与范德瓦耳斯方程

伦纳德-琼斯势为

$$
\phi(r)=\phi_0\left[
\left(\frac{r_0}{r}\right)^{12}
-2\left(\frac{r_0}{r}\right)^6
\right]. \tag{9.6.12}
$$

在 $r=r_0$ 处势能极小，$\phi(r_0)=-\phi_0$。为解析计算，采用硬核加吸引尾部模型：

$$
\begin{cases}
\phi(r)=+\infty,&r<r_0,\\[1mm]
\phi(r)=-\phi_0\left(\dfrac{r_0}{r}\right)^6,&r\geqslant r_0.
\end{cases} \tag{9.6.13}
$$

式 (9.6.11) 用球坐标写成

$$
B=-2\pi N_A
\int_0^\infty
\left[e^{-\beta\phi(r)}-1\right]r^2\mathrm dr. \tag{9.6.14}
$$

代入式 (9.6.13)：

$$
B=2\pi N_A\left\{
\int_0^{r_0}r^2\mathrm dr
-\int_{r_0}^\infty
\left[e^{-\beta\phi(r)}-1\right]r^2\mathrm dr
\right\}. \tag{9.6.15}
$$

高温 $\phi_0\ll kT$ 时，吸引尾项满足

$$
e^{-\beta\phi(r)}-1\approx-\beta\phi(r).
$$

因此

$$
\begin{aligned}
B
&=2\pi N_A\left[
\frac{r_0^3}{3}
-\int_{r_0}^\infty\beta\phi_0
\left(\frac{r_0}{r}\right)^6r^2\mathrm dr
\right]\\
&=2\pi N_A\left(
\frac{r_0^3}{3}-\frac{\phi_0r_0^3}{3kT}
\right).
\end{aligned}
$$

写成

$$
B=b-\frac a{N_AkT}, \tag{9.6.16}
$$

其中

$$
\begin{cases}
b=\dfrac{2\pi}{3}N_A r_0^3,\\[2mm]
a=\dfrac{2\pi}{3}N_A^2\phi_0r_0^3.
\end{cases} \tag{9.6.17}
$$

代回式 (9.6.10)：

$$
pV=NkT\left(1+\frac{nb}{V}\right)-\frac{n^2a}{V}.
$$

利用 $nb/V\ll1$，近似

$$
1+\frac{nb}{V}\approx\frac1{1-nb/V},
$$

得到

$$
p=\frac{nkT}{V-nb}-\frac{n^2a}{V^2},
$$

即

$$
\left(p+\frac{n^2a}{V^2}\right)(V-nb)=n kT. \tag{9.6.18}
$$

这就是范德瓦耳斯方程。$b$ 对应排除体积，$a$ 对应吸引作用。教材指出，由这个粗略模型得到的 $a,b$ 与温度无关；真实气体在通常范围可近似视作常量。

---

## §9.7　固体的热容

### 9.7.1 简正振动

固体有 $N$ 个原子，每个原子在平衡位置附近有三个方向的小振动，总自由度 $3N$。设第 $i$ 个自由度的位移为 $\xi_i$，动能为

$$
\sum_{i=1}^{3N}\frac{p_{\xi_i}^2}{2m}.
$$

势能在平衡位置展开到二阶：

$$
\Phi=\Phi_0
+\sum_i\left(\frac{\partial\Phi}{\partial\xi_i}\right)_0\xi_i
+\frac12\sum_{i,j}
\left(\frac{\partial^2\Phi}{\partial\xi_i\partial\xi_j}\right)_0
\xi_i\xi_j+\cdots.
$$

平衡位置无受力，一阶项为零。令

$$
a_{ij}=\left(\frac{\partial^2\Phi}{\partial\xi_i\partial\xi_j}\right)_0,
$$

则振动能量为

$$
E=\sum_{i=1}^{3N}\frac{p_{\xi_i}^2}{2m}
+\frac12\sum_{i,j}a_{ij}\xi_i\xi_j+\Phi_0. \tag{9.7.1}
$$

通过正交线性变换可把二次型化为平方和：

$$
E=\frac12\sum_{i=1}^{3N}
\left(p_i^2+\omega_i^2q_i^2\right)+\Phi_0. \tag{9.7.2}
$$

$q_i$ 称为**简正坐标**，$\omega_i$ 称为**简正频率**。原来强耦合的 $N$ 个原子被变换成 $3N$ 个近独立谐振子；展开的高阶项描述简正模式间的相互作用。

> 教材脚注：$3N$ 个简正模式中有 $3$ 个整体平动和 $3$ 个整体转动自由度，它们不感受原子间作用力，故至少有 $6$ 个 $\omega_i=0$。因 $6\ll3N$，以下忽略这一差别。

### 9.7.2 配分函数与内能

量子化后

$$
E=\Phi_0+\sum_{i=1}^{3N}
\hbar\omega_i\left(n_i+\frac12\right),
\qquad n_i=0,1,2,\cdots. \tag{9.7.3}
$$

正则配分函数为

$$
\begin{aligned}
Z
&=\sum_s e^{-\beta E_s}\\
&=e^{-\beta\Phi_0}
\prod_i\left[
\sum_{n_i=0}^{\infty}e^{-\beta\hbar\omega_i(n_i+1/2)}
\right]\\
&=e^{-\beta\Phi_0}
\prod_i\frac{e^{-\beta\hbar\omega_i/2}}
{1-e^{-\beta\hbar\omega_i}}.
\end{aligned}
\tag{9.7.4}
$$

因此

$$
U=-\frac{\partial}{\partial\beta}\ln Z
=U_0+\sum_{i=1}^{3N}
\frac{\hbar\omega_i}{e^{\beta\hbar\omega_i}-1}, \tag{9.7.5}
$$

其中

$$
U_0=\Phi_0+\sum_{i=1}^{3N}\frac{\hbar\omega_i}{2}.
$$

$\Phi_0$ 是负的结合能，其绝对值大于零点能，因此 $U_0<0$；式 (9.7.5) 第二项是温度 $T$ 时的热振动能。

### 9.7.3 德拜频谱

把固体看成连续弹性介质。声波分纵波和横波，色散关系为

$$
\omega=c_lk,
\qquad
\omega=c_tk. \tag{9.7.6}
$$

对给定 $k$，纵波有一个偏振，横波有两个偏振。仿照光子态密度：

$$
D(\omega)\mathrm d\omega
=\frac V{2\pi^2}
\left(\frac1{c_l^3}+\frac2{c_t^3}\right)
\omega^2\mathrm d\omega.
$$

引入常量

$$
B=\frac V{2\pi^2}
\left(\frac1{c_l^3}+\frac2{c_t^3}\right),
$$

则

$$
D(\omega)\mathrm d\omega=B\omega^2\mathrm d\omega. \tag{9.7.7}
$$

固体总简正模式数必须为 $3N$，故设截止频率 $\omega_D$：

$$
\int_0^{\omega_D}B\omega^2\mathrm d\omega=3N.
$$

因此

$$
\omega_D^3=\frac{9N}{B}. \tag{9.7.8}
$$

$\omega_D$ 为**德拜频率**。由式 (9.7.5)：

$$
\begin{aligned}
U
&=U_0+\int_0^{\omega_D}
D(\omega)\frac{\hbar\omega}{e^{\beta\hbar\omega}-1}\mathrm d\omega\\
&=U_0+B\int_0^{\omega_D}
\frac{\hbar\omega^3}{e^{\beta\hbar\omega}-1}\mathrm d\omega.
\end{aligned}
\tag{9.7.9}
$$

令

$$
y=\frac{\hbar\omega}{kT},
\qquad
x=\frac{\hbar\omega_D}{kT}
=\frac{\theta_D}{T}, \tag{9.7.10}
$$

其中 $\theta_D=\hbar\omega_D/k$ 称为**德拜特征温度**。

定义德拜函数

$$
\mathscr D(x)=\frac3{x^3}
\int_0^x\frac{y^3\mathrm dy}{e^y-1}. \tag{9.7.11}
$$

利用式 (9.7.8)，内能可写成

$$
U=U_0+3NkT\mathscr D(x). \tag{9.7.12}
$$

### 9.7.4 高温与低温极限

高温 $T\gg\theta_D$ 即 $x\ll1$。此时 $e^y-1\approx y$：

$$
\mathscr D(x)\approx
\frac3{x^3}\int_0^xy^2\mathrm dy=1.
$$

故

$$
\begin{cases}
U=U_0+3NkT,\\
C_V=3Nk.
\end{cases} \tag{9.7.13}
$$

这就是杜隆-珀蒂定律。

低温 $T\ll\theta_D$ 即 $x\gg1$。积分上限可取 $\infty$：

$$
\mathscr D(x)\approx
\frac3{x^3}\int_0^\infty
\frac{y^3\mathrm dy}{e^y-1}
=\frac3{x^3}\frac{\pi^4}{15}
=\frac{\pi^4}{5x^3}.
$$

因此

$$
\begin{cases}
U=U_0+3Nk\dfrac{\pi^4}{5}\dfrac{T^4}{\theta_D^3},\\[2mm]
C_V=3Nk\dfrac{4\pi^4}{5}
\left(\dfrac T{\theta_D}\right)^3.
\end{cases} \tag{9.7.14}
$$

这就是**德拜 $T^3$ 定律**。金属在约 $3\ \mathrm K$ 以下还要计入自由电子的正比于 $T$ 的热容。

德拜模型在 $\lambda\gg a$ 的低频长波极限准确；当波长与原子间距 $a$ 可比时，连续介质近似失效。不过热容是所有模式贡献之和，对频谱细节并不十分敏感，所以 $T^3$ 定律仍与实验很好符合。

### 9.7.5 声子图像

一个波矢为 $\boldsymbol k$、偏振确定的简正模式，其能量量子可看成准粒子，称为**声子**。声子的动量与能量为

$$
\begin{cases}
\boldsymbol p=\hbar\boldsymbol k,\\
\varepsilon=\hbar\omega.
\end{cases} \tag{9.7.15}
$$

声子可按任意非负整数数目产生，化学势为零，服从玻色分布。温度 $T$ 时平均声子数为

$$
\frac1{e^{\hbar\omega/(kT)}-1}.
$$

因此

$$
U=U_0+\sum_{i=1}^{3N}
\frac{\hbar\omega_i}{e^{\hbar\omega_i/(kT)}-1}. \tag{9.7.16}
$$

这与式 (9.7.5) 完全相同。

---

## *§9.8　液 $^4\mathrm{He}$ 的性质和朗道超流理论　补充材料

### 9.8.1 元激发与二流体模型

相互作用多粒子系统的低激发态可写成元激发能量之和：

$$
E=E_0+\sum_p\varepsilon(p)n(p). \tag{9.7.17}
$$

这里 $E_0$ 是基态能量，$\varepsilon(p)$ 是动量为 $p$ 的元激发能量，$n(p)$ 是占据数。若知道 $\varepsilon(p)$ 和元激发的统计类型，就可把它们当作元激发理想气体处理。固体声子、液 He II 的声子和旋子是玻色型，金属准电子和液 $^3\mathrm{He}$ 准粒子是费米型。

天然氦有 $^3\mathrm{He}$ 和 $^4\mathrm{He}$。$^3\mathrm{He}$ 自旋 $1/2$，是费米子；$^4\mathrm{He}$ 自旋为零，是玻色子。常压下液 $^4\mathrm{He}$ 在 $T_\lambda=2.17\ \mathrm K$ 从 He I 连续转变为 He II，没有潜热和体积突变，热容在两侧趋于对数发散。

He II 可沿 $0.1$ 至 $4\ \mu\mathrm m$ 的细管流动而几乎无黏性；但扭摆实验测得与 He I 同量级的黏度。二流体模型解释这一表面矛盾：He II 是正常成分与超流成分的混合，

$$
\rho_n+\rho_s=\rho,
$$

正常成分有黏度 $\rho_n$，超流成分无黏度 $\rho_s$。$T=0\ \mathrm K$ 时全部为超流成分，$T=T_\lambda$ 时全部为正常成分。二流体是宏观描述，并不表示某些原子永久“属于”某一成分。

### 9.8.2 力热效应与化学平衡

两容器由超漏连通，超漏只能通过超流成分。若温度分别为 $T$ 和 $T+\Delta T$，平衡时出现压强差

$$
\Delta p=\rho s\,\Delta T,
$$

其中 $s$ 是单位质量比熵，这称为**热力效应**。反过来的温度差由压强差产生，称为**力热效应**。喷泉效应是力热效应的一种表现。

化学平衡条件是

$$
\mu(T,p)=\mu(T+\Delta T,p+\Delta p). \tag{9.8.1}
$$

小量展开：

$$
\mu(T+\Delta T,p+\Delta p)
=\mu(T,p)-s\,\Delta T+v\,\Delta p.
$$

其中

$$
s=-\left(\frac{\partial\mu}{\partial T}\right)_p,
\qquad
v=\left(\frac{\partial\mu}{\partial p}\right)_T.
$$

代回式 (9.8.1)：

$$
\Delta p=\frac1v s\,\Delta T=\rho s\,\Delta T. \tag{9.8.2}
$$

He II 的超流成分可视为零熵载体，所以超流流动本身不携热，这是上述温度—压强耦合的微观图像。

### 9.8.3 声子与旋子谱

朗道把接近 $T_\lambda$ 的 He II 看成弱激发玻色系统：基态背景相当于超流成分，元激发气体相当于正常成分。系统的能量和动量可写为

$$
\begin{cases}
E=E_0+\displaystyle\sum_p n(p)\varepsilon(p),\\[2mm]
\boldsymbol P=\displaystyle\sum_p n(p)\boldsymbol p.
\end{cases} \tag{9.8.3}
$$

$T\ll T_\lambda$ 时热容呈 $T^3$，说明低能元激发是声子：

$$
\varepsilon=c_1p, \tag{9.8.4}
$$

其中 $c_1=238\ \mathrm{m\cdot s^{-1}}$ 是第一声速度。

在 $p_0$ 附近，能谱取极小值，展开到二阶：

$$
\varepsilon=\Delta+\frac{(p-p_0)^2}{2m^*}. \tag{9.8.5}
$$

这类元激发称为**旋子**，$\Delta$ 为能隙。实验值约为

$$
\frac{p_0}{\hbar}=19.2\ \mathrm{nm^{-1}},
\qquad
\frac{\Delta}{k}=8.65\ \mathrm K,
\qquad
m^*=0.16m_4.
$$

> 教材脚注：$m_4$ 是 $^4\mathrm{He}$ 原子质量。

由于元激发数不守恒，其化学势为零。能量为 $\varepsilon$ 的状态上的平均占据数为

$$
\bar n=\frac1{e^{\beta\varepsilon}-1}. \tag{9.8.6}
$$

### 9.8.4 声子与旋子对热容的贡献

He II 中声只有纵波，态数为

$$
\frac{4\pi Vp^2\mathrm dp}{h^3}.
$$

声子内能为

$$
\begin{aligned}
U_{\mathrm{ph}}
&=\frac{4\pi V}{h^3}
\int_0^\infty
\frac{c_1p\cdot p^2\mathrm dp}{e^{\beta c_1p}-1}\\
&=\frac{4\pi Vc_1}{h^3}
\left(\frac1{c_1\beta}\right)^4
\int_0^\infty\frac{y^3\mathrm dy}{e^y-1}\\
&=\frac{4\pi^5}{15}VkT
\left(\frac{kT}{hc_1}\right)^3.
\end{aligned}
$$

所以

$$
(C_V)_{\mathrm{ph}}
=\frac{16\pi^5}{15}Vk
\left(\frac{kT}{hc_1}\right)^3. \tag{9.8.7}
$$

旋子部分在 $T\lesssim2\ \mathrm K$ 时 $\Delta\gg kT$，故

$$
\bar n\approx e^{-\beta\varepsilon(p)}. \tag{9.8.8}
$$

可按玻耳兹曼分布处理。旋子单粒子配分函数为

$$
Z_1=\sum_l\omega_le^{-\beta\varepsilon_l}
=\frac{4\pi V}{h^3}e^{-\Delta/(kT)}
\int_0^\infty
e^{-(p-p_0)^2/(2m^*kT)}p^2\mathrm dp. \tag{9.8.9}
$$

积分主要由 $p\approx p_0$ 贡献，将 $p^2$ 换成 $p_0^2$，令

$$
p=p_0+(2m^*kT)^{1/2}y.
$$

因为 $p_0\gg(2m^*kT)^{1/2}$，积分下限可取 $-\infty$：

$$
Z_1=
\frac{4\pi^{3/2}Vp_0^2}{h^3}
(2m^*kT)^{1/2}e^{-\Delta/(kT)}. \tag{9.8.10}
$$

旋子气体非定域，自由能为

$$
F_r=-N_rkT\ln Z_1
+N_rkT(\ln N_r-1). \tag{9.8.11}
$$

由 $\partial F_r/\partial N_r=0$：

$$
\overline N_r=Z_1=
\frac{4\pi Vp_0^2}{h^3}
(2\pi m^*kT)^{1/2}e^{-\Delta/(kT)}. \tag{9.8.12}
$$

代回自由能：

$$
F_r=-kT\overline N_r
=-\frac{4\pi Vp_0^2}{h^3}
(2\pi m^*)^{1/2}(kT)^{3/2}e^{-\Delta/(kT)}. \tag{9.8.13}
$$

因此

$$
\begin{cases}
S_r=-\left(\dfrac{\partial F_r}{\partial T}\right)_V
=k\overline N_r\left(\dfrac32+\dfrac{\Delta}{kT}\right),\\[2mm]
(C_V)_r=T\left(\dfrac{\partial S_r}{\partial T}\right)_V
=\overline N_rk\left[
\dfrac34+\dfrac{\Delta}{kT}
+\left(\dfrac{\Delta}{kT}\right)^2
\right].
\end{cases} \tag{9.8.14}
$$

低温时 $\overline N_r\propto e^{-\Delta/(kT)}\to0$，旋子热容指数衰减。$T<0.6\ \mathrm K$ 时声子占主导，$T>1\ \mathrm K$ 时旋子逐渐占主导。

> 补充材料：1971 年发现 $^3\mathrm{He}$ 在约 $2.6\ \mathrm mK$ 以下超流。零场中有 A、B 两相；两个 $^3\mathrm{He}$ 原子组成自旋 1 的束缚对，表现为玻色子。教材还简述了 $^3\mathrm{He}$–$^4\mathrm{He}$ 稀释制冷：$^3\mathrm{He}$ 从正常相溶入超流相时吸收混合热，可在局部维持约 $2\ \mathrm{mK}$。

---

## *§9.9　伊辛模型的平均场理论　补充材料

### 9.9.1 模型与配分函数

伊辛模型把 $N$ 个自旋 $1/2$ 的磁原子固定在晶格上，只允许自旋沿一个晶轴取

$$
\sigma_i=\pm1.
$$

交换作用近似为近邻对的相互作用：

$$
-J\sum_{i,j}'\sigma_i\sigma_j. \tag{9.9.1}
$$

右上角撇号表示只对近邻对求和。若 $J>0$，平行自旋能量较低，低温可出现自发磁化。有外磁场 $B$ 时，磁矩 $-\mu B\sigma$ 给出塞曼能。位形 $\{\sigma_i\}$ 的总能量为

$$
E\{\sigma_i\}
=-J\sum_{i,j}'\sigma_i\sigma_j
-\mu B\sum_i\sigma_i. \tag{9.9.2}
$$

配分函数为

$$
Z=\sum_{\{\sigma_i\}}e^{-\beta E\{\sigma_i\}}
=\sum_{\sigma_1=\pm1}\cdots\sum_{\sigma_N=\pm1}
e^{\beta J\sum'_{i,j}\sigma_i\sigma_j
+\beta\mu B\sum_i\sigma_i}. \tag{9.9.3}
$$

共有 $2^N$ 个位形。严格求解一般困难：一维可严格解且有限温度无相变，二维昂萨格解存在相变，三维无一般严格解。

### 9.9.2 平均场近似

把能量写成

$$
E=-\mu B\sum_i\sigma_i
-\frac12\sum_{i,j}J_{ij}\sigma_i\sigma_j, \tag{9.9.4}
$$

其中近邻 $J_{ij}=J$，非近邻为零，并对所有 $i,j$ 独立求和，以 $1/2$ 消除重复计数。作用在第 $i$ 个自旋上的力为

$$
-\frac{\partial E}{\partial\sigma_i}
=\mu B+\sum_jJ_{ij}\sigma_j.
$$

因此等效场为

$$
B_i=B+\frac1\mu\sum_jJ_{ij}\sigma_j.
$$

把近邻自旋替换为平均值。设晶格平移不变且每个自旋有 $z$ 个近邻：

$$
\overline B_i
=B+\frac1\mu Jz\overline\sigma. \tag{9.9.5}
$$

平均场忽略局部涨落，使系统化为近独立自旋。配分函数写为

$$
Z=\prod_i\sum_{\sigma_i}e^{\beta\mu\overline B\sigma_i}
=(Z_1)^N, \tag{9.9.6}
$$

其中

$$
Z_1=e^{\beta\mu\overline B}
+e^{-\beta\mu\overline B}. \tag{9.9.7}
$$

### 9.9.3 自洽方程、临界温度与自发磁化

总磁矩为

$$
\begin{aligned}
m&=\frac1\beta\frac{\partial}{\partial B}\ln Z\\
&=\frac N\beta\frac{\partial}{\partial B}\ln Z_1\\
&=N\mu\frac{e^{\beta\mu\overline B}
-e^{-\beta\mu\overline B}}
{e^{\beta\mu\overline B}
+e^{-\beta\mu\overline B}}\\
&=N\mu\tanh\left(\frac{\mu\overline B}{kT}\right)
=N\mu\overline\sigma.
\end{aligned}
\tag{9.9.8}
$$

因此

$$
\overline\sigma
=\tanh\left(\frac{\mu\overline B}{kT}\right). \tag{9.9.9}
$$

无外场 $B=0$ 时，

$$
\overline\sigma
=\tanh(\beta Jz\overline\sigma). \tag{9.9.10}
$$

此超越方程可由 $y=\overline\sigma$ 与 $y=\tanh(\beta Jz\overline\sigma)$ 的交点求解。

- $\beta Jz<1$：只有 $\overline\sigma=0$，为顺磁相；
- $\beta Jz>1$：有 $\overline\sigma=0,\pm\overline\sigma_0$，非零解对应铁磁相。

零解在低温是自由能极大，应舍弃。临界条件为 $\beta Jz=1$：

$$
T_c=\frac{Jz}{k}. \tag{9.9.11}
$$

### 9.9.4 补充材料：平均场临界指数

将式 (9.9.9) 写成

$$
\overline\sigma
=\tanh\left(\frac{\mu B}{kT}
+\overline\sigma\tau\right),
\qquad
\tau=\frac{T_c}{T}. \tag{9.9.12}
$$

利用

$$
\tanh(x+y)=\frac{\tanh x+\tanh y}{1+\tanh x\tanh y},
$$

整理可得

$$
\tanh\left(\frac{\mu B}{kT}\right)
=\frac{\overline\sigma-\tanh(\overline\sigma\tau)}
{1-\overline\sigma\tanh(\overline\sigma\tau)}. \tag{9.9.13}
$$

临界区中 $\mu B/(kT)$ 和 $\overline\sigma\tau$ 都很小，用

$$
\tanh x\approx x-\frac{x^3}{3}.
$$

展开并保留到三阶：

$$
\begin{aligned}
\frac{\mu B}{kT}
&\approx
\left[
\overline\sigma(1-\tau)
+\frac{(\overline\sigma\tau)^3}{3}
\right]
\left[
1+\overline\sigma
\left(\overline\sigma\tau
-\frac{(\overline\sigma\tau)^3}{3}\right)
\right]\\
&\approx\overline\sigma\left[
(1-\tau)+\overline\sigma^2
\left(\tau-\tau^2+\frac{\tau^3}{3}\right)
\right].
\end{aligned}
\tag{9.9.14}
$$

令 $B=0$：

$$
\overline\sigma\left[
(1-\tau)+\overline\sigma^2
\left(\tau-\tau^2+\frac{\tau^3}{3}\right)
\right]=0. \tag{9.9.15}
$$

$T>T_c$ 即 $\tau<1$ 时，只有零解。$T<T_c$ 即 $\tau>1$ 时，非零解满足

$$
\overline\sigma^2
\approx-3(1-\tau)
=3\frac{T_c-T}{T}. \tag{9.9.16}
$$

因此

$$
\overline\sigma\propto(T_c-T)^\beta,
\qquad
\beta=\frac12.
$$

在临界温度 $\tau=1$ 处，式 (9.9.14) 给出

$$
\overline\sigma\approx B^{1/3},
$$

所以临界指数

$$
\delta=3. \tag{9.9.17}
$$

将式 (9.9.14) 对 $B$ 求导：

$$
\frac{\mu}{kT}
=\frac{\partial\overline\sigma}{\partial B}(1-\tau)
+3\overline\sigma^2
\left(\tau-\tau^2+\frac{\tau^3}{3}\right)
\frac{\partial\overline\sigma}{\partial B}. \tag{9.9.18}
$$

定义磁化率

$$
\chi=N\mu\frac{\partial\overline\sigma}{\partial B}.
$$

$T>T_c$ 时 $\overline\sigma=0$，由式 (9.9.18)：

$$
\chi=\frac{N\mu^2}{k}\frac1{T-T_c}. \tag{9.9.19}
$$

$T<T_c$ 时以式 (9.9.16) 代入，得

$$
\chi=\frac12\frac{N\mu^2}{k}\frac1{T_c-T}. \tag{9.9.20}
$$

所以 $\gamma=\gamma'=1$。

零场下配分函数

$$
Z_1=e^{\tau\overline\sigma}+e^{-\tau\overline\sigma}, \tag{9.9.21}
$$

内能为

$$
\begin{aligned}
U
&=-\frac{\partial}{\partial\beta}\ln Z\\
&=-N\frac{\partial}{\partial\beta}\ln Z_1\\
&=NkT^2\frac{\partial}{\partial T}
\ln\left(e^{\tau\overline\sigma}
+e^{-\tau\overline\sigma}\right)
=-NkT_c\overline\sigma^2.
\end{aligned}
\tag{9.9.22}
$$

$T>T_c$ 时 $\overline\sigma=0$，$U=0$；$T<T_c$ 时由式 (9.9.16)，

$$
U=-NkT_c\cdot3\frac{T_c-T}{T}
=3Nk(T-T_c).
$$

故零场热容在 $T_c$ 由零跳至 $3Nk$，临界指数 $\alpha=\alpha'=0$。

总括：

$$
\alpha=0,\qquad
\beta=\frac12,\qquad
\gamma=1,\qquad
\delta=3. \tag{9.9.23}
$$

这与朗道连续相变理论相同。平均场理论忽略涨落，所以对低维系统给出过高 $T_c$；维数越低，涨落破坏有序越明显。

---

## §9.10　巨正则系综

### 9.10.1 系统、热源和粒子源

系统有确定 $V,T,\mu$，与可交换能量和粒子的大源接触。复合系统孤立：

$$
\begin{cases}
E+E_r=E^{(0)},\\
N+N_r=N^{(0)}.
\end{cases} \tag{9.10.1}
$$

且 $E\ll E^{(0)},N\ll N^{(0)}$。

系统处于粒子数 $N$、微观状态 $s$、能量 $E_s$ 时，源有

$$
N^{(0)}-N,\qquad E^{(0)}-E_s.
$$

概率正比于源的微观状态数：

$$
\rho_{N,s}\propto
\Omega_r(N^{(0)}-N,E^{(0)}-E_s). \tag{9.10.2}
$$

将源的 $\ln\Omega_r$ 对 $N_r,E_r$ 展开，只留一阶：

$$
\begin{aligned}
\ln\Omega_r(N^{(0)}-N,E^{(0)}-E_s)
&=\ln\Omega_r(N^{(0)},E^{(0)})\\
&\quad-\left(\frac{\partial\ln\Omega_r}{\partial N_r}\right)N
-\left(\frac{\partial\ln\Omega_r}{\partial E_r}\right)E_s\\
&=\ln\Omega_r(N^{(0)},E^{(0)})
-\alpha N-\beta E_s.
\end{aligned}
\tag{9.10.3}
$$

源参数为

$$
\alpha=\left(\frac{\partial\ln\Omega_r}{\partial N_r}\right)
=-\frac{\mu}{kT},
\qquad
\beta=\left(\frac{\partial\ln\Omega_r}{\partial E_r}\right)
=\frac1{kT}.
$$

因此

$$
\rho_{N,s}\propto e^{-\alpha N-\beta E_s}. \tag{9.10.4}
$$

归一化为

$$
\rho_{N,s}
=\frac1\Xi e^{-\alpha N-\beta E_s}, \tag{9.10.5}
$$

其中巨配分函数

$$
\Xi=\sum_{N=0}^\infty\sum_s
e^{-\alpha N-\beta E_s}. \tag{9.10.6}
$$

量子巨正则分布由式 (9.10.5) 给出。经典表达式为

$$
\rho_N(q,p)\mathrm dq\,\mathrm dp
=\frac1{N!h^{Nr}}
\frac{e^{-\alpha N-\beta E(q,p)}}{\Xi}
\mathrm d\Omega, \tag{9.10.7}
$$

$$
\Xi=\sum_N\frac{e^{-\alpha N}}{N!h^{Nr}}
\int e^{-\beta E(q,p)}\mathrm d\Omega. \tag{9.10.8}
$$

---

## §9.11　巨正则系综理论的热力学公式

### 9.11.1 平均量与热力学函数

平均粒子数：

$$
\begin{aligned}
\overline N
&=\frac1\Xi\sum_N\sum_sN
e^{-\alpha N-\beta E_s}\\
&=\frac1\Xi\left(-\frac{\partial}{\partial\alpha}\right)\Xi
=-\frac{\partial\ln\Xi}{\partial\alpha}.
\end{aligned}
\tag{9.11.1}
$$

内能：

$$
\begin{aligned}
U=\overline E
&=\frac1\Xi\sum_N\sum_sE_s
e^{-\alpha N-\beta E_s}\\
&=-\frac{\partial\ln\Xi}{\partial\beta}.
\end{aligned}
\tag{9.11.2}
$$

广义力：

$$
\begin{aligned}
Y
&=\frac1\Xi\sum_N\sum_s
\frac{\partial E_s}{\partial y}
e^{-\alpha N-\beta E_s}\\
&=-\frac1\beta
\frac{\partial\ln\Xi}{\partial y}.
\end{aligned}
\tag{9.11.3}
$$

对体积功：

$$
p=\frac1\beta\frac{\partial\ln\Xi}{\partial V}. \tag{9.11.4}
$$

### 9.11.2 熵与巨热力势

由式 (9.11.1)–(9.11.3)：

$$
\begin{aligned}
\beta\left(\mathrm dU-Y\,\mathrm dy
+\frac{\alpha}{\beta}\mathrm d\overline N\right)
={}&-\beta\,\mathrm d\left(
\frac{\partial\ln\Xi}{\partial\beta}
\right)
+\frac{\partial\ln\Xi}{\partial y}\mathrm dy\\
&-\alpha\,\mathrm d\left(
\frac{\partial\ln\Xi}{\partial\alpha}
\right).
\end{aligned}
\tag{9.11.5}
$$

因为 $\ln\Xi=\ln\Xi(\alpha,\beta,y)$，

$$
\mathrm d\ln\Xi
=\frac{\partial\ln\Xi}{\partial\alpha}\mathrm d\alpha
+\frac{\partial\ln\Xi}{\partial\beta}\mathrm d\beta
+\frac{\partial\ln\Xi}{\partial y}\mathrm dy.
$$

故

$$
\beta\left(\mathrm dU-Y\,\mathrm dy
+\frac{\alpha}{\beta}\mathrm d\overline N\right)
=\mathrm d\left(
\ln\Xi-\alpha\frac{\partial\ln\Xi}{\partial\alpha}
-\beta\frac{\partial\ln\Xi}{\partial\beta}
\right). \tag{9.11.6}
$$

与开系关系

$$
\frac1T(\mathrm dU-Y\,\mathrm dy-\mu\,\mathrm dN)=\mathrm dS
$$

比较：

$$
\beta=\frac1{kT},
\qquad
\alpha=-\frac{\mu}{kT},
$$

$$
S=k\left(
\ln\Xi-\alpha\frac{\partial\ln\Xi}{\partial\alpha}
-\beta\frac{\partial\ln\Xi}{\partial\beta}
\right). \tag{9.11.7}
$$

以 $T,V,\mu$ 为自变量的特性函数为巨热力势

$$
J=U-TS-\mu\overline N=-kT\ln\Xi.
$$

### 9.11.3 粒子数涨落

$$
\overline{(N-\overline N)^2}
=\overline{N^2}-(\overline N)^2. \tag{9.11.8}
$$

由巨正则分布：

$$
\begin{aligned}
\frac{\partial\overline N}{\partial\alpha}
&=-\left[\overline{N^2}-(\overline N)^2\right],
\end{aligned}
$$

故

$$
\overline{(N-\overline N)^2}
=-\left(\frac{\partial\overline N}{\partial\alpha}\right)_{\beta,y}
=kT\left(\frac{\partial\overline N}{\partial\mu}\right)_{T,V}. \tag{9.11.9}
$$

相对涨落：

$$
\frac{\overline{(N-\overline N)^2}}{(\overline N)^2}
=\frac{kT}{(\overline N)^2}
\left(\frac{\partial\overline N}{\partial\mu}\right)_{T,V}. \tag{9.11.10}
$$

利用热力学恒等式，可把它写成

$$
\frac{\overline{(N-\overline N)^2}}{(\overline N)^2}
=\frac{kT}{V}\kappa_T, \tag{9.11.11}
$$

其中

$$
\kappa_T=-\frac1V\left(\frac{\partial V}{\partial p}\right)_T
$$

是等温压缩系数。正常热力学极限下 $V\propto\overline N$，故相对涨落 $\propto1/\overline N$，巨正则与正则系综等价。在气液临界点和两相共存区，$\kappa_T$ 很大甚至发散，粒子数涨落可宏观显著，这时巨正则描述尤其方便。

---

## *§9.12　巨正则系综理论的简单应用　补充材料

### 9.12.1 吸附现象

设有 $N_0$ 个独立吸附中心，每个中心最多吸附一个分子，被吸附分子的能量为 $-\varepsilon_0$。若已有 $N$ 个分子，能量为 $-N\varepsilon_0$，排列数为

$$
\frac{N_0!}{N!(N_0-N)!}.
$$

巨配分函数为

$$
\begin{aligned}
\Xi
&=\sum_{N=0}^{N_0}
e^{\beta(\mu+\varepsilon_0)N}
\frac{N_0!}{N!(N_0-N)!}\\
&=\left[1+e^{\beta(\mu+\varepsilon_0)}\right]^{N_0}.
\end{aligned}
\tag{9.12.1}
$$

平均吸附数为

$$
\begin{aligned}
\overline N
&=-\frac{\partial}{\partial\alpha}\ln\Xi
=kT\frac{\partial}{\partial\mu}\ln\Xi\\
&=\frac{N_0}{1+e^{-\beta(\varepsilon_0+\mu)}}.
\end{aligned}
\tag{9.12.2}
$$

被吸附分子与气体达到平衡，化学势相同。由理想气体化学势：

$$
e^{\mu/(kT)}
=\frac p{kT}
\left(\frac{h^2}{2\pi mkT}\right)^{3/2}.
$$

因此覆盖度

$$
\theta=\frac{\overline N}{N_0}
=\frac1{1+\dfrac{kT}{p}
\left(\dfrac{2\pi mkT}{h^2}\right)^{3/2}
e^{-\varepsilon_0/(kT)}}. \tag{9.12.3}
$$

### 9.12.2 严格重导玻色与费米分布

巨正则分布：

$$
\rho_{N,s}
=\frac1\Xi e^{-\alpha N-\beta E_s}. \tag{9.12.4}
$$

$$
\Xi=\sum_N\sum_s e^{-\alpha N-\beta E_s}. \tag{9.12.5}
$$

近独立粒子的分布为 $\{a_l\}$，满足

$$
N=\sum_la_l,
\qquad
E=\sum_l\varepsilon_la_l.
$$

巨正则对所有 $N,s$ 求和等价于对所有占据数组 $\{a_l\}$ 求和，故

$$
\begin{aligned}
\Xi
&=\sum_{\{a_l\}}
\exp\left[
-\sum_l(\alpha+\beta\varepsilon_l)a_l
\right]\\
&=\prod_l\sum_{a_l}
e^{-(\alpha+\beta\varepsilon_l)a_l}\\
&=\prod_l\Xi_l,
\end{aligned}
\tag{9.12.6}
$$

$$
\Xi_l=\sum_{a_l}
e^{-(\alpha+\beta\varepsilon_l)a_l}. \tag{9.12.7}
$$

第 $l$ 个能级的平均粒子数：

$$
\begin{aligned}
\overline{a_l}
&=\frac1\Xi\sum_N\sum_s a_l
e^{-\alpha N-\beta E_s}\\
&=\frac1{\Xi_l}
\sum_{a_l}a_l
e^{-(\alpha+\beta\varepsilon_l)a_l}\\
&=-\frac{\partial}{\partial\alpha}\ln\Xi_l.
\end{aligned}
\tag{9.12.8}
$$

玻色子 $a_l=0,1,2,\cdots$：

$$
\begin{aligned}
\Xi_l
&=\sum_{a_l=0}^{\infty}
e^{-(\alpha+\beta\varepsilon_l)a_l}
=\frac1{1-e^{-\alpha-\beta\varepsilon_l}},\\
\ln\Xi_l
&=-\ln\left(1-e^{-\alpha-\beta\varepsilon_l}\right),\\
\overline{a_l}
&=\frac1{e^{\alpha+\beta\varepsilon_l}-1}.
\end{aligned}
\tag{9.12.9}
$$

费米子 $a_l=0,1$：

$$
\begin{aligned}
\Xi_l
&=1+e^{-\alpha-\beta\varepsilon_l},\\
\ln\Xi_l
&=\ln\left(1+e^{-\alpha-\beta\varepsilon_l}\right),\\
\overline{a_l}
&=\frac1{e^{\alpha+\beta\varepsilon_l}+1}.
\end{aligned}
\tag{9.12.10}
$$

这两式不依赖第六章推导中 $a_l\gg1,\omega_l\gg1$ 的近似，因而逻辑上更严格。

### 9.12.3 玻色与费米分布的涨落

把一个能级上的粒子视为开系，由式 (9.11.9)：

$$
\overline{(a_l-\overline{a_l})^2}
=-\frac{\partial\overline{a_l}}{\partial\alpha}. \tag{9.12.11}
$$

代入式 (9.12.9)、(9.12.10)：

$$
\overline{(a_l-\overline{a_l})^2}
=\overline{a_l}
\left(1\pm\frac{\overline{a_l}}{\omega_l}\right), \tag{9.12.12}
$$

其中“$+$”对应玻色，“$-$”对应费米；式中的 $\omega_l$ 是能级简并度，等价于相加各单态涨落。

费米强简并时，$\overline{a_l}/\omega_l$ 只能为 $0$ 或 $1$，所以

$$
\overline{(a_l-\overline{a_l})^2}
=\overline{a_l}\left(1-\frac{\overline{a_l}}{\omega_l}\right),
$$

占据数涨落小。玻色气体没有此上限，尤其 $\mu\to0$ 时涨落可很大。

不同能级的涨落互不相关：

$$
\begin{aligned}
\overline{a_l a_m}
&=\frac1\Xi\sum_N\sum_s a_la_m
e^{-\alpha N-\beta E_s}\\
&=\frac1{\Xi_l}\sum_{a_l}a_l e^{-(\alpha+\beta\varepsilon_l)a_l}\\
&\quad\cdot\frac1{\Xi_m}\sum_{a_m}a_m
e^{-(\alpha+\beta\varepsilon_m)a_m}\\
&=\overline{a_l}\,\overline{a_m}.
\end{aligned}
\tag{9.12.13}
$$

因此

$$
\overline{(a_l-\overline{a_l})(a_m-\overline{a_m})}=0. \tag{9.12.14}
$$

### 9.12.4 补充材料：带简并度的严格推导

先列两个组合恒等式。对整数 $m>0$，

$$
(1+x)^m
=\sum_{n=0}^m\frac{m!}{n!(m-n)!}x^n. \tag{9.12.15}
$$

负二项式系数满足

$$
\binom{-m}{n}
=(-1)^n\binom{m+n-1}{n},
$$

故

$$
(1-x)^{-m}
=\sum_{n=0}^\infty
\binom{m+n-1}{n}x^n. \tag{9.12.16}
$$

对能级简并度 $\omega_l$，给定分布 $\{a_l\}$ 的微观状态数为

$$
\Omega_l=
\begin{cases}
\dfrac{(\omega_l+a_l-1)!}
{a_l!(\omega_l-1)!},&\text{玻色},\\[3mm]
\dfrac{\omega_l!}{a_l!(\omega_l-a_l)!},&\text{费米}.
\end{cases}
$$

因此

$$
\begin{aligned}
\Xi
&=\sum_{\{a_l\}}
\Omega e^{-\sum_l(\alpha+\beta\varepsilon_l)a_l}\\
&=\prod_l\sum_{a_l}
\Omega_l
e^{-(\alpha+\beta\varepsilon_l)a_l}\\
&=\prod_l\Xi_l,
\end{aligned}
\tag{9.12.17}
$$

$$
\Xi_l=
\sum_{a_l}\Omega_l
e^{-(\alpha+\beta\varepsilon_l)a_l}. \tag{9.12.18}
$$

玻色情形：

$$
\begin{align}
\Xi_l
&=\sum_{a_l=0}^\infty
\binom{\omega_l+a_l-1}{a_l}
e^{-(\alpha+\beta\varepsilon_l)a_l}\\
&=\left[1-e^{-\alpha-\beta\varepsilon_l}
\right]^{-\omega_l},\\
\ln\Xi_l
&=-\omega_l\ln\left(
1-e^{-\alpha-\beta\varepsilon_l}
\right), \tag{9.12.19}\\
\overline{a_l}
&=-\frac{\partial}{\partial\alpha}\ln\Xi_l
=\frac{\omega_l}{e^{\alpha+\beta\varepsilon_l}-1}. \tag{9.12.20}
\end{align}
$$

费米情形：

$$
\begin{align}
\Xi_l
&=\sum_{a_l=0}^{\omega_l}
\binom{\omega_l}{a_l}
e^{-(\alpha+\beta\varepsilon_l)a_l}\\
&=\left[1+e^{-\alpha-\beta\varepsilon_l}
\right]^{\omega_l},\\
\ln\Xi_l
&=\omega_l\ln\left(
1+e^{-\alpha-\beta\varepsilon_l}
\right), \tag{9.12.21}\\
\overline{a_l}
&=-\frac{\partial}{\partial\alpha}\ln\Xi_l
=\frac{\omega_l}{e^{\alpha+\beta\varepsilon_l}+1}. \tag{9.12.22}
\end{align}
$$

式 (9.12.19)、(9.12.21) 是第八章所用巨配分函数的严格推导。

---

## 本章小结

### 一、公式速查

| 内容 | 公式 |
| --- | --- |
| 哈密顿方程与刘维尔定理 | (9.1.1)；$\mathrm d\rho/\mathrm dt=0$ (9.1.6)；连续性方程 (9.1.9) |
| 微正则分布 | 经典 (9.2.7)；量子 $\rho_s=1/\Omega$ (9.2.8)；$\Omega$ 的相空间表示 (9.2.9)(9.2.10) |
| 微正则热力学 | $\beta=\partial_E\ln\Omega$ (9.3.5)；$S=k\ln\Omega$ (9.3.10)；$\gamma=p/(kT),\alpha=-\mu/(kT)$ (9.3.17) |
| 正则分布与配分函数 | (9.4.4)(9.4.5)(9.4.9) |
| 正则热力学 | $U=-\partial_\beta\ln Z$ (9.5.1)；$Y=-\beta^{-1}\partial_y\ln Z$ (9.5.2)；$S=k(\ln Z-\beta\partial_\beta\ln Z)$ (9.5.4)；$F=-kT\ln Z$ (9.5.5) |
| 能量涨落 | $\overline{(\Delta E)^2}=kT^2C_V$ (9.5.6) |
| 实际气体 | 位形积分 (9.6.4)；迈耶函数 (9.6.5)；$pV=NkT(1+nB/V)$ (9.6.10)；$B=b-a/(N_AkT)$ (9.6.16)；范德瓦耳斯方程 (9.6.18) |
| 德拜固体 | $Z$ 与 $U$ (9.7.4)(9.7.5)；$\omega_D^3=9N/B$ (9.7.8)；$U=U_0+3NkT\mathscr D(x)$ (9.7.12) |
| 德拜极限 | 高温 (9.7.13)；低温 $C_V\propto T^3$ (9.7.14) |
| 液 He II | $\Delta p=\rho s\Delta T$ (9.8.2)；声子谱 (9.8.4)；旋子谱 (9.8.5)；声子、旋子热容 (9.8.7)(9.8.14) |
| 伊辛平均场 | 自洽方程 (9.9.10)；$T_c=Jz/k$ (9.9.11)；平均场临界指数 (9.9.23) |
| 巨正则分布 | (9.10.5)；$\Xi$ (9.10.6) |
| 巨正则热力学 | $\overline N=-\partial_\alpha\ln\Xi$ (9.11.1)；$U=-\partial_\beta\ln\Xi$ (9.11.2)；$p=\beta^{-1}\partial_V\ln\Xi$ (9.11.4)；$S$ (9.11.7) |
| 粒子数涨落 | $\overline{(\Delta N)^2}=kT(\partial_\mu\overline N)_{T,V}$ (9.11.9)；相对涨落 (9.11.11) |
| 严格量子分布 | 玻色 (9.12.20)；费米 (9.12.22)；涨落 (9.12.12) |

### 二、概念脉络

```text
Γ 相空间 + Hamilton 方程
        │
        ▼
Liouville：相体积与代表点密度守恒
        │
        ├──► 孤立系 N,E,V
        │       └──► 微正则：等概率 Ω(N,E,V)
        │               S = k lnΩ
        │               β=∂lnΩ/∂E
        │
        ├──► 恒温闭系 N,V,T
        │       └──► 正则：ρs ∝ exp(-βEs)
        │               Z = Σ exp(-βEs)
        │               F = -kT lnZ
        │               ΔE² = kT² C_V
        │
        └──► 恒温开系 V,T,μ
                └──► 巨正则：ρ ∝ exp(-αN-βEs)
                        Ξ = Σ exp(-αN-βEs)
                        J = -kT lnΞ
                        ΔN² = kT ∂N/∂μ

相互作用系统
    ├──► 实际气体：位形积分 → Mayer 集团 → B(T)
    │                    └──► van der Waals 方程
    ├──► 固体：原子振动 → 简正模式 → 声子 → Debye T³
    ├──► 液 He II：声子 + 旋子 → 超流二流体
    └──► Ising：平均场 → 自洽磁化 → 临界指数
```

三种系综的对应关系：

| 系综 | 固定量 | 分布 | 特性函数 |
| --- | --- | --- | --- |
| 微正则 | $N,E,V$ | $\rho_s=1/\Omega$ | $S=k\ln\Omega$ |
| 正则 | $N,T,V$ | $\rho_s\propto e^{-\beta E_s}$ | $F=-kT\ln Z$ |
| 巨正则 | $\mu,T,V$ | $\rho_{N,s}\propto e^{-\alpha N-\beta E_s}$ | $J=-kT\ln\Xi$ |

### 三、易错点

1. **$\Gamma$ 空间不是单粒子 $\mu$ 空间**：$\Gamma$ 空间维度为 $2Nr$，描述整个系统；$\mu$ 空间只描述一个粒子的 $2r$ 个变量。
2. **刘维尔定理不含统计假设**：它只说明哈密顿流不可压缩；微正则分布仍依赖等概率原理。
3. **相空间态数要除 $h^{Nr}N!$**：$h$ 负责相格计数，$N!$ 负责全同粒子不可分辨。只除 $h^{Nr}$ 会把交换全同粒子误算成新状态。
4. **微正则能量壳层不能写成严格曲面上的连续均匀测度**：通常取有限厚度 $\Delta E$，再令其宏观小、微观大；物理量与 $\Delta E$ 无关。
5. **$S=k\ln\Omega$ 对相互作用系统仍然成立**：它不是近独立粒子的专门结论，只要 $\Omega$ 包含完整相互作用哈密顿量。
6. **正则配分函数 $Z$ 与单粒子 $Z_1$ 不同**：$Z$ 是整个 $N$ 粒子系统的配分函数；近独立系统才有 $Z=Z_1^N/N!$（经典极限）。
7. **$\beta$ 由积分因子确定，不能只凭量纲**：式 (9.5.4) 前的微分比较同时给出 $\beta=1/(kT)$。
8. **能量涨落与热容相连**：$\overline{(\Delta E)^2}=kT^2C_V$；热容非负也是正则分布的稳定性条件。
9. **正则与微正则等价是热力学极限意义下的等价**：宏观量相同，但涨落信息不同；正则允许能量涨落，微正则能量确定。
10. **第二位力系数不只来自硬核排斥**：排斥给出正贡献，吸引给出负贡献；$B$ 随温度可改变符号。
11. **Mayer 集团展开不是简单地把高阶项丢弃**：截断到二体项就是稀薄气体近似；密度增大后三体、多体集团不可忽略。
12. **德拜 $T^3$ 不是所有温度都成立**：它是 $T\ll\theta_D$ 的低频极限；高温恢复 $3Nk$。
13. **声子与真实原子不同**：声子是简正模式的量子，数目不守恒、化学势为零；不能把声子当作可在真空中自由传播的普通物质粒子。
14. **液 He II 的二流体不是两种原子**：它是同一量子液体的两种宏观响应成分；超流成分与正常成分可按温度和流动相互转化。
15. **超流不是零惯性**：超流成分可输运质量，但无黏性、无熵流；正常成分才携带熵和耗散。
16. **平均场低估临界区涨落**：临界指数在某些模型由平均场给出，但真实系统在临界点附近涨落长程关联，临界指数可能不同。
17. **平均场 $T_c$ 可能高估真实临界温度**：涨落破坏有序，尤其低维体系；平均场忽略这种破坏。
18. **巨正则分布中 $N$ 不固定**：$\overline N$ 只是一个平均值，不能把 $N$ 当作每个微观状态都相同的常数。
19. **巨正则粒子数相对涨落在临界点可异常增大**：$\kappa_T\to\infty$ 时，固定体积内的粒子数可明显涨落；这不是统计错误，而是临界涨落的物理特征。

### 四、名词中英对照

| 中文 | 英文 | 中文 | 英文 |
| --- | --- | --- | --- |
| 相空间 | phase space | $\Gamma$ 空间 | gamma space |
| 代表点 | representative point | 刘维尔定理 | Liouville theorem |
| 系综 | ensemble | 微正则系综 | microcanonical ensemble |
| 正则系综 | canonical ensemble | 巨正则系综 | grand canonical ensemble |
| 配分函数 | partition function | 巨配分函数 | grand partition function |
| 微观状态数 | number of microstates | 系综平均 | ensemble average |
| 等概率原理 | principle of equal a priori probabilities | 涨落 | fluctuation |
| 热源 | heat reservoir | 粒子源 | particle reservoir |
| 位形积分 | configuration integral | 迈耶集团展开 | Mayer cluster expansion |
| 第二位力系数 | second virial coefficient | 范德瓦耳斯方程 | van der Waals equation |
| 德拜模型 | Debye model | 德拜频率 | Debye frequency |
| 德拜温度 | Debye temperature | 声子 | phonon |
| 旋子 | roton | 元激发 | elementary excitation |
| 超流 | superfluidity | 二流体模型 | two-fluid model |
| 第二声 | second sound | 稀释制冷 | dilution refrigeration |
| 伊辛模型 | Ising model | 平均场近似 | mean-field approximation |
| 自发磁化 | spontaneous magnetization | 临界指数 | critical exponent |
| 巨热力势 | grand thermodynamic potential | 吸附 | adsorption |
