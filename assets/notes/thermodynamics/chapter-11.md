---
tags:
  - 热统
  - 统计物理
  - 第十一章
created: 2026-10-22
---

# 第十一章　非平衡态统计理论初步

> **教材**：汪志诚《热力学·统计物理》（第六版）
> **范围**：§11.1–§11.6（教材 pp. 273–289；本章无带星号正文节）

## 本章地图

| 节 | 标题 | 一句话 |
| --- | --- | --- |
| §11.1 | 玻耳兹曼方程的弛豫时间近似 | 用分布函数 $f(\boldsymbol r,\boldsymbol v,t)$ 描述非平衡气体，以偏离局域平衡的线性碰撞项闭合演化方程 |
| §11.2 | 气体的黏性现象 | 速度梯度使定向动量输运，弛豫时间近似给出黏度 $\eta=nkT\overline{\tau}_0$ |
| §11.3 | 金属的电导率 | 电场使费米面附近电子分布偏移，得 $\sigma=ne^2\tau_F/m$ 及高温 $\sigma\propto1/T$ |
| §11.4 | 玻耳兹曼积分微分方程 | 由二体弹性碰撞的元碰撞与反碰撞导出碰撞积分，得到玻耳兹曼积分微分方程 |
| §11.5 | $H$ 定理 | $H=\iint f\ln f\,\mathrm d\tau\mathrm d\omega$ 单调减少，直到细致平衡；统计层次解释不可逆性 |
| §11.6 | 细致平衡原理与平衡态的分布函数 | 细致平衡和运动项同时为零，确定平衡分布为局部麦克斯韦分布及相应的空间密度分布 |

**主线**：分布函数 → 相空间通量 → 弛豫时间近似 → 黏性与电导 → 严格二体碰撞积分 → $H$ 定理 → 细致平衡 → 平衡分布。

---

## §11.1　玻耳兹曼方程的弛豫时间近似

### 11.1.1 非平衡分布函数

稀薄气体满足

$$
\frac{h}{\sqrt{2\pi mkT}}
\left(\frac NV\right)^{1/3}\ll1,
$$

即分子平均热波长远小于平均间距，可把分子看成经典粒子。忽略内部结构时，用位置 $\boldsymbol r$ 和速度 $\boldsymbol v$ 描述单分子状态。

定义

$$
f(\boldsymbol r,\boldsymbol v,t)\mathrm d\tau\mathrm d\omega \tag{11.1.1}
$$

为时刻 $t$ 位于体积元 $\mathrm d\tau=\mathrm dx\mathrm dy\mathrm dz$、速度间隔 $\mathrm d\omega=\mathrm dv_x\mathrm dv_y\mathrm dv_z$ 内的平均分子数。$\mathrm d\tau\mathrm d\omega$ 微观上足够大以含大量分子，宏观上足够小可视为一点。

### 11.1.2 分布函数的时间变化率

经过 $\mathrm dt$ 后，同一体积元和速度间隔中的分子数为

$$
f(\boldsymbol r,\boldsymbol v,t+\mathrm dt)\mathrm d\tau\mathrm d\omega.
$$

一阶展开并与初值相减，得到分子数增加量

$$
\frac{\partial f}{\partial t}
\mathrm dt\,\mathrm d\tau\mathrm d\omega. \tag{11.1.2}
$$

变化来自两部分：

1. 分子以速度运动，在外力下速度改变，使分子穿过相空间边界；
2. 分子碰撞使速度突然改变。

### 11.1.3 运动引起的相空间通量

先看 $x$ 方向。通过平面 $x$ 进入相格的速度通量正比于 $f\dot x=v_xf$，通过平面 $x+\mathrm dx$ 流出的通量作一阶展开。两者相减，给出 $x$ 方向净流入

$$
-\frac{\partial}{\partial x}(fv_x)\mathrm dx\,\mathrm dt\mathrm dA
=-\frac{\partial}{\partial x}(fv_x)\mathrm dt\,\mathrm d\tau\mathrm d\omega.
$$

类似考察速度面的通量。六个面的贡献合起来为

$$
-\left[
\frac{\partial}{\partial x}(fv_x)
+\frac{\partial}{\partial y}(fv_y)
+\frac{\partial}{\partial z}(fv_z)
+\frac{\partial}{\partial v_x}(f\dot v_x)
+\frac{\partial}{\partial v_y}(f\dot v_y)
+\frac{\partial}{\partial v_z}(f\dot v_z)
\right]\mathrm dt\,\mathrm d\tau\mathrm d\omega. \tag{11.1.3}
$$

坐标和速度是独立变量，故

$$
\frac{\partial v_x}{\partial x}
=\frac{\partial v_y}{\partial y}
=\frac{\partial v_z}{\partial z}=0.
$$

设单位质量外力分量为 $X,Y,Z$，则

$$
\dot v_x=X,
\qquad
\dot v_y=Y,
\qquad
\dot v_z=Z. \tag{11.1.4}
$$

对质量 $m$、电荷 $q$ 的粒子，

$$
m\boldsymbol F
=q\left(\boldsymbol E+\boldsymbol v\times\boldsymbol B\right). \tag{11.1.5}
$$

重力与速度无关；洛伦兹力的 $x$ 分量不依赖 $v_x$。因此采用

$$
\frac{\partial X}{\partial v_x}
+\frac{\partial Y}{\partial v_y}
+\frac{\partial Z}{\partial v_z}=0. \tag{11.1.6}
$$

通量式 (11.1.3) 化为

$$
-\left(
v_x\frac{\partial f}{\partial x}
+v_y\frac{\partial f}{\partial y}
+v_z\frac{\partial f}{\partial z}
+X\frac{\partial f}{\partial v_x}
+Y\frac{\partial f}{\partial v_y}
+Z\frac{\partial f}{\partial v_z}
\right)\mathrm dt\,\mathrm d\tau\mathrm d\omega, \tag{11.1.7}
$$

故运动导致的分布函数变化率为

$$
-\left(
v_x\frac{\partial f}{\partial x}
+v_y\frac{\partial f}{\partial y}
+v_z\frac{\partial f}{\partial z}
+X\frac{\partial f}{\partial v_x}
+Y\frac{\partial f}{\partial v_y}
+Z\frac{\partial f}{\partial v_z}
\right). \tag{11.1.8}
$$

### 11.1.4 弛豫时间近似

分子碰撞很频繁，先在宏观小区域建立局域平衡，再通过扩散、热传导等缓慢过程达到整体平衡。局域平衡分布仍为麦克斯韦分布：

$$
f^{(0)}
=n\left(\frac{m}{2\pi kT}\right)^{3/2}
\exp\left[
-\frac{m}{2kT}
(\boldsymbol v-\boldsymbol v_0)^2
\right], \tag{11.1.9}
$$

只是 $n,T,\boldsymbol v_0$ 可以是 $\boldsymbol r,t$ 的缓变函数。

设碰撞使 $f-f^{(0)}$ 按比例衰减：

$$
\left[
\frac{\partial}{\partial t}(f-f^{(0)})
\right]_{\mathrm c}
=-\frac{f-f^{(0)}}{\tau_0}. \tag{11.1.10}
$$

积分得

$$
f(t)-f^{(0)}
=\left[f(0)-f^{(0)}\right]e^{-t/\tau_0}. \tag{11.1.11}
$$

$\tau_0$ 是**局域平衡弛豫时间**。进一步取其平均值 $\overline\tau_0$，它与两次碰撞之间的平均自由时间同量级。

把运动项与碰撞项相加：

$$
\frac{\partial f}{\partial t}
+v_x\frac{\partial f}{\partial x}
+v_y\frac{\partial f}{\partial y}
+v_z\frac{\partial f}{\partial z}
+X\frac{\partial f}{\partial v_x}
+Y\frac{\partial f}{\partial v_y}
+Z\frac{\partial f}{\partial v_z}
=-\frac{f-f^{(0)}}{\tau_0}. \tag{11.1.12}
$$

这就是**玻耳兹曼方程的弛豫时间近似**。

定态时 $\partial f/\partial t=0$：

$$
v_x\frac{\partial f}{\partial x}
+v_y\frac{\partial f}{\partial y}
+v_z\frac{\partial f}{\partial z}
+X\frac{\partial f}{\partial v_x}
+Y\frac{\partial f}{\partial v_y}
+Z\frac{\partial f}{\partial v_z}
=-\frac{f-f^{(0)}}{\tau_0}. \tag{11.1.13}
$$

---

## §11.2　气体的黏性现象

### 11.2.1 黏性定律

设气体沿 $y$ 方向流动，宏观速度 $v_0=v_0(x)$ 随 $x$ 变化。取平面 $x=x_0$，$x>x_0$ 一侧为正侧。正侧通过单位面积对负侧的作用力为 $p_{xy}$。牛顿黏性定律为

$$
p_{xy}
=\eta\frac{\mathrm dv_0}{\mathrm dx}, \tag{11.2.1}
$$

其中 $\eta$ 是黏度。

微观上，$x>x_0$ 的分子平均具有较大的 $v_y$，穿过平面后把定向动量输运到负侧；反向穿过的分子输运较小动量。单位时间单位面积输运的净 $y$ 动量为

$$
p_{xy}
=-\int_{-\infty}^{+\infty}
mv_xv_yf\,\mathrm d\omega. \tag{11.2.2}
$$

### 11.2.2 局部平衡分布与一阶修正

若 $v_0$ 与 $x$ 无关，系统处于平衡，分布为

$$
f^{(0)}
=n\left(\frac{m}{2\pi kT}\right)^{3/2}
\exp\left[
-\frac{m}{2kT}
\left(v_x^2+(v_y-v_0)^2+v_z^2\right)
\right]. \tag{11.2.3}
$$

代入式 (11.2.2) 时被积函数是 $v_x$ 的奇函数，故 $p_{xy}=0$。

若 $v_0=v_0(x)$，局域平衡分布 $f^{(0)}$ 只含缓变的 $v_0(x)$，并不满足定态玻耳兹曼方程。设速度梯度很小，令

$$
f=f^{(0)}+f^{(1)}, \tag{11.2.5}
$$

其中 $f^{(1)}\ll f^{(0)}$。无外力且 $f$ 只依赖 $x$ 时，式 (11.1.13) 化为

$$
v_x\frac{\partial f}{\partial x}
=-\frac{f-f^{(0)}}{\tau_0}. \tag{11.2.4}
$$

保留一阶小量：

$$
v_x\frac{\partial f^{(0)}}{\partial x}
=-\frac{f^{(1)}}{\tau_0}. \tag{11.2.6}
$$

由于 $f^{(0)}$ 通过 $v_0(x)$ 依赖 $x$，

$$
\frac{\partial f^{(0)}}{\partial x}
=\frac{\partial f^{(0)}}{\partial v_y}
\frac{\mathrm dv_0}{\mathrm dx}.
$$

因此

$$
f^{(1)}
=-\tau_0v_x
\frac{\partial f^{(0)}}{\partial x}
=\tau_0v_x
\frac{\partial f^{(0)}}{\partial v_y}
\frac{\mathrm dv_0}{\mathrm dx}, \tag{11.2.7}
$$

$$
f
=f^{(0)}
+\tau_0v_x
\frac{\partial f^{(0)}}{\partial v_y}
\frac{\mathrm dv_0}{\mathrm dx}. \tag{11.2.8}
$$

### 11.2.3 黏度的微观表达式

将式 (11.2.8) 代入式 (11.2.2)。$f^{(0)}$ 项积分为零，只剩修正项：

$$
p_{xy}
=-\int
mv_x^2v_y\tau_0
\frac{\partial f^{(0)}}{\partial v_y}
\frac{\mathrm dv_0}{\mathrm dx}
\mathrm d\omega.
$$

与式 (11.2.1) 比较：

$$
\eta
=-m\int
v_x^2v_y\tau_0
\frac{\partial f^{(0)}}{\partial v_y}
\mathrm d\omega. \tag{11.2.9}
$$

分部积分：

$$
\int_{-\infty}^{+\infty}
v_y\frac{\partial f^{(0)}}{\partial v_y}
\mathrm dv_y
=\left.f^{(0)}v_y\right|_{-\infty}^{+\infty}
-\int_{-\infty}^{+\infty}f^{(0)}\mathrm dv_y
=-\int f^{(0)}\mathrm dv_y.
$$

因此

$$
\eta
=m\overline{\tau}_0\int v_x^2f^{(0)}\mathrm d\omega
=nm\overline{\tau}_0\overline{v_x^2}. \tag{11.2.10}
$$

由能量均分：

$$
\frac12m\overline{v_x^2}
=\frac12kT,
$$

所以

$$
\eta=nkT\overline{\tau}_0. \tag{11.2.11}
$$

### 11.2.4 与初级动理论比较

以 $\bar l$ 表示平均自由程、$\bar v$ 表示平均速率：

$$
\bar l
=\bar v\overline{\tau}_0. \tag{11.2.12}
$$

于是

$$
\eta
=nkT\frac{\bar l}{\bar v}. \tag{11.2.13}
$$

$\bar l\propto1/n$，而 $\bar v\propto\sqrt T$，故

$$
\eta\propto\sqrt T. \tag{11.2.14}
$$

在温度一定时，黏度与压强无关。这一结论由麦克斯韦在 1860 年首先得到。

若再近似 $\overline{v_x^2}\approx\overline{v^2}/3$，忽略均方根速率与平均速率差别，并令 $\bar v\overline\tau_0=\bar l$，则

$$
\eta
=\frac13nm\bar v\bar l. \tag{11.2.15}
$$

这就是初级动理论的黏度公式。

---

## §11.3　金属的电导率

### 11.3.1 欧姆定律的微观出发点

金属内部有均匀恒定电场 $E_z$，电流密度与电场成正比：

$$
J_z=\sigma E_z. \tag{11.3.1}
$$

$\sigma$ 是电导率。单位体积内速度间隔 $\mathrm d\omega$ 的电子数写成

$$
f\frac{2m^3}{h^3}\mathrm d\omega,
$$

其中因子 $2$ 来自电子自旋。电流密度等于电子数密度乘速度再乘电荷 $-e$：

$$
J_z
=(-e)\int
fv_z\frac{2m^3}{h^3}\mathrm d\omega. \tag{11.3.2}
$$

无电场时，$f$ 是费米分布

$$
f^{(0)}
=\frac1{
\exp\left[
\beta\left(\frac{p^2}{2m}-\mu\right)
\right]+1
}. \tag{11.3.3}
$$

被积函数为 $v_z$ 的奇函数，故 $J_z=0$。

### 11.3.2 弱电场下的一阶修正

定态下，电子受电场力

$$
\dot v_z
=-\frac{eE_z}{m}.
$$

无空间不均匀时，式 (11.1.13) 化为

$$
-\frac{eE_z}{m}
\frac{\partial f}{\partial v_z}
=-\frac{f-f^{(0)}}{\tau_0}. \tag{11.3.4}
$$

令

$$
f=f^{(0)}+f^{(1)}, \tag{11.3.5}
$$

保留一阶项：

$$
\frac{eE_z}{m}
\frac{\partial f^{(0)}}{\partial v_z}
=\frac{f^{(1)}}{\tau_0}.
$$

因此

$$
f
=f^{(0)}
+\frac{eE_z}{m}\tau_0
\frac{\partial f^{(0)}}{\partial v_z}. \tag{11.3.6}
$$

### 11.3.3 电导率

代入电流密度：

$$
J_z
=-\frac{e^2E_z}{m}
\int\tau_0v_z
\frac{\partial f^{(0)}}{\partial v_z}
\frac{2m^3}{h^3}\mathrm d\omega.
$$

$\partial f^{(0)}/\partial v_z$ 仅在费米能 $\varepsilon\approx\mu$ 附近显著，因此取该处弛豫时间的值 $\tau_F$：

$$
J_z
=-\frac{e^2E_z}{m}\tau_F
\int v_z
\frac{\partial f^{(0)}}{\partial v_z}
\frac{2m^3}{h^3}\mathrm d\omega. \tag{11.3.7}
$$

对 $v_z$ 作分部积分：

$$
\int_{-\infty}^{+\infty}
v_z\frac{\partial f^{(0)}}{\partial v_z}\mathrm dv_z
=\left.f^{(0)}v_z\right|_{-\infty}^{+\infty}
-\int_{-\infty}^{+\infty}f^{(0)}\mathrm dv_z
=-\int f^{(0)}\mathrm dv_z.
$$

于是

$$
J_z
=\frac{e^2E_z}{m}\tau_F
\int f^{(0)}
\frac{2m^3}{h^3}\mathrm d\omega
=\frac{ne^2\tau_F}{m}E_z. \tag{11.3.8}
$$

与式 (11.3.1) 比较：

$$
\sigma
=\frac{ne^2\tau_F}{m}. \tag{11.3.9}
$$

### 11.3.4 高温下的温度依赖

高温时自由电子主要受离子振动（声子）散射。设 $l_F$、$v_F$ 分别为费米面附近电子的自由程和速率：

$$
l_F=\tau_Fv_F. \tag{11.3.10}
$$

爱因斯坦模型中，高温

$$
\frac{\hbar\omega}{kT}\ll1
$$

时声子数近似为

$$
n(\omega)
\sim\frac1{e^{\hbar\omega/(kT)}-1}
\sim\frac{kT}{\hbar\omega}.
$$

声子密度正比于 $T$，电子平均自由程反比于声子密度，因此 $l_F\propto1/T$。$v_F$ 随温度变化很弱，故

$$
\sigma\propto\frac1T. \tag{11.3.11}
$$

这与高温金属电导率的实验趋势一致。

---

## §11.4　玻耳兹曼积分微分方程

### 11.4.1 二体弹性碰撞

采用弹性刚球模型，只考虑平动能的交换，不讨论内部激发。稀薄气体中三体以上同时碰撞概率很小，只需考虑两体碰撞。

设两球质量为 $m_1,m_2$，直径为 $d_1,d_2$，碰前速度为 $\boldsymbol v_1,\boldsymbol v_2$，碰后速度为 $\boldsymbol v_1',\boldsymbol v_2'$。动量和能量守恒：

$$
\begin{cases}
m_1\boldsymbol v_1+m_2\boldsymbol v_2
=m_1\boldsymbol v_1'+m_2\boldsymbol v_2',\\[1mm]
\dfrac12m_1v_1^2+\dfrac12m_2v_2^2
=\dfrac12m_1v_1'^2+\dfrac12m_2v_2'^2.
\end{cases} \tag{11.4.1}
$$

共有四个标量方程，但碰后速度有六个分量，因此还需确定碰撞方向。令 $\boldsymbol n$ 为由第一个球心指向第二个球心的单位矢量。碰撞力沿 $\boldsymbol n$，故速度改变量也沿 $\boldsymbol n$：

$$
\boldsymbol v_1'-\boldsymbol v_1
=\lambda_1\boldsymbol n,
\qquad
\boldsymbol v_2'-\boldsymbol v_2
=\lambda_2\boldsymbol n. \tag{11.4.2}
$$

代入式 (11.4.1) 解得

$$
\lambda_1
=\frac{2m_2}{m_1+m_2}
(\boldsymbol v_2-\boldsymbol v_1)\cdot\boldsymbol n,
$$

$$
\lambda_2
=-\frac{2m_1}{m_1+m_2}
(\boldsymbol v_2-\boldsymbol v_1)\cdot\boldsymbol n.
$$

因此

$$
\begin{cases}
\boldsymbol v_1'
=\boldsymbol v_1
+\dfrac{2m_2}{m_1+m_2}
\left[(\boldsymbol v_2-\boldsymbol v_1)\cdot\boldsymbol n
\right]\boldsymbol n,\\[2mm]
\boldsymbol v_2'
=\boldsymbol v_2
-\dfrac{2m_1}{m_1+m_2}
\left[(\boldsymbol v_2-\boldsymbol v_1)\cdot\boldsymbol n
\right]\boldsymbol n.
\end{cases} \tag{11.4.3}
$$

两式相减：

$$
\boldsymbol v_2'-\boldsymbol v_1'
=\boldsymbol v_2-\boldsymbol v_1
-2\left[
(\boldsymbol v_2-\boldsymbol v_1)\cdot\boldsymbol n
\right]\boldsymbol n. \tag{11.4.4}
$$

平方后

$$
(\boldsymbol v_2'-\boldsymbol v_1')^2
=(\boldsymbol v_2-\boldsymbol v_1)^2, \tag{11.4.5}
$$

即相对速率不变。与 $\boldsymbol n$ 点乘：

$$
(\boldsymbol v_2'-\boldsymbol v_1')\cdot\boldsymbol n
=-(\boldsymbol v_2-\boldsymbol v_1)\cdot\boldsymbol n. \tag{11.4.6}
$$

相对速度沿碰撞法向的投影在碰撞中反号。

### 11.4.2 反碰撞

由式 (11.4.3) 可反解出碰前速度：

$$
\begin{cases}
\boldsymbol v_1
=\boldsymbol v_1'
+\dfrac{2m_2}{m_1+m_2}
\left[(\boldsymbol v_2'-\boldsymbol v_1')\cdot(-\boldsymbol n)
\right](-\boldsymbol n),\\[2mm]
\boldsymbol v_2
=\boldsymbol v_2'
-\dfrac{2m_1}{m_1+m_2}
\left[(\boldsymbol v_2'-\boldsymbol v_1')\cdot(-\boldsymbol n)
\right](-\boldsymbol n).
\end{cases} \tag{11.4.7}
$$

因此若以速度 $\boldsymbol v_1',\boldsymbol v_2'$ 和碰撞方向 $-\boldsymbol n$ 相碰，碰后恰为 $\boldsymbol v_1,\boldsymbol v_2$。这叫**反碰撞**。

### 11.4.3 碰撞圆柱与元碰撞数

以第一个分子球心为球心，以

$$
d_{12}=\frac12(d_1+d_2)
$$

为半径作虚球。第二个分子球心在碰撞瞬间必须落在虚球上。相对速度为 $\boldsymbol v_2-\boldsymbol v_1$，相对速率为 $v_r=|\boldsymbol v_2-\boldsymbol v_1|$。若相对速度与 $\boldsymbol n$ 的夹角为 $\theta$，则在 $\mathrm dt$ 内发生指定方向碰撞的碰撞圆柱体积为

$$
d_{12}^2v_r\cos\theta\,\mathrm d\Omega\,\mathrm dt.
$$

只有 $0\leqslant\theta\leqslant\pi/2$。一个速度在 $\mathrm d\omega_1$ 中的分子在 $\mathrm dt$ 内与速度在 $\mathrm d\omega_2$ 中分子发生碰撞的次数为

$$
f_2\mathrm d\omega_2\,
d_{12}^2v_r\cos\theta\,\mathrm d\Omega\,\mathrm dt. \tag{11.4.8}
$$

定义

$$
\Lambda\,\mathrm d\Omega
=d_{12}^2(\boldsymbol v_1-\boldsymbol v_2)\cdot\boldsymbol n
\,\mathrm d\Omega
=d_{12}^2v_r\cos\theta\,\mathrm d\Omega.
$$

则

$$
f_2\Lambda\,\mathrm d\omega_2\,\mathrm d\Omega\,\mathrm dt. \tag{11.4.9}
$$

再乘以 $\mathrm d\tau\mathrm d\omega_1$ 内分子数 $f_1\mathrm d\tau\mathrm d\omega_1$，元碰撞数为

$$
f_1f_2
\mathrm d\omega_1\mathrm d\omega_2
\Lambda\mathrm d\Omega\,\mathrm dt\,\mathrm d\tau. \tag{11.4.10}
$$

对应的反碰撞数，用碰后速度与 $-\boldsymbol n$ 表示，为

$$
f_1'f_2'
\mathrm d\omega_1'\mathrm d\omega_2'
\Lambda'\mathrm d\Omega\,\mathrm dt\,\mathrm d\tau, \tag{11.4.11}
$$

其中

$$
f_1'=f(\boldsymbol r,\boldsymbol v_1',t),
\qquad
f_2'=f(\boldsymbol r,\boldsymbol v_2',t),
\qquad
\Lambda'
=d_{12}^2(\boldsymbol v_1'-\boldsymbol v_2')\cdot\boldsymbol n'.
$$

### 11.4.4 速度体积元的变换

元碰撞使 $\mathrm d\omega_1$ 中分子数减少，反碰撞使其增加：

$$
\left(\frac{\partial f_1}{\partial t}\right)_{\mathrm c}
\mathrm dt\,\mathrm d\omega_1\mathrm d\tau. \tag{11.4.12}
$$

把式 (11.4.11) 中的 $\mathrm d\omega_1'\mathrm d\omega_2'$ 换回碰前变量：

$$
\mathrm d\omega_1'\mathrm d\omega_2'
=|J|\,\mathrm d\omega_1\mathrm d\omega_2, \tag{11.4.13}
$$

其中

$$
J=
\frac{\partial(v_{1x}',v_{1y}',v_{1z}',v_{2x}',v_{2y}',v_{2z}')}
{\partial(v_{1x},v_{1y},v_{1z},v_{2x},v_{2y},v_{2z})}.
$$

式 (11.4.3) 和式 (11.4.7) 给出互逆映射，故 $J'=J^{-1}$；行列式相乘法则又给出 $J'J=1$，于是 $J^2=1$，即

$$
|J|=1.
$$

同时 $\Lambda'=\Lambda$。反碰撞数变为

$$
f_1'f_2'
\mathrm d\omega_1\mathrm d\omega_2
\Lambda\mathrm d\Omega\,\mathrm dt\,\mathrm d\tau. \tag{11.4.14}
$$

净增加减去元碰撞损失，并对 $\boldsymbol v_2$ 与碰撞方向积分：

$$
\left(\frac{\partial f}{\partial t}\right)_{\mathrm c}
=\iint(f_1'f_2'-f_1f_2)
\mathrm d\omega_1\Lambda\,\mathrm d\Omega. \tag{11.4.15}
$$

### 11.4.5 玻耳兹曼积分微分方程

运动项与碰撞项相加：

$$
\begin{aligned}
&\frac{\partial f}{\partial t}
+v_x\frac{\partial f}{\partial x}
+v_y\frac{\partial f}{\partial y}
+v_z\frac{\partial f}{\partial z}
+X\frac{\partial f}{\partial v_x}
+Y\frac{\partial f}{\partial v_y}
+Z\frac{\partial f}{\partial v_z}\\
&\qquad
=\iint(f'f_1'-ff_1)
\mathrm d\omega_1\Lambda\,\mathrm d\Omega.
\end{aligned} \tag{11.4.16}
$$

积分限为

$$
\int\mathrm d\omega_1
=\int_{-\infty}^{+\infty}
\int_{-\infty}^{+\infty}
\int_{-\infty}^{+\infty}
\mathrm dv_{1x}\mathrm dv_{1y}\mathrm dv_{1z},
$$

$$
\int\mathrm d\Omega
=\int_0^{2\pi}\mathrm d\varphi
\int_0^{\pi/2}\sin\theta\,\mathrm d\theta.
$$

式 (11.4.16) 是**玻耳兹曼积分微分方程**，其中碰撞项对 $f$ 非线性。

### 11.4.6 分子混沌假设

严格地说，两分子分别处于 $\mathrm d\tau_1\mathrm d\omega_1$、$\mathrm d\tau_2\mathrm d\omega_2$ 的概率由双粒子分布 $f_2$ 给出：

$$
\frac1{N^2}
f(\boldsymbol r_1,\boldsymbol v_1,\boldsymbol r_2,\boldsymbol v_2,t)
\mathrm d\tau_1\mathrm d\omega_1
\mathrm d\tau_2\mathrm d\omega_2. \tag{11.4.17}
$$

若两分子概率相互独立：

$$
\frac{
f(\boldsymbol r_1,\boldsymbol v_1,t)\mathrm d\tau_1\mathrm d\omega_1
}{N}
\frac{
f(\boldsymbol r_2,\boldsymbol v_2,t)\mathrm d\tau_2\mathrm d\omega_2
}{N}. \tag{11.4.18}
$$

这就是**分子混沌假设**。碰撞时两分子显然存在力学关联，因此式 (11.4.18) 是近似。

> **补充说明**：从刘维尔方程出发可得到 $f_1,f_2,\cdots$ 的 BBGKY 方程链，每一级含下一级分布函数。要得到不可逆的玻耳兹曼型方程，必须引入额外统计假设。

---

## §11.5　$H$ 定理

### 11.5.1 $H$ 函数

玻耳兹曼定义

$$
H
=\iint f(\boldsymbol r,\boldsymbol v,t)
\ln f(\boldsymbol r,\boldsymbol v,t)
\mathrm d\tau\mathrm d\omega
=\iint f\ln f\,\mathrm d\tau\mathrm d\omega. \tag{11.5.1}
$$

其时间变化率为

$$
\frac{\mathrm dH}{\mathrm dt}
=\iint(1+\ln f)
\frac{\partial f}{\partial t}
\mathrm d\tau\mathrm d\omega.
$$

玻耳兹曼方程为

$$
\begin{aligned}
&\frac{\partial f}{\partial t}
+v_x\frac{\partial f}{\partial x}
+v_y\frac{\partial f}{\partial y}
+v_z\frac{\partial f}{\partial z}
+X\frac{\partial f}{\partial v_x}
+Y\frac{\partial f}{\partial v_y}
+Z\frac{\partial f}{\partial v_z}\\
&\qquad
=\iint(f'f_1'-ff_1)
\mathrm d\omega_1\Lambda\,\mathrm d\Omega.
\end{aligned} \tag{11.5.2}
$$

代入后

$$
\begin{aligned}
\frac{\mathrm dH}{\mathrm dt}
={}&-\iint(1+\ln f)
\left(
v_x\frac{\partial f}{\partial x}
+v_y\frac{\partial f}{\partial y}
+v_z\frac{\partial f}{\partial z}
\right)\mathrm d\tau\mathrm d\omega\\
&-\iint(1+\ln f)
\left(
X\frac{\partial f}{\partial v_x}
+Y\frac{\partial f}{\partial v_y}
+Z\frac{\partial f}{\partial v_z}
\right)\mathrm d\tau\mathrm d\omega\\
&-\iiint(1+\ln f)
(ff_1-f'f_1')
\mathrm d\tau\mathrm d\omega_1\Lambda\,\mathrm d\Omega\,\mathrm d\omega.
\end{aligned} \tag{11.5.3}
$$

### 11.5.2 运动项不改变 $H$

空间项可写成散度：

$$
-\int(1+\ln f)
(\boldsymbol v\cdot\nabla f)\mathrm d\tau
=-\int\nabla\cdot(\boldsymbol vf\ln f)\mathrm d\tau
=-\oint\mathrm d\boldsymbol\Sigma\cdot\boldsymbol vf\ln f.
$$

分子不能穿出器壁，边界上 $f=0$，所以该项为零。

速度项：

$$
-\int(1+\ln f)
\left[
\frac{\partial}{\partial v_x}(Xf)
+\frac{\partial}{\partial v_y}(Yf)
+\frac{\partial}{\partial v_z}(Zf)
\right]\mathrm d\omega.
$$

利用式 (11.1.6)，它等于

$$
-\int\left[
\frac{\partial}{\partial v_x}(Xf\ln f)
+\frac{\partial}{\partial v_y}(Yf\ln f)
+\frac{\partial}{\partial v_z}(Zf\ln f)
\right]\mathrm d\omega.
$$

当 $v_i\to\pm\infty$ 时 $f=0$，每一项的无限远边界值都为零。故运动项不改变 $H$。只剩碰撞项：

$$
\frac{\mathrm dH}{\mathrm dt}
=-\iiint(1+\ln f)
(ff_1-f'f_1')
\mathrm d\tau\mathrm d\omega_1\Lambda\,\mathrm d\Omega\,\mathrm d\omega. \tag{11.5.4}
$$

### 11.5.3 对称化

式 (11.5.4) 中把 $\boldsymbol v$ 与 $\boldsymbol v_1$ 互换，积分值不变：

$$
\frac{\mathrm dH}{\mathrm dt}
=-\iiint(1+\ln f_1)
(ff_1-f'f_1')
\mathrm d\tau\mathrm d\omega_1\Lambda\,\mathrm d\Omega\,\mathrm d\omega.
$$

两式相加除以 2：

$$
\frac{\mathrm dH}{\mathrm dt}
=-\frac12\iiint
\left[2+\ln(ff_1)\right]
(ff_1-f'f_1')
\mathrm d\tau\mathrm d\omega_1\Lambda\,\mathrm d\Omega\,\mathrm d\omega. \tag{11.5.5}
$$

再利用碰撞与反碰撞的对称性，在积分中替换

$$
\boldsymbol v\leftrightarrow\boldsymbol v',
\qquad
\boldsymbol v_1\leftrightarrow\boldsymbol v_1',
\qquad
\mathrm d\omega\,\mathrm d\omega_1
=\mathrm d\omega'\,\mathrm d\omega_1',
\qquad
\Lambda'=\Lambda.
$$

得到

$$
\frac{\mathrm dH}{\mathrm dt}
=-\frac12\iiint
\left[2+\ln(f'f_1')\right]
(f'f_1'-ff_1)
\mathrm d\tau\mathrm d\omega_1\Lambda\,\mathrm d\Omega\,\mathrm d\omega.
$$

与式 (11.5.5) 相加除以 2：

$$
\frac{\mathrm dH}{\mathrm dt}
=-\frac14\iiint
\left[\ln(ff_1)-\ln(f'f_1')\right]
(ff_1-f'f_1')
\mathrm d\tau\mathrm d\omega_1\Lambda\,\mathrm d\Omega\,\mathrm d\omega. \tag{11.5.6}
$$

### 11.5.4 单调性与平衡条件

在式 (11.5.6) 中令

$$
x=\ln(ff_1),
\qquad
y=\ln(f'f_1'),
$$

则被积函数的核心是

$$
F(x,y)
=(x-y)(e^x-e^y).
$$

- $x>y$ 时，$e^x>e^y$，所以 $F>0$；
- $x<y$ 时，$e^x<e^y$，仍有 $F>0$；
- 只有 $x=y$ 时 $F=0$。

所以积分非负，

$$
\frac{\mathrm dH}{\mathrm dt}\leqslant0. \tag{11.5.7}
$$

当且仅当

$$
ff_1=f'f_1' \tag{11.5.8}
$$

时等号成立。这就是 **$H$ 定理**：碰撞使 $H$ 单调减少；达到最小值且不再改变时，系统达到平衡。

$H$ 与熵的方向相反，$H$ 定理与熵增加原理等价。但它不是普遍力学定律，而是由玻耳兹曼方程得到的统计规律，适用于稀薄单原子经典气体。

### 11.5.5 可逆性与逆转向题

哈密顿正则方程：

$$
\dot q_i=\frac{\partial H}{\partial p_i},
\qquad
\dot p_i=-\frac{\partial H}{\partial q_i},
\qquad i=1,2,\cdots,s. \tag{11.5.9}
$$

动量反号后哈密顿量不变：

$$
H(q,-p)=H(q,p). \tag{11.5.10}
$$

若

$$
q_i=f_i(t),
\qquad
p_i=g_i(t) \tag{11.5.11}
$$

是解，则

$$
q_i'=f_i(-t),
\qquad
p_i'=-g_i(-t) \tag{11.5.12}
$$

也是解。微观力学运动可逆。

洛施密特据此提出**逆转向题**：若在某一时刻把所有速度反向，系统沿原相轨道逆行，似乎会使 $H$ 增加。玻耳兹曼的回答是：$H$ 不是单一力学量，而是统计平均量。由式 (11.1.1)，$f\mathrm d\tau\mathrm d\omega$ 是分子数的统计平均值，而 $H$ 是 $\ln f$ 的统计平均：

$$
H=\overline{\ln f}. \tag{11.5.13}
$$

$\mathrm dH/\mathrm dt$ 实质上是宏观小、微观长时间间隔内的统计平均变化率。因此 $H$ 定理说明系统对初始条件作统计平均后具有不可逆方向，并不否定微观力学方程的可逆性。

---

## §11.6　细致平衡原理与平衡态的分布函数

### 11.6.1 细致平衡

平衡时 $H$ 取极小，由式 (11.5.8)：

$$
f_1f_2=f_1'f_2'. \tag{11.6.1}
$$

这意味着每个元碰撞与对应反碰撞的速率恰好相等。凡是一个元过程与相应反过程相互抵消，就称为**细致平衡**。细致平衡保证总平衡，但总平衡不一定普遍保证细致平衡；在分子碰撞中，细致平衡由 $H$ 定理证明。

定态时 $\partial f/\partial t=0$，玻耳兹曼方程给出

$$
v_x\frac{\partial f}{\partial x}
+v_y\frac{\partial f}{\partial y}
+v_z\frac{\partial f}{\partial z}
+X\frac{\partial f}{\partial v_x}
+Y\frac{\partial f}{\partial v_y}
+Z\frac{\partial f}{\partial v_z}=0. \tag{11.6.2}
$$

### 11.6.2 平衡分布的五种守恒量

对式 (11.6.1) 取对数：

$$
\ln f_1+\ln f_2
=\ln f_1'+\ln f_2'. \tag{11.6.3}
$$

$\ln f$ 是碰撞前后的守恒量。二体碰撞中粒子数、三个动量分量和能量守恒，因此五个特解为

$$
1,\qquad
mv_x,\qquad
mv_y,\qquad
mv_z,\qquad
\frac12mv^2. \tag{11.6.4}
$$

通解是它们的线性组合：

$$
\ln f
=\alpha_0
+\alpha_1mv_x
+\alpha_2mv_y
+\alpha_3mv_z
+\alpha_4\frac12m(v_x^2+v_y^2+v_z^2). \tag{11.6.5}
$$

把五个常数改写成粒子数密度 $n$、温度 $T$、整体速度 $\boldsymbol v_0=(v_{0x},v_{0y},v_{0z})$：

$$
f
=n\left(\frac{m}{2\pi kT}\right)^{3/2}
\exp\left[
-\frac{m}{2kT}
\left(
(v_x-v_{0x})^2
+(v_y-v_{0y})^2
+(v_z-v_{0z})^2
\right)
\right]. \tag{11.6.6}
$$

五个参数的物理意义为

$$
\overline{v_x}=v_{0x},
\qquad
\overline{v_y}=v_{0y},
\qquad
\overline{v_z}=v_{0z}. \tag{11.6.7}
$$

一般来说，$n,T,\boldsymbol v_0$ 可以依赖坐标。

### 11.6.3 运动项为零的条件

把式 (11.6.6) 代入式 (11.6.2)，两边约去 $f$：

$$
\boldsymbol v\cdot\nabla
\left[
\ln n
+\frac32\ln\frac{m}{2\pi kT}
-\frac{m}{2kT}
(\boldsymbol v-\boldsymbol v_0)^2
\right]
-\frac{m}{kT}
\boldsymbol F\cdot(\boldsymbol v-\boldsymbol v_0)=0. \tag{11.6.8}
$$

其中

$$
\boldsymbol F=(X,Y,Z).
$$

该式对所有速度成立，所以各阶速度项系数分别相等。

### 11.6.4 温度与整体速度的限制

比较 $v$ 的三次项，得到

$$
\boldsymbol v\cdot\nabla T=0,
$$

即

$$
\frac{\partial T}{\partial x}
=\frac{\partial T}{\partial y}
=\frac{\partial T}{\partial z}=0. \tag{11.6.9}
$$

平衡系统的温度必须均匀。

比较速度的二次项：

$$
\boldsymbol v\cdot\nabla(\boldsymbol v\cdot\boldsymbol v_0)=0. \tag{11.6.10}
$$

等价于

$$
\begin{cases}
\dfrac{\partial v_{0x}}{\partial x}
=\dfrac{\partial v_{0y}}{\partial y}
=\dfrac{\partial v_{0z}}{\partial z}=0,\\[2mm]
\dfrac{\partial v_{0y}}{\partial z}
+\dfrac{\partial v_{0z}}{\partial y}=0,\quad
\dfrac{\partial v_{0z}}{\partial x}
+\dfrac{\partial v_{0x}}{\partial z}=0,\quad
\dfrac{\partial v_{0x}}{\partial y}
+\dfrac{\partial v_{0y}}{\partial x}=0.
\end{cases} \tag{11.6.11}
$$

满足这些条件的整体速度场为

$$
\boldsymbol v_0
=\boldsymbol a+\boldsymbol\omega\times\boldsymbol r, \tag{11.6.12}
$$

其中 $\boldsymbol a,\boldsymbol\omega$ 是常矢量。平衡气体只能作匀速平动或绕固定轴匀速转动。

### 11.6.5 密度分布

比较速度的一次项：

$$
\nabla\left(
\ln n-\frac{m}{2kT}v_0^2
\right)
-\frac{m}{kT}\boldsymbol F=0. \tag{11.6.13}
$$

若外力可写成势函数的负梯度，

$$
\boldsymbol F=-\nabla\varphi,
$$

则积分得

$$
n
=n_0
\exp\left[
\frac{m}{2kT}v_0^2
-\frac{m}{kT}\varphi
\right]. \tag{11.6.14}
$$

比较速度的零次项：

$$
\boldsymbol v_0\cdot\boldsymbol F=0. \tag{11.6.15}
$$

因此平衡整体速度场必须与外力垂直。例如在重力场中，整体运动只能发生在水平面内。

### 11.6.6 旋转气体

绕 $z$ 轴以角速度 $\omega$ 转动：

$$
v_{0x}=-\omega y,
\qquad
v_{0y}=\omega x,
\qquad
v_{0z}=0.
$$

取重力势

$$
\varphi=gz.
$$

则

$$
n
=n_0
\exp\left[
\frac{m\omega^2}{2kT}(x^2+y^2)
-\frac{mgz}{kT}
\right]. \tag{11.6.16}
$$

其中

$$
-\frac12m\omega^2(x^2+y^2)
$$

可理解为旋转参考系中的离心势能。

至此，玻耳兹曼积分微分方程既给出了趋向平衡的不可逆性，也确定了平衡态分布和空间分布。

---

## 本章小结

### 一、公式速查

**玻耳兹曼方程与输运**

| 内容 | 公式 |
| --- | --- |
| 弛豫时间近似 | (11.1.12)；定态 (11.1.13) |
| 黏度 | $\eta=nkT\overline\tau_0=nkT\bar l/\bar v$ (11.2.11)(11.2.13) |
| 黏度温度依赖 | $\eta\propto\sqrt T$ (11.2.14) |
| 金属电导率 | $\sigma=ne^2\tau_F/m$ (11.3.9) |
| 高温电导率 | $\sigma\propto1/T$ (11.3.11) |

**碰撞积分与 $H$ 定理**

| 内容 | 公式 |
| --- | --- |
| 弹性碰撞 | 式 (11.4.1)–(11.4.7) |
| 元碰撞与反碰撞 | 式 (11.4.10)(11.4.11)(11.4.14) |
| 碰撞变化率 | (11.4.15) |
| 玻耳兹曼积分微分方程 | (11.4.16) |
| $H$ 函数 | (11.5.1) |
| $H$ 定理 | $\mathrm dH/\mathrm dt\leqslant0$ (11.5.7) |
| 细致平衡 | $ff_1=f'f_1'$ (11.5.8)(11.6.1) |

**平衡分布**

| 内容 | 公式 |
| --- | --- |
| 平衡分布 | 麦克斯韦分布 (11.6.6) |
| 温度均匀 | (11.6.9) |
| 速度场限制 | (11.6.11) |
| 刚体转动速度场 | $\boldsymbol v_0=\boldsymbol a+\boldsymbol\omega\times\boldsymbol r$ (11.6.12) |
| 势场中的数密度 | (11.6.14) |
| 旋转气体数密度 | (11.6.16) |

### 二、概念脉络

```text
§11.1 非平衡分布函数
    f(r,v,t)
        ├── 运动项：相空间通量
        └── 碰撞项：局域平衡偏离的指数衰减
                    │
                    ▼
    弛豫时间近似
        f - f0 按 e^(-t/τ0) 衰减
                    │
        ┌───────────┴───────────┐
        ▼                       ▼
§11.2 黏性                   §11.3 电导
速度梯度 → 动量输运          电场 → 费米面偏移
η = nkT τ0                  σ = n e² τF/m
η ∝ √T                      σ ∝ 1/T（高温）
        └───────────┬───────────┘
                    ▼
§11.4 玻耳兹曼积分微分方程
    二体弹性碰撞
        ├── 动量、能量守恒
        ├── 碰撞圆柱 ΛdΩ
        ├── 元碰撞 - 反碰撞
        └── 分子混沌假设
                    │
                    ▼
§11.5 H 定理
    H = ∫ f ln f
        ├── 运动项只重排相空间，不改 H
        ├── 碰撞项对称化
        └── ∫(x-y)(e^x-e^y) ≥ 0
                    │
                    ▼
    dH/dt ≤ 0，平衡时 ff1 = f'f1'
                    │
                    ▼
§11.6 细致平衡
    ln f 是二体碰撞守恒量
        ├── 五个独立守恒量
        ├── 麦克斯韦分布
        ├── 温度、速度场限制
        └── 势场与旋转气体中的密度分布
```

### 三、易错点

1. **$f(\boldsymbol r,\boldsymbol v,t)$ 是统计平均分布，不是单个分子轨迹**：玻耳兹曼方程描述大量分子的平均行为。
2. **相空间通量包含位置和速度两类边界**：运动项不能只写成 $\boldsymbol v\cdot\nabla f$，外力对速度的作用同样重要。
3. **式 (11.1.6) 是外力条件**：重力和洛伦兹力满足速度散度为零；不是所有任意外力都满足。
4. **弛豫时间近似是唯象模型**：$\tau_0$ 并非由该近似本身确定，需要另行由碰撞机制计算。
5. **局域平衡分布 $f^{(0)}$ 不等于整体平衡分布**：$n,T,\boldsymbol v_0$ 可随位置缓慢变化。
6. **黏性来自定向动量输运**：不是分子随机热运动本身的净质量流；平衡时两种方向的动量输运互相抵消。
7. **一阶修正 $f^{(1)}$ 与速度梯度同阶**：只保留线性响应，梯度很大时该近似失效。
8. **$\eta\propto\sqrt T$ 与压强无关是稀薄气体结果**：稠密流体和液体不适用。
9. **电子电导率不由全部电子共同贡献**：只有费米面附近 $kT$ 壳层中的电子能被电场显著重排。
10. **$\tau_F$ 是费米面附近的弛豫时间**：不能把普通热运动弛豫时间直接代入金属电导公式。
11. **玻耳兹曼积分微分方程的碰撞项是非线性的**：碰撞项含 $ff_1$，完整方程不是线性微分方程。
12. **元碰撞数与反碰撞数不能随意相等**：平衡时才由细致平衡保证相等；非平衡时二者之差驱动分布演化。
13. **速度体积元变换的雅可比绝对值为 1**：碰撞映射保相空间体积，但速度本身会改变。
14. **分子混沌假设忽略了碰撞瞬间的关联**：它把双粒子分布拆成两个单粒子分布之积，是玻耳兹曼方程的关键统计假设。
15. **$H$ 定理不是微观可逆性的否定**：它描述统计平均行为的不可逆性；逆转向题针对的是单个微观轨道。
16. **$H$ 单调减少不意味着每个微观状态都减少**：在统计涨落意义上，最概然演化方向是 $H$ 减小。
17. **细致平衡比总平衡更强**：它要求每个元过程与对应逆过程分别抵消。
18. **平衡分布不只要求速度麦克斯韦分布**：还要求温度均匀、整体速度场满足刚体运动形式、密度满足势场分布。
19. **$\ln f$ 的五个特解不是任意五个函数**：它们分别对应粒子数、三个动量分量和能量的碰撞守恒。
20. **旋转气体中的离心项不是额外温度**：它是旋转参考系中的有效势能，决定密度分布，不改变局域温度。

### 四、名词中英对照

| 中文 | 英文 | 中文 | 英文 |
| --- | --- | --- | --- |
| 非平衡态统计理论 | nonequilibrium statistical theory | 分布函数 | distribution function |
| 相空间通量 | phase-space flux | 弛豫时间 | relaxation time |
| 局域平衡 | local equilibrium | 玻耳兹曼方程 | Boltzmann equation |
| 弛豫时间近似 | relaxation-time approximation | 黏性 | viscosity |
| 速度梯度 | velocity gradient | 平均自由程 | mean free path |
| 平均速率 | mean speed | 电导率 | electrical conductivity |
| 费米面 | Fermi surface | 费米分布 | Fermi distribution |
| 弹性碰撞 | elastic collision | 元碰撞 | direct collision |
| 反碰撞 | inverse collision | 碰撞圆柱 | collision cylinder |
| 碰撞截面 | collision cross-section | 碰撞积分 | collision integral |
| 分子混沌假设 | molecular chaos assumption | 玻耳兹曼积分微分方程 | Boltzmann integro-differential equation |
| BBGKY 方程链 | BBGKY hierarchy | $H$ 定理 | $H$-theorem |
| $H$ 函数 | $H$-function | 细致平衡 | detailed balance |
| 可逆性 | reversibility | 逆转向题 | reversibility paradox |
| 麦克斯韦分布 | Maxwell distribution | 宏观速度 | bulk velocity |
| 局域温度 | local temperature | 刚体转动 | rigid-body rotation |
| 离心势能 | centrifugal potential | 统计平均 | statistical average |
