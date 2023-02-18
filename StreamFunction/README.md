# 渦度流れ関数法
## 渦度流れ関数法
2次元非圧縮流れの支配方程式は以下に示す連続の式とNavier-Stokes方程式である.
$$
\begin{aligned}
    \nabla\cdot\bm{u} &= 0 \\
    \frac{\partial \bm{u}}{\partial t} + (\bm{u}\cdot\nabla)\bm{u} &= -\nabla p + \frac{1}{\mathrm{Re}} \nabla^2 \bm{u}
\end{aligned}
$$
$\bm{u} = (u,v)$を速度, $t$を時間, $p$を圧力, $\mathrm{Re}$をReynolds数とする. 
連続の式の成分を書き下すと
$$
\frac{\partial u}{\partial x} + 
\frac{\partial v}{\partial y} =0  
$$
を得る.
ここで
$$
u = \frac{\partial \Psi}{\partial y},\quad v = -\frac{\partial \Psi}{\partial x} 
$$
とする流れ関数$\Psi$を定義すると連続の式は自動的に満たされる.
次に渦度$\omega$を
$$
\omega = \frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}
$$
と定義する. 流れ関数のラプラシアンを取ると
以下のPoisson方程式を得る.
$$
-\nabla ^2 \Psi = \omega
$$

Navier-Stokes方程式のrotを取って以下の渦度輸送方程式を得る
$$
\frac{\partial \omega}{\partial t} + u\frac{\partial \omega}{\partial x} +v\frac{\partial \omega}{\partial y} = \frac{1}{\rm Re}\nabla^2\omega
$$
これを$u$, $v$を流れ関数の定義から流れ関数を用いて表記し, 整理すると
$$
\frac{\partial \omega}{\partial t} = -\frac{\partial \Psi}{\partial y}\frac{\partial \omega}{\partial x} +\frac{\partial \Psi}{\partial x}\frac{\partial \omega}{\partial y} + \frac{1}{\rm Re}\nabla^2\omega
$$
となる.

## 離散化
流れ関数の方程式を離散化する. 空間の離散化の幅を$h$, 時間の離散化の幅を$\Delta t$とする.
流れ関数の方程式の離散化は差分法で離散化を行い, 以下の式を得る.
$$
\frac{-\Psi_{i-1,j}-\Psi_{i,j-1}-\Psi_{i+1,j}-\Psi_{i,j+1}+4\Psi_{i,j}}{h^2}=\omega_{i,j}
$$
次に
渦度輸送方程式を離散化する.