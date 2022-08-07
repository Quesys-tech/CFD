# 1次元移流方程式

## 問題設定

$$
\frac{\partial u }{\partial t} + c\frac{\partial u}{\partial x} = 0 
$$
ここで, $c>0$とする.
初期条件は
$$
u(x,t) = 
\begin{cases}
    1 & (0 < x <\frac{1}{2})\\
    0 & otherwise
\end{cases}
$$
ここで境界条件は
$u(0,t)=u(1,t)=0$とする.

解析解は
$$
\begin{aligned}
    u(x,t) &= f(x-ct)
\end{aligned}
$$
ここで, 
$$
f(x) = 
\begin{cases}
    1 & (0 < x <\frac{1}{2})\\
    0 & otherwise
\end{cases}
$$
と与えられる.

## FTCSスキーム

時間の離散化幅を$\Delta t$, 空間の離散化幅を$\Delta x$とし,
$u(i\Delta x, j\Delta t) = u_i^j$とする.

時刻については前進差分,空間については中心差分をとるものとすると,以下の差分方程式が得られる.
$$
\frac{u_i^{j+1} - u_i^{j}}{\Delta t} + c\frac{u_{i+1}^{j}-u_{i-1}^j}{2\Delta x} =0
$$
未知数である$u_i^{j+1}$について解く.ここで$\nu=c\Delta t /\Delta x$とすると,
$$
\begin{aligned}
u_i^{j+1} = u_i^j - c\Delta t\frac{u_{i+1}^{j}-u_{i-1}^j}{2\Delta x }\\
u_i^{j+1} = u_i^j -\frac{1}{2}\nu(u_{i+1}^{j}-u_{i-1}^j)
\end{aligned}
$$
となる.