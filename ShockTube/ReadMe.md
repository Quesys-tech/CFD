# 衝撃波管問題 
## 問題設定
衝撃波管問題は一定の直径を持つ1次元とみなせるパイプに、比熱比が等しく熱力学的状態が異なる気体をが仕切られて入れている。時刻 $t=0$ で気体を仕切る壁を取り去って衝撃波、接触不連続、膨張波が発生する問題である。
Sodの問題[[Sod1978](https://doi.org/10.1016/0021-9991(78)90023-2)]が有名である。

本プログラムではSodの問題を解く。
支配方程式は
$$
\frac{\partial \bm{Q}}{\partial t} + 
\frac{\partial \bm{E}}{\partial x} = 0, \quad 
\bm{Q} = 
\begin{bmatrix}
\rho \\ 
\rho u \\
e
\end{bmatrix}, \quad 
\bm{E} = 
\begin{bmatrix}
\rho u \\ 
(\gamma -1)e + \frac{1-\gamma}{2}\rho u^2 \\
\gamma e u - \frac{\gamma - 1}{2}\rho u^3
\end{bmatrix}  
$$
となる。$x$は座標、$\rho$は密度、$u$は速度、$e$は内部エネルギーである。
$\bm{Q}$は保存量、$\bm{E}$は流束である。

$t=0$のときにSodの問題と同じく、
- $x< 0$では$p=1$、$\rho=1$、$u=0$
- $x \geq 1$では$p=0.1$、$\rho=0.125$、$u=0$
とする。
