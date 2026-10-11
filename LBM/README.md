# Lattice Boltzmann solvers

`cavity.jl` は既存の単相 density-based BGK solver、`fakhari.jl` は
Fakhari, Mitchell, Leonardi & Bolster (2017) の conservative phase-field / velocity-based
D2Q9 solver です。共通の lattice と平衡分布の velocity polynomial を `d2q9.jl` に
まとめています。配列は既存と同じ `(direction, x, y)`、速度は `(component, x, y)`。
CPU 用で、`--threads=auto` による行方向の並列実行に対応します。

## 実行

リポジトリのルートで既存の Julia project を使います。

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=auto LBM/test/runtests.jl
julia --project=. LBM/test/output_tests.jl
julia --project=. --threads=auto LBM/rayleigh_taylor.jl
```

以下はファイル生成・起動だけを確認する100ステップの smoke test です。
`tstar=100/16000=0.00625` までなので、界面の移動はほぼ見えません。
RTI の発展を見るには上の既定コマンドで48,000ステップまで実行してください。
smoke test は専用の出力ディレクトリに保存します。

```sh
julia --project=. --threads=4 LBM/rayleigh_taylor.jl --nx=32 --steps=100 --output-every=50 --output-dir=LBM/results/rti_smoke
```

`--help` で設定を表示します。既定の計算は `256 × 1024`、48,000 steps
（`tstar=3`）です。`L` を変更すると `Nx=L, Ny=4L` とし、無次元入力から
物性値を再計算します。`t0` は別の入力で、既定値は常に16,000です。

結果は `LBM/results/rti_fakhari/` に出力します。

- `rti.pvd`: ParaView の時系列。時刻は `tstar=t/t0`。
- `step_*.vtr`: phase field、密度、normalized / physical pressure、速度、total force、viscous correction。
- `diagnostics.csv`: phase mass と drift、phi の min/max、最大 Mach 数、有限性、最小 rho/tau、bubble / spike front と交点数。
- `parameters.txt`: 使用した無次元値と導出した lattice parameters、境界・初期化・時刻配置。

ParaView では `Phase` の `0.5` contour で界面、`Velocity` の glyph / streamline
で流れを確認できます。既定では `tstar=0,1,2,3` に必ず snapshot を保存します。
各 snapshot で PVD を更新するため、計算途中でも `rti.pvd` の Reload Files で
完了済みの時刻を追加できます。単独の `step_0.vtr` は初期状態だけなので、
時系列を見るときは **`rti.pvd` を開き、時間スライダーを進めてください**。
コンソールにも bubble / spike の位置と初期状態からの変位を表示します。
front は `x=0` の最高交点と `x=L/2` の最低交点を線形補間して測定します。
複数交点が生じた場合は交点数も保存します。

## モデルと時刻配置

- `phi=sum(h)`, `rho=rho_L+phi*(rho_H-rho_L)`。
- `pstar=sum(g)`。流体分布の zeroth moment を密度として使いません。
- `g_eq=pstar*w+(Gamma-w)`、`h_eq=phi*Gamma`。
- 両 collision は modified equilibrium `eq-force_i/2` を使います。
- `tau=nu/cs2`、`tau_phi=M/cs2`、collision denominator は `tau+1/2`。
- `Fs=mu_phi*grad(phi)`、`Fb=rho*gvec`、`Fp=-pstar*cs2*grad(rho)`。
- `Fmu=-tau*second_moment(K*(g-g_eq))*grad(rho)`。BGK は `K=I/(tau+1/2)`、
  MRT は `K=Mmat^-1*S*Mmat`。velocity の finite difference は使いません。
- `u=sum(g*e)+F/(2rho)`。
- phi の gradient / Laplacian は指定された D2Q9 isotropic stencil。
  `grad(rho)=(rho_H-rho_L)*grad(phi)` をすべての補正で共用します。

一 step は `collide! → stream! → macroscopic!` です。初期状態では場と分布を
準備済みなので、そのまま `step!` を呼びます。streaming 後に phi と pstar を
再構築し、phi derivatives、rho、tau、力、速度の順に更新します。

stress と half-force 速度の相互依存は explicit に処理します。
`g^(n+1), phi^(n+1), pstar^(n+1)` の stress equilibrium に **前時刻 `u^n`** を使い、
`F^(n+1)` から `u^(n+1)` を回復します。次の collision equilibrium には
この更新後の速度を使います。pressure predictor-corrector は行いません。
この速度の遅れは time-level convention であり、論文の離散アルゴリズムを完全に
再現したと断定するものではありません。

初期分布は `hbar_eq`, `gbar_eq`。odd half-force shift により初期の指定速度と
`sum(g*e)+F/(2rho)` が一致し、non-equilibrium second moment はゼロになります。
毎 step の phi clamp は行いません。overshoot、NaN/Inf、rho/tau の破綻は診断に
残し、runner は非有限値または非正の rho/tau を検出した時点で停止します。

## MRT と境界条件

MRT の具体的 matrix は `fakhari.jl` の `MRT_BASIS` に明記しています。
列は `(rest,E,N,W,S,NE,NW,SW,SE)`、行は
`(density,energy,energy²,jx,qx,jy,qy,pxx,pxy)`。
`S=diag(1,1,1,1,1,1,1,snu,snu)` とし、最後の2行が
`ex²-ey²` と `ex*ey` の stress modes です。実装ではその他の rate が1であることを
使って2つの stress projection だけを計算します。dense な式との一致をテストします。

水平は periodic、上下の壁は halfway hydrodynamic bounce-back。phase field は
壁外の ghost value を隣接 fluid node に反射し、h の壁 link は bounce-back
にして、neutral / zero-normal-gradient と phase no-flux を組み合わせています。
壁 link を含む streaming は population の置換なので質量を保存します。
この離散壁条件は実装上の選択で、論文が指定した境界式とはしていません。

## RTI baseline

`rti_parameters()` が以下の関係から物性を求めます。

```text
rho_H = rho_L*(1+At)/(1-At)
g = L/(At*t0²), U0 = sqrt(g*L)
mu_H = rho_H*U0*L/Re, mu_L = mu_H/mu_star
nu_H = mu_H/rho_H, nu_L = mu_L/rho_L
tau_H = nu_H/cs2, tau_L = nu_L/cs2
sigma = mu_H*U0/Ca, M = U0*L/Pe, tau_phi = M/cs2
```

既定値は `L=256, At=0.5, Re=3000, mu_star=1, Ca=0.26, Pe=1000,
t0=16000, xi=5, rho_L=1`。導出値は以下です。

| 量 | lattice value |
|---|---:|
| rho_H | 3 |
| g | 2.0e-6 |
| U0 | 0.02262741699796952 |
| mu_H = mu_L | 0.005792618751480197 |
| nu_H | 0.001930872917160066 |
| nu_L | 0.005792618751480197 |
| tau_H | 0.005792618751480198 |
| tau_L | 0.01737785625444059 |
| sigma | 0.0005041230769230768 |
| M | 0.005792618751480197 |
| tau_phi | 0.01737785625444059 |

RTI は MRT と `tau=tau_L+phi*(tau_H-tau_L)` の linear interpolation を既定とします。
`Parameters(viscosity_model=:harmonic)` と `:dynamic` も利用できます。

初期界面は `yI=0.1L*cos(2pi*x/L)`、heavy fluid は上側です。
phi は `0.5+0.5tanh(2d/xi)`、**初期化だけ `d≈y-yI(x)`** とします。
これは Eq.(40) の厳密な signed perpendicular distance ではありません。
格子は `x=0,...,L-1`, `y=-2L+0.5,...,2L-0.5`、上下壁は `y=±2L`。
physical pressure は各列で top wall をゼロとする hydrostatic 積分から初期化し、
`pstar=p/(rho*cs2)` に変換します。

## API

```julia
include("LBM/fakhari.jl")
using .Fakhari
setup = rti_parameters()
solver = initialize_rti(setup; collision=:mrt)
step!(solver, 100)
diagnostics(solver)
interface_fronts(solver)
```

一般の初期値は `PhaseFieldLBM(phi, Parameters(...); u, pstar, collision=:mrt)`。
`pstar` の既定値はゼロです。境界は RTI 用の水平 periodic / 上下 wall に固定されています。
テストは equilibrium moments、stencil、MRT、BGK/MRT stress、四つの力、modified
equilibria、half-force velocity、streaming、uniform state、RTI parameters、threading を
検証します。標準格子の1,000ステップ回帰テストでは、bubble の上昇と heavy spike
の下降も確認します。論文 Fig.5 の定量再現は、baseline の長時間計算と界面位置の比較で
別途評価する必要があります。

## 長時間の動作検証

標準 `256×1024`、MRT、linear tau interpolation、指定の Fig.5 パラメータで
48,000ステップ（`tstar=3`）まで実行しました。phi=0.5 の測定結果は以下です。

| tstar | bubble front y | heavy spike y |
|---:|---:|---:|
| 0 | 25.6038 | -25.6038 |
| 1 | 73.4804 | -89.9283 |
| 2 | 141.7309 | -215.8660 |
| 3 | 212.6157 | -356.4171 |

最終相質量の相対 drift は `-5.13e-12`、最大 Mach 数は `0.07190`、
phi の範囲は `[-0.0007953, 1.0007864]`。NaN/Inf はありませんでした。
VTK の実データから、spike の下降、bubble の上昇、逆向きの渦対とその発展を
画像でも確認しています。検証結果は `LBM/results/rti_debug_baseline/rti.pvd`、
診断は同ディレクトリの `diagnostics.csv`、画像は `rti_evolution.png` に保存しました。
結果ディレクトリは既存の `.gitignore` に従って Git 管理対象外です。

再実行するコマンド:

```sh
julia --project=. --threads=16 LBM/rayleigh_taylor.jl --steps=48000 --output-every=4000 --diagnostic-every=1000 --output-dir=LBM/results/rti_debug_baseline
```
