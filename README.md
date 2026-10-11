# 数値流体力学のコード

## ディレクトリ
- `LBM` 格子ボルツマン法のコード
  - `LBM/cavity.jl`: 単相 BGK cavity solver
  - `LBM/rayleigh_taylor.jl`: Fakhari et al. (2017) phase-field / velocity-based MRT Rayleigh-Taylor instability
- `ShockTube` 衝撃波管問題のコード
- `StreamFunction/2DCavity` 作者が2019年に作っていて、発掘した数値流体力学のコードの備忘録。流れ関数渦度法でキャビティ流れを解いている
