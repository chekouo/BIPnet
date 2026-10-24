# BIPnet (development version)

Changes made on 2026-10-24.

## Bug fixes (these change results)

* **Wrong model-vector length (`dim`) in `mainfunction`** (`src/BIP.c`).
  `dim` was computed as `sum(P[m]^2) * r + (Np - 1) * r` instead of
  `sum(P[m]) * r + r * #{m : IndVar[m] != 1}`, because `dim += P[m] * r` sat inside an
  unnecessary loop over `j`.
  - Memory: each stored MCMC model (`rhomodel`) used about `sum(P^2) * r` bytes instead of
    `sum(P) * r`. With two platforms of 500 markers, r = 4 and 5000 samples, that was about
    10 GB requested instead of about 20 MB.
  - Results: `UniqueModel()` compared the uninitialised bytes beyond the filled entries, so
    identical models could be counted as different. This affects the number of unique models,
    the posterior model probabilities (`PostGam`) and the BMA loadings (`EstLoadModel`).
  - `dim` now uses the same rule (`IndVar[m] != 1`) as the code that fills and reads
    `rhomodel`. Covariate platforms (`IndVar == 2`) are still stored. Their `rho` is always 1,
    so storing them does not affect which models count as unique.

* **Wrong prior variance in the loading updates** (`src/myfunction.c`, `LoadAOther()` and
  `EstimateLoad()`). The prior precision term used `Tau[s][j]`, where `s` is the position
  among the active components, instead of `Tau[NZ1[s]][j]`. Loadings therefore got the
  wrong prior variance whenever the active components were not exactly `0, ..., |S|-1`.
  `logGausQuadForm()` already used the correct index.

* **`Simulate()` was not reproducible** (`R/Simudata.R`). `set.seed(seed)` was only called
  after the true loadings (`runif`) and, in settings 3 and 4, the zeroed variables (`sample`)
  had been drawn. So `Simulate(seed = 1)` gave different `X1`/`X2` on every call, and the
  training (`seed = 1`) and test (`seed = 2`) sets in the examples came from different true
  models.
  - New argument `seedA = 1`. `set.seed(seedA)` is called before the loadings are drawn, and
    `set.seed(seed)` still controls the latent factors and noise. With a fixed `seedA`,
    datasets simulated with different `seed` values share the same true loadings.
  - For a given `seed`, the simulated data differ from earlier versions, which were not
    reproducible anyway.
  - `man/Simulate.Rd` is updated.

## Performance (no change in results)

The `BIPpredict()` example (n = 200, two platforms of 500 markers plus the response,
r = 4, K = 11 groups, 5000 + 1000 MCMC iterations) now runs in about 14 s instead of
111 s, roughly 8x faster.

| Version                                      | MCMC time |
|----------------------------------------------|-----------|
| Bug fixes only                               | 111 s     |
| + precomputed `U'U`, `U'X`                   | 26 s      |
| + cached group sums, restructured `SampleUU` | 14 s      |

* **Precomputed cross-products for the likelihood** (`src/myfunction.c`, `src/BIP.c`,
  `src/utils.c`, `src/header.h`).
  - `logGausQuadForm()` is called about 4·p·r times per platform per iteration. It used to
    rebuild `U_S'U_S` and `U_S'x_j` from the n × r matrix `U` on every call, at a cost of
    O(n·|S|^2), and allocated four GSL objects each time.
  - `U'U` (r × r), `U'X1[m]` (r × p) and the column sums of squares `||x_j||^2` are now
    computed once per iteration (`U'U`) and once per platform (`U'X`, `||x_j||^2`). Each call
    then costs O(|S|^3), independent of n, and makes no heap allocations.
  - `LoadAOther()` and `EstimateLoad()` use the same precomputed quantities.
  - After sampling, they are computed once from `meanU` and reused for every model evaluated
    by `logPostGam()` and `EstimateLoad()`. Centring the response by `InterceptMean`, which
    gave the same result for every model, was moved out of the per-model loop.
  - New helpers in `utils.c`: `crossprod()` and `colsumsq()`.
  - The signatures of `logGausQuadForm`, `logPostGam`, `SamplerhoGamma`, `LoadAOther`,
    `EstimateLoad`, and of the unused `logpost` and `rho1`, now take `UtU`/`UtX`/`xx`
    instead of `U`/`X`.
  - Results agree with the previous code to within floating-point rounding (relative
    difference of about 1e-15). With the same seed, `SamplerhoGamma` made identical draws.
    Summation order differs slightly, so long chains are equivalent but not bit-identical.

* **`GroupEffect()`: cached group sums** (`src/myfunction.c`).
  - The sum `kb[j] = sum_k Path[j][k]·B[k]` was recomputed for every variable, for every
    group, in each of the three moves, at a cost of O(p·K^2) per component.
  - It is now computed once per call for the active variables (`Gam[j][l] == 1`) and
    recomputed, in the same summation order, only for members of group k after a move is
    accepted.
  - The loops skip variables outside group k, for which the original terms were exactly
    zero, and the constant `lgamma`/`log` terms are computed once per call. The cost is now
    about O(p·K).
  - Results are bit-identical to the previous version (2000 repeated calls with K = 40 and
    p = 500).

* **`SampleUU()` restructured** (`src/myfunction.c`).
  - The loops over markers no longer scan all `Np` platforms to find which one each marker
    belongs to.
  - The GSL work vectors are allocated once instead of for each sample.
  - Results are bit-identical.

## Verification

* Old and new C functions were compiled side by side and compared on random inputs with
  the same RNG seed (see the results above).
* The `BIPpredict()` example gives identical output on repeated runs (test MSE = 2.029291)
  before and after the bit-identical performance changes.
