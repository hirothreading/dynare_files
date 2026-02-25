# Replication log files

## 2026-02-19: Created `hdnwr.mod` — HDNWR model implementation

Created the Dynare model file `sgu_model/hdnwr.mod` implementing the Heterogeneous Downward Nominal Wage Rigidity (HDNWR) model from Schmitt-Grohe and Uribe (2025).

**Model version:** Endogenous labor supply (Definition C1, Appendix C of the paper).

**Structure:**
- 10 endogenous variables: `jstar`, `y`, `h`, `u`, `w`, `ii`, `pi`, `piW`, `mu`, `a`
- 8 structural equations (B1–B5, eq 17, eq 27, eq 28) + 2 shock processes
- Calibration from Table 1 of the paper
- All integrals in equations (17) and (28) computed in closed form using model-local variables (`#` prefix), exploiting the linear specification of `gamma(j)`

**Steady state:**
- `steady_state_model` block with MATLAB `fzero` to solve the nonlinear wage aggregation equation (17) for `jstar`
- All other steady-state values computed analytically

**Verified steady-state values (Python):**
- `jstar = 0.6503` (fraction of unconstrained varieties)
- `u = 0.0596` (unemployment rate, above natural rate of 0.04)
- `h = 0.9056`, `y = 0.9283`, `w = 0.7688`
- `pi = piW = pistar = 0.0074` (quarterly, 3% annual)
- All 10 equation residuals at machine precision (~1e-15)

**Remaining:** Run in Dynare 6.5 / MATLAB R2025b to verify Blanchard-Kahn conditions and generate IRFs for comparison with Figures 8 and G2 of the paper.

---

## 2026-02-23: Nonlinearity visualisation — three-approach suite

Added five files to visualise the nonlinearities implied by the HDNWR model. The existing `hdnwr.mod` uses `order=1` (linear by construction), so the nonlinearities must be demonstrated via separate scripts/mod files.

### Why order=1 cannot show nonlinearities
First-order perturbation produces a linear policy function: all IRFs are symmetric and proportional to shock size. The nonlinearity is embedded in eqs (17) and (28), which define a convex `(u, piW)` locus parametrised by `jstar`. Visualising it requires either (i) tracing the static locus directly, or (ii) using order=2 perturbation.

### Files added

| File | Role |
|------|------|
| `plot_phillips_curve.m` | Approach 1: static Phillips curve (no Dynare) |
| `hdnwr_nl.mod` | Approaches 2 & 3: order=2 Dynare model |
| `hdnwr_nl_steadystate.m` | Steady-state companion for `hdnwr_nl.mod` (Dynare requires filename to match `.mod`) |
| `plot_asymmetric_irfs.m` | Approach 2: asymmetric ±shock IRFs via `simult_` |
| `plot_stochastic_scatter.m` | Approach 3: stochastic simulation scatter vs static PC |

---

### Approach 1 — Static Phillips curve (`plot_phillips_curve.m`)

Standalone MATLAB script, no Dynare needed. Traces eqs (17) and (28) parametrically over a grid of 2,000 `jstar` values to recover the exact nonlinear `(u, piW)` locus. The key algebraic step: dividing eq (17) by `(1+pistar)^(1-eta)` gives a clean expression with `F(jstar)` on the right, from which `piW = F(jstar)^{1/(1-eta)} * (1+pistar) - 1` directly. Steady-state `jstar_ss` is solved via `fzero` of the same normalised equation used in `hdnwr_nl_steadystate.m`. Overlays the linear tangent at SS for comparison. Outputs `fig_phillips_curve.pdf`.

---

### Approach 2 — Asymmetric IRFs (`plot_asymmetric_irfs.m`)

Runs `dynare hdnwr_nl`, then uses `simult_()` to compute IRFs for both `+1σ` and `−1σ` monetary policy shocks (`eps_mu`) from the deterministic steady state `oo_.dr.ys`. IRFs are computed as deviation from a zero-shock baseline (to correctly absorb any second-order constant correction). Outputs `fig_asymmetric_irfs.pdf`.

**Sign convention (important):** A positive `eps_mu` shock is **contractionary** (`mu↑ → ii↑ → u↑, piW↓, y↓`). A negative shock is **expansionary** (`u↓, piW↑, y↑`). Both IRFs are plotted with their **actual signs** (no normalisation). The dotted black line `irf_pos + irf_neg` is identically zero at order=1 and nonzero at order=2 — this is the asymmetry term. The expected direction: `irf_pos(u) + irf_neg(u) > 0` (unemployment rises more from a contractionary shock than it falls from an expansionary one), consistent with a convex Phillips curve.

**`simult_` signature in Dynare 6.x:**
```matlab
y_ = simult_(M_, options_, y0, oo_.dr, ex_, 2)
```
Note: older Dynare versions (≤ 5.x) used `simult_(y0, dr, ex_, iorder)` without `M_` and `options_`. Verify with `help simult_` if uncertain.

---

### Approach 3 — Stochastic simulation scatter (`plot_stochastic_scatter.m`)

Extracts `u_t` and `piW_t` from the 10,000-quarter stochastic simulation in `oo_.endo_simul` (size: `n_endo × (periods+2)`, columns 2:end-1 are the usable simulation). Scatter-plots simulated `(u_t, piW_t)` and overlays the theoretical static PC from Approach 1. The simulated cloud should track the nonlinear locus, confirming the second-order approximation respects the model's curvature. Outputs `fig_stochastic_scatter.pdf`.

---

### `hdnwr_nl.mod` vs `hdnwr.mod`

Equations, calibration, and steady-state logic are identical. Two changes:

1. **Shock calibration:** Changed from `var eps_mu = 0.01^2` (1% quarterly std dev = ~4% p.a.) to `stderr 0.0025` (25 bp quarterly = 1% p.a., a more conventional macro calibration). Technology shock: `stderr 0.01` (1% quarterly).

2. **`stoch_simul` options:** `order=1, irf=20` → `order=2, irf=20, periods=10000, pruning, nograph`. Pruning (Andreasen et al. 2013) prevents explosive paths in the second-order simulation.

---

### How to run

```
% Approach 1 (no Dynare):
run('plot_phillips_curve.m')

% Approaches 2 and 3 (run in same MATLAB session to avoid re-solving):
run('plot_asymmetric_irfs.m')    % calls dynare hdnwr_nl internally
run('plot_stochastic_scatter.m') % reuses M_, oo_ from above
```

All scripts should be run from within `sgu_model/`. Outputs are `fig_phillips_curve.pdf`, `fig_asymmetric_irfs.pdf`, `fig_stochastic_scatter.pdf`.
