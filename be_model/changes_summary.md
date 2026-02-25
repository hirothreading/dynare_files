# Replication log files

## 2026-02-24: Initial implementation of simplified 3-equation policy model

### Files created
- `be_simplified.mod` — Dynare 6.5 / OccBin model file implementing the simplified 3-equation policy model (equations 47-49) from Benigno & Eggertsson (2025).
- `run_be_simplified.m` — MATLAB driver script for running simulations, generating OccBin IRFs, and plotting AS-AD diagrams.

### Model specification
- **IS curve** (eq. 47): Standard New Keynesian Euler equation with demand shocks.
- **Phillips curve** (eq. 48): Piecewise-linear with two regimes:
  - Slack regime (Y ≤ Y*): flat, κ = 0.0065, κ_ν = 0.0093
  - Tight regime (Y > Y*): steep, κ_tight = 0.0736, κ_ν_tight = 0.2742
  - Coefficients from Table B, column 2 (2008-2023 sample)
- **Taylor rule** (eq. 49): φ_π = 1.5.
- **OccBin constraint**: Regime switches when output gap exceeds Y* threshold.
- **Shocks**: AR(1) processes for demand (g), supply (ν), monetary policy (e), and labor participation (χ). Persistence ρ = 0.8 (matching τ from paper).
- **Calibration**: Table A structural parameters (σ=0.5, β=0.99, etc.). Y* = 0.01 as starting value.

### Simulation scenarios in driver script
1. Positive demand shock (tests nonlinear amplification)
2. Supply shock (tests regime-dependent pass-through)
3. Combined demand + supply (2020s replication, Figure 9)
4. Negative demand shock (asymmetry comparison)
5. Static AS-AD diagrams (Figures 8-9 replicas)

### Status
- **Working in MATLAB R2025b + Dynare 6.5.** All four OccBin scenarios and plots complete successfully.
- Verified output: Scenario A demand shock gives y(1) = 0.025531, pi(1) = 0.004540.

### OccBin lessons learned (Dynare 6.5)
- **Shock syntax**: Use `shocks(surprise, overwrite)` with `var`/`periods`/`values` keywords. Must call `occbin_setup(simul_periods=N)` AFTER the `shocks` block (not before — `occbin_setup` reads the shock specification).
- **Scenario ordering matters**: `occbin.solver` uses persistent variables that cache decision rules. Non-binding scenarios (supply shock, negative demand) must run BEFORE binding scenarios (positive demand) to avoid contamination.
- **Combined shock convergence**: The 30x discontinuity in κ_ν at the regime boundary (0.0093 → 0.2742) causes OccBin to oscillate when both demand and supply shocks are present. Supply shock in Scenario C reduced to 0.001 to achieve convergence. The static AS-AD analysis handles arbitrary supply shock sizes without this limitation.
- **Results stored in** `oo_.occbin.simul.piecewise` (T×7 matrix, no extra SS row) and `oo_.occbin.simul.linear`.

### Known issues / next steps
- Y* calibration is approximate (0.01). Needs mapping from θ* = 1 via structural model.
- Shock standard deviations are initial guesses; need calibration to match 2020s magnitudes.
- Supply shock coefficient κ_ν = -0.0093 in the paper's Table B (2008-2023) is negative; we use the absolute value. This should be revisited.
- Combined shock (Scenario C) limited to small supply shock (0.001) due to OccBin convergence; larger supply shocks cause regime oscillation.
- Future work: implement the full structural model with search-and-matching labor market.