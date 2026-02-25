# Benigno and Eggertsson (2025): Paper Summary

## 1. Central Thesis

The paper proposes a **nonlinear (Inverse-L shaped) New Keynesian Phillips Curve** to explain the 2020s inflation surge. The key insight is that inflation responds very differently to labor market conditions depending on whether the labor market is in a state of **labor shortage** ($\theta_t > \theta^*$, where $\theta = V/U$ is the vacancy-to-unemployment ratio) versus normal times ($\theta_t \leq \theta^*$).

- **Flat region** ($\theta \leq \theta^* \approx 1$): The Phillips curve is nearly flat. Tightening the labor market has minimal inflationary effects.
- **Steep region** ($\theta > \theta^*$): The Phillips curve becomes very steep. Small increases in tightness generate large inflation.

This shape resembles the "crude Keynesian" backward-L and Phillips' (1958) original nonlinear curve, but is embedded in a modern New Keynesian framework with forward-looking expectations.

---

## 2. How the Nonlinearity is Modelled

### 2.1 The Source: Asymmetric Wage Setting (The "Kink")

The nonlinearity comes from a **max operator** in the wage-setting equation for newly hired workers:

$$w_t^{new} = \max \left\{w_t^{ex}, w_t^{flex}\right\}$$

where:

- $w_t^{ex}$ = wage of workers in existing employment relationships (a "wage norm" that adjusts slowly)
- $w_t^{flex}$ = the flexible (market-clearing) wage, which is an increasing function of labor market tightness $\theta_t$

**Economic intuition (directly from Phillips 1958):**

- When unemployment is high ($\theta$ low), the flexible wage is below the existing wage. Workers **refuse to accept wages below the prevailing rate**, so $w_t^{new} = w_t^{ex}$. Wages fall only slowly.
- When the labor market is very tight ($\theta$ high), the flexible wage exceeds the existing wage. Employers **bid wages up rapidly** to attract scarce workers, so $w_t^{new} = w_t^{flex}$. Wages become fully flexible.

### 2.2 The Flexible Wage

The flexible wage is determined by the zero-profit condition of employment agencies:

$$w_t^{flex} = \frac{\gamma_t^c}{\gamma_t^b} \frac{1}{m_t} \theta_t^{\eta}$$

This is increasing and concave in $\theta_t$ (labor market tightness), derived from a standard search-and-matching framework where employment agencies post vacancies at cost $\gamma^c$ and earn a fee $\gamma^b$ proportional to the new hire's wage.

### 2.3 The Wage Norm (Existing Wages)

Existing wages evolve according to:

$$W_t^{ex} = \left(W_{t-1}^{ex} (\Pi_{t+1}^e)^{\delta}\right)^{\lambda} \left(P_t w_t^{flex}\right)^{1-\lambda} \phi_t$$

- $\lambda = 1$: wages are fully rigid at previous nominal level (Keynesian extreme)
- $\lambda = 0$: wages are fully flexible
- $0 < \lambda < 1$: existing wages are **gradually pulled toward flexible wages** but with inertia
- $\delta > 0$: inflation expectations feed into wage setting (wage-price spiral channel)

### 2.4 The Threshold $\theta^*$

The threshold is endogenously determined by the point where $w_t^{flex} = w_t^{ex}$:

$$\theta_t^* = \frac{\gamma_t^b}{\gamma_t^c} m_t \left(w_{t-1}^{ex} \frac{(\Pi_{t+1}^e)^{\delta}}{\Pi_t}\right)^{1/\eta} (\phi_t)^{1/(\lambda \eta)}$$

In principle $\theta^*$ varies over time, but empirically the authors approximate it as $\theta^* \approx 1$ (the "Beveridge threshold").

---

## 3. The Resulting Phillips Curve (Log-Linearized)

After log-linearizing the price-setting equation (Rotemberg adjustment costs) and substituting in the wage-setting equations, the Inv-L NK Phillips Curve is:

$$\pi_t = \begin{cases} -c + \kappa^{tight} \hat{\theta}_t + \kappa_v^{tight} \hat{v}*t + \beta E_t \pi*{t+1} & \text{if } \hat{\theta}_t > \hat{\theta}*t^* \text{ (labor shortage)}  \kappa_w \hat{w}*{t-1} + \kappa \hat{\theta}*t + \kappa_v \hat{v}t + \kappa\beta E_t \pi*{t+1} & \text{if } \hat{\theta}_t \leq \hat{\theta}_t^* \text{ (normal)} \end{cases}$$

**Three key theoretical predictions confirmed empirically:**

1. **$\kappa^{tight} > \kappa > 0$**: The slope w.r.t. tightness is much steeper during labor shortages. Empirically, the slope changes from ~0.52 to ~5.88 (column 4, Table 1).
2. **$\kappa_v^{tight} > \kappa_v$**: Supply shocks have a larger effect when the labor market is tight. This means supply shocks are "supercharged" during labor shortages.
3. **Lagged wages matter only in the normal regime**: During labor shortages, the Phillips curve is purely forward-looking (no $\hat{w}*{t-1}$ term). In normal times, wage persistence ($\kappa_w \hat{w}*{t-1}$) enters, creating inflation inertia.

---

## 4. Model Components

### Households

- GHH preferences (no wealth effects on labor supply)
- Extensive margin labor supply: household chooses how many members participate in labor force
- Members ordered by disutility of working (Galí 2009 device)

### Firms

- Continuum of monopolistically competitive firms (Dixit-Stiglitz)
- Production: $y_t(i) = A_t N_t(i)^\alpha O_t(i)^{1-\alpha}$ (labor + intermediate input, e.g. oil)
- Rotemberg price adjustment costs
- **Crucial feature**: firms first hire existing workers at $W^{ex}$, then new workers at $(1+\gamma^b) W^{new}$
- Marginal cost of production depends on **new hire wages**, not average wages

### Labor Market

- Search and matching: $M_t = m_t U_t^{\eta} V_t^{1-\eta}$
- Fraction $s$ of labor force must search; $(1-s)$ are attached to existing jobs
- Employment agencies post vacancies, charge firms a fee

### Closing

- Goods market clearing
- Taylor rule: $i_t = \rho \hat{r}*t^e + \phi*\pi (\pi_t - \pi^*) + e_t$

---

## 5. Empirical Specification

The benchmark regression:

$$\pi_t = \beta_c + \beta_\pi \pi_{t-1} + (\beta_\theta + \beta_{\theta_d} D_t) \ln \theta_t + (\beta_v + \beta_{v_d} D_t) \nu_t + \beta_{\pi^e} \pi_t^e + \varepsilon_t$$

where $D_t = 1$ if $\theta_t \geq 1$. Results (Table 1, 2008-2024):

- $\beta_\theta = 0.52$ (not significant alone), but $\beta_{\theta_d} = 5.36$ (significant at 1%)
- Combined slope when $\theta > 1$: about 5.88
- Supply shock coefficient doubles+ when $\theta > 1$

---

## 6. Policy Implications

1. **2020s inflation surge**: ~2/3 demand shocks, ~1/3 supply shocks. The nonlinearity amplified both. A "soft landing" is feasible because the steep Phillips curve means small output reductions can achieve large inflation reductions.
2. **1970s Great Inflation**: Primarily driven by unanchored inflation expectations + oil shocks. Labor market was *not* in shortage territory. The flat Phillips curve meant that reducing inflation required a deep recession (Volcker disinflation).
3. **2008 missing disinflation**: The flat Phillips curve when $\theta < 1$ explains why the large output drop post-2008 did not produce significant deflation.

---

## 7. Implementation in Dynare with OccBin: Feasibility Assessment

### The core question: Can this model be implemented in Dynare using OccBin?

**Yes, this is feasible.** The model is well-suited for OccBin because:

### Why OccBin is appropriate:

1. **Piecewise-linear structure**: The model is log-linearized and features exactly two regimes defined by whether $\hat{\theta}_t > \hat{\theta}_t^*$ or not. This is precisely what OccBin is designed for — it solves models with occasionally binding constraints that create piecewise-linear dynamics.
2. **The constraint is on a single variable**: The regime switch depends on $\theta_t$ relative to $\theta_t^*$. In OccBin, this can be formulated as an occasionally binding constraint on the wage-setting equation.
3. **Standard NK backbone**: The rest of the model (Euler equation, Taylor rule, goods market clearing) is entirely standard and easily coded in Dynare.

### How to implement:

The key is to formulate the max operator $w_t^{new} = \max \left\{w_t^{ex}, w_t^{flex}\right\}$ as an occasionally binding constraint in OccBin:

- **Reference regime** (normal times, $\theta \leq \theta^*$): New wages are constrained by existing wages. $w_t^{new} = w_t^{ex}$.
- **Alternative regime** (labor shortage, $\theta > \theta^*$): The constraint is slack, and $w_t^{new} = w_t^{flex}$.

In Dynare/OccBin syntax, this would be set up as:

- Define a slack variable (e.g., the gap between flexible and existing wages)
- The constraint binds when $w_t^{flex} \geq w_t^{ex}$ (or equivalently $\theta_t \geq \theta_t^*$)

### Potential challenges:

1. **Threshold determination**: $\theta_t^*$ is in principle time-varying (eq. 41). For a first pass, treating it as a constant ($\theta^* = 1$) as the authors do empirically is reasonable and simplifies implementation considerably.
2. **Two wage types**: The model distinguishes between $w^{new}$ and $w^{ex}$, which requires tracking both as state/jump variables.
3. **Search and matching block**: The matching function, employment agency problem, and Beveridge curve add equations but are all standard.
4. **Rotemberg vs Calvo**: The paper uses Rotemberg pricing, which is convenient for nonlinear work but equivalent to Calvo at first order. Either works in Dynare.
5. **IRFs and simulation**: OccBin naturally produces asymmetric impulse responses — shocks hitting during a labor shortage will produce different dynamics than the same shock in normal times. This is exactly the paper's key insight.

### Recommended implementation strategy:

1. Start with the simplified policy model (Section 5, equations 47-49 + the piecewise Phillips curve 48) — this is a 3-equation system.
2. Implement the full structural model later if needed.
3. Use OccBin's `occbin_setup` and constraint syntax.
4. Validate against the paper's AS-AD diagrams (Figures 9-11).

