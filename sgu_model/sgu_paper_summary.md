# Summary: Schmitt-Grohe and Uribe (2025) — "Heterogeneous Downward Nominal Wage Rigidity: Foundations of a Nonlinear Phillips Curve"

## 1. Main Contribution

The paper builds a model with **heterogeneous downward nominal wage rigidity (HDNWR)** across a continuum of labour varieties. The key innovation is that the degree of wage rigidity varies across workers/occupations, rather than being uniform. This delivers:

1. A **smooth, convex wage Phillips curve** — steep at high inflation, flat at low inflation.
2. A Phillips curve that is **contemporaneous** (no expected future inflation term), unlike the New Keynesian forward-looking Phillips curve.
3. A model that is **amenable to perturbation methods** despite featuring occasionally binding constraints at the micro level.

---

## 2. Model Setup

### 2.1 Firms

- Firms produce output using a CES composite of labour varieties $h_{jt}$ for $j \in [0,1]$, with elasticity of substitution $\eta > 0$:

$$h_t = \left[\int_0^1 h_{jt}^{1-1/\eta} dj\right]^{1/(1-1/\eta)}$$

- Cost minimisation yields the standard demand for variety $j$:

$$h_{jt} = \left(\frac{W_{jt}}{W_t}\right)^{-\eta} h_t$$

- The aggregate wage index is:

$$W_t = \left[\int_0^1 W_{jt}^{1-\eta} dj\right]^{1/(1-\eta)}$$

### 2.2 Households

- Representative household with standard preferences over consumption. Labour supply is **inelastic** in the baseline (endogenous supply is treated as an extension).
- Effective labour supply per variety is $\bar{h}(1 - u_t^n)$, where $u_t^n$ is the natural rate of unemployment (capturing non-nominal frictions).
- Employment is bounded: $h_{jt} \leq \bar{h}(1 - u_t^n)$.
- Standard Euler equation for consumption/bond holdings.

### 2.3 Heterogeneous Downward Nominal Wage Rigidity

This is the paper's key modelling innovation. Each variety $j$ faces a **variety-specific** wage floor:

$$W_{jt} \geq \gamma(j) W_{t-1}$$

where $\gamma(j)$ is **positive and increasing** in $j$. Crucially:

- The lower bound depends on the **past aggregate wage** $W_{t-1}$, not the past variety-specific wage $W_{jt-1}$. This is what enables clean aggregation.
- The function $\gamma(\cdot)$ captures heterogeneity in fairness standards, regulations, or productivity across workers — grounded in the empirical literature (Bewley, 1999; Fehr and Goette, 2005; Davis and Krolikowski, 2024).

The labour market closes with a **complementary slackness condition** for each variety:

$$[\bar{h}(1 - u_t^n) - h_{jt}][W_{jt} - \gamma(j) W_{t-1}] = 0$$

This says: if a variety has unemployment above the natural rate, its wage must be at its floor; if its wage is above the floor, it must be at full employment.

### 2.4 No Market Power

Unlike the New Keynesian tradition (Erceg, Henderson, and Levin, 2000), **workers do not have monopoly power**. Households and firms are price takers. This is motivated by the low unionisation rate in the U.S.

The consequence: the resulting Phillips curve is **contemporaneous**, not forward-looking. In the NK model, workers with market power make forward-looking wage-setting decisions, generating an expectations-augmented Phillips curve. Here, wages are constrained by floors, not set strategically.

---

## 3. The Cutoff Variety $j_t^*$ and Aggregation

### 3.1 How Individual Constraints Vanish in Aggregate

This is the central mechanism that makes perturbation feasible. The key insight:

**In equilibrium, there exists a cutoff variety $j_t^*$ that partitions all labour varieties into two groups:**

| | Varieties $j \leq j_t^*$ | Varieties $j > j_t^*$ |
|---|---|---|
| **Employment** | Full employment: $h_{jt} = \bar{h}(1 - u_t^n)$ | Involuntary unemployment: $h_{jt} < \bar{h}(1 - u_t^n)$ |
| **Wage** | Unconstrained (all paid $\gamma(j_t^*) W_{t-1}$) | Stuck at floor: $W_{jt} = \gamma(j) W_{t-1}$ |

**Why does this partition arise?**

- Consider variety $j_t^*$: it is at full employment AND its wage floor binds with equality.
- For any $j < j_t^*$: the wage floor $\gamma(j)W_{t-1}$ is lower than $\gamma(j_t^*)W_{t-1}$ (since $\gamma$ is increasing). If these varieties paid less than $\gamma(j_t^*)W_{t-1}$, their employment would exceed full employment (by the downward-sloping demand curve), which is impossible. So they must also pay $\gamma(j_t^*)W_{t-1}$ and be at full employment.
- For any $j > j_t^*$: the wage floor $\gamma(j)W_{t-1}$ is higher. The demand at these wages falls below full employment. By slackness, these varieties are stuck at their (binding) wage floors.

### 3.2 Aggregation into Two Smooth Equations

Given the cutoff structure, the aggregate equilibrium conditions reduce to two key equations that are **smooth functions of $j_t^*$** (no kinks, no occasionally binding constraints):

**Wage inflation equation** (from the wage aggregation):

$$(1 + \pi_t^W)^{1-\eta} = j_t^* \gamma(j_t^*)^{1-\eta} + \int_{j_t^*}^1 \gamma(j)^{1-\eta} dj \tag{17}$$

**Unemployment equation:**

$$u_t = u_t^n + (1 - u_t^n)\left[(1 - j_t^*) - \int_{j_t^*}^1 \left(\frac{\gamma(j)}{\gamma(j_t^*)}\right)^{-\eta} dj\right] \tag{18}$$

**These are the building blocks of the Phillips curve.** Together, they parametrically define a relationship $u_t = \mathcal{PC}(\pi_t^W; u_t^n)$ with $j_t^*$ as the parameter.

### 3.3 Why Perturbation Works

The critical point: **equations (17) and (18) are smooth, differentiable functions of $j_t^*$, $\pi_t^W$, and $u_t$.** The occasionally binding constraints at the individual variety level have been "integrated out" — the cutoff $j_t^*$ is a continuous variable, and the integrals over constrained/unconstrained varieties are smooth.

Compare this with the **homogeneous DNWR** model (Schmitt-Grohe and Uribe, 2016), where $\gamma(j) = \gamma$ for all $j$. In that case:
- Either ALL varieties are constrained or NONE are.
- The aggregate equilibrium features a single occasionally binding constraint: $W_t \geq \gamma W_{t-1}$.
- This produces an **L-shaped Phillips curve** (horizontal for positive unemployment, vertical at zero unemployment) — a non-differentiable kink.
- Perturbation around steady state is not possible because the equilibrium conditions are not differentiable.

With heterogeneity, by contrast, a **marginal change** in aggregate conditions shifts $j_t^*$ continuously. A small decrease in demand moves $j_t^*$ down slightly — a few more varieties become constrained, unemployment rises smoothly, and wage inflation falls smoothly. There is no discrete jump or kink.

**In short:** micro-level occasionally binding constraints $\to$ smooth macro-level equations through the continuous cutoff $j_t^*$.

---

## 4. Properties of the Wage Phillips Curve

### 4.1 Shape: Smooth and Convex

The Phillips curve is:
- **Downward-sloping:** higher $j_t^*$ (tighter labour market) $\to$ higher wage inflation (eq. 17) and lower unemployment (eq. 18).
- **Convex:** steep at high inflation/low unemployment, flat at low inflation/high unemployment.

**Intuition for convexity:**
- At **high inflation** (large $j_t^*$): most varieties are unconstrained. An additional increase in inflation only releases a *small* mass of workers from their constraints $\to$ small reduction in unemployment $\to$ steep curve.
- At **low inflation** (small $j_t^*$): many varieties are constrained. An increase in inflation releases a *large* mass from constraints $\to$ large reduction in unemployment $\to$ flat curve.

### 4.2 Quantitative Illustration

With the baseline calibration (linear $\gamma(j) = (1+\pi^*)(\Gamma_0 + \Gamma_1 j)$, $\Gamma_0 = 0.978$, $\Gamma_1 = 0.031$, $\eta = 11$, $u^n = 0.04$):

| Inflation change | Unemployment increase |
|---|---|
| 6% $\to$ 5% | 0.3 pp |
| 2% $\to$ 1% | 3.0 pp |

This asymmetry is large: reducing inflation is 10 times cheaper (in unemployment terms) at high inflation than at low inflation.

### 4.3 Shifters

- **Supply shocks ($u_t^n$):** shift the Phillips curve right/up. The natural rate is the only shifter; demand shocks, monetary shocks, and productivity shocks cause movements *along* the curve.
- **Long-run inflation ($\pi^*$):** higher target shifts the curve up/right (more grease needed).
- **Elasticity of substitution ($\eta$):** higher $\eta$ flattens the curve.

### 4.4 Pandemic-Era Decomposition

Using the model to back out supply shocks from observed $(u_t, \pi_t^W)$ pairs:
- **2020-2021:** Large negative supply shocks (consistent with lockdowns).
- **2022:** Supply shocks had dissipated $\to$ the inflation surge was primarily **demand-driven**.
- This aligns with Bergholt et al. (2024) and Giannone and Primiceri (2024).

---

## 5. Empirical Evidence from Regional Data

Using a panel of 13 large U.S. Combined Statistical Areas (CSAs), 2006:Q4–2025:Q2:

- The authors estimate: $\pi_{it}^W = \alpha_0 + \alpha_1 u_{it} + \alpha_2 u_{it}^{-1} + \text{controls} + \epsilon_{it}$
- The coefficient $\alpha_2$ on the reciprocal of unemployment (the nonlinear term) is **significantly positive** across all four specifications, rejecting a linear Phillips curve.
- The nonlinearity also holds in the **pre-Covid subsample** (2006:Q4–2019:Q4).
- The model, calibrated only to match position and slope (not curvature), captures the estimated nonlinearity well.

Key innovation in estimation strategy:
- No expected inflation term (consistent with the model's contemporaneous Phillips curve).
- Smooth nonlinear specification (not piecewise linear).
- ECI-based wage inflation (mitigates composition bias).

---

## 6. Extensions

### 6.1 Long-Run Phillips Curve

- Under **full indexation** ($\delta = 1$): the long-run Phillips curve is **vertical** — inflation cannot permanently grease the labour market.
- Under **partial indexation** ($\delta < 1$): the long-run curve is downward-sloping but steeper than the short-run curve.
- The NAIRU is the unemployment rate at which the short-run and long-run curves intersect at the target inflation rate.

### 6.2 Endogenous Labour Supply

- Adding elastic labour supply (disutility $V(h) = h^{1+\theta}/(1+\theta)$) barely changes the Phillips curve for empirically relevant labour supply elasticities ($1/\theta = 0.2$).
- The key aggregation result survives: the cutoff structure and smooth aggregate equations remain.

### 6.3 Heterogeneous Labour Productivity

- An alternative microfoundation: workers differ in productivity $z$, and the wage floor is $W_t(z) \geq z^\xi \gamma W_{t-1}$.
- Less productive workers are more likely to be constrained — consistent with evidence that low-skill employment is more procyclical.
- The resulting Phillips curve is again smooth and convex.

---

## 7. Local Equivalence with the New Keynesian Model

For **regular fluctuations** around the inflation target, the HDNWR model produces dynamics that are **quantitatively similar** to the standard NK model with Calvo wage staggering.

**Why?** The linearised HDNWR Phillips curve is $\hat{\pi}_t^W = \kappa_1 \hat{u}_t$, while the NK Phillips curve is $\hat{\pi}_t^W = \beta E_t \hat{\pi}_{t+1}^W + \kappa_2 \hat{u}_t$. Iterating the NK version forward gives $\hat{\pi}_t^W = \frac{\kappa_2}{1 - \lambda\beta} \hat{u}_t$, where $\lambda$ captures the endogenous persistence of unemployment. For standard calibrations, $\kappa_1 \approx \frac{\kappa_2}{1-\lambda\beta}$, so the two models generate similar impulse responses.

This means the HDNWR model:
- **Locally** (near steady state): behaves like the workhorse NK model — good for understanding the Great Moderation.
- **Globally** (far from steady state): delivers a nonlinear Phillips curve — good for understanding the pandemic inflation and the Great Recession.

---

## 8. Key Takeaways

1. **Heterogeneity in DNWR is the crucial innovation.** It transforms a kink (homogeneous DNWR) or an L-shaped curve into a smooth, convex curve.

2. **The cutoff $j_t^*$ is the bridge between micro and macro.** Individual varieties face occasionally binding constraints, but the aggregate is characterised by the continuous variable $j_t^*$, yielding smooth, differentiable equilibrium conditions.

3. **Perturbation works because:** the aggregate equilibrium conditions (17) and (18) are smooth functions — the integrals over the distribution of constrained/unconstrained workers are differentiable. No aggregate occasionally binding constraint appears.

4. **No market power needed.** The model delivers nominal rigidity with competitive markets, producing a contemporaneous (not forward-looking) Phillips curve.

5. **Policy implications:**
   - Disinflation from high inflation is cheap in unemployment terms (steep part of the curve).
   - Reducing high unemployment is cheap in inflation terms (flat part of the curve).
   - Both the "missing inflation" (post-2008) and "missing unemployment" (post-Covid disinflation) can be rationalised by a single nonlinear curve.

---

## 9. Comparison with Benigno and Eggertsson (2023)

| Feature | SGU (2025) | BE (2023) |
|---|---|---|
| **Nonlinearity** | Smooth, convex | Piecewise linear with a kink |
| **Labour market** | Competitive (price takers) | Search frictions |
| **DNWR** | Heterogeneous across varieties | Wages flexible when market tight, rigid when slack |
| **Phillips curve** | Contemporaneous | Forward-looking (expectations-augmented) |
| **Solution method** | Perturbation (standard) | Piecewise linear / global methods |
| **Captures flat part** | Yes (low inflation) | No (single kink only captures steepening) |
| **DNWR when market tight** | Yes (some varieties always constrained) | No (wages flexible when tight) |
