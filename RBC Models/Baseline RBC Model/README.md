---
author:
- David Murakami
date: 2023-04-12
title: "**A Baseline RBC Model**"
---

# Ramsey social planner problem {#sec: baseline RBC model .unnumbered}

Consider a Hansen-style RBC model.Final good firms use the production
function $Y_t = A_tK_{t-1}^\alpha L_t^{1-\alpha}$. The representative
household utility function is
$u(C_t,L_t) = \ln C_t - \chi L^{1+\nu^{-1}}/(1+\nu^{-1})$. Set
$\alpha=0.3$, $\nu=2$, $\chi=4.5$, $\beta=0.99$, $\delta=0.025$,
$\rho_a=0.95$, and $\sigma_a=0.01$.

As there are no distortions, we can solve the model from the perspective
of the Ramsey social planner who aims to
$$\max_{\{C_t,L_t,K_t\}} \, \mathbb{E}_t \sum_{i=0}^\infty \beta^i \left(\ln C_{t+i} - \chi \frac{L_{t+i}^{1+\frac{1}{\nu}}}{1+\frac{1}{\nu}} \right),$$
subject to $$\begin{aligned}
        \label{eq: aggregate resource constraint RBC}
        Y_t &= C_t + I_t, \\
        \label{eq: law of motion of capital RBC}
        K_t &= I_t + (1-\delta)K_{t-1}, \\
        \label{eq: production function RBC}
        Y_t &= A_tK_{t-1}^\alpha L_t^{1-\alpha}, \\
        \label{eq: TFP process RBC}
        \ln A_t &= \rho_a \ln A_{t-1} + \varepsilon_t^a, \quad \varepsilon_t^a \overset{IID}{\sim} \mathcal{N}(0,\sigma_a^2).
    
\end{aligned}$$ Solving the problem yields the following first-order
conditions: $$\begin{aligned}
        \frac{1}{C_t} &= \beta \mathbb{E}_t\left[\frac{1}{C_{t+1}}\left(\alpha\frac{Y_{t+1}}{K_t} + 1 - \delta\right)\right], \\
        \chi L_t^{\frac{1}{\nu}} &= \frac{(1-\alpha)Y_t}{C_tL_t}.
    
\end{aligned}$$

Define the marginal value of an additional unit of capital in period
$t+1$ as $$\label{eq: capital return RBC}
        R_{t+1} = \alpha\frac{Y_{t+1}}{K_t} + 1 - \delta,$$ and define
the marginal product of labour as the real wage, $$\label{eq: wages RBC}
        W_t = (1-\alpha)\frac{Y_t}{L_t}$$ Use these to write the FOCs
as: $$\begin{aligned}
        \label{eq: consumption euler equation RBC}
        \frac{1}{C_t} &= \beta \mathbb{E}_t\left[\frac{R_{t+1}}{C_{t+1}}\right], \\
        \label{eq: intratemporal euler equation RBC}
        W_t &= \chi L_t^\frac{1}{\nu}C_t.
    
\end{aligned}$$

#### Centralised equilibrium.

The equilibrium is a set of prices, $R_t$ and $W_t$; allocations, $Y_t$,
$C_t$, $L_t$, $I_t$, and $K_t$; and productivity, $A_t$, which satisfy
eight equations,
[\[eq: aggregate resource constraint RBC\]](#eq: aggregate resource constraint RBC){reference-type="eqref"
reference="eq: aggregate resource constraint RBC"}-[\[eq: intratemporal euler equation RBC\]](#eq: intratemporal euler equation RBC){reference-type="eqref"
reference="eq: intratemporal euler equation RBC"}.

#### Steady state.

Steady state quantities are found by first using the fact that $A=1$,
and from
[\[eq: capital return RBC\]](#eq: capital return RBC){reference-type="eqref"
reference="eq: capital return RBC"} and
[\[eq: consumption euler equation RBC\]](#eq: consumption euler equation RBC){reference-type="eqref"
reference="eq: consumption euler equation RBC"} we have: $$\begin{split}
            \frac{1}{\beta} &= \alpha\left(\frac{K}{L}\right)^{\alpha-1} + 1 - \delta \\
            \implies \frac{K}{L} &= \left(\frac{\alpha}{\beta^{-1}-1+\delta}\right)^{\frac{1}{1-\alpha}}.
        \end{split}$$ Then express $W$ in
[\[eq: wages RBC\]](#eq: wages RBC){reference-type="eqref"
reference="eq: wages RBC"} as:
$$W = (1-\alpha)\left(\frac{K}{L}\right)^\alpha.$$

From
[\[eq: law of motion of capital RBC\]](#eq: law of motion of capital RBC){reference-type="eqref"
reference="eq: law of motion of capital RBC"} we know: $$I = \delta K.$$
Then, combine this with
[\[eq: aggregate resource constraint RBC\]](#eq: aggregate resource constraint RBC){reference-type="eqref"
reference="eq: aggregate resource constraint RBC"} and
[\[eq: production function RBC\]](#eq: production function RBC){reference-type="eqref"
reference="eq: production function RBC"} to get the consumption-output
ratio:
$$\frac{C}{L} = \left(\frac{K}{L}\right)^\alpha - \delta\left(\frac{K}{L}\right).$$

Next, turn to the intratemporal Euler equation
[\[eq: intratemporal euler equation RBC\]](#eq: intratemporal euler equation RBC){reference-type="eqref"
reference="eq: intratemporal euler equation RBC"}, and multiply both
sides by $L$:
$$\chi L^{\frac{1+\nu}{\nu}} = \frac{(1-\alpha)L}{C}\left(\frac{K}{L}\right)^\alpha,$$
and then rearrange to get $C/L$:
$$\frac{C}{L} = \frac{(1-\alpha)\left(\frac{K}{L}\right)^\alpha}{\chi L^{\frac{1+\nu}{\nu}}}.$$
Set this equal to the previous expression for $C/L$ to solve for $L$:
$$L = \left\{\frac{(1-\alpha)\left(\frac{K}{L}\right)^\alpha}{\chi\left[\left(\frac{K}{L}\right)^\alpha - \delta\frac{K}{L}\right]}\right\}^{\frac{\nu}{1+\nu}}.$$

With $L$ in hand, turn back to the capital-labour ratio to obtain $K$:
$$K = \left(\frac{\alpha}{\beta^{-1}-1+\delta}\right)^{\frac{1}{1-\alpha}}L.$$
Thus, $I$ and $Y$ can easily be computed, after which $C$ can obtained
easily using the aggregate resource constraint
[\[eq: aggregate resource constraint RBC\]](#eq: aggregate resource constraint RBC){reference-type="eqref"
reference="eq: aggregate resource constraint RBC"}.
