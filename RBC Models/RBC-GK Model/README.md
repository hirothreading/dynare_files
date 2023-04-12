---
author:
- David Murakami
date: 2023-04-12
title: "**RBC Model with Financial Frictions**"
---

Note: You will need to download and place the s_find.m in the same directory as your Dynare mod file. It assists Dynare's steady state solver to compute the steady state.  

# Competitive equilibrium {#sec: RBC model with financial frictions .unnumbered}

Consider a simple Hansen-style RBC model but with financial frictions in
line with Gertler-Kiyotaki/Gertler-Karadi.

#### Households.

The representative household faces the following problem:
$$\max_{\{C_t,L_t,D_t,K_t^h\}} \, \mathbb{E}_t\sum_{i=0}^\infty\beta^i\left(\ln C_{t+i} - \chi\frac{L_{t+i}^{1+\frac{1}{\nu}}}{1+\frac{1}{\nu}}\right),$$
subject to their period budget constraint:
$$%       \label{eq: household period BC, RBC-GK}
        C_t + D_t + Q_tK_t^h + \chi_t^h = W_tL_t + (Z_t + \lambda Q_t)K_{t-1}^h + R_{t-1}D_{t-1} + \Pi_t^p,$$
where $\lambda=1-\delta$; $Q_t$ is the equity price in terms of final
goods; $\chi_t^h$ are portfolio management costs of the workers in the
household; $\Pi_t^p$ are profits earned by the household from the
production of intermediate goods, production of investment goods, and
banking; $Z_t$ is the net rental rate of capital; and $R_t$ is the gross
nominal interest rate on deposits. Denote $\sigma$ as the continuation
probability of banker maintaining its franchise, and $\gamma$ as the
fraction of total household assets given to new bank franchises as
initial funds. The efficiency cost of workers directly purchasing equity
is given by:
$$\chi_{t}^{h} = \frac{\varkappa^h}{2}\left(\frac{K_{t}^{h}}{K_{t}}\right)^2 K_{t}.$$

The FOCs for labour, savings in equity, and deposits which emerge from
the representative household problem, respectively, are:
$$\begin{aligned}
        \label{eq: intratemporal Euler equation, RBC-GK}
        W_t &= \chi L_t^{\frac{1}{\nu}}C_t, \\
        \label{eq: household FOC wrt equity}
        1 &= \mathbb{E}_{t}\left[ \Lambda_{t,t+1}\frac{Z_{t+1} + \lambda Q_{t+1}}{Q_{t} + \varkappa^{h}\frac{K_{t}^{h}}{K_{t}}}\right], \\
        \label{eq:household FOC wrt deposits}
        1 &= \mathbb{E}_{t}\left[ \Lambda_{t,t+1}R_t \right],
    
\end{aligned}$$ and where the household stochastic discount factor (SDF)
is defined as:
$$\Lambda_{t,\tau} = \beta^{\tau-t}\frac{C_t}{C_{\tau}}.$$

#### Firms.

There are two representative, perfectly competitive firms in the
economy: final good firms and investment good firms. Final good firms
produce according to the following constant returns to scale technology:
$$\label{eq: production, RBC-GK}
        Y_t = A_tK_{t-1}^\alpha L_t^{1-\alpha}.$$ Wages and return to
capital are then given by the following FOCs: $$\begin{aligned}
        \label{eq: wages, RBC-GK}
        W_t &= (1-\alpha)\frac{Y_t}{L_t}, \\
        \label{eq: rental rate, RBC-GK}
        Z_t &= \alpha\frac{Y_t}{K_{t-1}}.
    
\end{aligned}$$

Investment goods are produced by firms such that the aggregate capital
stock grows according to the following law of motion:
$$\label{eq: law of motion of capital, RBC-GK}
        K_t = \lambda K_{t-1} + I_t.$$ Total investment costs are given
by $$I_t\left[1 + \Phi\left(\frac{I_t}{I}\right)\right],$$ where
$\Phi(\cdot)$ are investment adjustment costs:
$$\Phi\left(\frac{I_t}{I}\right) = \frac{\kappa_I}{2}\left(\frac{I_t}{I} - 1\right)^2,$$
with $\Phi(1) = \Phi'(1) = 0$ and $\Phi''(\cdot) > 0$. Thus, the
representative investment good firm aims to maximises its profits:
$$\max_{I_t} \, \left\{Q_tI_t - I_t - \Phi\left(\frac{I_t}{I}\right)I_t\right\}.$$
Differentiating with respect to $I_t$ gives the following FOC:
$$\label{eq: equity price, RBC-GK}
        Q_t = 1 + \Phi\left(\frac{I_t}{I}\right) + \left(\frac{I_t}{I}\right)\Phi'\left(\frac{I_t}{I}\right).$$

#### Banks.

A banker will seek to maximise their franchise value, $\mathbb{V}_t^b$,
which is the expected present discount value of future dividends:
$$%\label{eq:franchise value of a bank}
        \mathbb{V}_t^b = \mathbb{E}_t\left[ \sum_{s=1}^\infty \Lambda_{t,t+s}\sigma^{s-1}(1-\sigma)n_{t+s} \right],$$
by choosing quantities $k_t^b$, $d_t$, and $d_t^*$.

A financial friction is used to limit a banker's ability to raise funds,
whereby a banker faces a moral hazard problem: they can either abscond
with the funds they have raised from depositors, or they can operate
honestly and pay out their obligations. Absconding is costly, however,
and so the banker can only divert a fraction $\theta > 0$ of assets they
have accumulated.

The caveat to absconding, in addition to only being able to take a
fraction of assets away, is that it takes time -- i.e. it take a full
period for the banker to abscond. Thus, the banker must decide to
abscond in period $t$, in addition to announcing what value of amount of
deposits they will choose and prior to realising next period's rental
rate of capital. If a banker chooses to abscond in period $t$, their
creditors will force the bank to shutdown in period $t+1$, causing the
banker's franchise value to become zero.

Therefore, the banker will choose to abscond in period $t$ if and only
if the return to absconding is greater than the franchise value of the
bank at the end of period $t$, $\mathbb{V}_t^b$. It is assumed that the
depositors act rationally, and that no rational depositor will supply
funds to the bank if they clearly have an incentive to abscond. In other
words, the bankers face the following incentive constraint:
$$%       \label{eq:banker incentive constraint}
        \mathbb{V}_t^b \geq \theta Q_tk_t^b,$$ where I assume that the
banker will not abscond in the case of the constraint holding with
equality.

The balance sheet constraint that the banker faces is:
$$%       \label{eq:banker balance sheet constraint}
        Q_tk_t^b = d_t + n_t,$$ the flow of funds constraint for a
banker is: $$%       \label{eq:banker flow of funds constraint}
        n_t = (Z_t + \lambda Q_t)k_{t-1}^b - R_{t-1}d_{t-1},$$ Note that
for the case of a new banker, the net worth is the startup fund given by
the household (fraction $\gamma$ of the household's assets):
$$n_t = \gamma\left(Z_t + \lambda Q_t\right)k_{t-1}.$$

With the constraints of the banker established, one can proceed to write
the banker's problem as: $$%\label{eq:banker's problem}
        \max_{k_t^b,d_t} \mathbb{V}_t^b = \mathbb{E}_t \left[ \Lambda_{t,t+1} \left\{ (1-\sigma)n_{t+1}+\sigma \mathbb{V}_{t+1}^b \right\} \right],$$
subject to the incentive constraint and the balance sheet constraint.

Since $\mathbb{V}_t^b$ is the franchise value of the bank, which can be
interpreted as a \"market value\", divide $\mathbb{V}_t^b$ by the bank's
net worth to obtain a Tobin's Q ratio for the bank denoted by $\psi_t$:
$$%       \label{eq:Tobin's Q for bank}
        \psi_t \equiv \frac{\mathbb{V}_t^b}{n_t} = \mathbb{E}_t \left[ \Lambda_{t,t+1}(1-\sigma+\sigma\psi_{t+1})\frac{n_{t+1}}{n_t} \right].$$

Define $\phi_t$ as the maximum feasible asset to net worth ratio, or,
rather, the leverage ratio of a bank:
$$%       \label{eq:banker leverage ratio}
        \phi_t = \frac{Q_t k_t^b}{n_t}.$$ Additionally, define
$\Omega_{t,t+1}$ as the stochastic discount factor of the banker,
$\mu_t$ as the excess return on capital over deposits, and $\upsilon_t$
as the marginal cost of deposits. The banker's problem can then be
written as the following:
$$%       \label{eq:banker's problem rewritten}
        \psi_t = \max_{\phi_t}\left\{ \mu_t\phi_t + \upsilon_t \right\},$$
subject to $$\psi_t \geq \theta\phi_t.$$ Solving this problem yields:
$$\begin{aligned}
        \label{eq:psi solution}
        \psi_t &= \theta\phi_t, \\
        \label{eq:phi solution Tobin Q}
        \phi_t &= \frac{\upsilon_t}{\theta - \mu_t},
    
\end{aligned}$$ where $$\begin{aligned}
        \label{eq:mu excess return on capital}
        \mu_t &= \mathbb{E}_t \left[ \Omega_{t,t+1} \left(\frac{Z_{t+1}+\lambda Q_{t+1}}{Q_t} - R_t \right) \right], \\
        \label{eq:upsilon marginal cost of deposits}
        \upsilon_t &= \mathbb{E}_t \left[ \Omega_{t,t+1}R_t \right],
    
\end{aligned}$$ with $$%       \label{eq:Omega SDF of banker}
        \Omega_{t,t+1} = \Lambda_{t,t+1}(1-\sigma+\sigma\psi_{t+1}).$$

#### Market clearing.

The aggregate resource constraint of the economy is:
$$\label{eq: aggregate resource constraint RBC-GK}
        Y_t = C_t + I_t + \chi_t^h.$$ Aggregate capital is the sum of
capital owned by workers and bankers: $$\label{eq: sum of capital}
        K_t = K_t^h + K_t^b.$$ Aggregate net worth of banks are given
by: $$\label{eq: aggregate net worth}
        N_t = \sigma\left[(Z_t + \lambda Q_t)K_{t-1}^b - R_{t-1}D_{t-1}\right] + \gamma(Z_t + \lambda Q_t)K_{t-1},$$
and the aggregate balance sheet of the banking sector is given by:
$$\begin{aligned}
        \label{eq: aggregate leverage rate}
        Q_tK_t^b &= \phi_t N_t, \\
        \label{eq: aggregate bank balance sheet}
        Q_tK_t^b &= D_t + N_t.
    
\end{aligned}$$

#### Competitive equilibrium.

A competitive equilibrium is a set of prices, $Q_t$, $R_t$, $W_t$, and
$Z_t$; quantities, $C_t$, $D_t$, $I_t$, $K_t$, $K_t^b$, $K_t^h$, $L_t$,
$N_t$, and $Y_t$; productivity, $A_t$; and bank variables, $\psi_t$,
$\phi_t$, $\mu_t$, and $\upsilon_t$, which satisfy 18 equations
[\[eq: intratemporal Euler equation, RBC-GK\]](#eq: intratemporal Euler equation, RBC-GK){reference-type="eqref"
reference="eq: intratemporal Euler equation, RBC-GK"}-[\[eq: aggregate bank balance sheet\]](#eq: aggregate bank balance sheet){reference-type="eqref"
reference="eq: aggregate bank balance sheet"} (including an additional
equation for the law of motion for $A_t$).

#### Steady state.

From
[\[eq: equity price, RBC-GK\]](#eq: equity price, RBC-GK){reference-type="eqref"
reference="eq: equity price, RBC-GK"} when $I_t = I$, we have $$Q = 1,$$
and from
[\[eq:household FOC wrt deposits\]](#eq:household FOC wrt deposits){reference-type="eqref"
reference="eq:household FOC wrt deposits"} we have:
$$R = \frac{1}{\beta}.$$

In equilibrium, the incentive compatibility constraint of the banker is
binding. Define the discounted spread between equity and deposits, $s$,
as: $$\label{eq: discounted spread}
        s = \beta(Z + \lambda) - 1,$$ which is considered to be
endogenous.

From the household's FOC wrt equity
[\[eq: household FOC wrt equity\]](#eq: household FOC wrt equity){reference-type="eqref"
reference="eq: household FOC wrt equity"} we have:
$$\frac{K^h}{K} = \frac{s}{\varkappa^h}.$$ Additionally, in steady
steady we have: $$\begin{aligned}
        \Omega &= \beta(1-\sigma+\sigma\psi), \\
        \upsilon &= \frac{\Omega}{\beta}, \\
        \mu &= \Omega\left(Z + \lambda - \frac{1}{\beta}\right).
    
\end{aligned}$$ So, using
[\[eq: discounted spread\]](#eq: discounted spread){reference-type="eqref"
reference="eq: discounted spread"}, we can write:
$$\frac{\mu}{\upsilon} = s.$$

Next, define $J$ as:
$$J = \frac{n_{t+1}}{n} = (Z + \lambda)\frac{K^b}{N} - \frac{D}{N}R,$$
and use $$\begin{aligned}
        \frac{D}{N} &= \phi - 1, \\
        \phi &= \frac{K^b}{N},
    
\end{aligned}$$ to write $J$ as $$\begin{split}
            J &= \left(Z + \lambda - R\right)\phi + R \\
            &= \frac{1}{\beta}\left(s\phi + 1\right).
        \end{split}$$

Then, from
[\[eq: aggregate net worth\]](#eq: aggregate net worth){reference-type="eqref"
reference="eq: aggregate net worth"} we have: $$\begin{aligned}
        N &= \sigma\left[(Z + \lambda)K^b - RD\right] + \gamma\left(Z+\lambda\right) \\
        \frac{N}{N} &= \sigma\left[(Z+\lambda)\frac{K^b}{N} - \frac{D}{N}R\right] + \frac{\gamma}{N}(Z + \lambda)K \\
        \beta &= \sigma\beta J + \frac{\gamma\beta}{N}(Z + \lambda)K \\
        & = \sigma\beta J + \frac{\gamma K^b}{N}\left(1 + \varkappa^h\frac{K^h}{K}\right)\frac{K}{K^b} \\
        &= \sigma\beta J + \gamma(1+s)\phi\frac{1}{\frac{K^b}{K}} \\
        &= \sigma\left(s\phi + 1\right) + \gamma(1+s)\phi\frac{1}{1-\frac{s}{\varkappa^h}} \\
        \beta &= \sigma + \left(\sigma s + \gamma\frac{1+s}{1-\frac{s}{\varkappa^h}}\right)\phi,
    
\end{aligned}$$ or
$$\phi = \frac{\beta - \sigma}{\sigma s + \gamma \frac{1+s}{1-\frac{s}{\varkappa^h}}}.$$

The Tobin's Q ratio for the bank,
$\psi_t \equiv \frac{\mathbb{V}_t^b}{n_t}$, in steady state is
$$\begin{split}
            \psi &= \beta(1-\sigma+\sigma\psi)J \\
            &= \frac{(1-\sigma)(s\phi + 1)}{1-\sigma-\sigma s\phi},
        \end{split}$$ and from
[\[eq:psi solution\]](#eq:psi solution){reference-type="eqref"
reference="eq:psi solution"} we have $\psi = \theta\phi$. Combine the
expressions for $\phi$ and $\psi$ to get
$$\frac{\theta(\beta-\sigma)}{\sigma s + \gamma\frac{1+s}{1-\frac{s}{\varkappa^h}}} = \frac{(1-\sigma)\left[\frac{s(\beta-\sigma)}{\sigma s + \gamma\frac{1+s}{1-\frac{s}{\varkappa^h}}} + 1\right]}{1-\sigma - \sigma\left[\frac{s(\beta-\sigma)}{\sigma s + \gamma\frac{1+s}{1-\frac{s}{\varkappa^h}}}\right]},$$
then rearrange: $$\begin{split}
            0 &= H(s,\gamma) \\
            &= (1-\sigma)\left[\beta s + \gamma\frac{1+s}{1-\frac{s}{\varkappa^h}}\right]\left[\sigma s + \gamma\frac{1+s}{1-\frac{s}{\varkappa^h}}\right] - \theta(\beta-\sigma)\left[\sigma(1-\beta)s + (1-\sigma)\gamma\frac{1+s}{1-\frac{s}{\varkappa^h}}\right].
        \end{split}$$ We can observe that as $\gamma \rightarrow 0$,
$$\begin{split}
            H(s,0) &= (1-\sigma)s^2\beta\sigma - \theta(\beta-\sigma)\left(\sigma(1-\beta)s\right) \\
            \implies s &\rightarrow \theta\frac{(\beta-\sigma)(1-\beta)}{(1-\sigma)\beta}.
        \end{split}$$ Thus, there exists a unique steady state
equilibrium with positive spread. Due to the constant returns to scale
property of the banker, bank variables depend only on parameters of the
banker ($s,\theta,\gamma,\beta,\sigma$) and not on the firm/supply side
parameters in steady state.

Given $s$, we can then get $$Z = \frac{1}{\beta}(1+s) - \lambda.$$ From
[\[eq: household FOC wrt equity\]](#eq: household FOC wrt equity){reference-type="eqref"
reference="eq: household FOC wrt equity"} we also have:
$$\frac{K^h}{K} = \frac{s}{\varkappa^h}.$$

We know that in this economy marginal cost is equal to unity, it then
follows that
$$1 = \frac{1}{A}\left(\frac{Z}{\alpha}\right)^\alpha\left(\frac{W}{1-\alpha}\right)^{1-\alpha}.$$
As we have $Z$, we then get $W$:
$$W = (1-\alpha)\left(\frac{\alpha}{Z}\right)^{\frac{\alpha}{1-\alpha}}.$$

The capital-output ratio is standard:
$$\frac{K}{Y} = \frac{\alpha}{Z}.$$ Then from
[\[eq: intratemporal Euler equation, RBC-GK\]](#eq: intratemporal Euler equation, RBC-GK){reference-type="eqref"
reference="eq: intratemporal Euler equation, RBC-GK"} and
[\[eq: wages, RBC-GK\]](#eq: wages, RBC-GK){reference-type="eqref"
reference="eq: wages, RBC-GK"}, rearrange to get $L$:
$$L = \left[\frac{(1-\alpha)}{\chi}\frac{Y}{C}\right]^{\frac{\nu}{1+\nu}},$$
and then combine with
[\[eq: aggregate resource constraint RBC-GK\]](#eq: aggregate resource constraint RBC-GK){reference-type="eqref"
reference="eq: aggregate resource constraint RBC-GK"} and
[\[eq: law of motion of capital, RBC-GK\]](#eq: law of motion of capital, RBC-GK){reference-type="eqref"
reference="eq: law of motion of capital, RBC-GK"} to write:
$$L = \left[\frac{(1-\alpha)}{\chi}\frac{1}{1-\delta\frac{K}{Y}-\frac{s^2}{2\varkappa^h}\frac{K}{Y}}\right]^{\frac{\nu}{1+\nu}}.$$
Then, get $K$ by dividing
[\[eq: rental rate, RBC-GK\]](#eq: rental rate, RBC-GK){reference-type="eqref"
reference="eq: rental rate, RBC-GK"} by
[\[eq: wages, RBC-GK\]](#eq: wages, RBC-GK){reference-type="eqref"
reference="eq: wages, RBC-GK"} and rearranging:
$$K = \frac{\alpha}{1-\alpha}\frac{WL}{Z}.$$ With $K$ and $L$, get $Y$
through
[\[eq: production, RBC-GK\]](#eq: production, RBC-GK){reference-type="eqref"
reference="eq: production, RBC-GK"} and then back out $C$ from
[\[eq: intratemporal Euler equation, RBC-GK\]](#eq: intratemporal Euler equation, RBC-GK){reference-type="eqref"
reference="eq: intratemporal Euler equation, RBC-GK"}.
