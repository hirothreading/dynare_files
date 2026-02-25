// =========================================================================
// Benigno and Eggertsson (2025): Simplified 3-Equation Policy Model
// with Nonlinear (Inverse-L) Phillips Curve
//
// Implementation using Dynare 6.5 OccBin toolkit
// Quarterly frequency
//
// Reference: Benigno, P. and Eggertsson, G. (2025)
// "It's Baaack: The Surge in Inflation"
// =========================================================================

// -------------------------------------------------------------------------
// Variable declarations
// -------------------------------------------------------------------------

var
    y       (long_name='Output gap')
    pi      (long_name='Inflation deviation from target')
    ii      (long_name='Nominal interest rate deviation')
    g       (long_name='Demand shock state')
    nu      (long_name='Supply/cost-push shock state')
    e       (long_name='Monetary policy shock state')
    chi     (long_name='Labor participation shock state')
;

varexo
    eps_d   (long_name='Demand shock innovation')
    eps_s   (long_name='Supply shock innovation')
    eps_m   (long_name='Monetary policy shock innovation')
    eps_chi (long_name='Labor participation shock innovation')
;

// -------------------------------------------------------------------------
// Parameters
// -------------------------------------------------------------------------

parameters
    // Structural parameters (Table A)
    sigma       $\sigma$        (long_name='Inverse intertemporal elasticity of substitution')
    betta       $\beta$         (long_name='Discount factor')
    phi_pi      $\phi_\pi$      (long_name='Taylor rule coefficient on inflation')

    // Phillips curve coefficients - slack regime (Table B, 2008-2023)
    kappa       $\bar{\kappa}$              (long_name='PC slope on output gap, slack regime')
    kappa_v     $\bar{\kappa}_\nu$          (long_name='PC slope on supply shock, slack regime')

    // Phillips curve coefficients - tight regime (Table B, 2008-2023)
    kappa_tight $\tilde{\kappa}^{tight}$    (long_name='PC slope on output gap, tight regime')
    kappa_v_tight $\tilde{\kappa}_\nu^{tight}$ (long_name='PC slope on supply shock, tight regime')

    // Regime switch parameters
    c_tilde     $\tilde{c}$     (long_name='Constant shift in tight regime')
    y_star      $\hat{Y}^*$     (long_name='Output gap threshold for regime switch')

    // Composite parameter
    alpha_omega $\alpha/\omega$  (long_name='Labor share over inverse Frisch elasticity')

    // Shock persistence
    rho_g       $\rho_g$        (long_name='Demand shock persistence')
    rho_nu      $\rho_\nu$      (long_name='Supply shock persistence')
    rho_e       $\rho_e$        (long_name='Monetary policy shock persistence')
    rho_chi     $\rho_\chi$     (long_name='Labor participation shock persistence')
;

// -------------------------------------------------------------------------
// Calibration
// -------------------------------------------------------------------------

// Structural parameters (Table A of BE 2025)
sigma       = 0.5;
betta       = 0.99;
phi_pi      = 1.5;

// Phillips curve coefficients (Table B, column 2: 2008-2023 sample)
kappa       = 0.0065;       // Slack regime slope
kappa_v     = 0.0093;       // Slack regime supply shock coefficient (abs value)
kappa_tight = 0.0736;       // Tight regime slope
kappa_v_tight = 0.2742;     // Tight regime supply shock coefficient

// Composite parameter: alpha/omega = 0.9/1.0
alpha_omega = 0.9;

// Regime switch threshold
// Y* is the output gap at which theta = theta* = 1
// Starting value; to be refined via structural model mapping
y_star      = 0.01;         // 1% output gap threshold

// Continuity condition: at y = y_star, both regimes give same inflation
// -c_tilde + kappa_tight * y_star = kappa * y_star
// => c_tilde = (kappa_tight - kappa) * y_star
c_tilde     = (kappa_tight - kappa) * y_star;

// Shock persistence (based on tau = 0.8 from the paper's Markov structure)
rho_g       = 0.8;
rho_nu      = 0.8;
rho_e       = 0.5;          // Monetary policy shocks less persistent
rho_chi     = 0.8;

// -------------------------------------------------------------------------
// Model equations
// -------------------------------------------------------------------------

model;

    // ----- Equation 1: IS Curve (Euler equation, eq. 47) -----
    // Y_t - G_t = E_t(Y_{t+1} - G_{t+1}) - sigma^{-1}(i_t - E_t(pi_{t+1}))
    // r_t^epsilon absorbed into demand shock g
    [name='IS curve']
    y - g = y(+1) - g(+1) - (1/sigma) * (ii - pi(+1));

    // ----- Equation 2: Phillips Curve (eq. 48) -----
    // Reference (slack) regime: flat Phillips curve
    // pi_t = kappa*(Y_t + alpha/omega * chi_t) + kappa_v * nu_t + beta * E_t(pi_{t+1})
    [name='Phillips curve', relax='tight_labor']
    pi = kappa * (y + alpha_omega * chi) + kappa_v * nu + betta * pi(+1);

    // Alternative (tight) regime: steep Phillips curve
    // pi_t = -c_tilde + kappa_tight*(Y_t + alpha/omega * chi_t) + kappa_v_tight * nu_t + beta * E_t(pi_{t+1})
    [name='Phillips curve', bind='tight_labor']
    pi = -c_tilde + kappa_tight * (y + alpha_omega * chi) + kappa_v_tight * nu + betta * pi(+1);

    // ----- Equation 3: Taylor Rule (eq. 49) -----
    // i_t = phi_pi * (pi_t - pi*) + e_t
    // Since pi is already deviation from target, and rho*r_t^e absorbed into e
    [name='Taylor rule']
    ii = phi_pi * pi + e;

    // ----- Shock processes -----
    [name='Demand shock']
    g = rho_g * g(-1) + eps_d;

    [name='Supply shock']
    nu = rho_nu * nu(-1) + eps_s;

    [name='Monetary policy shock']
    e = rho_e * e(-1) + eps_m;

    [name='Labor participation shock']
    chi = rho_chi * chi(-1) + eps_chi;

end;

// -------------------------------------------------------------------------
// OccBin constraint: regime switch at Y* threshold
// -------------------------------------------------------------------------
// When y > y_star, the labor market is tight (theta > theta* approx 1)
// and the steep Phillips curve applies (wages become flexible)

occbin_constraints;
    name 'tight_labor';
    bind y > y_star;
    relax y <= y_star;
end;

// -------------------------------------------------------------------------
// Steady state (all zeros - model is in deviations)
// -------------------------------------------------------------------------

steady_state_model;
    y   = 0;
    pi  = 0;
    ii  = 0;
    g   = 0;
    nu  = 0;
    e   = 0;
    chi = 0;
end;

// -------------------------------------------------------------------------
// Check steady state and Blanchard-Kahn conditions
// -------------------------------------------------------------------------

steady;
check;

// -------------------------------------------------------------------------
// OccBin deterministic scenarios using occbin_solver
// Uses shocks(surprise, overwrite) blocks with var/periods/values syntax.
// Each scenario overwrites the previous shock specification.
// -------------------------------------------------------------------------

// Run non-binding scenarios first (before regime-switching contaminates
// persistent variables inside occbin.solver)

// --- Scenario B: Supply shock (no regime switch expected) ---
shocks(surprise, overwrite);
    var eps_s;
    periods 1;
    values 0.01;
end;

occbin_setup(simul_periods=20, simul_check_ahead_periods=20);
occbin_solver;

verbatim;
    results_supply.piecewise = oo_.occbin.simul.piecewise;
    results_supply.linear    = oo_.occbin.simul.linear;
    fprintf('\n=== Scenario B (supply shock) completed ===\n');
end;

// --- Scenario D: Negative demand shock (no regime switch expected) ---
shocks(surprise, overwrite);
    var eps_d;
    periods 1;
    values -0.05;
end;

occbin_setup(simul_periods=20, simul_check_ahead_periods=20);
occbin_solver;

verbatim;
    results_neg_demand.piecewise = oo_.occbin.simul.piecewise;
    results_neg_demand.linear    = oo_.occbin.simul.linear;
    fprintf('\n=== Scenario D (negative demand shock) completed ===\n');
end;

// --- Scenario A: Large positive demand shock (triggers tight regime) ---
shocks(surprise, overwrite);
    var eps_d;
    periods 1;
    values 0.05;
end;

occbin_setup(simul_periods=20, simul_check_ahead_periods=20);
occbin_solver;

verbatim;
    results_demand.piecewise = oo_.occbin.simul.piecewise;
    results_demand.linear    = oo_.occbin.simul.linear;
    fprintf('\n=== Scenario A (demand shock) completed ===\n');
    fprintf('y(1) = %.6f, pi(1) = %.6f\n', ...
            results_demand.piecewise(1,1), results_demand.piecewise(1,2));
end;

// --- Scenario C: Combined demand + supply (triggers tight regime) ---
// NOTE: Supply shock reduced to 0.001 to avoid OccBin convergence failure.
// The 30x discontinuity in kappa_v at the regime boundary (0.0093 vs 0.2742)
// causes regime oscillation with larger supply shocks. The AS-AD static
// analysis captures the full supply shock amplification without this issue.
shocks(surprise, overwrite);
    var eps_d;
    periods 1;
    values 0.05;
    var eps_s;
    periods 1;
    values 0.001;
end;

occbin_setup(simul_periods=20, simul_check_ahead_periods=20);
occbin_solver;

verbatim;
    results_combined.piecewise = oo_.occbin.simul.piecewise;
    results_combined.linear    = oo_.occbin.simul.linear;
    fprintf('\n=== Scenario C (combined shocks) completed ===\n');
end;
