% Bayesian estimation for Blanchard, L'Huillier, and Lorenzoni (2013).
% Modern scalar port for Dynare 7.1. The authors' pre-Dynare 4.5 files used
% matrix expressions inside the model block; those matrices are now solved
% in bll_signal_extraction.m and exposed to Dynare as scalar parameters.

@#ifndef RUN_ESTIMATION
    @#define RUN_ESTIMATION = 1
@#endif
@#ifndef MODE_COMPUTE
    @#define MODE_COMPUTE = 6
@#endif
@#ifndef MH_REPLIC
    @#define MH_REPLIC = 0
@#endif
@#ifndef MH_JSCALE
    @#define MH_JSCALE = 0.01
@#endif
@#ifndef USE_MODE_FILE
    @#define USE_MODE_FILE = 0
@#endif
@#ifndef RUN_MODE_CHECK
    @#define RUN_MODE_CHECK = 0
@#endif
@#ifndef COMPUTE_POSTERIOR_IRFS
    @#define COMPUTE_POSTERIOR_IRFS = 1
@#endif
@#ifndef PAPER_SPECIFICATION
    @#define PAPER_SPECIFICATION = 0
@#endif

addpath('../Model');
@#include "../Model/bll_declarations.inc"

% Calibrated values and prior means. Estimated values are never assigned in
% the steady-state file, so Dynare's candidate values remain authoritative.
rho = 0.6;
rho_d = 0.6;
rho_q = 0.4;
rho_p = 0.6;
rho_w = 0.6;
rho_g = 0.6;
sig_u = 0.5;
sig_nu = 1;
sig_d = 0.15;
sig_q = 0.15;
sig_p = 0.15;
sig_w = 0.15;
sig_g = 0.5;
zet = 2;
h = 0.5;
del = 0.025;
bet = 0.99;
alp = 0.3;
chi = 4;
xi = 5;
iot = 0;
iot_w = 0;
psi = 0.2;
cy_ratio = 1.28;
cal = 0.66;
cal_w = 0.66;
mu_p = 0.30;
psi_p = 0.5;
mu_w = 0.05;
psi_w = 0.5;
rho_r = 0.5;
gam_pi = 1.5;
gam_y = 0.005;

% Initial placeholders for coefficients updated by estimate_steadystate.m.
BM_11 = 0; BM_12 = 0; BM_21 = 0; BM_22 = 0; BM_31 = 0; BM_32 = 0;
FA_11 = 0; FA_12 = 0; FA_13 = 0; FA_21 = 0; FA_22 = 0; FA_23 = 0;
CM_11 = 0; CM_12 = 0; CM_21 = 0; CM_22 = 0;

@#include "../Model/bll_model.inc"

steady;
check;

shocks;
var e_1; stderr 1;
var e_2; stderr 1;
var eps_d; stderr 1;
var eps_q; stderr 1;
var eps_p; stderr 1;
var eps_w; stderr 1;
var eps_g; stderr 1;
end;

estimated_params;
rho,     0.9203, 0.20, 0.9999, BETA_PDF,      0.6,   0.2;
rho_d,   0.3468, 0.01, 0.99,   BETA_PDF,      0.6,   0.2;
rho_q,   0.1500, 0.01, 0.99,   BETA_PDF,      0.4,   0.2;
rho_p,   0.9288, 0.01, 0.99,   BETA_PDF,      0.6,   0.2;
rho_w,   0.4815, 0.01, 0.99,   BETA_PDF,      0.6,   0.2;
rho_g,   0.9933, 0.01, 0.9999, BETA_PDF,      0.6,   0.2;

sig_u,   1.1588, 0.01, 10,     INV_GAMMA_PDF, 0.5,   1;
sig_nu,  0.9938, 0.01, 10,     INV_GAMMA_PDF, 1,     1;
@#if PAPER_SPECIFICATION == 1
sig_d,  15.4434, 0.01, 30,     INV_GAMMA_PDF, 0.15,  1.5;
@#else
sig_d,  15.4434, 0.01, 30,     INV_GAMMA_PDF, 5,     1.5;
@#endif
sig_q,   0.3821, 0.01, 10,     INV_GAMMA_PDF, 0.15,  1;
sig_p,   0.1878, 0.01, 10,     INV_GAMMA_PDF, 0.15,  1;
sig_w,   0.3593, 0.01, 10,     INV_GAMMA_PDF, 0.15,  1;
sig_g,   0.2988, 0.01, 10,     INV_GAMMA_PDF, 0.5,   1;

h,       0.6105, 0,    0.999,  BETA_PDF,      0.5,   0.1;
alp,     0.1783, 0,    0.5,    NORMAL_PDF,    0.3,   0.05;
zet,     4.0000, 1.10, 6,      GAMMA_PDF,     2,     0.75;
@#if PAPER_SPECIFICATION == 1
chi,     5.5000, 1.10, 15,     GAMMA_PDF,     4,     1;
xi,      2.0000, 0.01, 10,     NORMAL_PDF,    5,     1;
@#else
chi,     5.5000, 1.10, 15,     NORMAL_PDF,    4,     1;
xi,      2.0000, 0.01, 10,     GAMMA_PDF,     5,     1;
@#endif
cal,     0.8500, 0.001,0.9999, BETA_PDF,      0.66,  0.1;
cal_w,   0.8825, 0.001,0.9999, BETA_PDF,      0.66,  0.1;
psi_p,   0.7700, 0,    0.99,   BETA_PDF,      0.5,   0.2;
psi_w,   0.9100, 0,    0.99,   BETA_PDF,      0.5,   0.2;
rho_r,   0.4861, 0.10, 0.99,   BETA_PDF,      0.5,   0.2;
gam_pi,  1.5000, 1.0001,2,     NORMAL_PDF,    1.5,   0.3;
@#if PAPER_SPECIFICATION == 1
gam_y,   0.0044, 0.0001,1,     NORMAL_PDF,    0.005, 0.05;
@#else
gam_y,   0.0044, 0.0001,1,     NORMAL_PDF,    0.005, 0.005;
@#endif
end;

varobs dca diad dn r pi dwa dya;

@#if RUN_ESTIMATION == 1
    @#if MH_REPLIC > 0
estimation(
    @#if COMPUTE_POSTERIOR_IRFS == 1
    bayesian_irf,
    consider_all_endogenous,
    @#endif
    datafile='data.mat',
    graph_format=pdf,
    irf=20,
    lik_init=2,
    mh_replic=@{MH_REPLIC},
    mh_nblocks=2,
    mh_jscale=@{MH_JSCALE},
    mh_drop=0.2,
    @#if RUN_MODE_CHECK == 1
    mode_check,
    @#endif
    mode_compute=@{MODE_COMPUTE},
    @#if USE_MODE_FILE == 1
    mode_file='mode/estimate_mode.mat',
    @#endif
    nodisplay,
    optim=('MaxIter',200),
    posterior_max_subsample_draws=1200,
    posterior_sampling_method='random_walk_metropolis_hastings',
    prefilter=1,
    presample=4,
    silent_optimizer,
    use_penalized_objective_for_hessian);
    @#else
estimation(datafile='data.mat',
    lik_init=2,
    @#if RUN_MODE_CHECK == 1
    mode_check,
    @#endif
    mode_compute=@{MODE_COMPUTE},
    @#if USE_MODE_FILE == 1
    mode_file='mode/estimate_mode.mat',
    @#endif
    nodisplay,
    optim=('MaxIter',200),
    prefilter=1,
    presample=4,
    silent_optimizer,
    use_penalized_objective_for_hessian);
    @#endif
save('estimate_results.mat', 'M_', 'oo_', 'options_', 'estim_params_');
@#else
model_diagnostics;
stoch_simul(order=1, irf=0, nograph, noprint);
@#endif
