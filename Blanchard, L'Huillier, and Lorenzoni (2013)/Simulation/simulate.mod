% BLL (2013) model evaluated at the posterior means reported in Table 5.
% The equations and shock normalizations are shared with estimate.mod.

@#ifndef PAPER_SPECIFICATION
    @#define PAPER_SPECIFICATION = 0
@#endif

addpath('../Model');
addpath('../Estimation');
@#include "../Model/bll_declarations.inc"

rho = 0.9426;
rho_d = 0.4641;
rho_q = 0.0413;
rho_p = 0.7722;
rho_w = 0.9530;
rho_g = 0.9972;
sig_u = 1.1977;
sig_nu = 1.4738;
sig_d = 11.0982;
sig_q = 0.3500;
sig_p = 0.1778;
sig_w = 0.3057;
sig_g = 0.2877;
zet = 2.0871;
h = 0.5262;
del = 0.025;
bet = 0.99;
alp = 0.1859;
chi = 4.3311;
xi = 3.4919;
iot = 0;
iot_w = 0;
psi = 0.2;
cy_ratio = 1.28;
cal = 0.8770;
cal_w = 0.8690;
mu_p = 0.30;
psi_p = 0.4953;
mu_w = 0.05;
psi_w = 0.9683;
rho_r = 0.5583;
gam_pi = 1.0137;
gam_y = 0.0050;

BM_11 = 0; BM_12 = 0; BM_21 = 0; BM_22 = 0; BM_31 = 0; BM_32 = 0;
FA_11 = 0; FA_12 = 0; FA_13 = 0; FA_21 = 0; FA_22 = 0; FA_23 = 0;
CM_11 = 0; CM_12 = 0; CM_21 = 0; CM_22 = 0;

@#include "../Model/bll_model.inc"

steady;
check;
model_diagnostics;

shocks;
var e_1; stderr 1;
var e_2; stderr 1;
var eps_d; stderr 1;
var eps_q; stderr 1;
var eps_p; stderr 1;
var eps_w; stderr 1;
var eps_g; stderr 1;
end;

stoch_simul(order=1, irf=100, nograph, noprint);
