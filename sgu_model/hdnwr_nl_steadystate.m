function [ys, params, info] = hdnwr_nl_steadystate(ys, exo, M_, options_)
% Steady-state file for hdnwr_nl.mod
% Identical to hdnwr_steadystate.m — required because Dynare automatically
% looks for a file matching the .mod filename (hdnwr_nl_steadystate.m).
%
% Schmitt-Grohe and Uribe (2025) -- Endogenous Labor Supply version
%
% Inputs (Dynare convention):
%   ys       -- initial guess for steady-state vector
%   exo      -- exogenous steady-state values (unused)
%   M_       -- model structure
%   options_ -- options structure
%
% Outputs:
%   ys       -- steady-state vector (ordered as M_.endo_names)
%   params   -- (possibly updated) parameter vector
%   info     -- 0 if successful, >0 signals an error to Dynare

info   = 0;
params = M_.params;

% --- Unpack parameters ------------------------------------------------
idx_p  = @(name) find(strcmp(name, cellstr(M_.param_names)));

beta     = params(idx_p('beta'));
sig      = params(idx_p('sig'));
theta    = params(idx_p('theta'));
alpha    = params(idx_p('alpha'));
eta      = params(idx_p('eta'));
Gam0     = params(idx_p('Gam0'));
Gam1     = params(idx_p('Gam1'));
del      = params(idx_p('del'));
pistar   = params(idx_p('pistar'));
un       = params(idx_p('un'));

% --- Steady-state inflation -------------------------------------------
pi_ss  = pistar;
piW_ss = pistar;

% --- Solve eq (17) for jstar -----------------------------------------
% With piW = pistar and del = 1: gamma(j) = (1+pistar)*(Gam0+Gam1*j).
% Dividing (1+piW)^{1-eta} = jstar*gamjstar^{1-eta} + INT17
% by (1+pistar)^{1-eta} gives:
%   1 = j*(Gam0+Gam1*j)^{1-eta}
%       + [(Gam0+Gam1)^{2-eta} - (Gam0+Gam1*j)^{2-eta}] / [Gam1*(2-eta)]

f_eq17 = @(j) j*(Gam0 + Gam1*j)^(1-eta) ...
         + ((Gam0+Gam1)^(2-eta) - (Gam0+Gam1*j)^(2-eta)) / (Gam1*(2-eta)) ...
         - 1;

jstar_ss = fzero(f_eq17, 0.65);

if jstar_ss <= 0 || jstar_ss >= 1
    info = 1;   % signal failure to Dynare
    return;
end

% --- Inner gamma at cutoff -------------------------------------------
gstar = Gam0 + Gam1*jstar_ss;

% --- Aggregate hours from eq (27) in steady state --------------------
% (h/(1-un))^theta * gstar^{-eta*theta} * h^{alpha*sig}
%     = gstar * alpha * h^{alpha-1}
% => h^{theta + alpha*sig - alpha + 1} = (1-un)^theta * gstar^{1+eta*theta} * alpha

power_h = theta + alpha*sig - alpha + 1;
h_ss    = ((1-un)^theta * gstar^(1+eta*theta) * alpha)^(1/power_h);

% --- Remaining steady-state values -----------------------------------
a_ss  = 1;
mu_ss = 1;
y_ss  = a_ss * h_ss^alpha;
w_ss  = alpha * h_ss^(alpha-1);
ii_ss = (1+pistar)/beta - 1;

% --- Unemployment (eq 28) --------------------------------------------
ratio_ss = (Gam0+Gam1) / gstar;
INT_s_ss = gstar / (Gam1*(1/theta+1)) * (ratio_ss^(1/theta+1) - 1);
INT_d_ss = gstar / (Gam1*(1-eta))     * (ratio_ss^(1-eta)     - 1);
u_ss     = un + (1-un) * (INT_s_ss - INT_d_ss) / (jstar_ss + INT_s_ss);

% --- Pack into ys (order must match M_.endo_names) -------------------
idx_e = @(name) find(strcmp(name, cellstr(M_.endo_names)));

ys(idx_e('jstar')) = jstar_ss;
ys(idx_e('y'))     = y_ss;
ys(idx_e('h'))     = h_ss;
ys(idx_e('u'))     = u_ss;
ys(idx_e('w'))     = w_ss;
ys(idx_e('ii'))    = ii_ss;
ys(idx_e('pi'))    = pi_ss;
ys(idx_e('piW'))   = piW_ss;
ys(idx_e('mu'))    = mu_ss;
ys(idx_e('a'))     = a_ss;

end
