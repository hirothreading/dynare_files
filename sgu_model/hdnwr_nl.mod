%--------------------------------------------------------------------------
% Heterogeneous Downward Nominal Wage Rigidity (HDNWR) Model
% Schmitt-Grohe and Uribe (2025)
% Endogenous Labor Supply Version (Definition C1, Appendix C)
%
% Second-order perturbation version — used for nonlinearity visualisation:
%   Approach 2: Asymmetric IRFs (±monetary shock)
%   Approach 3: Stochastic simulation scatter vs static Phillips curve
%
% Differs from hdnwr.mod only in the stoch_simul options:
%   order=2, irf=20, periods=10000, pruning
%
% Dynare 6.5 / MATLAB R2025b
%--------------------------------------------------------------------------

var jstar   $j^{*}$        (long_name='cutoff variety')
    y       $y$            (long_name='output')
    h       $h$            (long_name='aggregate labor')
    u       $u$            (long_name='unemployment rate')
    w       $w$            (long_name='real wage')
    ii      $i$            (long_name='nominal interest rate')
    pi      $\pi$          (long_name='price inflation')
    piW     $\pi^{W}$      (long_name='wage inflation')
    mu      $\mu$          (long_name='monetary policy shock')
    a       $a$            (long_name='technology shock');

varexo eps_mu $\varepsilon^{\mu}$ (long_name='monetary policy innovation')
       eps_a  $\varepsilon^{a}$   (long_name='technology innovation');

parameters beta sig theta alpha eta Gam0 Gam1 del pistar un
           alpha_pi alpha_y rho_mu rho_a;

%--------------------------------------------------------------------------
% Calibration (Table 1)
%--------------------------------------------------------------------------
beta     = 0.99;                % Discount factor
sig      = 1;                   % CRRA parameter (log utility)
theta    = 5;                   % Inverse Frisch elasticity
alpha    = 0.75;                % Labor elasticity of output
eta      = 11;                  % Elasticity of substitution across labor varieties
Gam0     = 0.978;               % Wage lower bound: intercept
Gam1     = 0.031;               % Wage lower bound: slope
del      = 1;                   % Wage indexation to steady-state inflation
pistar   = 1.03^(1/4) - 1;      % Quarterly inflation target (3% annual)
un       = 0.04;                % Natural rate of unemployment
alpha_pi = 1.5;                 % Taylor rule: inflation response
alpha_y  = 0.125;               % Taylor rule: output gap response (= 0.5/4)
rho_mu   = 0.5;                 % Monetary shock persistence
rho_a    = 0.9;                 % Technology shock persistence

%--------------------------------------------------------------------------
% Model
%--------------------------------------------------------------------------
model;

  %-- Model-local variables: closed-form integrals --%

  % gamma(j) = (1+pistar)^del * (Gam0 + Gam1*j)
  # gamjstar = (1+pistar)^del * (Gam0 + Gam1*jstar);

  % Inner linear terms
  # g1  = Gam0 + Gam1;          % at j = 1
  # gjs = Gam0 + Gam1*jstar;    % at j = j*

  % Wage aggregation integral (eq 17):
  %   INT17 = integral_{j*}^{1} gamma(j)^{1-eta} dj
  %         = (1+pistar)^{del*(1-eta)} / [Gam1*(2-eta)] * [g1^{2-eta} - gjs^{2-eta}]
  # INT17 = (1+pistar)^(del*(1-eta)) / (Gam1*(2-eta))
            * (g1^(2-eta) - gjs^(2-eta));

  % Ratio for unemployment integrals
  # ratio = g1 / gjs;

  % Labor supply integral (eq 28):
  %   INT_s = integral_{j*}^{1} (gamma(j)/gamma(j*))^{1/theta} dj
  %         = gjs / [Gam1*(1/theta+1)] * [ratio^{1/theta+1} - 1]
  # INT_s = gjs / (Gam1*(1/theta+1)) * (ratio^(1/theta+1) - 1);

  % Labor demand integral (eq 28):
  %   INT_d = integral_{j*}^{1} (gamma(j)/gamma(j*))^{-eta} dj
  %         = gjs / [Gam1*(1-eta)] * [ratio^{1-eta} - 1]
  # INT_d = gjs / (Gam1*(1-eta)) * (ratio^(1-eta) - 1);

  %-- Structural Equations --%

  % (B1) Production function
  y = a * h^alpha;

  % (B2) Euler equation
  y^(-sig) = beta * (1+ii) * y(+1)^(-sig) / (1+pi(+1));

  % (B3) Labor demand (MPL = real wage)
  a * alpha * h^(alpha-1) = w;

  % (B4) Taylor rule
  1 + ii = (1+pistar)/beta * ((1+pi)/(1+pistar))^alpha_pi
           * (y/steady_state(y))^alpha_y * mu;

  % (B5) Wage-price inflation identity
  1 + piW = (w/w(-1)) * (1+pi);

  % (27) Cutoff variety condition (endogenous labor supply)
  %   V'(h_j*) / U'(y) = gamma(j*) * w(-1) / (1+pi)
  %   where h_j* = [h/(1-un)] * [gamma(j*)/(1+piW)]^{-eta}
  (h/(1-un) * (gamjstar/(1+piW))^(-eta))^theta * y^sig
      = gamjstar * w(-1) / (1+pi);

  % (17) Wage aggregation
  (1+piW)^(1-eta) = jstar * gamjstar^(1-eta) + INT17;

  % (28) Unemployment (endogenous labor supply version)
  u = un + (1-un) * (INT_s - INT_d) / (jstar + INT_s);

  %-- Shock Processes --%

  % (32) Monetary policy shock
  log(mu) = rho_mu * log(mu(-1)) + eps_mu;

  % Technology shock
  log(a) = rho_a * log(a(-1)) + eps_a;

end;

%--------------------------------------------------------------------------
% Steady State  (computed in hdnwr_nl_steadystate.m)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Shocks
%--------------------------------------------------------------------------
shocks;
  var eps_mu; stderr 0.0025;    % 1% p.a. std dev monetary shock
  var eps_a;  stderr 0.01;      % 0.5% technology shock
end;

%--------------------------------------------------------------------------
% Computation
%--------------------------------------------------------------------------

% Verify steady state
steady;

% Check Blanchard-Kahn conditions
check;

% Second-order perturbation
%   irf=20    : IRFs for 20 quarters (positive shocks stored in oo_.irfs)
%   periods=10000 : long stochastic simulation for scatter plot (Approach 3)
%   pruning   : Andreasen et al. (2013) pruning to prevent explosive paths
stoch_simul(order=2, irf=20, periods=10000, pruning, nograph); %y h u w pi piW ii jstar;

%--------------------------------------------------------------------------
% Steady-State Verification
%--------------------------------------------------------------------------
fprintf('\n=== Steady-State Verification ===\n');
fprintf('jstar   = %.4f  (expected ~0.65)\n', oo_.steady_state(strmatch('jstar',M_.endo_names,'exact')));
fprintf('u       = %.4f  (expected ~0.06)\n', oo_.steady_state(strmatch('u',M_.endo_names,'exact')));
fprintf('h       = %.4f\n', oo_.steady_state(strmatch('h',M_.endo_names,'exact')));
fprintf('y       = %.4f\n', oo_.steady_state(strmatch('y',M_.endo_names,'exact')));
fprintf('w       = %.4f\n', oo_.steady_state(strmatch('w',M_.endo_names,'exact')));
fprintf('pi      = %.6f  (pistar = %.6f)\n', oo_.steady_state(strmatch('pi',M_.endo_names,'exact')), pistar);
fprintf('piW     = %.6f  (pistar = %.6f)\n', oo_.steady_state(strmatch('piW',M_.endo_names,'exact')), pistar);
fprintf('ii      = %.6f\n', oo_.steady_state(strmatch('ii',M_.endo_names,'exact')));
fprintf('1-jstar = %.4f  (expected ~0.35)\n', 1-oo_.steady_state(strmatch('jstar',M_.endo_names,'exact')));
fprintf('================================\n\n');
