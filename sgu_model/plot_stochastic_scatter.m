%--------------------------------------------------------------------------
% Approach 3: Stochastic Simulation Scatter — SGU (2025) HDNWR Model
%
% Generates a long stochastic simulation (10,000 quarters) at order=2 and
% scatter-plots the simulated (u_t, piW_t) pairs against the theoretical
% nonlinear Phillips curve from Approach 1. If the model is working
% correctly, the simulated cloud should lie along the nonlinear locus.
%
% Run from: sgu_model/
% Can be run in the same MATLAB session as plot_asymmetric_irfs.m —
% if M_ is already in the workspace the dynare call is skipped.
%--------------------------------------------------------------------------

clc;

%--- Run Dynare (skip if already done in this session) ------------------%
if ~exist('M_', 'var')
    dynare hdnwr_nl
end

%--- Extract simulated time series from oo_.endo_simul ------------------%
% oo_.endo_simul is n_endo x (periods+2): first and last columns are
% the initial and terminal conditions added by Dynare.
% The usable simulation is columns 2:end-1.
get_idx = @(name) find(strcmp(name, cellstr(M_.endo_names)));
idx_u   = get_idx('u');
idx_piW = get_idx('piW');

u_sim   = oo_.endo_simul(idx_u,   2:end-1) * 100;       % percent
piW_sim = oo_.endo_simul(idx_piW, 2:end-1) * 4 * 100;   % annualised percent

fprintf('Simulation: %d periods\n', length(u_sim));
fprintf('  u   : mean = %.3f%%, std = %.3f%%\n', mean(u_sim),   std(u_sim));
fprintf('  piW : mean = %.3f%%, std = %.3f%%\n', mean(piW_sim), std(piW_sim));

%--- Static Phillips curve (inline from Approach 1 logic) ---------------%
% Calibration
eta    = M_.params(find(strcmp('eta',   cellstr(M_.param_names))));
Gam0   = M_.params(find(strcmp('Gam0',  cellstr(M_.param_names))));
Gam1   = M_.params(find(strcmp('Gam1',  cellstr(M_.param_names))));
del    = M_.params(find(strcmp('del',   cellstr(M_.param_names))));
pistar = M_.params(find(strcmp('pistar',cellstr(M_.param_names))));
un     = M_.params(find(strcmp('un',    cellstr(M_.param_names))));
theta  = M_.params(find(strcmp('theta', cellstr(M_.param_names))));

g1     = Gam0 + Gam1;
c_17   = Gam1 * (2 - eta);
g1_pow = g1^(2 - eta);

jstar_vec = linspace(0.001, 0.999, 2000)';
gjs       = Gam0 + Gam1 * jstar_vec;
F_jstar   = jstar_vec .* gjs.^(1 - eta) + (g1_pow - gjs.^(2 - eta)) / c_17;
piW_pc    = (F_jstar.^(1/(1-eta)) .* (1+pistar) - 1) * 4 * 100;  % annualised %
ratio     = g1 ./ gjs;
INT_s     = gjs ./ (Gam1*(1/theta+1)) .* (ratio.^(1/theta+1) - 1);
INT_d     = gjs ./ (Gam1*(1-eta))     .* (ratio.^(1-eta)     - 1);
u_pc      = (un + (1-un) .* (INT_s - INT_d) ./ (jstar_vec + INT_s)) * 100;

% Steady-state point
f_eq17   = @(j) (Gam0+Gam1*j)^(1-eta)*j + (g1_pow-(Gam0+Gam1*j)^(2-eta))/c_17 - 1;
jstar_ss = fzero(f_eq17, 0.65);
gjs_ss   = Gam0 + Gam1*jstar_ss;
rat_ss   = g1/gjs_ss;
IS_ss    = gjs_ss/(Gam1*(1/theta+1))*(rat_ss^(1/theta+1)-1);
ID_ss    = gjs_ss/(Gam1*(1-eta))    *(rat_ss^(1-eta)    -1);
u_ss_pct = (un + (1-un)*(IS_ss-ID_ss)/(jstar_ss+IS_ss)) * 100;
piW_ss_ann = pistar * 4 * 100;

%--- Figure --------------------------------------------------------------%
% Trim PC to the range covered by the simulation (with some margin)
u_lo   = max(0,   min(u_sim)   - 1);
u_hi   =          max(u_sim)   + 1;
piW_lo =          min(piW_sim) - 0.5;
piW_hi =          max(piW_sim) + 0.5;
mask_pc = u_pc >= u_lo & u_pc <= 15 & piW_pc >= piW_lo & piW_pc <= piW_hi;

fig = figure('Units', 'inches', 'Position', [1 1 6 4.5]);

% Simulated scatter
scatter(u_sim, piW_sim, 2, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.4);
hold on;

% Theoretical PC
plot(u_pc(mask_pc), piW_pc(mask_pc), 'b-', 'LineWidth', 2.5);

% Steady-state point
plot(u_ss_pct, piW_ss_ann, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');

xlim([u_lo u_hi]);
ylim([piW_lo piW_hi]);
xlabel('Unemployment rate (%)');
ylabel('Wage inflation (%, annualised)');
title('Simulated (u, \pi^W) vs Theoretical Phillips Curve — SGU (2025) order=2');
legend({'Simulated data (n=10,000)', 'Nonlinear PC (exact)', 'Steady state'}, ...
       'Location', 'northeast', 'FontSize', 10);
grid on; box on;
set(gca, 'FontSize', 11);

% Save
exportgraphics(fig, 'fig_stochastic_scatter.pdf');
fprintf('Figure saved to fig_stochastic_scatter.pdf\n');
