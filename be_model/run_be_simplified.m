%% run_be_simplified.m
% =========================================================================
% Driver script for the Benigno & Eggertsson (2025) simplified policy model
% with nonlinear (Inverse-L) Phillips curve implemented in Dynare/OccBin.
%
% This script:
%   1. Runs the Dynare model (incl. OccBin scenarios)
%   2. Plots OccBin impulse responses for different shock scenarios
%   3. Plots the Inverse-L Phillips curve (replicating Figure 8)
%   4. Plots AS-AD diagrams (replicating Figures 9-11)
%
% Requirements: MATLAB R2025b, Dynare 6.5
% =========================================================================

clear; close all; clc;

%% 1. Run Dynare model (includes OccBin scenarios via occbin_solver)
% =========================================================================
dynare be_simplified noclearall;

% M_, oo_, options_ are now global variables in the workspace
global M_ oo_ options_;

%% 2. Extract variable indices
% =========================================================================

var_names = M_.endo_names;
y_idx  = strmatch('y',  var_names, 'exact');
pi_idx = strmatch('pi', var_names, 'exact');
ii_idx = strmatch('ii', var_names, 'exact');

% Retrieve key parameters
ystar  = M_.params(strmatch('y_star', M_.param_names, 'exact'));
kap       = M_.params(strmatch('kappa', M_.param_names, 'exact'));
kap_tight = M_.params(strmatch('kappa_tight', M_.param_names, 'exact'));
kap_v     = M_.params(strmatch('kappa_v', M_.param_names, 'exact'));
kap_v_tight = M_.params(strmatch('kappa_v_tight', M_.param_names, 'exact'));
sig       = M_.params(strmatch('sigma', M_.param_names, 'exact'));
phipi     = M_.params(strmatch('phi_pi', M_.param_names, 'exact'));
ctilde    = M_.params(strmatch('c_tilde', M_.param_names, 'exact'));

%% 3. Plot OccBin IRFs
% =========================================================================
% The occbin_solver results are stored in results_demand, results_supply,
% results_combined, results_neg_demand (saved in verbatim blocks of .mod)
%
% oo_.occbin.piecewise is a matrix: rows = periods (including initial SS),
% columns = endogenous variables in declaration order.
% Row 1 = initial steady state, rows 2:end = simulation periods.

T_simul = 20;
quarters = 1:T_simul;

% Extract piecewise IRFs (all rows = simulation periods, 40x7 matrix)
y_demand  = results_demand.piecewise(:, y_idx);
pi_demand = results_demand.piecewise(:, pi_idx);
ii_demand = results_demand.piecewise(:, ii_idx);

y_neg     = results_neg_demand.piecewise(:, y_idx);
pi_neg    = results_neg_demand.piecewise(:, pi_idx);
ii_neg    = results_neg_demand.piecewise(:, ii_idx);

y_combined  = results_combined.piecewise(:, y_idx);
pi_combined = results_combined.piecewise(:, pi_idx);
ii_combined = results_combined.piecewise(:, ii_idx);

% Also extract linear (no regime switching) for comparison
y_demand_lin  = results_demand.linear(:, y_idx);
pi_demand_lin = results_demand.linear(:, pi_idx);

% --- Figure 1: Demand shock IRFs (positive vs negative) ---
figure('Name', 'Demand Shock: Asymmetric IRFs', 'Position', [100 100 900 700]);

subplot(2,2,1);
plot(quarters, y_demand*100, 'b-', 'LineWidth', 2); hold on;
plot(quarters, y_neg*100, 'r--', 'LineWidth', 2);
plot(quarters, y_demand_lin*100, 'b:', 'LineWidth', 1.5);
yline(ystar*100, 'k:', 'LineWidth', 1);
title('Output Gap');
ylabel('Percent');
xlabel('Quarters');
legend('Positive (piecewise)', 'Negative (piecewise)', ...
       'Positive (linear)', 'Y^*', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(quarters, pi_demand*400, 'b-', 'LineWidth', 2); hold on;
plot(quarters, pi_neg*400, 'r--', 'LineWidth', 2);
plot(quarters, pi_demand_lin*400, 'b:', 'LineWidth', 1.5);
title('Inflation (annualized)');
ylabel('Percentage points');
xlabel('Quarters');
legend('Positive (piecewise)', 'Negative (piecewise)', ...
       'Positive (linear)', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(quarters, ii_demand*400, 'b-', 'LineWidth', 2); hold on;
plot(quarters, ii_neg*400, 'r--', 'LineWidth', 2);
title('Nominal Interest Rate (annualized)');
ylabel('Percentage points');
xlabel('Quarters');
legend('Positive shock', 'Negative shock', 'Location', 'best');
grid on;

subplot(2,2,4);
% Dynamic Phillips curve paths
plot(y_demand*100, pi_demand*400, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10); hold on;
plot(y_neg*100, pi_neg*400, 'r.--', 'LineWidth', 1.5, 'MarkerSize', 10);
xline(ystar*100, 'k:', 'LineWidth', 1);
title('Phillips Curve (dynamic path)');
xlabel('Output Gap (%)');
ylabel('Inflation (ann., pp)');
legend('Positive shock', 'Negative shock', 'Location', 'best');
grid on;

sgtitle('Benigno-Eggertsson (2025): Asymmetric Demand Shock Responses');
exportgraphics(gcf, 'fig_demand_asymmetry.pdf');

% --- Figure 2: 2020s scenario (combined shocks) ---
figure('Name', '2020s Inflation Surge', 'Position', [100 100 900 400]);

subplot(1,3,1);
plot(quarters, y_combined*100, 'b-', 'LineWidth', 2);
yline(ystar*100, 'k:', 'LineWidth', 1);
title('Output Gap');
ylabel('Percent');
xlabel('Quarters');
grid on;

subplot(1,3,2);
plot(quarters, pi_combined*400, 'b-', 'LineWidth', 2);
title('Inflation (annualized)');
ylabel('Percentage points');
xlabel('Quarters');
grid on;

subplot(1,3,3);
plot(quarters, ii_combined*400, 'b-', 'LineWidth', 2);
title('Nominal Interest Rate (ann.)');
ylabel('Percentage points');
xlabel('Quarters');
grid on;

sgtitle('2020s Scenario: Combined Demand + Supply Shocks');
exportgraphics(gcf, 'fig_2020s_scenario.pdf');

%% 4. Static AS-AD Diagram (replicating Figure 9)
% =========================================================================
% Uses the paper's two-period Markov structure for the AS-AD representation
%
% AS curve (eq. 55):
%   Tight:  pi_S - pi* = -c_tilde/(1-tau) + kappa_tight/(1-tau) * Y_S
%                         + kappa_v_tight/(1-tau) * (1-alpha)*q_S + (pi_L^e - pi*)
%   Slack:  pi_S - pi* = kappa/(1-tau) * Y_S
%                         + kappa_v/(1-tau) * (1-alpha)*q_S + (pi_L^e - pi*)
%
% AD curve (eq. 53):
%   Y_S = D_S - sigma^{-1} * (phi_pi - tau)/(1-tau) * (pi_S - pi*)
%              + sigma^{-1} * (pi_L^e - pi*)

% Parameters for AS-AD analysis
tau = 0.8;          % Shock persistence probability
alpha = 0.9;        % Labor share

% Output gap grid
y_grid = linspace(-0.08, 0.08, 500);

% --- Baseline AS curve (no supply shock) ---
as_slack = (kap / (1-tau)) .* y_grid;
as_tight = -ctilde/(1-tau) + (kap_tight / (1-tau)) .* y_grid;
as_baseline = as_slack;
as_baseline(y_grid > ystar) = as_tight(y_grid > ystar);

% --- Shifted AS curve (with adverse supply shock) ---
q_shock = 0.10;    % 10% oil price increase
as_shift_slack = (kap / (1-tau)) .* y_grid + (kap_v / (1-tau)) * (1-alpha) * q_shock;
as_shift_tight = -ctilde/(1-tau) + (kap_tight / (1-tau)) .* y_grid ...
                  + (kap_v_tight / (1-tau)) * (1-alpha) * q_shock;
as_shifted = as_shift_slack;
as_shifted(y_grid > ystar) = as_shift_tight(y_grid > ystar);

% --- AD curves (eq. 53): Y = D - (1/sigma)*(phi_pi - tau)/(1-tau) * pi ---
% Rearranging for pi as function of Y:
% pi = -(1-tau)*sigma/(phi_pi - tau) * (Y - D)
D_baseline = 0;        % No demand shock
D_shifted  = 0.025;    % Positive demand shock

ad_coeff = -(1-tau) * sig / (phipi - tau);  % slope: d(pi)/d(Y)
pi_ad_baseline = ad_coeff .* (y_grid - D_baseline);
pi_ad_shifted  = ad_coeff .* (y_grid - D_shifted);

% --- Figure 3: AS-AD diagram (Figure 9 replica) ---
figure('Name', 'AS-AD Diagram: 2020s', 'Position', [100 100 800 600]);

% Convert to annualized percentage points (multiply by 400) and add 2% target
plot(y_grid*100, as_baseline*400 + 2, 'b-', 'LineWidth', 2); hold on;
plot(y_grid*100, as_shifted*400 + 2, 'b--', 'LineWidth', 2);
plot(y_grid*100, pi_ad_baseline*400 + 2, 'r-', 'LineWidth', 2);
plot(y_grid*100, pi_ad_shifted*400 + 2, 'r--', 'LineWidth', 2);

% Mark threshold
xline(ystar*100, 'k:', 'LineWidth', 1);

xlabel('Output Gap (%)');
ylabel('Inflation (annualized, %)');
title('AS-AD Diagram: 2020s Inflation Surge (BE 2025, Figure 9)');
legend('AS', 'AS'' (supply shock)', 'AD', 'AD'' (demand shock)', ...
       'Y^* threshold', 'Location', 'northwest');
ylim([-2 10]);
xlim([-8 8]);
grid on;

exportgraphics(gcf, 'fig_as_ad_2020s.pdf');

%% 5. Plot the Inverse-L Phillips Curve (Figure 8 replica)
% =========================================================================

figure('Name', 'Inverse-L Phillips Curve', 'Position', [100 100 700 500]);

% Static Phillips curve shape (no shocks, expectations at target)
pc_slack = (kap / (1-tau)) .* y_grid;
pc_tight = -ctilde/(1-tau) + (kap_tight / (1-tau)) .* y_grid;
pc_combined_static = pc_slack;
pc_combined_static(y_grid > ystar) = pc_tight(y_grid > ystar);

% Show the flat and steep regions separately for clarity
plot(y_grid(y_grid <= ystar)*100, pc_slack(y_grid <= ystar)*400 + 2, ...
    'b-', 'LineWidth', 2.5); hold on;
plot(y_grid(y_grid > ystar)*100, pc_tight(y_grid > ystar)*400 + 2, ...
    'r-', 'LineWidth', 2.5);

xline(ystar*100, 'k--', 'LineWidth', 1);
yline(2, 'k:', '\pi^* = 2%', 'LineWidth', 0.5);

xlabel('Output Gap (%)');
ylabel('Inflation (annualized, %)');
title('The Inverse-L Phillips Curve (BE 2025, Figure 8)');
legend('Slack regime (\theta \leq \theta^*)', ...
       'Tight regime (\theta > \theta^*)', ...
       'Y^* threshold', 'Location', 'northwest');
grid on;
ylim([-1 10]);

exportgraphics(gcf, 'fig_inverse_L_phillips_curve.pdf');

fprintf('\n=== All simulations and plots completed ===\n');
fprintf('Figures saved to current directory.\n');
