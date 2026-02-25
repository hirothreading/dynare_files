%--------------------------------------------------------------------------
% Static Wage Phillips Curve — SGU (2025) HDNWR Model
% Endogenous Labour Supply Version (Definition C1)
%
% Traces equations (17) and (28) parametrically over a grid of jstar values
% to produce the exact nonlinear (u, piW) locus. No Dynare required.
%
% Output: Figure showing the nonlinear PC, its linear tangent at SS, and
%         the steady-state point.
%
% Run from: sgu_model/
%--------------------------------------------------------------------------

clear; clc;

%--- Calibration (Table 1, identical to hdnwr.mod) -----------------------%
eta    = 11;
Gam0   = 0.978;
Gam1   = 0.031;
del    = 1;
pistar = 1.03^(1/4) - 1;    % quarterly ≈ 0.0074
un     = 0.04;
theta  = 5;

%--- Derived scalar quantities -------------------------------------------%
g1     = Gam0 + Gam1;       % Gam0 + Gam1*j at j = 1
c_inf  = (1 + pistar)^(del * (1 - eta));   % scalar prefactor in INT17
c_17   = Gam1 * (2 - eta);  % denominator coefficient in INT17 (< 0 since eta=11)
g1_pow = g1^(2 - eta);      % g1^{2-eta}, scalar

%--- jstar grid ----------------------------------------------------------%
% Avoid exact endpoints to prevent numerical issues at boundaries.
% At jstar -> 0: ratio -> g1/Gam0 (finite); at jstar -> 1: ratio -> 1,
% INT_s and INT_d -> 0, piW determined by eq (17) alone.
jstar_vec = linspace(0.001, 0.999, 2000)';

%--- Model-local variables (vectorised) ----------------------------------%
gjs  = Gam0 + Gam1 * jstar_vec;          % Gam0 + Gam1*j  at j = jstar

%-- Eq (17): Wage aggregation --
% (1+piW)^{1-eta} = jstar * gamjstar^{1-eta} + INT17
% gamjstar = (1+pistar)^del * gjs
% Factor out (1+pistar)^{1-eta}:
% [(1+piW)/(1+pistar)]^{1-eta} = jstar * gjs^{1-eta}
%                               + [g1^{2-eta} - gjs^{2-eta}] / [Gam1*(2-eta)]
%                           =: F(jstar)
%
% => piW = (F(jstar))^{1/(1-eta)} * (1+pistar) - 1
F_jstar = jstar_vec .* gjs.^(1 - eta) + (g1_pow - gjs.^(2 - eta)) / c_17;
piW_vec = F_jstar.^(1 / (1 - eta)) .* (1 + pistar) - 1;

%-- Eq (28): Unemployment (endogenous labour supply) --
ratio = g1 ./ gjs;                                          % g1 / gjs
INT_s = gjs ./ (Gam1 * (1/theta + 1)) .* (ratio.^(1/theta + 1) - 1);
INT_d = gjs ./ (Gam1 * (1 - eta))     .* (ratio.^(1 - eta)     - 1);
u_vec = un + (1 - un) .* (INT_s - INT_d) ./ (jstar_vec + INT_s);

%--- Convert to plot units -----------------------------------------------%
u_pct   = u_vec   * 100;                 % percent
piW_ann = piW_vec * 4 * 100;            % annualised percent

%--- Steady-state point --------------------------------------------------%
% Solve for jstar_ss using the normalised eq (17) (del=1 case):
%   F(j) - 1 = 0
f_eq17  = @(j) (Gam0 + Gam1*j)^(1 - eta) * j ...
               + (g1_pow - (Gam0 + Gam1*j)^(2 - eta)) / c_17 - 1;
jstar_ss = fzero(f_eq17, 0.65);

gjs_ss  = Gam0 + Gam1 * jstar_ss;
ratio_ss = g1 / gjs_ss;
INT_s_ss = gjs_ss / (Gam1 * (1/theta + 1)) * (ratio_ss^(1/theta + 1) - 1);
INT_d_ss = gjs_ss / (Gam1 * (1 - eta))     * (ratio_ss^(1 - eta)     - 1);
u_ss     = un + (1 - un) * (INT_s_ss - INT_d_ss) / (jstar_ss + INT_s_ss);
piW_ss   = pistar;

u_ss_pct   = u_ss   * 100;
piW_ss_ann = piW_ss * 4 * 100;

fprintf('Steady state: jstar = %.4f, u = %.4f%%, piW = %.4f%% (annual)\n', ...
        jstar_ss, u_ss_pct, piW_ss_ann);

%--- Plot range (defined here so the linear tangent can use it) ----------%
piW_lo = -2; piW_hi = 14;
u_lo   =  0; u_hi   = 22;

%--- Linear approximation (tangent at steady state) ----------------------%
% Numerical slope dpiW_ann/du_pct from the parametric curve
[~, idx_ss] = min(abs(jstar_vec - jstar_ss));
% Use a window of +/- 5 points for the finite-difference slope
w = 5;
lo = max(1, idx_ss - w);
hi = min(length(jstar_vec), idx_ss + w);
slope_lin = (piW_ann(hi) - piW_ann(lo)) / (u_pct(hi) - u_pct(lo));

% Tangent line over the full u range of the plot
u_lin   = linspace(u_lo, u_hi, 500)';
piW_lin = piW_ss_ann + slope_lin * (u_lin - u_ss_pct);

%--- Figure --------------------------------------------------------------%
fig = figure('Units', 'inches', 'Position', [1 1 6 4.5]);

% Trim extremes that go outside a sensible plot range
mask = piW_ann >= piW_lo & piW_ann <= piW_hi & u_pct >= u_lo & u_pct <= 15;

% Nonlinear PC
plot(u_pct(mask), piW_ann(mask), 'b-', 'LineWidth', 2); hold on;

% Linear approximation (clip to plot range)
lin_mask = piW_lin >= piW_lo & piW_lin <= piW_hi & u_lin >= u_lo & u_lin <= u_hi;
plot(u_lin(lin_mask), piW_lin(lin_mask), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

% Steady-state point
plot(u_ss_pct, piW_ss_ann, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');

% Formatting
xlim([u_lo u_hi]);
ylim([piW_lo piW_hi]);
xlabel('Unemployment rate (%)');
ylabel('Wage inflation (%, annualised)');
title('Wage Phillips Curve — SGU (2025) HDNWR Model');
legend({'Nonlinear PC (model)', 'Linear approximation at SS', 'Steady state'}, ...
       'Location', 'northeast');
grid on; box on;
set(gca, 'FontSize', 11);

% Save
exportgraphics(fig, 'fig_phillips_curve.pdf');
fprintf('Figure saved to fig_phillips_curve.pdf\n');
