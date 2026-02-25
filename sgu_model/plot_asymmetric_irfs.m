%--------------------------------------------------------------------------
% Approach 2: Asymmetric Impulse Responses — SGU (2025) HDNWR Model
%
% Compares the model's response to a +1 std dev vs -1 std dev monetary
% policy shock (eps_mu) at second order. At order=1, the responses are
% exactly antisymmetric. At order=2, they differ — revealing the curvature
% of the wage Phillips curve.
%
% A positive eps_mu shock is CONTRACTIONARY (mu rises -> ii rises above
% rule): u rises, piW falls, y falls.
% A negative eps_mu shock is EXPANSIONARY (mu falls -> ii falls below
% rule): u falls, piW rises, y rises.
%
% Both IRFs are plotted with their correct signs. At order=1 they would be
% exact mirror images. At order=2 the magnitudes differ — e.g. unemployment
% rises more from the contractionary shock than it falls from the
% expansionary shock, consistent with the convex Phillips curve.
%
% Method:
%   1. Run dynare hdnwr_nl (order=2, pruning).
%   2. Use simult_() to simulate paths for +shock, -shock, and zero-shock
%      baseline starting from the deterministic steady state.
%   3. Plot actual IRFs for both shocks. Overlay irf_pos + irf_neg (dashed
%      black): exactly zero at order=1, nonzero at order=2.
%
% Run from: sgu_model/
%--------------------------------------------------------------------------

clear; clc;

%--- Run Dynare ----------------------------------------------------------%
% Populates M_, oo_, options_ in the workspace.
% If hdnwr_nl has already been run in this MATLAB session, comment out
% the next line to skip re-solving.
dynare hdnwr_nl

%--- Shock setup ---------------------------------------------------------%
idx_mu  = find(strcmp('eps_mu', cellstr(M_.exo_names)));
sig_mu  = sqrt(M_.Sigma_e(idx_mu, idx_mu));   % = 0.01 (1% std dev)

T       = 20;       % IRF horizon (quarters)
ex_zero = zeros(T, M_.exo_nbr);
ex_pos  = ex_zero;  ex_pos(1, idx_mu) = +sig_mu;
ex_neg  = ex_zero;  ex_neg(1, idx_mu) = -sig_mu;

%--- Simulate paths using simult_ ----------------------------------------%
% simult_ signature (Dynare >= 6.x):  y_ = simult_(M_, options_, y0, dr, ex_, iorder)
% Output y_ is n_endo x (T+1); columns are t=0,1,...,T.
y0     = oo_.dr.ys;
y_base = simult_(M_, options_, y0, oo_.dr, ex_zero, 2);
y_pos  = simult_(M_, options_, y0, oo_.dr, ex_pos,  2);
y_neg  = simult_(M_, options_, y0, oo_.dr, ex_neg,  2);

% IRFs = deviation from zero-shock baseline (correctly accounts for
% any second-order constant correction in the baseline)
irf_pos = y_pos - y_base;   % n_endo x (T+1)
irf_neg = y_neg - y_base;

% Drop t=0 (initial condition); IRF starts at t=1
irf_pos = irf_pos(:, 2:end);   % n_endo x T
irf_neg = irf_neg(:, 2:end);
t_vec   = 1:T;

%--- Variable indices ----------------------------------------------------%
% simult_ output ordering matches declaration order (M_.endo_names).
% Note: if results look inconsistent, check oo_.dr.order_var for the
% mapping between DR ordering and declaration ordering.
get_idx = @(name) find(strcmp(name, cellstr(M_.endo_names)));
idx_u   = get_idx('u');
idx_piW = get_idx('piW');
idx_y   = get_idx('y');

%--- Convert to convenient units ----------------------------------------%
% u, piW in percentage points (pp); y in percent deviation from SS
u_pos   =  irf_pos(idx_u,   :) * 100;
u_neg   =  irf_neg(idx_u,   :) * 100;
piW_pos =  irf_pos(idx_piW, :) * 4 * 100;   % annualised pp
piW_neg =  irf_neg(idx_piW, :) * 4 * 100;
y_pos   =  irf_pos(idx_y,   :) / oo_.dr.ys(idx_y) * 100;
y_neg   =  irf_neg(idx_y,   :) / oo_.dr.ys(idx_y) * 100;

%--- Report asymmetry ---------------------------------------------------%
% irf_pos + irf_neg = 0 at order=1 (antisymmetry); nonzero at order=2.
% The sign of the asymmetry for u tells us about Phillips curve curvature:
%   irf_pos(u) + irf_neg(u) > 0  =>  unemployment rises more from a
%   contractionary shock than it falls from an equally-sized expansionary
%   shock, consistent with a convex (concave toward origin) Phillips curve.
asym_u   = max(abs(u_pos   + u_neg));
asym_piW = max(abs(piW_pos + piW_neg));
asym_y   = max(abs(y_pos   + y_neg));
fprintf('\n=== Asymmetry: max |IRF(+shock) + IRF(-shock)| ===\n');
fprintf('  (Zero at order=1; nonzero reveals Phillips curve curvature)\n');
fprintf('  u   : %.6f pp\n',             asym_u);
fprintf('  piW : %.6f pp (annualised)\n', asym_piW);
fprintf('  y   : %.6f %%\n',             asym_y);
fprintf('  Peak u asymmetry at quarter %d\n', find(abs(u_pos + u_neg) == asym_u, 1));
if max([asym_u, asym_piW, asym_y]) < 1e-10
    warning('Asymmetry is at machine precision — check that order=2 is active.');
end
fprintf('===================================================\n\n');

%--- Figure --------------------------------------------------------------%
% Each panel plots the ACTUAL IRF for both shocks without sign normalisation.
% +1sigma shock (contractionary): u rises, piW falls, y falls.
% -1sigma shock (expansionary):   u falls, piW rises, y rises.
% At order=1 the IRFs would be exact mirror images; at order=2 the
% magnitudes differ, revealing the curvature of the Phillips curve.
% The dotted line (irf_pos + irf_neg) is exactly zero at order=1.
col_pos = [0.15 0.45 0.75];   % blue  — contractionary (+sigma)
col_neg = [0.80 0.20 0.20];   % red   — expansionary   (-sigma)
col_sum = [0.00 0.00 0.00];   % black — asymmetry term

var_names   = {'u', 'piW', 'y'};
irf_pos_all = {u_pos, piW_pos, y_pos};
irf_neg_all = {u_neg, piW_neg, y_neg};
ylabels     = {'Unemployment (pp)', 'Wage inflation (pp, ann.)', 'Output (%)'};

fig = figure('Units', 'inches', 'Position', [1 1 11 3.8]);
for k = 1:3
    subplot(1, 3, k);
    ip = irf_pos_all{k};
    in = irf_neg_all{k};

    plot(t_vec, ip,      '-',  'Color', col_pos, 'LineWidth', 2);   hold on;
    plot(t_vec, in,      '--', 'Color', col_neg, 'LineWidth', 2);
    plot(t_vec, ip + in, ':',  'Color', col_sum, 'LineWidth', 1.5); % asymmetry: 0 iff order=1
    yline(0, 'k-', 'LineWidth', 0.5);

    xlim([1 T]);
    xlabel('Quarters');
    ylabel(ylabels{k});
    title(var_names{k});
    grid on; box on;
    set(gca, 'FontSize', 10);

    if k == 1
        legend({'+\sigma shock (contractionary)', ...
                '-\sigma shock (expansionary)',   ...
                'Sum (0 \Leftrightarrow linear)'}, ...
               'Location', 'best', 'FontSize', 9);
    end
end
sgtitle('Asymmetric IRFs to Monetary Policy Shock — SGU (2025) order=2', 'FontSize', 12);

% Save
exportgraphics(fig, 'fig_asymmetric_irfs.pdf');
fprintf('Figure saved to fig_asymmetric_irfs.pdf\n');
