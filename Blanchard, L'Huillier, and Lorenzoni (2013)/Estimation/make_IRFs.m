%% Make the IRFs of the "News, Noise, and Fluctuations"
% Replicate Figures 5 and 6 in the AER paper

clear all;
clc; 

% run the estimation
dynare estimate.mod %parallel conffile=[your conf file here] % optional but highly recommended

% save results
save('estimate_results.mat', 'M_', 'oo_', 'options_');


%% Make IRFs
% ensure you have the results saved

load('estimate_results.mat');

% Define the number of periods, shocks, and variables
T = 20; % Number of periods for IRFs (this needs to match the irf option in estimation()
shock_names = cellstr(M_.exo_names); % Get the shock names from M_.exo_names
var_names = cellstr(M_.endo_names); % Get the variable names from M_.endo_names
n_shock = length(shock_names); % Number of shocks
n_var = length(var_names); % Number of variables

% Initialise the IRF_matrix for mean and HPD bands
IRF_matrix_mean = zeros(T, n_shock, n_var);
IRF_matrix_low  = zeros(T, n_shock, n_var);
IRF_matrix_high = zeros(T, n_shock, n_var);

% Fill the IRF_matrix with the IRFs from oo_.PosteriorIRF.dsge
for j = 1:n_var
    for i = 1:n_shock
        var_name = var_names{j};
        shock_name = shock_names{i};
        field_name = [var_name '_' shock_name];
        if isfield(oo_.PosteriorIRF.dsge.Mean, field_name)
            IRF_matrix_mean(:, i, j) = oo_.PosteriorIRF.dsge.Mean.(field_name);
        else
            warning(['Field ' field_name ' not found in oo_.irfs']);
        end
    end
end

for j = 1:n_var
    for i = 1:n_shock
        var_name = var_names{j};
        shock_name = shock_names{i};
        field_name = [var_name '_' shock_name];
        if isfield(oo_.PosteriorIRF.dsge.HPDinf, field_name)
            IRF_matrix_low(:, i, j) = oo_.PosteriorIRF.dsge.HPDinf.(field_name);
        else
            warning(['Field ' field_name ' not found in oo_.irfs']);
        end
    end
end

for j = 1:n_var
    for i = 1:n_shock
        var_name = var_names{j};
        shock_name = shock_names{i};
        field_name = [var_name '_' shock_name];
        if isfield(oo_.PosteriorIRF.dsge.HPDsup, field_name)
            IRF_matrix_high(:, i, j) = oo_.PosteriorIRF.dsge.HPDsup.(field_name);
        else
            warning(['Field ' field_name ' not found in oo_.irfs']);
        end
    end
end

% Save the IRF matrices
save('IRF_matrices.mat', 'IRF_matrix_mean', 'IRF_matrix_low', 'IRF_matrix_high');

% Load the IRF matrices
load('IRF_matrices.mat');

% Initialize IRF objects
IRF_mean = zeros(T, n_shock+1, n_var);
IRF_low = zeros(T, n_shock+1, n_var);
IRF_high = zeros(T, n_shock+1, n_var);

% Process each IRF matrix
IRF_mean(:, 4:n_shock+1, 1:n_var) = IRF_matrix_mean(:, 3:n_shock, 1:n_var);
IRF_low(:, 4:n_shock+1, 1:n_var) = IRF_matrix_low(:, 3:n_shock, 1:n_var);
IRF_high(:, 4:n_shock+1, 1:n_var) = IRF_matrix_high(:, 3:n_shock, 1:n_var);

% Load parameter values from M_.params and assign them to variables
param_names = cellstr(M_.param_names); % Get the parameter names from M_.param_names
param_values = M_.params; % Get the parameter values from M_.params

% Assign the parameters their values
for k = 1:length(param_names)
    assignin('base', param_names{k}, param_values(k));
end

% Form transition matrices and solve filtering problem
A = [ 1+rho -rho 0 ; 1 0 0 ; 0 0 rho ];
S_1 = [ (1-rho)*sig_u 0 0 ; 0 0 0 ; 0 0 (rho^.5)*sig_u ]^2;
F = [ 1 0 1 ; 1 0 0 ];
S_2 = [0 0; 0 sig_nu]^2;

P = .01*eye(3,3);
K = P*F'/(F*P*F' + S_2);
for iter = 1:10000
    P = A*(P-K*F*P)*A' + S_1;
    Ki = P*F'/(F*P*F' + S_2);  
    dif = max(max(abs(Ki-K))); K = Ki;
    if dif < 1e-15
        break
    end
end

% Compute transitions for expected values
IKFA = (eye(3,3) - K*F)*A;
FA = F*A;
C = real((F*P*F'+S_2)^.5);
B = K*C;

% Compute IRFs for mean, low, and high matrices
IRF_mean = compute_irf(IRF_matrix_mean, T, n_var, A, F, K, C, rho, sig_u, sig_nu);
IRF_low = compute_irf(IRF_matrix_low, T, n_var, A, F, K, C, rho, sig_u, sig_nu);
IRF_high = compute_irf(IRF_matrix_high, T, n_var, A, F, K, C, rho, sig_u, sig_nu);

% Save the computed IRFs
save('computed_IRFs.mat', 'IRF_mean', 'IRF_low', 'IRF_high');

%%
% Load the IRF matrices
load('IRF_matrices.mat');

% Define the variables to plot and the number of shocks
variables_to_plot = {'ca', 'ia', 'ya', 'n'};
n_shocks = 3; % Number of shocks (e_1, e_2, news shock)
figure;
plot_index = 1;
T_plot = 20; % Ensure this matches the IRF length in estimation

% Define row labels and their vertical positions
row_labels = {'Permanent tech.', 'Transitory tech.', 'Noise'};
row_label_offsets = [-0.25, 0.2, 0.35]; % Adjust these values to set the vertical position of each label

% Loop through each shock and each variable to plot the IRFs
for shock = 1:n_shocks
    for i = 1:length(variables_to_plot)
        var_name = variables_to_plot{i};
        var_index = find(strcmp(var_names, var_name));
        if ~isempty(var_index)
            subplot(n_shocks, length(variables_to_plot), plot_index);
            % Plot the lower bound
            plot(1:T_plot, IRF_low(1:T_plot, shock, var_index), 'r--', 'LineWidth', 1);
            hold on;
            % Plot the point estimate
            plot(1:T_plot, IRF_mean(1:T_plot, shock, var_index), 'b-', 'LineWidth', 1);
            % Plot the upper bound
            plot(1:T_plot, IRF_high(1:T_plot, shock, var_index), 'r--', 'LineWidth', 1);
            if shock == 1
                title(M_.endo_names_long(strmatch(var_name, M_.endo_names, 'exact'), :));
            end
            hold on;
            plot_index = plot_index + 1;
        end
    end
    % Add row label
    % Select the first subplot in the current row
    subplot(n_shocks, length(variables_to_plot), (shock-1)*length(variables_to_plot) + 1);
    % Add the text label outside the plot area
    ax = gca;
    yl = ylim;
    text(ax.Position(1) - 0.4, mean(yl) + row_label_offsets(shock), row_labels{shock}, 'FontSize', 10, 'HorizontalAlignment', 'center', 'Rotation', 90, 'Units', 'normalized');
end

set(gcf, 'Color', 'w');


% Choose variables for Figure 6:
% Define the variables to plot and the number of shocks for the second figure
variables_to_plot = {'r', 'pi', 'wa'};
figure;
plot_index = 1;

row_label_offsets = [0.45, 0.5, 0.45]; % Adjust these values to set the vertical position of each label

% Loop through each shock and each variable to plot the IRFs
for shock = 1:n_shocks
    for i = 1:length(variables_to_plot)
        var_name = variables_to_plot{i};
        var_index = find(strcmp(var_names, var_name));
        if ~isempty(var_index)
            subplot(n_shocks, length(variables_to_plot), plot_index);
            % Plot the lower bound
            plot(1:T_plot, IRF_low(1:T_plot, shock, var_index), 'r--', 'LineWidth', 1);
            hold on;
            % Plot the point estimate
            plot(1:T_plot, IRF_mean(1:T_plot, shock, var_index), 'b-', 'LineWidth', 1);
            % Plot the upper bound
            plot(1:T_plot, IRF_high(1:T_plot, shock, var_index), 'r--', 'LineWidth', 1);
            if shock == 1
                title(M_.endo_names_long(strmatch(var_name, M_.endo_names, 'exact'), :));
            end
            hold on;
            plot_index = plot_index + 1;
        end
    end
    % Add row label
    % Select the first subplot in the current row
    subplot(n_shocks, length(variables_to_plot), (shock-1)*length(variables_to_plot) + 1);
    % Add the text label outside the plot area
    ax = gca;
    yl = ylim;
    text(ax.Position(1) - 0.4, mean(yl) + row_label_offsets(shock), row_labels{shock}, 'FontSize', 10, 'HorizontalAlignment', 'center', 'Rotation', 90, 'Units', 'normalized');
end

set(gcf, 'Color', 'w');


%% Function to compute IRF
function IRF = compute_irf(IRF_matrix, T, n_var, A, F, K, C, rho, sig_u, sig_nu)
    IRF = zeros(T, 3, n_var);
    for i = 1:n_var
        Z_e = IRF_matrix(:, 1:2, i); % select responses to e1-e2
        Z_u = zeros(T, 2);

        X = zeros(3, T); Xh = X; 
        Y = zeros(2, T); V = Y;
        e = zeros(2, T);

        for shock = 1:3
            if shock == 1    
                U1 = [(1-rho)*sig_u 0 0]'; U2 = [0 0]'; 
            elseif shock == 2
                U1 = [0 0 (rho^.5)*sig_u ]'; U2 = [0 0]'; 
            elseif shock == 3
                U1 = [0 0 0]'; U2 = [0 sig_nu]'; 
            end

            % Compute impulse responses 
            X(:, 1) = U1;
            Y(:, 1) = F*X(:, 1) + U2;
            Xh(:, 1) = K*Y(:, 1);
            % Get innovations
            V(:, 1) = Y(:, 1);
            e(:, 1) = C\V(:, 1);
            Z_u(1, shock) = Z_e(1, :)*e(:, 1);

            for t = 2:T
                % Update
                X(:, t) = A*X(:, t-1);
                Y(:, t) = F*X(:, t);
                Xh(:, t) = A*Xh(:, t-1) + K*(Y(:, t)-F*A*Xh(:, t-1));
                V(:, t) = Y(:, t)-F*A*Xh(:, t-1);
                e(:, t) = C\V(:, t);
                yy = 0;
                for l = 1:t % This is the convolution step
                    yy = yy + Z_e(l, :)*e(:, t+1-l);
                end
                Z_u(t, shock) = yy;
            end
        end
        IRF(:, 1:3, i) = Z_u;
    end
end