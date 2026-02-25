%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info
options_ = [];
M_.fname = 'be_simplified';
M_.dynare_version = '6.5';
oo_.dynare_version = '6.5';
options_.dynare_version = '6.5';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(4,1);
M_.exo_names_tex = cell(4,1);
M_.exo_names_long = cell(4,1);
M_.exo_names(1) = {'eps_d'};
M_.exo_names_tex(1) = {'eps\_d'};
M_.exo_names_long(1) = {'Demand shock innovation'};
M_.exo_names(2) = {'eps_s'};
M_.exo_names_tex(2) = {'eps\_s'};
M_.exo_names_long(2) = {'Supply shock innovation'};
M_.exo_names(3) = {'eps_m'};
M_.exo_names_tex(3) = {'eps\_m'};
M_.exo_names_long(3) = {'Monetary policy shock innovation'};
M_.exo_names(4) = {'eps_chi'};
M_.exo_names_tex(4) = {'eps\_chi'};
M_.exo_names_long(4) = {'Labor participation shock innovation'};
M_.endo_names = cell(7,1);
M_.endo_names_tex = cell(7,1);
M_.endo_names_long = cell(7,1);
M_.endo_names(1) = {'y'};
M_.endo_names_tex(1) = {'y'};
M_.endo_names_long(1) = {'Output gap'};
M_.endo_names(2) = {'pi'};
M_.endo_names_tex(2) = {'pi'};
M_.endo_names_long(2) = {'Inflation deviation from target'};
M_.endo_names(3) = {'ii'};
M_.endo_names_tex(3) = {'ii'};
M_.endo_names_long(3) = {'Nominal interest rate deviation'};
M_.endo_names(4) = {'g'};
M_.endo_names_tex(4) = {'g'};
M_.endo_names_long(4) = {'Demand shock state'};
M_.endo_names(5) = {'nu'};
M_.endo_names_tex(5) = {'nu'};
M_.endo_names_long(5) = {'Supply/cost-push shock state'};
M_.endo_names(6) = {'e'};
M_.endo_names_tex(6) = {'e'};
M_.endo_names_long(6) = {'Monetary policy shock state'};
M_.endo_names(7) = {'chi'};
M_.endo_names_tex(7) = {'chi'};
M_.endo_names_long(7) = {'Labor participation shock state'};
M_.endo_partitions = struct();
M_.param_names = cell(15,1);
M_.param_names_tex = cell(15,1);
M_.param_names_long = cell(15,1);
M_.param_names(1) = {'sigma'};
M_.param_names_tex(1) = {'\sigma'};
M_.param_names_long(1) = {'Inverse intertemporal elasticity of substitution'};
M_.param_names(2) = {'betta'};
M_.param_names_tex(2) = {'\beta'};
M_.param_names_long(2) = {'Discount factor'};
M_.param_names(3) = {'phi_pi'};
M_.param_names_tex(3) = {'\phi_\pi'};
M_.param_names_long(3) = {'Taylor rule coefficient on inflation'};
M_.param_names(4) = {'kappa'};
M_.param_names_tex(4) = {'\bar{\kappa}'};
M_.param_names_long(4) = {'PC slope on output gap, slack regime'};
M_.param_names(5) = {'kappa_v'};
M_.param_names_tex(5) = {'\bar{\kappa}_\nu'};
M_.param_names_long(5) = {'PC slope on supply shock, slack regime'};
M_.param_names(6) = {'kappa_tight'};
M_.param_names_tex(6) = {'\tilde{\kappa}^{tight}'};
M_.param_names_long(6) = {'PC slope on output gap, tight regime'};
M_.param_names(7) = {'kappa_v_tight'};
M_.param_names_tex(7) = {'\tilde{\kappa}_\nu^{tight}'};
M_.param_names_long(7) = {'PC slope on supply shock, tight regime'};
M_.param_names(8) = {'c_tilde'};
M_.param_names_tex(8) = {'\tilde{c}'};
M_.param_names_long(8) = {'Constant shift in tight regime'};
M_.param_names(9) = {'y_star'};
M_.param_names_tex(9) = {'\hat{Y}^*'};
M_.param_names_long(9) = {'Output gap threshold for regime switch'};
M_.param_names(10) = {'alpha_omega'};
M_.param_names_tex(10) = {'\alpha/\omega'};
M_.param_names_long(10) = {'Labor share over inverse Frisch elasticity'};
M_.param_names(11) = {'rho_g'};
M_.param_names_tex(11) = {'\rho_g'};
M_.param_names_long(11) = {'Demand shock persistence'};
M_.param_names(12) = {'rho_nu'};
M_.param_names_tex(12) = {'\rho_\nu'};
M_.param_names_long(12) = {'Supply shock persistence'};
M_.param_names(13) = {'rho_e'};
M_.param_names_tex(13) = {'\rho_e'};
M_.param_names_long(13) = {'Monetary policy shock persistence'};
M_.param_names(14) = {'rho_chi'};
M_.param_names_tex(14) = {'\rho_\chi'};
M_.param_names_long(14) = {'Labor participation shock persistence'};
M_.param_names(15) = {'occbin_tight_labor_bind'};
M_.param_names_tex(15) = {'occbin\_tight\_labor\_bind'};
M_.param_names_long(15) = {'occbin_tight_labor_bind'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 4;
M_.endo_nbr = 7;
M_.param_nbr = 15;
M_.orig_endo_nbr = 7;
M_.aux_vars = [];
M_.Sigma_e = zeros(4, 4);
M_.Correlation_matrix = eye(4, 4);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.learnt_shocks = [];
M_.learnt_endval = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
M_.matched_irfs = {};
M_.matched_irfs_weights = {};
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.ramsey_policy = false;
options_.discretionary_policy = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.eq_nbr = 7;
M_.ramsey_orig_eq_nbr = 0;
M_.ramsey_orig_endo_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 5 12;
 0 6 13;
 0 7 0;
 1 8 14;
 2 9 0;
 3 10 0;
 4 11 0;]';
M_.nstatic = 1;
M_.nfwrd   = 2;
M_.npred   = 3;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 4;
M_.ndynamic   = 6;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.equations_tags = {
  1 , 'name' , 'IS curve' ;
  2 , 'dynamic' , '' ;
  2 , 'name' , 'Phillips curve' ;
  3 , 'name' , 'Taylor rule' ;
  4 , 'name' , 'Demand shock' ;
  5 , 'name' , 'Supply shock' ;
  6 , 'name' , 'Monetary policy shock' ;
  7 , 'name' , 'Labor participation shock' ;
};
M_.mapping.y.eqidx = [1 2 ];
M_.mapping.pi.eqidx = [1 2 3 ];
M_.mapping.ii.eqidx = [1 3 ];
M_.mapping.g.eqidx = [1 4 ];
M_.mapping.nu.eqidx = [2 5 ];
M_.mapping.e.eqidx = [3 6 ];
M_.mapping.chi.eqidx = [2 7 ];
M_.mapping.eps_d.eqidx = [4 ];
M_.mapping.eps_s.eqidx = [5 ];
M_.mapping.eps_m.eqidx = [6 ];
M_.mapping.eps_chi.eqidx = [7 ];
M_.static_and_dynamic_models_differ = true;
M_.has_external_function = false;
M_.block_structure.time_recursive = false;
M_.block_structure.block(1).Simulation_Type = 1;
M_.block_structure.block(1).endo_nbr = 4;
M_.block_structure.block(1).mfs = 4;
M_.block_structure.block(1).equation = [ 4 5 6 7];
M_.block_structure.block(1).variable = [ 4 5 6 7];
M_.block_structure.block(1).is_linear = true;
M_.block_structure.block(1).NNZDerivatives = 8;
M_.block_structure.block(1).bytecode_jacob_cols_to_sparse = [1 2 3 4 5 6 7 8 ];
M_.block_structure.block(2).Simulation_Type = 7;
M_.block_structure.block(2).endo_nbr = 3;
M_.block_structure.block(2).mfs = 2;
M_.block_structure.block(2).equation = [ 3 2 1];
M_.block_structure.block(2).variable = [ 3 2 1];
M_.block_structure.block(2).is_linear = true;
M_.block_structure.block(2).NNZDerivatives = 8;
M_.block_structure.block(2).bytecode_jacob_cols_to_sparse = [0 1 2 0 0 ];
M_.block_structure.block(1).g1_sparse_rowval = int32([]);
M_.block_structure.block(1).g1_sparse_colval = int32([]);
M_.block_structure.block(1).g1_sparse_colptr = int32([]);
M_.block_structure.block(2).g1_sparse_rowval = int32([1 2 1 2 ]);
M_.block_structure.block(2).g1_sparse_colval = int32([1 1 2 2 ]);
M_.block_structure.block(2).g1_sparse_colptr = int32([1 3 5 ]);
M_.block_structure.variable_reordered = [ 4 5 6 7 3 2 1];
M_.block_structure.equation_reordered = [ 4 5 6 7 3 2 1];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [
 4 4;
 5 5;
 6 6;
 7 7;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [
 1 1;
 1 3;
 1 4;
 2 1;
 2 2;
 2 5;
 2 7;
 3 2;
 3 3;
 3 6;
 4 4;
 5 5;
 6 6;
 7 7;
];
M_.block_structure.incidence(3).lead_lag = 1;
M_.block_structure.incidence(3).sparse_IM = [
 1 1;
 1 2;
 1 4;
 2 2;
];
M_.block_structure.dyn_tmp_nbr = 0;
M_.state_var = [4 5 6 7 ];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(7, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(4, 1);
M_.params = NaN(15, 1);
M_.endo_trends = struct('deflator', cell(7, 1), 'log_deflator', cell(7, 1), 'growth_factor', cell(7, 1), 'log_growth_factor', cell(7, 1));
M_.NNZDerivatives = [26; 0; -1; ];
M_.dynamic_g1_sparse_rowval = int32([4 5 6 7 1 2 2 3 1 3 1 4 2 5 3 6 2 7 1 1 2 1 4 5 6 7 ]);
M_.dynamic_g1_sparse_colval = int32([4 5 6 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 16 16 18 22 23 24 25 ]);
M_.dynamic_g1_sparse_colptr = int32([1 1 1 1 2 3 4 5 7 9 11 13 15 17 19 20 22 22 23 23 23 23 24 25 26 27 ]);
M_.dynamic_g2_sparse_indices = int32([]);
M_.lhs = {
'y-g'; 
'(pi-(kappa*(y+alpha_omega*chi)+kappa_v*nu+pi(1)*betta))*(1-occbin_tight_labor_bind)+occbin_tight_labor_bind*(pi-(pi(1)*betta+(y+alpha_omega*chi)*kappa_tight-c_tilde+nu*kappa_v_tight))'; 
'ii'; 
'g'; 
'nu'; 
'e'; 
'chi'; 
};
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.block_structure_stat.block(1).Simulation_Type = 3;
M_.block_structure_stat.block(1).endo_nbr = 1;
M_.block_structure_stat.block(1).mfs = 1;
M_.block_structure_stat.block(1).equation = [ 4];
M_.block_structure_stat.block(1).variable = [ 4];
M_.block_structure_stat.block(2).Simulation_Type = 3;
M_.block_structure_stat.block(2).endo_nbr = 1;
M_.block_structure_stat.block(2).mfs = 1;
M_.block_structure_stat.block(2).equation = [ 5];
M_.block_structure_stat.block(2).variable = [ 5];
M_.block_structure_stat.block(3).Simulation_Type = 3;
M_.block_structure_stat.block(3).endo_nbr = 1;
M_.block_structure_stat.block(3).mfs = 1;
M_.block_structure_stat.block(3).equation = [ 6];
M_.block_structure_stat.block(3).variable = [ 6];
M_.block_structure_stat.block(4).Simulation_Type = 3;
M_.block_structure_stat.block(4).endo_nbr = 1;
M_.block_structure_stat.block(4).mfs = 1;
M_.block_structure_stat.block(4).equation = [ 7];
M_.block_structure_stat.block(4).variable = [ 7];
M_.block_structure_stat.block(5).Simulation_Type = 6;
M_.block_structure_stat.block(5).endo_nbr = 3;
M_.block_structure_stat.block(5).mfs = 3;
M_.block_structure_stat.block(5).equation = [ 2 3 1];
M_.block_structure_stat.block(5).variable = [ 1 2 3];
M_.block_structure_stat.variable_reordered = [ 4 5 6 7 1 2 3];
M_.block_structure_stat.equation_reordered = [ 4 5 6 7 2 3 1];
M_.block_structure_stat.incidence.sparse_IM = [
 1 2;
 1 3;
 2 1;
 2 2;
 2 5;
 2 7;
 3 2;
 3 3;
 3 6;
 4 4;
 5 5;
 6 6;
 7 7;
];
M_.block_structure_stat.tmp_nbr = 0;
M_.block_structure_stat.block(1).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(2).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(2).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(3).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(3).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(3).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(4).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(4).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(4).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(5).g1_sparse_rowval = int32([1 1 2 3 2 3 ]);
M_.block_structure_stat.block(5).g1_sparse_colval = int32([1 2 2 2 3 3 ]);
M_.block_structure_stat.block(5).g1_sparse_colptr = int32([1 2 5 7 ]);
M_.static_g1_sparse_rowval = int32([2 1 2 3 1 3 4 2 5 3 6 2 7 ]);
M_.static_g1_sparse_colval = int32([1 2 2 2 3 3 4 5 5 6 6 7 7 ]);
M_.static_g1_sparse_colptr = int32([1 2 5 7 8 10 12 14 ]);
M_.params(1) = 0.5;
sigma = M_.params(1);
M_.params(2) = 0.99;
betta = M_.params(2);
M_.params(3) = 1.5;
phi_pi = M_.params(3);
M_.params(4) = 0.0065;
kappa = M_.params(4);
M_.params(5) = 0.0093;
kappa_v = M_.params(5);
M_.params(6) = 0.0736;
kappa_tight = M_.params(6);
M_.params(7) = 0.2742;
kappa_v_tight = M_.params(7);
M_.params(10) = 0.9;
alpha_omega = M_.params(10);
M_.params(9) = 0.01;
y_star = M_.params(9);
M_.params(8) = (M_.params(6)-M_.params(4))*M_.params(9);
c_tilde = M_.params(8);
M_.params(11) = 0.8;
rho_g = M_.params(11);
M_.params(12) = 0.8;
rho_nu = M_.params(12);
M_.params(13) = 0.5;
rho_e = M_.params(13);
M_.params(14) = 0.8;
rho_chi = M_.params(14);
M_.params(15) = 0;
occbin_tight_labor_bind = M_.params(15);
M_.occbin.constraint_nbr = 1;
M_.occbin.pswitch = [
15 ];
options_.occbin = struct();
options_.occbin = occbin.set_default_options(options_.occbin, M_);
oo_.dr=set_state_space(oo_.dr,M_);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
M_.surprise_shocks = [
struct('exo_id',2,'periods',1:1,'value',0.01);
];
options_occbin_ = struct();
options_occbin_.simul.check_ahead_periods = 20;
options_occbin_.simul.periods = 20;
[M_, options_] = occbin.setup(M_, options_, options_occbin_);
if ~isfield(options_,'occbin')
    options_.occbin = struct();
end
[oo_.dr, oo_.occbin.simul]= occbin.solver(M_, options_, oo_.dr , oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    results_supply.piecewise = oo_.occbin.simul.piecewise;
    results_supply.linear    = oo_.occbin.simul.linear;
    fprintf('\n=== Scenario B (supply shock) completed ===\n');
M_.surprise_shocks = [
struct('exo_id',1,'periods',1:1,'value',(-0.05));
];
options_occbin_ = struct();
options_occbin_.simul.check_ahead_periods = 20;
options_occbin_.simul.periods = 20;
[M_, options_] = occbin.setup(M_, options_, options_occbin_);
if ~isfield(options_,'occbin')
    options_.occbin = struct();
end
[oo_.dr, oo_.occbin.simul]= occbin.solver(M_, options_, oo_.dr , oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    results_neg_demand.piecewise = oo_.occbin.simul.piecewise;
    results_neg_demand.linear    = oo_.occbin.simul.linear;
    fprintf('\n=== Scenario D (negative demand shock) completed ===\n');
M_.surprise_shocks = [
struct('exo_id',1,'periods',1:1,'value',0.05);
];
options_occbin_ = struct();
options_occbin_.simul.check_ahead_periods = 20;
options_occbin_.simul.periods = 20;
[M_, options_] = occbin.setup(M_, options_, options_occbin_);
if ~isfield(options_,'occbin')
    options_.occbin = struct();
end
[oo_.dr, oo_.occbin.simul]= occbin.solver(M_, options_, oo_.dr , oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    results_demand.piecewise = oo_.occbin.simul.piecewise;
    results_demand.linear    = oo_.occbin.simul.linear;
    fprintf('\n=== Scenario A (demand shock) completed ===\n');
    fprintf('y(1) = %.6f, pi(1) = %.6f\n', ...
            results_demand.piecewise(1,1), results_demand.piecewise(1,2));
M_.surprise_shocks = [
struct('exo_id',1,'periods',1:1,'value',0.05);
struct('exo_id',2,'periods',1:1,'value',0.001);
];
options_occbin_ = struct();
options_occbin_.simul.check_ahead_periods = 20;
options_occbin_.simul.periods = 20;
[M_, options_] = occbin.setup(M_, options_, options_occbin_);
if ~isfield(options_,'occbin')
    options_.occbin = struct();
end
[oo_.dr, oo_.occbin.simul]= occbin.solver(M_, options_, oo_.dr , oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    results_combined.piecewise = oo_.occbin.simul.piecewise;
    results_combined.linear    = oo_.occbin.simul.linear;
    fprintf('\n=== Scenario C (combined shocks) completed ===\n');


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'be_simplified_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'be_simplified_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'be_simplified_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'be_simplified_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'be_simplified_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'be_simplified_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'be_simplified_results.mat'], 'oo_recursive_', '-append');
end
if exist('options_mom_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'be_simplified_results.mat'], 'options_mom_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
