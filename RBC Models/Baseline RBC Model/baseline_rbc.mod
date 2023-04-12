// Baseline RBC Model 
// A simple Hansen-style RBC model
// By David Murakami

// DECLARE VARIABLES
var
// QUANTITIES
Y           $Y$                 (long_name='Output')
C           $C$                 (long_name='Consumption')
I           $I$                 (long_name='Investment')
L           $L$                 (long_name='Labour Supply')
K           $K$                 (long_name='Capital')
// PRICES
R           $R$                 (long_name='Interest Rate')
W           $w$                 (long_name='Wages')
// TECHNOLOGY
A           $A$                 (long_name='TFP')
;

// EXOGENOUS SHOCKS
varexo
eps_a       $\varepsilon_a$     (long_name='TFP Shock')
;

// PARAMETERS
parameters
betta        $\beta$            (long_name='Discount factor')
alphha       $\alpha$           (long_name='Capital share')
nu           $\nu$              (long_name='Frisch labour supply elasticity')
delta        $\delta$           (long_name='Depreciation rate')
chi          $\chi$             (long_name='Labour disutility')
rho_a        $\rho_a$           (long_name='TFP shock persistence')
sigma_a      $\sigma_a$         (long_name='TFP shock standard deviation')
;
// PARAMETERISE
betta = 0.99;
alphha = 0.3;
nu = 2;
delta = 0.025;
chi = 4.5;
rho_a = 0.95;
sigma_a = 0.01;

// DECLARE MODEL
model;
// Note: Period utility is: ln(C) - chi*L^(1+1/nu)/(1+1/nu)
[name='Consumption Euler equation']
1/C = betta*1/C(+1)*R(+1);

[name='Intratemporal Euler equation']
chi*L^(1/nu) = (1-alphha)*Y/(C*L);

[name='Resource constraint']
Y = C + I;

[name='Law of motion for capital']
K = (1-delta)*K(-1) + I;

[name='Production function']
Y = A*K(-1)^alphha*L^(1-alphha);

[name='Return on capital']
R = alphha*(Y/K(-1)) + 1 - delta;

[name='Wage rate']
W = (1-alphha)*(Y/L);

[name='TFP process']
log(A) = rho_a*log(A(-1)) + eps_a;
end;

// PROVIDE ANALYTICAL STEADY STATE EXPRESSIONS
steady_state_model;
A = 1;
KL = (alphha/(1/betta-1+delta))^(1/(1-alphha));
L = (((1-alphha)*KL^alphha)/(chi*KL^alphha - chi*delta*KL))^(nu/(1+nu));
K = KL*L;
I = delta*K;
Y = K^alphha*L^(1-alphha);
C = Y - delta*K;
R = 1/betta;
W = (1-alphha)*KL^alphha;
end;

steady;
check;
model_diagnostics;

// DECLARE SHOCKS
shocks;
var eps_a = sigma_a^2;
end;

stoch_simul(order=1,irf=40);

write_latex_original_model;
write_latex_static_model;
write_latex_dynamic_model(write_equation_tags);
write_latex_parameter_table;
write_latex_definitions;
write_latex_steady_state_model;
