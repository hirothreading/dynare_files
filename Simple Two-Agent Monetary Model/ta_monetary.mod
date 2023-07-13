// Baseline Two-Agent Endowment Economy
// Written by David Murakami, Ivan Shchapov, and Ganesh Viswanath-Natraj
// Version: 31 May 2022
// Implemented in Dynare 5.1


// Declare variables
// Quantity variables
var
    C_h       $C^h$        (long_name='BHH Consumption')
    C_u       $C^u$        (long_name='UHH Consumption')
    C         $C$          (long_name='Agg. Consumption')
    T_h       $T^h$        (long_name='BHH Transfers')
    T_u       $T^u$        (long_name='UHH Transfers')
    M         $M$          (long_name='Real Money')
// Prices
    R         $R$          (long_name='Nominal Interest Rate')
    Pi        $\pi$        (long_name='Inflation')
// Exogenous processes
    xi_M      $\xi^M$      (long_name='Money Growth Process')
// Other variables
    chi_M     $\chi^M$     (long_name='Money Adjustment Cost')
    lambda_u  $\lambda^u$  (long_name='BC Multiplier')
    mu_u      $\mu^u$      (long_name='CIA Multiplier')
    omega     $\omega$     (long_name='Consumption Inequality')
    V_h       $V^h$        (long_name='BHH Welfare')
    V_u       $V^u$        (long_name='UHH Welfare')
    V         $V$          (long_name='Agg. Welfare')
;

// Shocks
varexo
       eps_T        $\varepsilon^{T}$       (long_name='Transfer Shock')
       eps_M        $\varepsilon^M$         (long_name='Money Shock')
;

// Parameters
parameters
           beta          $\beta$            (long_name='HH Discount Factor')
           rho_A         $\rho_A$           (long_name='TFP Shock Persistence')
           rho_M         $\rho_M$           (long_name='Money Growth Persistence')
           rho_T         $\rho_{T}$         (long_name='Transfer Persistence')
           phi_M         $\phi_M$           (long_name='Money Adjustment Cost Parameter')
           Gamma         $\Gamma$           (long_name='Share of Banked Households')
           T_h_ss        $\bar{T}^h$        (long_name='Steady State BHH Transfers')
           T_u_ss        $\bar{T}^u$        (long_name='Steady State UHH Transfers')
;

beta = 0.99;
rho_A = 0.85;
rho_M = 0.85;
rho_T = 0.85;
phi_M = 2;
Gamma = 0.5;
T_h_ss = 1;
T_u_ss = 1;


model;
//Households
// Banked
[name='Euler Equation of BHH']
beta*R/Pi(+1)*(1/C_h(+1)) = 1/C_h;

// Unbanked
[name='Marginal Utility of Consumption of UHH']
1/C_u = lambda_u + mu_u;

[name='Euler Equation of UHH']
beta*lambda_u(+1)/Pi(+1) - lambda_u + beta*mu_u(+1)/Pi(+1)
  - phi_M*(M - steady_state(M)) = 0;

[name='Money Adjustment Cost']
chi_M = phi_M/2*(M - steady_state(M))^2;

[name='Budget Constraint of UHH']
C_u + M + chi_M = T_u + M(-1)/Pi;

[name='CIA Constraint for UHH']
C_u = M(-1)/Pi;

// Market clearing
[name='Economy Resource Constraint']
Gamma*T_h + (1-Gamma)*T_u = C;

[name='Aggregate Consumption']
C = Gamma*C_h + (1-Gamma)*C_u;

[name='Money Growth']
M = M(-1)*xi_M/Pi;

[name='Consumption Inequality']
omega = 1 - C_u/C_h;

[name='Welfare of BHH']
V_h = log(C_h) + beta*V_h(+1);

[name='Welfare of UHH']
V_u = log(C_u) + beta*V_u(+1);

[name='Welfare of Agg. HH']
V = Gamma*V_h + (1-Gamma)*V_u;

// Exogenous processes
[name='Real Money Shock']
log(xi_M) = rho_M*log(xi_M(-1)) + eps_M;

[name='Banked HH Transfer Shock']
log(T_h) = rho_T*log(T_h(-1)) + eps_T;

[name='Unbanked HH Transfer Shock']
log(T_u) = rho_T*log(T_u(-1)) + eps_T;


end;

steady_state_model;

R = 1/beta;
Pi = beta*R;
T_u = T_u_ss;
C_u = T_u;
M = T_u;
xi_M = 1;
chi_M = 0;
lambda_u = beta;
mu_u = 1 - beta;
T_h = T_h_ss;
C_h = T_h;
C = Gamma*C_h + (1-Gamma)*C_u;
omega = 1 - C_u/C_h;
V_h = log(C_h)/(1-beta);
V_u = log(C_u)/(1-beta);
V = Gamma*V_h + (1-Gamma)*V_u;

end;

steady;
model_diagnostics;
check;

write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;

shocks;
var eps_M; stderr 0.01;
var eps_T; stderr 0.01;
end;

stoch_simul(order=1,periods=200,irf=40);
