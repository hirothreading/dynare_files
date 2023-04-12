// RBC Model with Financial Frictions
// A simple RBC model with banks as in Gertler-Kiyotaki/Gertler-Karadi
// By David Murakami

// DECLARE VARIABLES
var
// QUANTITIES
Y           $Y$                 (long_name='Output')
C           $C$                 (long_name='Consumption')
D           $D$                 (long_name='Deposits')
I           $I$                 (long_name='Investment')
L           $L$                 (long_name='Labour Supply')
K           $K$                 (long_name='Capital')
K_b         $K^b$               (long_name='Bank Equity')
K_h         $K^h$               (long_name='Household Equity')
N           $N$                 (long_name='Bank Net Worth')
// PRICES
Q           $Q$                 (long_name='Equity Price')
R           $R$                 (long_name='Interest Rate')
W           $W$                 (long_name='Wages')
Z           $Z$                 (long_name='Capital Rental Rate')
LAMBDA      $\Lambda$           (long_name='Household SDF')
// BANK VARIABLES
PSI         $\psi$              (long_name='Tobin Q')
PHI         $\phi$              (long_name='Leverage Ratio')
MU          $\mu$               (long_name='Excess Return of Equity')
UPSILON     $\upsilon$          (long_name='Marginal Cost of Deposits')
OMEGA       $\Omega$            (long_name='Banker SDF')
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
varkappa_h   $\varkappa^h$      (long_name='Equity purchase cost')
lambda       $\lambda$          (long_name='Undepreciated capital')
kappa_I      $\kappa_I$         (long_name='Investment adjustment cost')
sigma        $\sigma$           (long_name='Banker survival probability')
theta        $\theta$           (long_name='Absconding proportion')
gamma        $\gamma$           (long_name='Startup funds proportion')
rho_a        $\rho_a$           (long_name='TFP shock persistence')
sigma_a      $\sigma_a$         (long_name='TFP shock standard deviation')
;
// PARAMETERISE
betta = 0.99;
alphha = 0.3;
nu = 2;
delta = 0.025;
chi = 2;
varkappa_h = 0.0197;
lambda = 1 - delta;
kappa_I = 2/3;
sigma = 0.98;
theta = 1/2;
gamma = 0.0025; 
rho_a = 0.95;
sigma_a = 0.01;

// DECLARE MODEL
model;
// Note: Period utility is: ln(C) - chi*L^(1+1/nu)/(1+1/nu)
[name='Intertemporal Euler equation']
W = chi*L^(1/nu)*C;

[name='Euler equation wrt equity']
1 = LAMBDA*(Z(+1) + lambda*Q(+1))/(Q + varkappa_h*K_h/K);

[name='Euler equation on deposits']
1 = LAMBDA*R;

[name='Household SDF']
LAMBDA = betta*C/C(+1);

[name='Production']
Y = A*K(-1)^alphha*L^(1-alphha);

[name='Wages']
W = (1-alphha)*Y/L;

[name='Capital rental rate']
Z = alphha*Y/K(-1);

[name='Law of motion of capital']
K = lambda*K(-1) + I;

#I_ss = steady_state(I);
[name='FOC wrt investment goods']
Q = 1 + kappa_I/2*(I/I_ss-1)^2 + (I/I_ss)*kappa_I*(I/I_ss - 1);

[name='Tobin Q']
PSI = theta*PHI;

[name='Max leverage ratio']
PHI = UPSILON/(theta - MU);

[name='Excess return of equity']
MU = OMEGA*((Z(+1) + lambda*Q(+1))/Q - R);

[name='Marginal cost of deposits']
UPSILON = OMEGA*R;

[name='Banker SDF']
OMEGA = LAMBDA*(1 - sigma + sigma*PSI(+1));

[name='Aggregate resource constraint']
Y = C + I + (varkappa_h/2)*(K_h/K)^2*K;

[name='Aggregate capital']
K = K_h + K_b;

[name='Aggregate net worth of banks']
N = sigma*((Z + lambda*Q)*K_b(-1) - R(-1)*D(-1)) + gamma*(Z + lambda*Q)*K(-1);

[name='Aggregate leverage']
PHI = Q*K_b/N;

[name='Aggregate balance sheet constraint']
Q*K_b = D + N;

[name='TFP process']
ln(A) = rho_a*log(A(-1)) + eps_a;

end;

// PROVIDE ANALYTICAL STEADY STATE EXPRESSIONS
steady_state_model;
A = 1;
R = 1/betta;
Q = 1;
LAMBDA = betta;
s = s_find(sigma,betta,gamma,theta,varkappa_h);
PHI = (betta - sigma)/(sigma*s + gamma*(1+s)/(1-s/varkappa_h));
PSI = theta*PHI;
OMEGA = LAMBDA*(1 - sigma + sigma*PSI);
KhonK = s/varkappa_h;
Z = R*(1+s) - lambda;
W = (1-alphha)*(alphha/Z)^(alphha/(1-alphha));
KonY = alphha/Z;
L = ((1-alphha)/(chi*(1 - delta*KonY - s^2/(2*varkappa_h)*KonY)))^(nu/(1+nu));
K = alphha/(1-alphha)*W*L/Z;
C = W/(chi*L^(1/nu));
Y = K^alphha*L^(1-alphha);
K_h = KhonK*K;
K_b = K - K_h;
I = delta*K;
UPSILON = OMEGA*R;
MU = s*UPSILON;
N = Q*K_b/PHI;
D = Q*K_b - N;
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
