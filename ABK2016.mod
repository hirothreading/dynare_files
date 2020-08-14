// Small Open Economy New Keynesian Model
// Replication of Aoki, Benigno, and Kiyotaki (2016, working paper)
// David Murakami

// Declare variables
// Price variables
var mc $mc$                 (long_name='marginal cost')
    Pi $\Pi$                (long_name='gross inflation')
    Z $Z$                   (long_name='rental price of capital')
    w $w$                   (long_name='real wages')
    n_int $i$               (long_name='net nominal interest rate')
    R $R$                   (long_name='gross real interest rate')
    epsilon $\epsilon$      (long_name='real exchange rate')
    Q $Q$                   (long_name='equity price')
//Quantity variables
    Y $Y$                   (long_name='output')
    M $M$                   (long_name='imports')
    L $L$                   (long_name='labour supply')
    C $C$                   (long_name='consumption')
    I $I$                   (long_name='investment')
    K $K$                   (long_name='capital stock')
    EX $EX$                 (long_name='exports')
    N $N$                   (long_name='net worth')
    K_b $K^{b}$             (long_name='bank capital')
    K_h $K^{h}$             (long_name='household capital')
    D $D$                   (long_name='domestic deposits')
    D_star $D^{*}$          (long_name='foreign deposits')
//Bank variables
    x $x$                   (long_name='fraction of assets financed by foreign borrowing')
    psi $\psi$              (long_name='Tobin Q ratio of the bank')
    phi $\phi$              (long_name='leverage multiple (Qk/n)')
    upsilon $\upsilon$      (long_name='marginal cost of deposits')
    mu $\mu$                (long_name='excess return on capital over home deposits')
    mu_Dstar $\mu^{D*}$     (long_name='cost advantage of foreign currency debt over home deposits')
//Exogenous processes
    A $A$                   (long_name='total factor productivity')
    R_star $R^{*}$          (long_name='foreign gross interest rate')
    Y_star $Y^{*}$          (long_name='foreign income')
//Other
    Lambda $\Lambda$        (long_name='stochastic discount factor between t and t+1')
    Omega $\Omega$          (long_name='stochastic discount factor of banker')
    Theta $\Theta(x_{t})$   (long_name='fraction of banker assets diverted')
    chi $\chi^{h}$          (long_name='worker extra management cost of buying equity')
    mu_star $\mu^{*}$       (long_name='fraction of mu_Dstar to mu')
;

//Specify shock process
//Policy rules
varexo  varepsilon_R $\varepsilon^{R}$        (long_name='shock process for domestic interest rate')
        varepsilon_A $\varepsilon^{A}$        (long_name='shock process for TFP')
        varepsilon_Rstar $\varepsilon^{R*}$   (long_name='shock process for foreign interest rate')
        varepsilon_Ystar $\varepsilon^{Y*}$   (long_name='shock process for foreign income')
;

//Specify parameters
//Banks
parameters  theta $\theta$ (long_name='divertable proportion of assets')
            gammma $\gamma$ (long_name='home bias in funding')
            siggma $\sigma$ (long_name='survival probability')
            xi $\xi$ (long_name='fraction of total assets brought by new banks')
//Households
            betta $\beta$ (long_name='discount rate')
            zeta $\zeta$ (long_name='inverse of Frisch elasticity of labour supply')
            zeta_0 $\zeta_{0}$ (long_name='inverse of labour supply capacity')
            varkappa $\varkappa$ (long_name='cost parameter of direct finance')
//Producers
            alphha_K $\alpha_{K}$ (long_name='cost share of capital')
            alphha_M $\alpha_{M}$ (long_name='cost share of imported intermediate goods')
            lambda $\lambda$ (long_name='one minus the depreciation rate')
            eta $\eta$ (long_name='elasticity of demand')
            omega $\omega$ (long_name='pins down kappa (slope of NKPC)')
            kappa_I $\kappa_{I}$ (long_name='cost of adjusting investment goods production')
            varphi $\varphi$ (long_name='price elasticity of export demand')
//Policy and exogenous processes
            rho_i $\rho_{i}$ (long_name='Taylor rule persistence')
            omega_pi $\omega_{\pi}$ (long_name='Taylor rule response to inflation')
            rho_A $\rho_{A}$ (long_name='TFP persistence')
            rho_Rstar $\rho_{R*}$ (long_name='foreign interest rate persistence')
            rho_Ystar $\rho_{Y*}$ (long_name='foreign income persistence')
            siggma_i $\sigma_{i}$ (long_name='standard deviation of interest rate shock')
            siggma_istar $\sigma_{i*}$ (long_name='standard deviation of foreign interest rate shock')
            siggma_A $\sigma_{A}$ (long_name='standard deviation of TFP shock')
            siggma_Ystar $\sigma_{Y*}$ (long_name='standard deviation of foreign income shock')
            omega_tauDstar $\omega_{\tau^{D*}}$ (long_name='coefficient for cyclical tax policy')
//Steady state values
            I_ss $\bar{I}$ (long_name='ss investment')
            A_ss $\bar{A}$ (long_name='ss productivity')
            R_star_ss $\bar{R}^{*}$ (long_name='ss foreign gross interest rate')
            Y_star_ss $\bar{Y}^{*}$ (long_name='ss foreign output')
;

//Parameterise
            theta = 0.475;
            gammma = 6.4;
            siggma = 0.94;
            xi = 0.000588;
            betta = 0.985;
            zeta = 0.2;
            zeta_0 = 5.89;
            varkappa = 0.000985;
            alphha_K = 0.3;
            alphha_M = 0.15;
            lambda = 0.98;
            eta = 9;
            omega = 0.66;
            kappa_I = 1;
            varphi = 1;
            rho_i = 0.85;
            omega_pi = 1.5;
            rho_A = 0.95;
            rho_Rstar = 0.95;
            rho_Ystar = 0.95;
            siggma_i = 0.01;
            siggma_istar = 0.01;
            siggma_A = 0.013;
            siggma_Ystar = 0.03;
            omega_tauDstar = 0.1;
            I_ss = 1.6;
            A_ss = 1;
            R_star_ss = (1.02)^(1/4);
            Y_star_ss = 0.412053;


model;
//Additional parameters and variables
//kappa $\kappa$ (long_name='slope of NKPC')
#kappa = (eta-1)*omega/((1-omega)*(1-betta*omega));
//R_ss $\bar{R}$ (long_name='ss gross interest rate')
#R_ss = 1/betta;
//i_ss $\bar{i}$ (long_name='ss net interest rate')
#i_ss = 1-1/betta;

[name='stochastic discount factor']
Lambda = betta*(C(-1)-zeta_0*L(-1)^(1+zeta)/(1+zeta))/(C-zeta_0*L^(1+zeta)/(1+zeta));

[name='Fisher equation']
R = (1+n_int(-1))/Pi;

[name='extra management cost of buying equity']
chi = varkappa/2*K_h^2;

[name='Stochastic discount factor of banker']
Omega = Lambda*(1-siggma+siggma*psi);

[name='fraction of assets diverted']
Theta = theta*(1+gammma/2*x^2);

[name='ratio of mu_Dstar to mu']
mu_star = mu_Dstar/mu;


//Production
[name='marginal cost, eq. (2)']
mc = 1/A*Z^alphha_K*epsilon^alphha_M*w^(1-alphha_K-alphha_M);

[name='FOC wrt P_{i,t}, eq. (3)']
(Pi-1)*Pi = 1/kappa*(eta*mc+1-eta)+Lambda(+1)*Y(+1)/Y*Pi(+1)*(Pi(+1)-1);

[name='domestic output, eq. (4)']
Y = A*(K(-1)/alphha_K)^alphha_K*(M/alphha_M)^alphha_M*(L/(1-alphha_K-alphha_M))^(1-alphha_K-alphha_M);

[name='import to capital share ratio, eq. (5)']
epsilon*M/(Z*K(-1)) = alphha_M/alphha_K;

[name='wage to capital share ratio, eq. (6)']
w*L/(Z*K(-1))=(1-alphha_K-alphha_M)/alphha_K;

[name='law of motion of capital, eq. (7)']
K = I + lambda*K(-1);

[name='exports, eq. (8)']
EX = epsilon^(varphi)*Y_star;

//Households
[name='FOC wrt labour, eq. (9)']
w = zeta_0*L^zeta;

[name='FOC wrt savings in equity,eq. (10)']
1 = Lambda(+1)*(Z(+1)+lambda*Q(+1))/(Q+varkappa*K_h);

[name='FOC wrt savings in deposits, eq. (11)']
1 = Lambda(+1)*R(+1);

[name='FOC wrt investment goods, eq. (12)']
Q = 1 + kappa_I/2*(I/I_ss-1)^2 + (I/I_ss)*kappa_I*(I/I_ss-1);


//Banks
[name='excess return on capital over home deposits, eq. (16)']
mu(-1) = Omega*((Z+lambda*Q)/Q(-1)-R);

[name='cost advantage of foreign currency debt over home deposits, eq. (17)']
mu_Dstar = Omega(+1)*(R(+1)-epsilon(+1)/epsilon*R_star(-1));

[name='marginal cost of deposit, eq. (18)']
upsilon = Omega(+1)*R(+1);

[name='bank leverage multiple, eq. (19)']
phi = upsilon/(Theta - (mu+mu_Dstar*x));

[name='Tobin Q ratio of the bank, eq. (20)']
psi = Theta*phi;

[name='fraction of assets financed by foreign borrowing, eq. (21)']
x = 1/mu_star*(-1+sqrt(1+2/gammma*mu_star^2));

//Market equilibrium
[name='output, eq. (22)']
Y = C+(1+kappa_I/2*(I/I_ss-1)^2)*I+EX+kappa/2*(Pi-1)^2*Y+varkappa/2*K_h^2;

[name='law of motion of net foreign debt, eq. (23)']
D_star = R_star(-1)*D_star(-1) + M - EX/epsilon;

[name='aggregate net worth of banks, eq. (24)']
N = (siggma+xi)*(Z+lambda*Q)*K_b(-1)-siggma*R*D(-1)-siggma*epsilon*R_star(-1)*D_star(-1);

[name='aggregate balance sheet of the bank, eq. (25)']
Q*K_b = phi*N;

[name='aggregate balance sheet of the bank, eq. (26)']
Q*K_b = N+D+epsilon*D_star;

[name='aggregate balance sheet of the bank, eq. (27)']
x = epsilon*D_star/(Q*K_b);

[name='aggregate capital, eq. (28)']
K = K_b + K_h;

[name='Taylor rule, eq. (30)']
n_int-i_ss = (1-rho_i)*omega_pi*(Pi-1) + rho_i*(n_int(-1)-i_ss) + varepsilon_R;

// Laws of motion and shock processes
[name='productivity']
log(A/A_ss) = rho_A*log(A(-1)/A_ss) + varepsilon_A;

[name='foreign income']
log(Y_star/Y_star_ss) = rho_Ystar*log(Y_star(-1)/Y_star_ss) + varepsilon_Ystar;

[name='foreign interest rate']
log(R_star/R_star_ss) = rho_Rstar*log(R_star(-1)/R_star_ss) + varepsilon_Rstar;

end;


initval;

mc =   0.888889;
Pi =   1.000000;
Z =   0.041627;
w =   4.571565;
n_int =   0.015228;
epsilon =   1.000000;
Q =   1.000000;
Y =   2.633920;
M =   0.351189;
L =   0.281675;
C =   1.913927;
I =   0.337465;
K =  16.873246;
EX =   0.362365;
N =   1.987162;
K_b =  10.474908;
K_h =   6.398339;
D =   6.252539;
D_star =   2.235207;
x =   0.213387;
psi =   3.458969;
phi =   5.271291;
upsilon =   3.311431;
mu =   0.020870;
mu_Dstar =   0.033363;
A =   1.000000;
Y_star =   0.362365;
R_star =   1.005000;
R =   1.015228;
Omega =   3.261760;
Lambda =   0.985000;
Theta =   0.544212;
mu_star =   1.598607;
chi =   0.020162;

end;


steady;
check;

shocks ;
//var varepsilon_R ; stderr siggma_i;
var varepsilon_Rstar ; stderr siggma_istar ;
//var varepsilon_A ; stderr siggma_A;
//var varepsilon_Ystar ; stderr siggma_Ystar;
end ;

stoch_simul(order=1,irf=80,periods=2000,nodisplay,loglinear);
write_latex_dynamic_model;

set(0,'DefaultAxesTitleFontWeight','normal');
figure('Name','Baseline impulse responses to foreign interest rate shock without policy');
subplot(4,4,1);
plot(Y_varepsilon_Rstar,'b');
axis tight;
title('Output');

subplot(4,4,2);
plot(C_varepsilon_Rstar,'b');
axis tight;
title('Consumption');

subplot(4,4,3);
plot(I_varepsilon_Rstar,'b');
axis tight;
title('Investment');

subplot(4,4,4);
plot(EX_varepsilon_Rstar,'b');
axis tight;
title('Exports');

subplot(4,4,5);
plot(M_varepsilon_Rstar,'b');
axis tight;
title('Imports');

subplot(4,4,6);
plot(D_star_varepsilon_Rstar,'b');
axis tight;
title('Net foreign debt');

subplot(4,4,7);
plot(epsilon_varepsilon_Rstar,'b');
axis tight;
title('Real exchange rate');

subplot(4,4,8);
plot(Q_varepsilon_Rstar,'b');
axis tight;
title('Capital price');

subplot(4,4,9);
plot(N_varepsilon_Rstar,'b');
axis tight;
title('Net worth');

subplot(4,4,10);
plot(Pi_varepsilon_Rstar,'b');
axis tight;
title('Inflation (gross)');

subplot(4,4,11);
plot(n_int_varepsilon_Rstar,'b');
axis tight;
title('Nominal interest (net)');

subplot(4,4,12);
plot(R_star_varepsilon_Rstar,'b');
axis tight;
title('Foreign interest (gross)');
supersizeme(0.6);
