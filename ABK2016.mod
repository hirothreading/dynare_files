%Small Open Economy New Keynesian Model
%Replication of Aoki, Benigno, and Kiyotaki (2016, Oct 2018 wp)
%Author: David Murakami (University of Oxford)

//Choose to run baseline, permanent policy, or cyclical policy simulation
@#define baseline_tax = 1
@#define permanent_policy_rule = 0
@#define cyclical_policy_rule = 0

// Declare variables
// Price variables
var mc        $mc$              (long_name='marginal cost')
    Pi        $\Pi$             (long_name='gross inflation')
    Z         $Z$               (long_name='rental price of capital')
    w         $w$               (long_name='real wages')
    R         $R$               (long_name='gross real interest rate')
    epsilon   $\epsilon$        (long_name='real exchange rate')
    Q         $Q$               (long_name='equity price')
//Quantity variables
    Y         $Y$               (long_name='output')
    M         $M$               (long_name='imports')
    L         $L$               (long_name='labour supply')
    C         $C$               (long_name='consumption')
    I         $I$               (long_name='investment')
    K         $K$               (long_name='capital stock')
    EX        $EX$              (long_name='exports')
    N         $N$               (long_name='net worth')
    K_b       $K^{b}$           (long_name='bank capital')
    K_h       $K^{h}$           (long_name='household capital')
    D         $D$               (long_name='domestic deposits')
    D_star    $D^{*}$           (long_name='foreign deposits')
//Bank variables
    x         $x$               (long_name='fraction of assets financed by foreign borrowing')
    psi       $\psi$            (long_name='Tobin Q ratio of the bank')
    phi       $\phi$            (long_name='leverage multiple (Qk/n)')
    upsilon   $\upsilon$        (long_name='marginal cost of deposits')
    mu        $\mu$             (long_name='excess return on capital over home deposits')
    mu_star   $\mu^{*}$         (long_name='cost advantage of foreign currency debt over home deposits')
//Exogenous processes
    A         $A$               (long_name='total factor productivity')
    R_star    $R^{*}$           (long_name='foreign gross interest rate')
    Y_star    $Y^{*}$           (long_name='foreign income')
//Other
    //Lambda    $\Lambda$         (long_name='stochastic discount factor between t and t+1')
    Phi       $\Phi$            (long_name='investment cost of adjustment')
    //Omega     $\Omega$          (long_name='stochastic discount factor of banker')
    Theta     $\Theta(x_{t})$   (long_name='fraction of banker assets diverted')
    chi_h     $\chi^{h}$        (long_name='worker extra management cost of buying equity')
    chi_b     $\chi^{b}$        (long_name='cost of borrowing from foreigners')
    Ynet      $Y^{net}$         (long_name='net output')
@#if baseline_tax==0
    tau_N     $\tau^{N}$        (long_name='subsidy rate on net worth')
    tau_K     $\tau^{K}$        (long_name='tax rate on risky asset holdings')
    tau_Dstar $\tau^{D*}$       (long_name='tax rates on foreign debt')
@#endif
@#if permanent_policy_rule==1
    W         $W$               (long_name='welfare')
@#endif
;

//Specify shock process
//Policy rules
varexo  varepsilon_R     $\varepsilon^{R}$    (long_name='shock process for domestic interest rate')
        varepsilon_A     $\varepsilon^{A}$    (long_name='shock process for TFP')
        varepsilon_Rstar $\varepsilon^{R*}$   (long_name='shock process for foreign interest rate')
        varepsilon_Ystar $\varepsilon^{Y*}$   (long_name='shock process for foreign income')
@#if permanent_policy_rule==1
        tau_K_bar        $\bar{\tau}^{K}$     (long_name='tax rate on risky asset holdings')
        tau_Dstar_bar    $\bar{\tau}^{D*}$    (long_name='tax rate on foreign borrowings')
@#endif
;

//Specify parameters
//Banks
parameters  theta $\theta$                      (long_name='elasticity of leverage wrt foreign borrowing')
            theta_0 $\theta_{0}$                (long_name='home bias in funding')
            sigma $\sigma$                      (long_name='survival probability')
            xi $\xi$                            (long_name='fraction of total assets brought by new banks')
            varkappa_b $\varkappa^{b}$          (long_name='management cost for foreign borrowing')
//Households
            betta $\beta$                       (long_name='discount rate')
            zeta $\zeta$                        (long_name='inverse of Frisch elasticity of labour supply')
            zeta_0 $\zeta_{0}$                  (long_name='inverse of labour supply capacity')
            varkappa_h $\varkappa_{h}$          (long_name='cost parameter of direct finance')
//Producers
            alphha_K $\alpha_{K}$               (long_name='cost share of capital')
            alphha_M $\alpha_{M}$               (long_name='cost share of imported intermediate goods')
            lambda $\lambda$                    (long_name='one minus the depreciation rate')
            eta $\eta$                          (long_name='elasticity of demand')
            omega $\omega$                      (long_name='pins down kappa (slope of NKPC)')
            kappa_I $\kappa_{I}$                (long_name='cost of adjusting investment goods production')
            varphi $\varphi$                    (long_name='price elasticity of export demand')
            kappa $\kappa$                      (long_name='slope of NKPC')
//Policy and exogenous processes
            A_ss $\bar{A}$                      (long_name='steady state productivity')
            rho_i $\rho_{i}$                    (long_name='Taylor rule persistence')
            omega_pi $\omega_{\pi}$             (long_name='Taylor rule response to inflation')
            rho_A $\rho_{A}$                    (long_name='TFP persistence')
            rho_Rstar $\rho_{R*}$               (long_name='foreign interest rate persistence')
            rho_Ystar $\rho_{Y*}$               (long_name='foreign income persistence')
            sigma_i $\sigma_{i}$                (long_name='standard deviation of interest rate shock')
            sigma_istar $\sigma_{i*}$           (long_name='standard deviation of foreign interest rate shock')
            sigma_A $\sigma_{A}$                (long_name='standard deviation of TFP shock')
            sigma_Ystar $\sigma_{Y*}$           (long_name='standard deviation of foreign income shock')
            omega_tauDstar $\omega_{\tau^{D*}}$ (long_name='coefficient for cyclical tax policy')
;

//Parameterise
            theta = 0.1;
            theta_0 = 0.401; //0.399;
            sigma = 0.94;
            xi = 0.0045; //0.0046;
            varkappa_b = 0.0197; //0.0219;
            betta = 1/1.015;
            zeta = 1/3;
            zeta_0 = 7.883;
            varkappa_h = 0.0197;
            alphha_K = 0.3;
            alphha_M = 0.18;
            lambda = 0.98;
            eta = 9;
            omega = 2/3;
            kappa = (eta-1)*omega/((1-omega)*(1-betta*omega));
            kappa_I = 2/3;
            varphi = 1;
            rho_i = 0.8;
            omega_pi = 1.5;
            rho_A = 0.9;
            rho_Rstar = 0.9;
            rho_Ystar = 0.9;
            sigma_i = 0.01; //need to check the scale of these
@#if permanent_policy_rule==1
            sigma_i = 0.005;
@#endif
            sigma_istar = 0.01;
            sigma_A = 0.013;
            sigma_Ystar = 0.03;
            omega_tauDstar = 0.05;
            A_ss = 1;


model;
//Additional parameters and variables
//R_ss $\bar{R}$ (long_name='ss gross interest rate')
#R_ss = 1/betta;

//R_star_ss $\bar{R}^{*}$ (long_name='ss foreign interest rate')
#R_star_ss = (steady_state(R_star));

//Y_star_ss $\bar{Y}^{*}$ (long_name='ss foreign GDP')
#Y_star_ss = (steady_state(Y_star));

//I_ss $\bar{I}$ (long_name='ss investment')
#I_ss = (steady_state(I));

//[name='stochastic discount factor']
#Lambda = betta*(C-zeta_0*L^(1+zeta)/(1+zeta))/(C(+1)-zeta_0*L(+1)^(1+zeta)/(1+zeta));

//[name='Stochastic discount factor of banker']
#Omega = Lambda*(1-sigma+sigma*psi(+1));

//Welfare measure for second-order approximation
@#if permanent_policy_rule==1
[name='Welfare']
W = log(C - zeta_0*(1+zeta)*L^(1+zeta)) + betta*W(+1);
@#endif

[name='investment cost of adjustment']
Phi = kappa_I/2*(I/I_ss-1)^2;

[name='extra management cost of buying equity']
chi_h = (varkappa_h/2)*K_h^2/K;

[name='cost of borrowing from foreigners']
chi_b = (varkappa_b/2)*Q*K_b*x^2;

[name='fraction of assets diverted']
Theta = theta_0*exp(-theta*x);


//Production
[name='marginal cost, eq. (2)']
mc = (1/A)*Z^alphha_K*epsilon^alphha_M*w^(1-alphha_K-alphha_M);

[name='FOC wrt P_{i,t}, eq. (3)']
(Pi-1)*Pi = 1/kappa*(eta*mc+1-eta)+Lambda*Y(+1)/Y*Pi(+1)*(Pi(+1)-1);

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

[name='FOC wrt savings in equity, eq. (10)']
1 = Lambda*(Z(+1)+lambda*Q(+1))/(Q+varkappa_h*K_h/K);

[name='FOC wrt savings in deposits, eq. (11)']
1 = Lambda*R/Pi(+1);

[name='FOC wrt investment goods, eq. (12)']
Q = 1 + Phi + (I/I_ss)*kappa_I*(I/I_ss-1);


//Banks
@#if baseline_tax==0
[name='tax on net worth, eq. (13)']
tau_N*N = tau_K*Q*K_b + tau_Dstar*epsilon*D_star;
//In the baseline case, there are no taxes
@#endif

@#if baseline_tax==1
[name='excess return on capital over home deposits, eq. (17)']
mu = Omega*((Z(+1)+lambda*Q(+1))/Q-R/Pi(+1));
@#else
[name='excess return on capital over home deposits, eq. (17)']
mu = Omega*((Z(+1)+lambda*Q(+1))/Q-(1+tau_K)*R/Pi(+1));
@#endif

@#if baseline_tax==1
[name='cost advantage of foreign currency debt over home deposits, eq. (18)']
mu_star = Omega*(R/Pi(+1)-epsilon(+1)/epsilon*R_star);
@#else
[name='cost advantage of foreign currency debt over home deposits, eq. (18)']
mu_star = Omega*((1-tau_Dstar)*R(+1)-epsilon(+1)/epsilon*R_star);
@#endif
//There is a typo in the ABK paper: the timing for R_star in the paper is -1.

[name='marginal cost of deposit, eq. (19)']
upsilon = Omega*R/Pi(+1);

@#if baseline_tax==1
[name='bank leverage multiple, eq. (20)']
phi = upsilon/(Theta + (varkappa_b/2)*(x^2)*upsilon - (mu + mu_star*x));
@#else
[name='bank leverage multiple, eq. (20)']
phi = (1+tau_N)*upsilon/(Theta + (varkappa_b/2)*(x^2)*upsilon - (mu + mu_star*x));
@#endif

[name='Tobin Q ratio of the bank, eq. (21)']
psi = Theta*phi;

[name='fraction of assets financed by foreign borrowing, eq. (22)']
x = mu_star/(varkappa_b*upsilon) - 1/theta + sqrt((mu_star/(varkappa_b*upsilon))^2 + (1/theta)^2 + 2*mu/(varkappa_b*upsilon));


//Market equilibrium
[name='output, eq. (23)']
Y = C + (1+Phi)*I + EX + (kappa/2)*(Pi-1)^2*Y + chi_h + chi_b;

[name='law of motion of net foreign debt, eq. (24)']
D_star = R_star(-1)*D_star(-1) + M - EX/epsilon;

[name='aggregate net worth of banks, eq. (25)']
N = sigma*((Z+lambda*Q)*K_b(-1)-D(-1)*R(-1)/Pi-epsilon*R_star(-1)*D_star(-1))+xi*(Z+lambda*Q)*K(-1);
//Using this equation seems to produce a smoother IRF in net foreign debt, closer to the paper
//N = (sigma+xi)*(Z+lambda*Q)*K_b(-1) - sigma*D(-1)*R(-1)/Pi - sigma*epsilon*R_star(-1)*D_star(-1);
//ABK use this expression in set of more recent presentation slides. The main difference is K_b

[name='aggregate balance sheet of the bank, eq. (26)']
Q*K_b*(1+varkappa_b/2*x^2) = (1+varkappa_b/2*x^2)*phi*N;

[name='aggregate balance sheet of the bank, eq. (27)']
Q*K_b*(1+varkappa_b/2*x^2) = N + D + epsilon*D_star;

[name='aggregate balance sheet of the bank, eq. (28)']
x = epsilon*D_star/(Q*K_b);

[name='aggregate capital, eq. (29)']
K = K_b + K_h;

[name='Taylor rule, eq. (30)']
R-1-(R_ss-1)=(1-rho_i)*omega_pi*(Pi-1) + rho_i*(R(-1)-1-(R_ss-1)) + varepsilon_R;

@#if permanent_policy_rule==1
[name='tax rate on foreign debt']
tau_Dstar = tau_Dstar(-1) + tau_Dstar_bar;
@#endif
@#if cyclical_policy_rule==1
[name='tax rate on foreign debt, eq. (31)']
tau_Dstar = omega_tauDstar*(log(K_b(-1))-log(K_b_ss));
@#endif

@#if permanent_policy_rule==1
[name='tax rate on risky assets of banks']
tau_K = tau_K(-1) + tau_K_bar;
@#endif
@#if cyclical_policy_rule==1
[name='tax rate on risky assets of banks']
tau_K = 0;
@#endif


// Laws of motion and shock processes
[name='productivity']
log(A/A_ss) = rho_A*log(A(-1)/A_ss) + varepsilon_A;

[name='foreign income']
log(Y_star/Y_star_ss) = rho_Ystar*log(Y_star(-1)/Y_star_ss) + varepsilon_Ystar;

[name='foreign interest rate']
log(R_star/R_star_ss) = rho_Rstar*log(R_star(-1)/R_star_ss) + varepsilon_Rstar;


//Observed variables
[name='Net output, pg. 21']
Ynet = Y - epsilon*M - kappa/2*(Pi-1)^2*Y - chi_h - chi_b;
end;


initval;
//Steady state values
//Price variables
mc = 0.888888889;
Pi = 1.000000000;
Z = 0.045331866;
w = 4.750961764;
R = 1.015000000;
epsilon = 1.000000000;
Q = 1.000000000;
//Quantity variables:
Y = 2.250094407;
M = 0.360015105;
L = 0.218912230;
C = 1.565202470;
I = 0.264725557;
K = 13.236277844;
EX = 0.379506937;
N = 2.762062537;
K_b = 6.396966472;
K_h = 6.839311372;
D = 1.691570865;
D_star = 1.949183226;
//Bank variables:
x = 0.304704305;
psi = 1.681936952;
phi = 2.316010730;
upsilon = 1.641020735;
mu = 0.016704243;
mu_star = 0.008083846;
//Exogenous process variables:
A = 1.000000000;
R_star = 1.010000000;
Y_star = 0.379506937;
//Other variables:
Phi = 0.000000000;
Theta = 0.388965635;
chi_h = 0.034809285;
chi_b = 0.005850157;
Ynet = 1.849419860;
@#if permanent_policy_rule==1
    tau_N = 0.000000000;
    tau_K = 0.000000000;
    tau_Dstar = 0.000000000;
    W = -116.609429144;
@#endif
end ;


//BASELINE CASE
@#if baseline_tax==1

write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;

steady ;
check ;


shocks ;
//var varepsilon_R ; stderr sigma_i;
var varepsilon_Rstar ; stderr sigma_istar ;
//var varepsilon_A ; stderr sigma_A;
//var varepsilon_Ystar ; stderr sigma_Ystar;
end ;

options_.TeX = 1;
write_latex_dynamic_model;
write_latex_parameter_table;

stoch_simul(order=1,irf=80,nodisplay) ; //Ynet C I EX M D_star epsilon Q N Pi R R_star ;


//Make plots
%%%Plots require supersizeme.m. Disable supersizeme() if package not installed.
%%%Plots below are only for the foreign interest rate shock.
set(0,'DefaultAxesTitleFontWeight','normal');
figure('Name','Baseline impulse responses to foreign interest rate shock without policy');
subplot(4,4,1);
plot(Ynet_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Net output');

subplot(4,4,2);
plot(C_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Consumption');

subplot(4,4,3);
plot(I_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Investment');

subplot(4,4,4);
plot(EX_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Exports');

subplot(4,4,5);
plot(M_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Imports');

subplot(4,4,6);
plot(D_star_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Net foreign debt');

subplot(4,4,7);
plot(epsilon_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Real exchange rate');

subplot(4,4,8);
plot(Q_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Capital price');

subplot(4,4,9);
plot(N_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Net worth');

subplot(4,4,10);
plot(Pi_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Inflation');

subplot(4,4,11);
plot(R_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Nominal interest (gross)');

subplot(4,4,12);
plot(R_star_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Foreign interest rate');
supersizeme(0.6);

@#endif


//PERMANENT POLICY EXPERIMENT
@#if permanent_policy_rule==1

steady ;
check ;

shocks ;
var varepsilon_R ; stderr sigma_i;
var varepsilon_Rstar ; stderr sigma_istar ;
var varepsilon_A ; stderr sigma_A;
var varepsilon_Ystar ; stderr sigma_Ystar;
//TAX RATES
//var tau_Dstar_bar;
//var tau_K_bar;
end ;

stoch_simul(order=2,pruning,irf=0) W Q epsilon N;

mc_pos    =strmatch('mc',M_.endo_names,'exact');
pi_pos    =strmatch('Pi',M_.endo_names,'exact');
z_pos     =strmatch('Z',M_.endo_names,'exact');
w_pos     =strmatch('w',M_.endo_names,'exact');
r_pos     =strmatch('R',M_.endo_names,'exact');
eps_pos   =strmatch('epsilon',M_.endo_names,'exact');
q_pos     =strmatch('Q',M_.endo_names,'exact');
taun_pos  =strmatch('tau_N',M_.endo_names,'exact');
y_pos     =strmatch('Y',M_.endo_names,'exact');
m_pos     =strmatch('M',M_.endo_names,'exact');
l_pos     =strmatch('L',M_.endo_names,'exact');
c_pos     =strmatch('C',M_.endo_names,'exact');
i_pos     =strmatch('I',M_.endo_names,'exact');
k_pos     =strmatch('K',M_.endo_names,'exact');
ex_pos    =strmatch('EX',M_.endo_names,'exact');
n_pos     =strmatch('N',M_.endo_names,'exact');
kb_pos    =strmatch('K_b',M_.endo_names,'exact');
kh_pos    =strmatch('K_h',M_.endo_names,'exact');
d_pos     =strmatch('D',M_.endo_names,'exact');
dstar_pos =strmatch('D_star',M_.endo_names,'exact');
x_pos     =strmatch('x',M_.endo_names,'exact');
psi_pos   =strmatch('psi',M_.endo_names,'exact');
phi_pos   =strmatch('phi',M_.endo_names,'exact');
ups_pos   =strmatch('upsilon',M_.endo_names,'exact');
mu_pos    =strmatch('mu',M_.endo_names,'exact');
mustar_pos=strmatch('mu_Dstar',M_.endo_names,'exact');
a_pos     =strmatch('A',M_.endo_names,'exact');
rstar_pos =strmatch('R_star',M_.endo_names,'exact');
ystar_pos =strmatch('Y_star',M_.endo_names,'exact');
lambda_pos=strmatch('Lambda',M_.endo_names,'exact');
phix_pos  =strmatch('Phi',M_.endo_names,'exact');
omega_pos =strmatch('Omega',M_.endo_names,'exact');
theta_pos =strmatch('Theta',M_.endo_names,'exact');
chi_h_pos =strmatch('chi_h',M_.endo_names,'exact');
chi_b_pos =strmatch('chi_b',M_.endo_names,'exact');
tauk_pos  =strmatch('tau_K',M_.endo_names,'exact');
taud_pos  =strmatch('tau_Dstar',M_.endo_names,'exact');
ynet_pos  =strmatch('Ynet',M_.endo_names,'exact');
wel_pos   =strmatch('W',M_.endo_names,'exact');


IRF_periods=200;
burnin=5000; %periods for convergence

shock_mat_with_zeros=zeros(burnin,M_.exo_nbr); %shocks set to 0 to simulate without uncertainty
out_noshock=simult_(oo_.dr.ys,oo_.dr,shock_mat_with_zeros,options_.order); %simulate series
log_deviations_SS_noshock=out_noshock-oo_.dr.ys*ones(1,burnin+M_.maximum_lag); %subtract steady state to get deviations from steady state
ergodicmean_no_shocks=out_noshock(:,end); %ergodic mean absent of shocks (EMAS) is the final product

%%% TAX POLICY SHOCKS
shock_mat = zeros(IRF_periods,M_.exo_nbr);
shock_mat(25,strmatch('tau_Dstar_bar',M_.exo_det_names,'exact'))= 0.0001;

sim_mat = simult_(ergodicmean_no_shocks,oo_.dr,shock_mat,options_.order);

%sim_mat_percent_from_SSS = (sim_mat(25+burnin+1:25+burnin+IRF_periods,:)-IRF_no_shock_mat(25+burnin+1:25+burnin+IRF_periods,:))./repmat(stochastic_steady_state,IRF_periods,1); %only valid for variables not yet logged
%IRF_mat_percent_from_SSS = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,:)-IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,:))./repmat(stochastic_steady_state,IRF_periods,1); %only valid for variables not yet logged

/*
//scale IRFs as required
ynet_vola_IRF   =100*sim_mat_percent_from_SSS(:,ynet_pos);
c_vola_IRF      =100*sim_mat_percent_from_SSS(:,c_pos);
i_vola_IRF      =100*sim_mat_percent_from_SSS(:,i_pos);
ex_vola_IRF     =100*sim_mat_percent_from_SSS(:,ex_pos);
m_vola_IRF      =100*sim_mat_percent_from_SSS(:,m_pos);
dstar_vola_IRF  =100*sim_mat_percent_from_SSS(:,dstar_pos);
eps_vola_IRF    =100*sim_mat_percent_from_SSS(:,eps_pos);
q_vola_IRF      =100*sim_mat_percent_from_SSS(:,q_pos);
n_vola_IRF      =100*sim_mat_percent_from_SSS(:,n_pos);
pi_vola_pos     =400*sim_mat_percent_from_SSS(:,pi_pos);
r_vola_IRF      =400*sim_mat_percent_from_SSS(:,r_pos);
rstar_vola_IRF  =400*sim_mat_percent_from_SSS(:,rstar_pos);
wel_vola_IRF    =100*sim_mat_percent_from_SSS(:,wel_pos);

hh=figure;
figure(hh)
subplot(2,2,1)
hold on
plot(wel_vola_IRF,'b-','LineWidth',1)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off');
xlim([1 IRF_periods]);
set(gca,'XTick',[50:50:IRF_periods],'FontSize',12);
title('Welfare','FontSize',14);
ylabel('Percent','FontSize',12);

figure(hh)
subplot(2,2,2)
hold on
plot(q_vola_IRF,'b-','LineWidth',1)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off');
xlim([1 IRF_periods]);
set(gca,'XTick',[50:50:IRF_periods],'FontSize',12);
title('Capital price','FontSize',14);
ylabel('Percent','FontSize',12);
%ylim([-0.3 0.1]);set(gca,'YTick',[-0.3:0.1:0.1],'FontSize',12);

figure(hh)
subplot(2,2,3)
hold on
plot(eps_vola_IRF,'b-','LineWidth',1)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off');
xlim([1 IRF_periods]);
set(gca,'XTick',[50:50:IRF_periods],'FontSize',12);
title('Exchange rate','FontSize',14);
ylabel('Percent','FontSize',12);
%ylim([-0.6 0.4]);set(gca,'YTick',[-0.6:0.2:0.4],'FontSize',12);

figure(hh)
subplot(2,2,4)
hold on
plot(n_vola_IRF,'b-','LineWidth',1)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off');
xlim([1 IRF_periods]);
set(gca,'XTick',[50:50:IRF_periods],'FontSize',12);
title('Net worth','FontSize',14);
ylabel('Percent','FontSize',12);

*/
@#endif


//CYCLICAL TAX POLICY
@#if cyclical_policy_rule==1

steady ;
check ;

shocks ;
//var varepsilon_R ; stderr sigma_i;
var varepsilon_Rstar ; stderr sigma_istar ;
//var varepsilon_A ; stderr sigma_A;
//var varepsilon_Ystar ; stderr sigma_Ystar;
end ;

options_.TeX = 1;
write_latex_dynamic_model;
write_latex_parameter_table;

stoch_simul(order=1,nodisplay) ; //Ynet C I EX M D_star epsilon Q N Pi R R_star ;


//Make plots
%%%Plots require supersizeme.m. Disable supersizeme() if package not installed.
%%%Plots below are only for the foreign interest rate shock.
set(0,'DefaultAxesTitleFontWeight','normal');
figure('Name','Impulse responses to foreign interest rate shock with cyclical tax policy');
subplot(4,4,1);
plot(Ynet_varepsilon_Rstar,'b-');
axis tight;
title('Net output');

subplot(4,4,2);
plot(C_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Consumption');

subplot(4,4,3);
plot(I_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Investment');

subplot(4,4,4);
plot(EX_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Exports');

subplot(4,4,5);
plot(M_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Imports');

subplot(4,4,6);
plot(D_star_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Net foreign debt');

subplot(4,4,7);
plot(epsilon_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Real exchange rate');

subplot(4,4,8);
plot(Q_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Capital price');

subplot(4,4,9);
plot(N_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Net worth');

subplot(4,4,10);
plot(Pi_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Inflation');

subplot(4,4,11);
plot(R_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Nominal interest (gross)');

subplot(4,4,12);
plot(R_star_varepsilon_Rstar,'b-','LineWidth',1);
axis tight;
title('Foreign interest rate');
supersizeme(0.6);

@#endif
