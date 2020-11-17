% New Keynesian Life-Cycle Model
% with social security, variable labour, and small open economy features
% Author: David Murakami (Keio University)

//Set experiment to run
// Experiment 1: A decline in the population growth rate combined with an increase in life expectancy.
// Experiment 2: An increase in the effective retirement age from 68 to 75.
// Experiment 3: A decrease in the productivity gap between young and old workers.
// Experiment 4: A combination of experiments 1 and 2.
// Experiment 5: A combination of experiments 1, 2, and 3.
// Experiment 6: A combination of experiments 1 and 3.

@#define EXP1 = 1
@#define EXP2 = 0
@#define EXP3 = 0
@#define EXP4 = 0
@#define EXP5 = 0
@#define EXP6 = 0

//Declare variables
//Production and output
var Y $Y$ (long_name='GDP')
    W $W$ (long_name='wages')
    RK $R^{K}$ (long_name='capital rental rate')
    MC $\varphi$ (long_name='marginal cost')
    I $I$ (long_name='investment')
    K $K$ (long_name='capital')
    L $L$ (long_name='aggregate labour supply')
//Intermediate goods manufactures
    PINT $P^{I}$ (long_name='real market value of intermediate firms')
    DI $D^{I}$ (long_name='profits from intermediate goods producers')
//Households and aggregation
    PSI $\psi$ (long_name='dependency ratio')
    LAMBDA $\lambda$ (long_name='ratio of financial assets held by the old')
    OMEGA $\Omega$ (long_name='young-old adjustment factor')
    XI $\Xi$ (long_name='Ratio of MPCs')
    V $V$ (long_name='aggregate welfare')
    C $C$ (long_name='aggregate consumption')
    H $H$ (long_name='human assets')
    A $A$ (long_name='aggregate asset holdings')
//The Elderly
    VR $V^{r}$ (long_name='welfare of the old')
    MPCR $\xi^{r}$ (long_name='MPC of old')
    CR $C^{r}$ (long_name='aggregate consumption by the old')
    AO $A^{r}$ (long_name='financial assets held by the old')
    SR $S^{r}$ (long_name='social security benefits to the old')
    LR $L^{r}$ (long_name='labour supply by the old')
    HR $H^{r}$ (long_name='human assets of the old')
//The Young
    VW $V^{w}$ (long_name='welfare of the young')
    MPCW $\xi^{w}$ (long_name='MPC of young')
    CW $C^{w}$ (long_name='aggregate consumption by the young')
    AW $A^{w}$ (long_name='financial assets held by young')
    SW $S^{w}$ (long_name='present discounted social security benefits to workers')
    LW $L^{w}$ (long_name='labour supply by the young')
    HW $H^{w}$ (long_name='human assets of the young')
//Fiscal and Monetary Policy
    B $B$ (long_name='domestic government bonds')
    G $G$ (long_name='government expenditure')
    E $E$ (long_name='government pension expenses')
    T $T$ (long_name='total government lump sum tax revenues')
    R $R$ (long_name='nominal interest rate')
    PI $\pi$ (long_name='domestic inflation rate')
//Open Economy
    NX $NX$ (long_name='net exports')
    PHI $\Phi$ (long_name='interest rate risk premium')
    FSTAR $F^{*}$ (long_name='net foreign assets')
//Shock process for ZLB disaster shock
    ZK $z^{K}$ (long_name='capital/investment shock process')
;

//Specify shock process
varexo    EK $e_{K}$ (long_name='capital disaster shock')
          omega $\omega$ (long_name='probability of remaining as worker')
          gamma $\gamma$ (long_name='probability of death')
          n ${n}$ (long_name='population growth rate')
          varsigma $\varsigma$ (long_name='relative productivity of the old')
;

//Specify parameters
parameters x $x$ (long_name='technology growth rate')
           alpha $\alpha$ (long_name='labour share of output')
           delta $\delta$ (long_name='depreciation rate')
           kappa $\kappa$ (long_name='elasticity of substitution amongst goods')
           phiI $\phi_{I}$ (long_name='cost adjustment parameter for intermediate firms')
           rhoK $\rho_{K}$ (long_name='persistence of AR(1) capital/investment shock')
           sigma $\sigma$ (long_name='IES')
           varsigma_inital $\varsigma_{0}$ (long_name='initial relative productivity')
           varsigma_final $\varsigma_{T}$ (long_name='final relative productivity')
           upsilon $\upsilon$ (long_name='utility weight on consumption')
           g $g$ (long_name='government spending as percentage of GDP')
           b $b$ (long_name='government debt-GDP ratio')
           phiPI $\phi_{\pi}$ (long_name='Taylor rule response to inflation')
           phiR $\phi_{R}$ (long_name='Taylor rule relative rate to natural rate of interest')
           beta $\beta$ (long_name='discount factor')
           varrho $\varrho$ (long_name='replacement rate')
           rbar $\bar{R}$ (long_name='natural rate of interest')
           phiNX $\phi_{NX}$ (long_name='net export demand')
           ystar $Y^{*}$ (long_name='foreign output')
           Rstar $R^{*}$ (long_name='foreign interest rate')
           phiFSTAR $\phi_{\bar{F^{*}}}$ (long_name='elasticity of foreign debt')
           fstarbar $\bar{F}^{*}$ (long_name='natural level of foreign debt')
           omega_initial $\omega_{0}$ (long_name='initial probability of remaining young')
           omega_final $\omega_{{T}}$ (long_name='final probability of remaining young')
           gamma_initial $\gamma_{0}$ (long_name='initial probability of survival')
           gamma_final $\gamma_{{T}}$ (long_name='final probability of survival')
           n_initial ${n}_{0}$ (long_name='initial population growth rate')
           n_final ${n}_{{T}}$ (long_name='final population growth rate')
 ;
//Parameterise
           x = 0.0014 ; %C&F (2018) and Hayashi & Prescott (2002)
           alpha = 0.377; %Braun et al. (2006)
           delta = 0.025; %Braun et al. (2006)
           kappa = 12; %C&F (2018)
           phiI = 132; %Higo and Saita (2007)
           rhoK = 0.895; %C&F (2018) (6 quarters)
           sigma = 0.5; %Hall (1988) and Yogo (2004)
           upsilon = 0.5 ; %Gertler (1999)
           g = 0.17; %C&F (2018)
           b = 4; %C&F (2018)
           phiPI = 2; %Bernanke and Gertler (1999)
           phiR = 0.45; %C&F (2018)
           beta = 0.99; %F&T (2008)
           varrho = 0.4; %C&F (2018)
           rbar = 1/beta;
           phiNX = 0.08 ; %Aoki et al.
           ystar = 1 ;
           Rstar = 1/beta ;
           phiFSTAR = 0.000742 ; %S&U
           fstarbar = 0.7442 ; %S&U
           varsigma_initial = 0.88 ; %Gertler (1999): 0.6
           varsigma_final = 0.95 ;
           omega_initial = 0.9948 ; //retire at 68 as at 1990
           omega_final = 0.9954 ; //retire at 75 as at 2020
           gamma_initial = 0.9750 ; //retire at 68, die at 78
           gamma_final = 0.9875 ; //retire at 68, die at 88 by 2050
           n_initial = 0.001075 ; //population growth as at 1990
           n_final = -0.001425 ; //population growth as at 2050

model;
//PRODUCTION AND OUTPUT
[name='Law of motion of capital with adjustment costs']
(1+n+x)*K = ((1-delta)*K(-1) + (1-((I/I(-1))-1)^2)*I)*exp(-ZK) ;

[name='Real wage']
W = MC*(1-alpha)*Y ;

[name='Rental rate']
RK = MC*alpha*Y/K(-1) ;

[name='Marginal cost']
MC = (RK^alpha*W^(1-alpha))/(alpha^alpha*(1-alpha)^(1-alpha)) ;

[name='Profits of intermediate goods producers']
DI = (1 - MC - (phiI/2)*(PI-1)^2)*Y ;

[name='Phillips Curve']
(PI-1)*PI = ((kappa-1)/phiI)*((kappa/(kappa-1))*MC-1) + (Y(+1)/Y)*(1+n+x)*(PI(+1)-1)*PI(+1)/(R/PI(+1)) ;


//HOUSEHOLDS AND AGGREGATION
[name='Dependency ratio']
(1+n)*PSI = (1-omega) + gamma*PSI(-1) ;

[name='Distribution of wealth']
(1+n+x)*(LAMBDA - (1-omega))*A = omega*((1-MPCR)*LAMBDA(-1)*R(-1)*A(-1)/PI + E - MPCR*SR) ;

[name='Ratio of MPCs']
XI = MPCR/MPCW ;

[name='Adjustment factor']
OMEGA = omega + (1-omega)*(1/varsigma)^(1-upsilon)*XI^(1/(1-sigma)) ;

[name='Aggregate welfare']
V = VR + VW ;

[name='Aggregate consumption']
C = CR + CW ;

[name='Aggregate value of human wealth']
H = HR + HW ;


//THE ELDERLY
[name='Old welfare']
VR = MPCR^(sigma/(1-sigma))*CR*(((1-upsilon)/upsilon)*1/(varsigma*W))^(1-upsilon) ;

[name='MPC for old']
MPCR = 1 - (W/W(+1))^((1-upsilon)*sigma/(1-sigma))*gamma*beta^sigma*(R/PI(+1))^(sigma-1)*(MPCR/MPCR(+1)) ;

[name='Aggregate consumption of the old']
CR = MPCR*((R(-1)/PI)*AO(-1) + SR + HR) ;

[name='Old assets']
LAMBDA = AO/A ;

[name='Social security benefits to retirees']
SR = E + PSI*SR(+1)/((1+n)*PSI(+1)*R/(gamma*PI(+1))) ;

[name='Human assets of the old']
HR = W*varsigma*LR + (1+x)*HR(+1)*(gamma*PI(+1))/((1+n)*R(+1)) ;

[name='Labour supply by the old']
LR = 1 - ((1-upsilon)/upsilon)*CR/(varsigma*W) ;


//THE YOUNG
[name='Young welfare']
VW = MPCW^(sigma/(1-sigma))*CW*(((1-upsilon)/upsilon)*1/W)^(1-upsilon) ;

[name='MPC for young']
MPCW = 1 - (W/W(+1))^((1-upsilon)*sigma/(1-sigma))*beta^sigma*(OMEGA(+1)*R/PI(+1))^(sigma-1)*(MPCW/MPCW(+1)) ;

[name='Aggregate consumption for workers']
CW = MPCW*(R(-1)*AW(-1)/PI + SW + HW) ;

[name='Young assets']
1 - LAMBDA = AW/A ;

[name='Present discounted social security benefits to workers']
SW = ((OMEGA(+1)-omega)*SR(+1)/PSI(+1))/(OMEGA(+1)*R/PI(+1)) + (omega*SW(+1))/(OMEGA(+1)*R/PI(+1)) ;

[name='Human assets of the young']
HW = W*LW - T + omega*(1+x)*HW(+1)*PI(+1)/((1+n)*R(+1)*OMEGA(+1))+(1-omega)*(1+x)*(XI(+1))^(1/(sigma-1))*(1/varsigma)^(1-upsilon)*((HR(+1)*PI(+1))/((1+n)*R(+1)*OMEGA(+1))) ;

[name='Labour supply by the young']
LW = 1 - ((1-upsilon)/upsilon)*CW/W ;


//FISCAL AND MONETARY POLICY
[name='Government budget constraint']
(1+n+x)*B = (R(-1)/PI)*B(-1) + G + E - T ;

[name='Fiscal rule']
(1+n+x)*B = b*Y ;

[name='Government consumption']
G = g*Y ;

[name='Government pension expenses']
E = varrho*PSI*((1-alpha)*MC*Y - T) ;

[name='Monetary policy rule']
R = rbar^phiR*R(-1)^(1-phiR)*PI^phiPI ;


//MARKET CLEARING CONDITIONS
[name='Resource constraint']
(1 - (phiI/2)*(PI-1)^2)*Y = C + I + G + NX ;

[name='No arbitrage in intermediate markets']
R/PI(+1) = (1+n+x)*(PINT(+1)+DI(+1))/PINT ;

[name='No arbitrage in capital markets']
R/PI(+1) = RK(+1) + 1 - delta ;

[name='Aggregate assets']
A = K + B + PINT + FSTAR ;

[name='Aggregate labour supply']
L = varsigma*LR + LW ;

//OPEN ECONOMY
[name='Net exports']
NX = phiNX*ystar ;

[name='Law of motion for net foreign assets']
FSTAR = (Rstar + PHI(-1) - 1)*FSTAR(-1) + NX ;

[name='Risk premium on foreign debt']
PHI = sign(FSTAR(-1))*phiFSTAR*(exp(FSTAR(-1)-fstarbar)-1) ;

//SHOCK PROCESSES
[name='Shock process for capital/investment']
ZK = rhoK*ZK(-1) - EK ;

end ;

initval ;
omega      =omega_initial ;
gamma      =gamma_initial ;
n          =n_initial ;
varsigma   =varsigma_initial ;
Y      		 =3.24148;
W      		 =1.85314;
RK     		 =0.0495455;
MC     		 =0.91765;
I      		 =0.621864;
K      		 =22.6338;
L      		 =0.808891;
PINT   		 =11.9589;
DI     		 =0.263287;
PSI    		 =0.199425;
LAMBDA 		 =0.213026;
OMEGA  		 =1.01266;
XI     		 =1.79482;
V      		 =0.049662;
C      		 =1.98492;
H      		 =12.6605;
A      		 =47.6074;
VR     		 =0.0332761;
MPCR   		 =0.0415784;
CR     		 =1.02202;
AO     		 =10.1416;
SR     		 =1.52096;
LR     		 =0.373287;
HR     		 =12.6691;
VW     		 =0.016386;
MPCW   		 =0.0231658;
CW     		 =0.962895;
AW     		 =37.4658;
SW     		 =3.18846;
LW     		 =0.480398;
HW     		 =-0.00856248;
B      		 =12.9339;
G      		 =0.551052;
E      		 =0.0751054;
T      		 =0.911615;
R      		 =1.02878;
PI     		 =1.00413;
NX     		 =0.08;
PHI    		 =-0.000359803;
FSTAR  		 =0.080787;
ZK     		 =0;
end ;

steady(solve_algo=0) ;

// Experiment 1: A decline in the population growth rate combined with an increase in life expectancy.
// Experiment 2: An increase in the effective retirement age from 68 to 75.
// Experiment 3: A decrease in the productivity gap between young and old workers.
// Experiment 4: A combination of experiments 1 and 2.
// Experiment 5: A combination of experiments 1, 2, and 3.
// Experiment 6: A combination of experiments 1 and 3.

//******************************************************************************
@#if EXP1 == 1
endval ;
omega      =omega_initial;
gamma      =gamma_final;
n          =n_final;
varsigma   =varsigma_initial;
end ;

steady(solve_algo=0) ;

load nsmoothed.mat ; //MATLAB file
nshocks = nsmoothed ;
load gammasmoothed.mat ; //MATLAB file: need to update if omega changes
gammashocks = gammasmoothed ;

shocks ;
var gamma ;
periods 1:240 ;
values (gammasmoothed) ;

var n ;
periods 1:240 ;
values (nsmoothed) ;
end ;

perfect_foresight_setup(periods=250);
perfect_foresight_solver;
//model does not converge to new SS by t=240, so need to run more periods

t = datetime(1990,01,01):calquarters(1):datetime(2050,12,31) ;

figure('Name','Change in population growth and survival rate')
subplot(4,4,1)
plot(t(1:241),Y(1:241),'r')
datetick('x','yyyy','keeplimits')
axis tight
title('y')

subplot(4,4,2)
plot(t(1:241),W(1:241),'r')
datetick('x','yyyy')
axis tight
title('W')

subplot(4,4,3)
plot(t(1:241),I(1:241),'r')
datetick('x','yyyy')
axis tight
title ('i')

subplot(4,4,4)
plot(t(1:241),K(1:241),'r')
datetick('x','yyyy')
axis tight
title('k')

subplot(4,4,5)
plot(t(1:241),L(1:241),'r')
datetick('x','yyyy')
axis tight
title('l')

subplot(4,4,6)
plot(t(1:241),PSI(1:241),'r')
datetick('x','yyyy')
axis tight
title('\Psi')

subplot(4,4,7)
plot(t(1:241),LAMBDA(1:241),'r')
datetick('x','yyyy')
axis tight
title('\lambda')

subplot(4,4,8)
plot(t(1:241),XI(1:241),'r')
datetick('x','yyyy')
axis tight
title('\Xi')

subplot(4,4,9)
plot(t(1:241),H(1:241),'r')
datetick('x','yyyy')
axis tight
title('h')

subplot(4,4,10)
plot(t(1:241),A(1:241),'r')
datetick('x','yyyy')
axis tight
title('a')

subplot(4,4,11)
plot(t(1:241),B(1:241),'r')
datetick('x','yyyy')
axis tight
title('b')

subplot(4,4,12)
plot(t(1:241),G(1:241),'r')
datetick('x','yyyy')
axis tight
title('g')

subplot(4,4,13)
plot(t(1:241),E(1:241),'r')
datetick('x','yyyy')
axis tight
title('e')

subplot(4,4,14)
plot(t(1:241),T(1:241),'r')
datetick('x','yyyy')
axis tight
title('t')

subplot(4,4,15)
plot(t(1:241),(R(1:241)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('R')

subplot(4,4,16)
plot(t(1:241),(PI(1:241)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('\pi')
@# endif

//******************************************************************************
@#if EXP2 == 1
endval ;
omega      =omega_final ;
gamma      =0.9167; //retire at 75 die at roughly 78
n          =n_initial;
varsigma   =varsigma_initial;
end ;

steady(solve_algo=0) ;

shocks ;
var gamma ;
periods 1:120 ;
values (gamma_initial) ;

var omega ;
periods 1:120 ;
values (omega_initial) ;
end ;

perfect_foresight_setup(periods=250);
perfect_foresight_solver;

t = datetime(1990,01,01):calquarters(1):datetime(2050,12,31) ;

figure('Name','Change in retirement age in 2020')
subplot(4,4,1)
plot(t(1:241),Y(1:241),'r')
datetick('x','yyyy','keeplimits')
axis tight
title('y')

subplot(4,4,2)
plot(t(1:241),W(1:241),'r')
datetick('x','yyyy')
axis tight
title('W')

subplot(4,4,3)
plot(t(1:241),I(1:241),'r')
datetick('x','yyyy')
axis tight
title ('i')

subplot(4,4,4)
plot(t(1:241),K(1:241),'r')
datetick('x','yyyy')
axis tight
title('k')

subplot(4,4,5)
plot(t(1:241),L(1:241),'r')
datetick('x','yyyy')
axis tight
title('l')

subplot(4,4,6)
plot(t(1:241),PSI(1:241),'r')
datetick('x','yyyy')
axis tight
title('\Psi')

subplot(4,4,7)
plot(t(1:241),LAMBDA(1:241),'r')
datetick('x','yyyy')
axis tight
title('\lambda')

subplot(4,4,8)
plot(t(1:241),XI(1:241),'r')
datetick('x','yyyy')
axis tight
title('\Xi')

subplot(4,4,9)
plot(t(1:241),H(1:241),'r')
datetick('x','yyyy')
axis tight
title('h')

subplot(4,4,10)
plot(t(1:241),A(1:241),'r')
datetick('x','yyyy')
axis tight
title('a')

subplot(4,4,11)
plot(t(1:241),B(1:241),'r')
datetick('x','yyyy')
axis tight
title('b')

subplot(4,4,12)
plot(t(1:241),G(1:241),'r')
datetick('x','yyyy')
axis tight
title('g')

subplot(4,4,13)
plot(t(1:241),E(1:241),'r')
datetick('x','yyyy')
axis tight
title('e')

subplot(4,4,14)
plot(t(1:241),T(1:241),'r')
datetick('x','yyyy')
axis tight
title('t')

subplot(4,4,15)
plot(t(1:241),(R(1:241)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('R')

subplot(4,4,16)
plot(t(1:241),(PI(1:241)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('\pi')
@# endif

//******************************************************************************
@#if EXP3 == 1
endval ;
omega      =omega_initial;
gamma      =gamma_initial;
n          =n_initial;
varsigma   =varsigma_final;
end ;

steady(solve_algo=0) ;

varsigmashocks = varsigma_initial:(varsigma_final-varsigma_initial)/(120-1):varsigma_final ;

shocks ;
var varsigma ;
periods 1:120 ;
values (varsigmashocks) ;
end ;

perfect_foresight_setup(periods=250);
perfect_foresight_solver;

t = datetime(1990,01,01):calquarters(1):datetime(2050,12,31) ;

figure('Name','Change in relative productivity of elderly')
subplot(4,4,1)
plot(t(1:121),Y(1:121),'r')
datetick('x','yyyy','keeplimits')
axis tight
title('y')

subplot(4,4,2)
plot(t(1:121),W(1:121),'r')
datetick('x','yyyy')
axis tight
title('W')

subplot(4,4,3)
plot(t(1:121),I(1:121),'r')
datetick('x','yyyy')
axis tight
title ('i')

subplot(4,4,4)
plot(t(1:121),K(1:121),'r')
datetick('x','yyyy')
axis tight
title('k')

subplot(4,4,5)
plot(t(1:121),L(1:121),'r')
datetick('x','yyyy')
axis tight
title('l')

subplot(4,4,6)
plot(t(1:121),PSI(1:121),'r')
datetick('x','yyyy')
axis tight
title('\Psi')

subplot(4,4,7)
plot(t(1:121),LAMBDA(1:121),'r')
datetick('x','yyyy')
axis tight
title('\lambda')

subplot(4,4,8)
plot(t(1:121),XI(1:121),'r')
datetick('x','yyyy')
axis tight
title('\Xi')

subplot(4,4,9)
plot(t(1:121),H(1:121),'r')
datetick('x','yyyy')
axis tight
title('h')

subplot(4,4,10)
plot(t(1:121),A(1:121),'r')
datetick('x','yyyy')
axis tight
title('a')

subplot(4,4,11)
plot(t(1:121),B(1:121),'r')
datetick('x','yyyy')
axis tight
title('b')

subplot(4,4,12)
plot(t(1:121),G(1:121),'r')
datetick('x','yyyy')
axis tight
title('g')

subplot(4,4,13)
plot(t(1:121),E(1:121),'r')
datetick('x','yyyy')
axis tight
title('e')

subplot(4,4,14)
plot(t(1:121),T(1:121),'r')
datetick('x','yyyy')
axis tight
title('t')

subplot(4,4,15)
plot(t(1:121),(R(1:121)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('R')

subplot(4,4,16)
plot(t(1:121),(PI(1:121)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('\pi')
@# endif

//******************************************************************************
@#if EXP4 == 1
endval ;
omega      =omega_final;
gamma      =0.9808; //retire at 75 die at 88
n          =n_final;
varsigma   =varsigma_initial;
end ;

steady(solve_algo=0) ;

load nsmoothed.mat ; //MATLAB file
nshocks = nsmoothed ;
load gammaomegasmoothed.mat ; //MATLAB file. Gamma updated for omega shock
gammashocks = gammaomegasmoothed ;

shocks ;
var omega ;
periods 1:120 ;
values (omega_initial) ;

var gamma ;
periods 1:240 ;
values (gammaomegasmoothed) ;

var n ;
periods 1:240 ;
values (nsmoothed) ;
end ;

perfect_foresight_setup(periods=250);
perfect_foresight_solver;

t = datetime(1990,01,01):calquarters(1):datetime(2050,12,31) ;

figure('Name','Change in population growth, survival rate and retirement age')
subplot(4,4,1)
plot(t(1:241),Y(1:241),'r')
datetick('x','yyyy','keeplimits')
axis tight
title('y')

subplot(4,4,2)
plot(t(1:241),W(1:241),'r')
datetick('x','yyyy')
axis tight
title('W')

subplot(4,4,3)
plot(t(1:241),I(1:241),'r')
datetick('x','yyyy')
axis tight
title ('i')

subplot(4,4,4)
plot(t(1:241),K(1:241),'r')
datetick('x','yyyy')
axis tight
title('k')

subplot(4,4,5)
plot(t(1:241),L(1:241),'r')
datetick('x','yyyy')
axis tight
title('l')

subplot(4,4,6)
plot(t(1:241),PSI(1:241),'r')
datetick('x','yyyy')
axis tight
title('\Psi')

subplot(4,4,7)
plot(t(1:241),LAMBDA(1:241),'r')
datetick('x','yyyy')
axis tight
title('\lambda')

subplot(4,4,8)
plot(t(1:241),XI(1:241),'r')
datetick('x','yyyy')
axis tight
title('\Xi')

subplot(4,4,9)
plot(t(1:241),H(1:241),'r')
datetick('x','yyyy')
axis tight
title('h')

subplot(4,4,10)
plot(t(1:241),A(1:241),'r')
datetick('x','yyyy')
axis tight
title('a')

subplot(4,4,11)
plot(t(1:241),B(1:241),'r')
datetick('x','yyyy')
axis tight
title('b')

subplot(4,4,12)
plot(t(1:241),G(1:241),'r')
datetick('x','yyyy')
axis tight
title('g')

subplot(4,4,13)
plot(t(1:241),E(1:241),'r')
datetick('x','yyyy')
axis tight
title('e')

subplot(4,4,14)
plot(t(1:241),T(1:241),'r')
datetick('x','yyyy')
axis tight
title('t')

subplot(4,4,15)
plot(t(1:241),(R(1:241)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('R')

subplot(4,4,16)
plot(t(1:241),(PI(1:241)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('\pi')
@# endif

//******************************************************************************
@#if EXP5 == 1
endval ;
omega      =omega_final;
gamma      =0.9808; //retire at 75, die at 88
n          =n_final;
varsigma   =varsigma_final;
end ;

steady(solve_algo=0) ;

load nsmoothed.mat ; //MATLAB file
nshocks = nsmoothed ;
load gammaomegasmoothed.mat ; //MATLAB file. Gamma updated for omega shock
gammashocks = gammaomegasmoothed ;

varsigmashocks = varsigma_initial:(varsigma_final-varsigma_initial)/(120-1):varsigma_final ;

shocks ;
var omega ;
periods 1:120 ;
values (omega_initial) ;

var gamma ;
periods 1:240 ;
values (gammaomegasmoothed) ;

var n ;
periods 1:240 ;
values (nsmoothed) ;

var varsigma ;
periods 1:120 ;
values (varsigmashocks) ;
end ;

perfect_foresight_setup(periods=250);
perfect_foresight_solver;

t = datetime(1990,01,01):calquarters(1):datetime(2050,12,31) ;

figure('Name','Change in population growth, survival rate, retirement age, and old productivity')
subplot(4,4,1)
plot(t(1:241),Y(1:241),'r')
datetick('x','yyyy','keeplimits')
axis tight
title('y')

subplot(4,4,2)
plot(t(1:241),W(1:241),'r')
datetick('x','yyyy')
axis tight
title('W')

subplot(4,4,3)
plot(t(1:241),I(1:241),'r')
datetick('x','yyyy')
axis tight
title ('i')

subplot(4,4,4)
plot(t(1:241),K(1:241),'r')
datetick('x','yyyy')
axis tight
title('k')

subplot(4,4,5)
plot(t(1:241),L(1:241),'r')
datetick('x','yyyy')
axis tight
title('l')

subplot(4,4,6)
plot(t(1:241),PSI(1:241),'r')
datetick('x','yyyy')
axis tight
title('\Psi')

subplot(4,4,7)
plot(t(1:241),LAMBDA(1:241),'r')
datetick('x','yyyy')
axis tight
title('\lambda')

subplot(4,4,8)
plot(t(1:241),XI(1:241),'r')
datetick('x','yyyy')
axis tight
title('\Xi')

subplot(4,4,9)
plot(t(1:241),H(1:241),'r')
datetick('x','yyyy')
axis tight
title('h')

subplot(4,4,10)
plot(t(1:241),A(1:241),'r')
datetick('x','yyyy')
axis tight
title('a')

subplot(4,4,11)
plot(t(1:241),B(1:241),'r')
datetick('x','yyyy')
axis tight
title('b')

subplot(4,4,12)
plot(t(1:241),G(1:241),'r')
datetick('x','yyyy')
axis tight
title('g')

subplot(4,4,13)
plot(t(1:241),E(1:241),'r')
datetick('x','yyyy')
axis tight
title('e')

subplot(4,4,14)
plot(t(1:241),T(1:241),'r')
datetick('x','yyyy')
axis tight
title('t')

subplot(4,4,15)
plot(t(1:241),(R(1:241)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('R')

subplot(4,4,16)
plot(t(1:241),(PI(1:241)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('\pi')
@# endif

//******************************************************************************
@#if EXP6 == 1
endval ;
omega      =omega_initial;
gamma      =gamma_final;
n          =n_final;
varsigma   =varsigma_final;
end ;

steady(solve_algo=0) ;

load nsmoothed.mat ; //MATLAB file
nshocks = nsmoothed ;
load gammasmoothed.mat ; //MATLAB file
gammashocks = gammasmoothed ;

varsigmashocks = varsigma_initial:(varsigma_final-varsigma_initial)/(120-1):varsigma_final ;

shocks ;
var gamma ;
periods 1:240 ;
values (gammasmoothed) ;

var n ;
periods 1:240 ;
values (nsmoothed) ;

var varsigma ;
periods 1:120 ;
values (varsigmashocks) ;
end ;

perfect_foresight_setup(periods=250);
perfect_foresight_solver;

t = datetime(1990,01,01):calquarters(1):datetime(2050,12,31) ;

figure('Name','Change in population growth, survival rate, and old productivity')
subplot(4,4,1)
plot(t(1:121),Y(1:121),'r')
datetick('x','yyyy','keeplimits')
axis tight
title('y')

subplot(4,4,2)
plot(t(1:121),W(1:121),'r')
datetick('x','yyyy')
axis tight
title('W')

subplot(4,4,3)
plot(t(1:121),I(1:121),'r')
datetick('x','yyyy')
axis tight
title ('i')

subplot(4,4,4)
plot(t(1:121),K(1:121),'r')
datetick('x','yyyy')
axis tight
title('k')

subplot(4,4,5)
plot(t(1:121),L(1:121),'r')
datetick('x','yyyy')
axis tight
title('l')

subplot(4,4,6)
plot(t(1:121),PSI(1:121),'r')
datetick('x','yyyy')
axis tight
title('\Psi')

subplot(4,4,7)
plot(t(1:121),LAMBDA(1:121),'r')
datetick('x','yyyy')
axis tight
title('\lambda')

subplot(4,4,8)
plot(t(1:121),XI(1:121),'r')
datetick('x','yyyy')
axis tight
title('\Xi')

subplot(4,4,9)
plot(t(1:121),H(1:121),'r')
datetick('x','yyyy')
axis tight
title('h')

subplot(4,4,10)
plot(t(1:121),A(1:121),'r')
datetick('x','yyyy')
axis tight
title('a')

subplot(4,4,11)
plot(t(1:121),B(1:121),'r')
datetick('x','yyyy')
axis tight
title('b')

subplot(4,4,12)
plot(t(1:121),G(1:121),'r')
datetick('x','yyyy')
axis tight
title('g')

subplot(4,4,13)
plot(t(1:121),E(1:121),'r')
datetick('x','yyyy')
axis tight
title('e')

subplot(4,4,14)
plot(t(1:121),T(1:121),'r')
datetick('x','yyyy')
axis tight
title('t')

subplot(4,4,15)
plot(t(1:121),(R(1:121)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('R')

subplot(4,4,16)
plot(t(1:121),(PI(1:121)-1)*400,'r')
datetick('x','yyyy')
axis tight
title('\pi')
@# endif

/*
//Code for when change in parameters is unexpected.
yy = oo_.steady_state ;
//First period
perfect_foresight_setup(periods=120) ;
perfect_foresight_solver ;
yy = [yy,oo_.endo_simul(:,2)] ;

//Following periods
for i=2:length(nshocks)
  oo_.exo_simul(4) = nshocks(i) ;
  oo_.exo_simul(3) = gammashocks(i) ;
  oo_.endo_simul(:,1) = yy(:,end) ;
  perfect_foresight_solver ;
  yy = [yy, oo_.endo_simul(:,2)] ;
end ;

yy = [yy, oo_.endo_simul(:,3:end)] ;
*/
