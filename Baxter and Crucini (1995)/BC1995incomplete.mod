// Business Cycles, Foreign Trade, and Incomplete Markets
// Replication of Baxter and Crucini (1995, IER)
// Code written by David Murakami (Oxford, MPhil Economics)
// For use with Dynare 4.6.3

// Incomplete markets model
// Define variables
var c           $\hat{c}$         (long_name='Domestic consumption')
    cstar       $\hat{c}^*$       (long_name='Foreign consumption')
    L           $\hat{L}$         (long_name='Domestic leisure time')
    Lstar       $\hat{L}^*$       (long_name='Foreign leisure time')
    N           $\hat{N}$         (long_name='Domestic labour supply')
    Nstar       $\hat{N}^*$       (long_name='Foreign labour supply')
    inv         $\hat{i}$         (long_name='Domestic investment')
    istar       $\hat{i}^*$       (long_name='Foreign investment')
    y           $\hat{y}$         (long_name='Domestic output')
    ystar       $\hat{y}^*$       (long_name='Foreign output')
    A           $\hat{A}$         (long_name='Domestic TFP')
    Astar       $\hat{A}^*$       (long_name='Foreign TFP')
    b           $\hat{b}$         (long_name='Domestic net debt')
    bstar       $\hat{b}^*$       (long_name='Foreign net debt')
    Q           $\hat{Q}$         (long_name='Price of investment goods')
    w           $\hat{w}$         (long_name='Domestic real wages')
    wstar       $\hat{w}^*$       (long_name='Foreign real wages')
    Rk          $R^k$             (long_name='Domestic MPK')
    Rkstar      $R^{k^*}$         (long_name='Foreign MPK')
    k           $\hat{k}$         (long_name='Domestic capital')
    kstar       $\hat{k}^*$       (long_name='Foreign capital')
    b           $\hat{b}$         (long_name='Domestic net debt')
    bstar       $\hat{b}^*$       (long_name='Foreign net debt')
    R           $\hat{R}$         (long_name='International interest rate')
    ;

//Exogenous processes
varexo epsa         $\varepsilon^{a}$       (long_name='Domestic TFP shock')
       epsastar     $\varepsilon^{a^*}$     (long_name='Foreign TFP shock')
       ;

//Declare parameters
parameters BETA         $\beta$                     (long_name='Discount factor')
           THETA        $\theta$                    (long_name='Consumption preference')
           SIGMA        $\sigma$                    (long_name='Inverse of intertemporal elasticity of substitution')
           GAMMA        $\gamma$                    (long_name='Trend growth factor')
           PI           $\pi$                       (long_name='Relative size of home economy')
           DELTA        $\delta$                    (long_name='Depreciation rate')
           KAPPA        $\kappa_{I}$                (long_name='Investment adjustment cost')
           ALPHA        $\alpha$                    (long_name='Capital share of income')
           RHO          $\rho$                      (long_name='Domestic TFP shock persistence')
           RHOSTAR      $\rho^*$                    (long_name='Foreign-Domestic TFP shock persistence')
           NU           $\nu$                       (long_name='Domestic-Foreign TFP shock persistence')
           NUSTAR       $\nu^*$                     (long_name='Foreign TFP shock persistence')
           SIGE         $\sigma_{\epsilon}$         (long_name='Domestic SD of TFP shock')
           SIGESTAR     $\sigma_{\epsilon^*}$       (long_name='Foreign SD of TFP shock')
           PSI          $\psi$                      (long_name='Covariance of TFP shock')
           PSIB         $\psi_b$                    (long_name='Portfolio adjustment cost')
           KNRATSS      $\frac{\bar{c}}{\bar{N}}$   (long_name='SS consumption-labour ratio')
           ;

// Set parameters
BETA = 1/(1.065^(1/4)); //BC (1995)
THETA = 0.24; //BC (1993) Chosen so that NSS=0.2
SIGMA = 0.5; //BC (1995) ***inverse of IES***
GAMMA = 1.004; //BC (1995)
PI = 0.5;
DELTA = 0.025; //BC (1995)
KAPPA = 0.67; //ABK (2016)
ALPHA = 0.42; //BC (1995). ***This is 1-ALPHA in the original paper. See notes.***
RHO = 0.906; //BC (1995). ***Baseline calibration in Section 3.3 of BC (1995) paper***
RHOSTAR = 0.906; //See above
NU = 0.088; //See above
NUSTAR = 0.088; //See above
SIGE = 1; //See above
SIGESTAR = 1; //See above
PSI = 0.258; //See above
PSIB = 0.00074;
KNRATSS = (ALPHA/(GAMMA^(1+ALPHA)/BETA - GAMMA^ALPHA*(1-DELTA)))^(1/(1-ALPHA));

model;
#ISS = steady_state(inv);

[name='Domestic labour and leisure constraint']
1 = L + N;

[name='Foreign labour and leisure constraint']
1 = Lstar + Nstar;

[name='Domestic law of motion of capital']
k = (1-DELTA)/GAMMA*k(-1) + inv;

[name='Foreign law of motion of capital']
kstar = (1-DELTA)/GAMMA*kstar(-1) + istar;

[name='Investment adjustment costs']
Q = PI*(1+0.5*KAPPA*(inv/ISS-1)^2 + inv/ISS*KAPPA*(inv/ISS-1)*(1/ISS)) + (1-PI)*(1+0.5*KAPPA*(istar/ISS-1)^2 + istar/ISS*KAPPA*(istar/ISS-1)*(1/ISS));

[name='Domestic production technology']
y = A*(k(-1)/GAMMA)^ALPHA*N^(1-ALPHA);

[name='Foreign production technology']
ystar = Astar*(kstar(-1)/GAMMA)^ALPHA*Nstar^(1-ALPHA);

[name='Domestic marginal product of capital']
Rk = GAMMA*ALPHA*y/k(-1);

[name='Foreign marginal product of capital']
Rkstar = GAMMA*ALPHA*ystar/kstar(-1);

[name='Domestic real wage']
w = (1-ALPHA)*y/N;

[name='Foreign real wage']
wstar = (1-ALPHA)*ystar/Nstar;

[name='Domestic consumption-leisure choice']
c = THETA/(1-THETA)*w*L;

[name='Foreign consumption-leisure choice']
cstar = THETA/(1-THETA)*wstar*Lstar;

[name='Domestic TFP process']
log(A/1) = RHO*log(A(-1)/1) + NU*log(Astar(-1)/1) + epsa/100;

[name='Foreign TFP process']
log(Astar/1) = RHOSTAR*log(Astar(-1)/1) + NUSTAR*log(A(-1)/1) + epsastar/100;

[name='Domestic investment-saving equation']
c^(THETA-SIGMA*THETA-1)*L^(1-THETA-SIGMA*(1-THETA)) = BETA*c(+1)^(THETA-SIGMA*THETA-1)*L(+1)^(1-THETA-SIGMA*(1-THETA))*(Rk(+1)/GAMMA + (1-DELTA)/GAMMA);

[name='Foreign investment-saving equation']
cstar^(THETA-SIGMA*THETA-1)*Lstar^(1-THETA-SIGMA*(1-THETA)) = BETA*cstar(+1)^(THETA-SIGMA*THETA-1)*Lstar(+1)^(1-THETA-SIGMA*(1-THETA))*(Rkstar(+1)/GAMMA + (1-DELTA)/GAMMA);

[name='Domestic Euler equation']
c^(THETA-SIGMA*THETA-1)*L^(1-THETA-SIGMA*(1-THETA))*(1+PSIB*b) = BETA*c(+1)^(THETA-SIGMA*THETA-1)*L(+1)^(1-THETA-SIGMA*(1-THETA))*(R/GAMMA);

[name='Foreign Euler equation']
cstar^(THETA-SIGMA*THETA-1)*Lstar^(1-THETA-SIGMA*(1-THETA))*(1+PSIB*bstar) = BETA*cstar(+1)^(THETA-SIGMA*THETA-1)*Lstar(+1)^(1-THETA-SIGMA*(1-THETA))*(R/GAMMA);

[name='Domestic resource constraint']
w*N + Rk*k(-1)/GAMMA + R(-1)*b(-1)/GAMMA = c + inv + b + 0.5*PSIB*(b)^2;

[name='Foreign resource constraint']
wstar*Nstar + Rkstar*kstar(-1)/GAMMA + R(-1)*bstar(-1)/GAMMA = cstar + istar + bstar + 0.5*PSIB*(bstar)^2;

[name='International bond market clearing']
0 = b + bstar;

end;

initval;
Q = 1;
R = GAMMA/BETA;
A = 1;
Rk = GAMMA^(1-ALPHA)*ALPHA*KNRATSS^(ALPHA-1);
w = (1-ALPHA)*GAMMA^(-ALPHA)*KNRATSS^ALPHA;
N = (THETA*(1-ALPHA)/(GAMMA^ALPHA*(1-THETA))*KNRATSS^ALPHA)/(GAMMA^(-ALPHA)*KNRATSS^ALPHA - KNRATSS*(1-(1-DELTA)/GAMMA) + THETA*(1-ALPHA)/(GAMMA^ALPHA*(1-THETA))*KNRATSS^ALPHA);
L = 1-N;
y = GAMMA^(-ALPHA)*KNRATSS^ALPHA*N;
k = KNRATSS*N;
c = THETA*(1-ALPHA)/(GAMMA^ALPHA*(1-THETA))*KNRATSS^ALPHA*(1-N);
inv = y - c;
b = 0;
Astar = 1;
Rkstar = GAMMA^(1-ALPHA)*ALPHA*KNRATSS^(ALPHA-1);
wstar = (1-ALPHA)*GAMMA^(-ALPHA)*KNRATSS^ALPHA;
Nstar = (THETA*(1-ALPHA)/(GAMMA^ALPHA*(1-THETA))*KNRATSS^ALPHA)/(GAMMA^(-ALPHA)*KNRATSS^ALPHA - KNRATSS*(1-(1-DELTA)/GAMMA) + THETA*(1-ALPHA)/(GAMMA^ALPHA*(1-THETA))*KNRATSS^ALPHA);
Lstar = 1-N;
ystar = GAMMA^(-ALPHA)*KNRATSS^ALPHA*N;
kstar = KNRATSS*N;
cstar = THETA*(1-ALPHA)/(GAMMA^ALPHA*(1-THETA))*KNRATSS^ALPHA*(1-N);
istar = ystar - cstar;
bstar = 0;
end;


steady;
check;
model_diagnostics;

shocks;
var epsa = SIGE^2;
var epsastar = SIGESTAR^2;
var epsa,epsastar = PSI;
end;

write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;

stoch_simul(order=1,irf=40,periods=200);

cpos=strmatch('c',M_.endo_names,'exact');
cstarpos=strmatch('cstar',M_.endo_names,'exact');
Lpos=strmatch('L',M_.endo_names,'exact');
Lstarpos=strmatch('Lstar',M_.endo_names,'exact');
Npos=strmatch('N',M_.endo_names,'exact');
Nstarpos=strmatch('Nstar',M_.endo_names,'exact');
kpos=strmatch('k',M_.endo_names,'exact');
kstarpos=strmatch('kstar',M_.endo_names,'exact');
ipos=strmatch('inv',M_.endo_names,'exact');
istarpos=strmatch('istar',M_.endo_names,'exact');
ypos=strmatch('y',M_.endo_names,'exact');
ystarpos=strmatch('ystar',M_.endo_names,'exact');
Apos=strmatch('A',M_.endo_names,'exact');
Astarpos=strmatch('Astar',M_.endo_names,'exact');
Qpos=strmatch('Q',M_.endo_names,'exact');
wpos=strmatch('w',M_.endo_names,'exact');
wstarpos=strmatch('wstar',M_.endo_names,'exact');
Rkpos=strmatch('Rk',M_.endo_names,'exact');
Rkstarpos=strmatch('Rkstar',M_.endo_names,'exact');
bpos=strmatch('b',M_.endo_names,'exact');
bstarpos=strmatch('bpos',M_.endo_names,'exact');
Rpos=strmatch('R',M_.endo_names,'exact');


save('BC1995_incompletemarkets','M_','oo_','options_','cpos','cstarpos','Lpos','Lstarpos','Npos','Nstarpos','kpos','kstarpos','ipos','istarpos','ypos','ystarpos','Apos','Astarpos','Qpos','wpos','wstarpos','Rpos','Rstarpos');
