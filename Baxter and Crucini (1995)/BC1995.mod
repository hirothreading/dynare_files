// Business Cycles, Foreign Trade, and Incomplete Markets
// Replication of Baxter and Crucini (1995, IER)
// Code written by David Murakami (Oxford, MPhil Economics)
// For use with Dynare 4.6.3

// Define variables
var c           $\hat{c}$         (long_name='Domestic consumption')
    cstar       $\hat{c}^*$       (long_name='Foreign consumption')
    L           $\hat{L}$         (long_name='Domestic leisure time')
    Lstar       $\hat{L}^*$       (long_name='Foreign leisure time')
    N           $\hat{N}$         (long_name='Domestic labour supply')
    Nstar       $\hat{N}^*$       (long_name='Foreign labour supply')
    k           $\hat{k}$         (long_name='Domestic capital')
    kstar       $\hat{k}^*$       (long_name='Foreign capital')
    inv         $\hat{i}$         (long_name='Domestic investment')
    istar       $\hat{i}^*$       (long_name='Foreign investment')
    y           $\hat{y}$         (long_name='Domestic output')
    ystar       $\hat{y}^*$       (long_name='Foreign output')
    A           $\hat{A}$         (long_name='Domestic TFP')
    Astar       $\hat{A}^*$       (long_name='Foreign TFP')
    Q           $\hat{Q}$         (long_name='Price of investment goods')
    w           $\hat{w}$         (long_name='Domestic real wages')
    wstar       $\hat{w}^*$       (long_name='Foreign real wages')
    R           $\hat{R}$         (long_name='Domestic interes rate')
    Rstar       $\hat{R}^*$       (long_name='Foreign rental rate of capital')
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
           IKRATSS      $\frac{\bar{i}}{\bar{k}}$   (long_name='SS investment-capital ratio')
           YKRATSS      $\frac{\bar{y}}{\bar{k}}$   (long_name='SS output-capital ratio')
           IYRATSS      $\frac{\bar{i}}{\bar{y}}$   (long_name='SS investment-output ratio')
           KNRATSS      $\frac{\bar{c}}{\bar{N}}$   (long_name='SS consumption-labour ratio')
           NSS          $\bar{N}$                   (long_name='SS labour supply')
           YSS          $\bar{y}$                   (long_name='SS output')
           ISS          $\bar{i}$                   (long_name='SS investment')
           CSS          $\bar{c}$                   (long_name='SS consumption')
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
IKRATSS = (GAMMA + DELTA - 1)/GAMMA;
YKRATSS = (BETA^(-1)*GAMMA - 1 + DELTA)/(ALPHA*GAMMA);
IYRATSS = (ALPHA*(GAMMA + DELTA -1))/(BETA^(-1)*GAMMA - 1 + DELTA);
KNRATSS = (ALPHA/(GAMMA^ALPHA*(BETA^(-1)-(1-DELTA)/GAMMA)))^(1/(1-ALPHA));
NSS = (GAMMA^(1-ALPHA)*THETA*(1-ALPHA)/(1-THETA))/(GAMMA^(1-ALPHA) - (GAMMA+DELTA-1)*KNRATSS^(1-ALPHA) + GAMMA^(1-ALPHA)*THETA*(1-ALPHA)/(1-THETA));
YSS = 1/(GAMMA^ALPHA)*KNRATSS^ALPHA*NSS;
ISS = (ALPHA*(GAMMA + DELTA - 1))*YSS/(BETA^(-1)*GAMMA - 1 + DELTA);
CSS = YSS - ISS;

model(linear);
[name='Domestic labour and leisure constraint']
N = -(1-NSS)/NSS*L;

[name='Foreign labour and leisure constraint']
Nstar = -(1-NSS)/NSS*Lstar;

[name='Domestic law of motion of capital']
k = ((1-DELTA)/GAMMA)*k(-1) + IKRATSS*inv;

[name='Foreign law of motion of capital']
kstar = ((1-DELTA)/GAMMA)*kstar(-1) + IKRATSS*istar;

[name='Investment adjustment costs']
Q = (PI*KAPPA/ISS)*inv + ((1-PI)*KAPPA/ISS)*istar;

[name='Global resource constraint']
y = (CSS/YSS)*c + (CSS/YSS)*cstar + IYRATSS*inv + IYRATSS*istar - ystar;

[name='Domestic production technology']
y = A + ALPHA*k(-1) + (1-ALPHA)*N;

[name='Foreign production technology']
ystar = Astar + ALPHA*kstar(-1) + (1-ALPHA)*Nstar;

[name='Domestic real wage']
w = y - N;

[name='Foreign real wage']
wstar = ystar - Nstar;

[name='Domestic interest rate']
R = ALPHA*BETA*YKRATSS*(y - k(-1));

[name='Foreign interest rate']
Rstar = ALPHA*BETA*YKRATSS*(ystar - kstar(-1));

[name='Domestic consumption-leisure choice']
c = w + L;

[name='Foreign consumption-leisure choice']
cstar = wstar + Lstar;

[name='Domestic Euler equation']
c = c(+1) + ((1+THETA*(SIGMA-1)-SIGMA)/(THETA*(1-SIGMA)-1))*(L(+1) - L) + (1/(THETA*(1-SIGMA)-1))*R(+1);

[name='Foreign Euler equation']
cstar = cstar(+1) + ((1+THETA*(SIGMA-1)-SIGMA)/(THETA*(1-SIGMA)-1))*(Lstar(+1) - Lstar) + (1/(THETA*(1-SIGMA)-1))*Rstar(+1);

[name='Marginal utility of consumption parity condition']
c = cstar + ((1+THETA*(SIGMA-1)-SIGMA)/(THETA*(1-SIGMA)-1))*(L - Lstar);

[name='Domestic TFP process']
A = RHO*A(-1) + NU*Astar(-1) + epsa;

[name='Foreign TFP process']
Astar = RHOSTAR*Astar(-1) + NUSTAR*A(-1) + epsastar;

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

stoch_simul(order=1,irf=40,nodisplay);

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
Rpos=strmatch('R',M_.endo_names,'exact');
Rstarpos=strmatch('Rstar',M_.endo_names,'exact');

save('BC1995_completemarkets','M_','oo_','options_','cpos','cstarpos','Lpos','Lstarpos','Npos','Nstarpos','kpos','kstarpos','ipos','istarpos','ypos','ystarpos','Apos','Astarpos','Qpos','wpos','wstarpos','Rpos','Rstarpos');
