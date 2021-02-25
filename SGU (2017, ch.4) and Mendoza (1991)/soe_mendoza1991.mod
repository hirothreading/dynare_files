% Small Open Economy based on Mendoza (1991) and SGU (2017, ch4)
% Calibrated to Canadian business cycle moments (as in SGU (2017))
% Author: David Murakami (MPhil Economics, Oxford)
% Implemented on Dynare 4.6.3

// Declare variables
var ld      $\ln d$         (long_name='foreign liabilities')
    lc      $\ln c$         (long_name='consumption')
    lR      $\ln R$         (long_name='gross interest rate')
    h       $h$             (long_name='household labour hours')
    lw      $\ln w$         (long_name='wages')
    llam    $\ln\lambda$    (long_name='shadow price')
    ly      $\ln y$         (long_name='output')
//    rk      $r_t$           (long_name='capital rental rate')
    lk      $\ln k$         (long_name='capital')
    li      $\ln i$         (long_name='investment')
    lA      $\ln A$         (long_name='total factor productivity')
    lnx     $\ln nx$        (long_name='trade balance')
    ca     $\ln ca$        (long_name='current account balance')
    ;


// Exogenous variables
varexo epsA $\varepsilon^{A}$  (long_name='shock process for productivity')
    ;


// Declare parameters
parameters betta    $\beta$     (long_name='household discount factor')
           sig      $\sigma$    (long_name='degree of relative risk aversion')
           del      $\delta$    (long_name='capital depreciation rate')
           omega    $\omega$    (long_name='governs wage elasticity of labour supply')
           alphha   $\alpha$    (long_name='capital share of income')
           rhoA     $\rho_A$    (long_name='persistence of TFP shock')
           sigA     $\sigma_A$  (long_name='standard deviation of shock process of TFP')
           phii     $\phi$      (long_name='governs magnitude of capital adjustment costs')
           psii     $\psi$      (long_name='debt sensitivity of the interest rate')
           kap      $\kappa$    (long_name='capital output ratio in steady state')
    ;
// Set values
betta = 0.96;
sig = 2;
del = 0.1;
omega = 1.455;
alphha = 0.32;
rhoA = 0.42; //initial guess
sigA = 1; //standard 1% shock
phii = 0.028; //initial guess
psii = 0.000742; //initial guess
kap = ((betta^-1 - 1 + del)/alphha)^(1/(alphha-1));

model;
// Additional parameters and variables
// Steady state labour hours
#barh = ((1-alphha)*kap^alphha)^(1/(omega-1));

// Steady state capital
#bark = kap*barh;

// Steady state R^*
#barRstar = 1/betta;

// Steady state y
#bary = ((1-alphha)*kap^(alphha*omega))^(1/(omega-1));

// Steady state trade balance to GDP ratio
#barnxy = 0.02;

// Steady state external debt
#bard = barnxy*bary/(barRstar-1);

// Steady state consumption
#barc = (1-barRstar)*bard + kap^alphha*barh - del*bark;

[name='eq.(3) process for TFP']
lA = (1-rhoA)*steady_state(lA) + rhoA*lA(-1) + epsA/100;

[name='eq.(4) gross interest rate/country interest premium']
exp(lR) = barRstar + psii*(exp(exp(ld)-bard)-1);

[name='eq.(5) output']
exp(ly) = lA*exp(lk(-1))^alphha*h^(1-alphha);

[name='eq.(6) investment']
exp(lk) = (1-del)*exp(lk(-1)) + exp(li);

[name='eq.(9) FOC wrt consumption']
exp(llam) = (exp(lc)-h^omega/omega)^(-sig);

[name='eq.(10) FOC wrt h']
exp(llam)*h^(omega-1) = (1-alphha)*exp(llam)*exp(ly)/h;

[name='eq.(11) FOC wrt k'] //replaces eq.(24)
exp(llam)*(1+phii*(exp(lk)-exp(lk(-1)))) = betta*exp(llam(+1))*(alphha*exp(ly)/exp(lk)+1-del+phii*(exp(lk(+1))-exp(lk)));

[name='eq.(12) FOC wrt external debt']
exp(llam) = betta*exp(lR)*exp(llam(+1));

// rewrote this equation to get rid of s (ls)
[name='eq.(18) period budget constraint']
exp(ld) = exp(lR(-1))*exp(ld(-1)) + exp(lc) + exp(li) + phii/2*(lk-lk(-1))^2 - exp(ly);

[name='eq.(21) marginal product of labour']
exp(lw) = (1-alphha)*lA*exp(lk(-1))^alphha*h^(-alphha);

// This is a redundant equation
//[name='eq.(22) marginal product of capital']
//rk(-1) = alphha*lA*exp(lk(-1))^(alphha-1)*h^(1-alphha);

[name='eq.(26) trade balance']
exp(lnx) = exp(ly) - exp(lc) - exp(li) - phii/2*(exp(lk)-exp(lk(-1)))^2;

[name='eq.(27.1) current account balance in terms of external debt']
ca = exp(ld(-1)) - exp(ld);

end;


initval;
ld = log(0.700919377);
lc = log(1.101199316);
lR = log(1.041666667);
h = 0.995162450;
lw = log(0.997796006);
llam = log(5.702929541);
ly = log(1.460248703);
//rk = 0.141666667;
lk = log(3.298444128);
li = log(0.329844413);
lA = 1.000000000;
lnx = log(0.029204974);
ca = 0;
end;


steady;
check;
model_diagnostics;

shocks;
var epsA = sigA^2;

end;

write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;

stoch_simul(order=1,irf=10,nodisplay);

h_pos=strmatch('h',M_.endo_names,'exact');

%Baseline impulse responses to TFP shock (log-deviations)
set(0,'DefaultAxesTitleFontWeight','normal');
figure();
subplot(3,2,1);
plot(100*ly_epsA,'b-','LineWidth',1);
axis tight;
title('Output');

subplot(3,2,2);
plot(100*lc_epsA,'b-','LineWidth',1);
axis tight;
title('Consumption');

subplot(3,2,3);
plot(100*li_epsA,'b-','LineWidth',1);
axis tight;
title('Investment');

subplot(3,2,4);
plot(100*(h_epsA/oo_.dr.ys(h_pos)),'b-','LineWidth',1);
axis tight;
title('Hours');

subplot(3,2,5);
plot(lnx_epsA-ly_epsA,'b-','LineWidth',1);
axis tight;
title('Trade Balance/Output');

subplot(3,2,6);
plot(100*lA_epsA,'b-','LineWidth',1);
axis tight;
title('TFP');
