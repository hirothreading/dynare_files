% Medium Size Model a la JPT with stationary and nonstationary tech shocks and with Imperfect Information
% This is a replication of Blanchard, L'Huillier, and Lorenzoni  (2013, AER)
% by David Murakami for modern Dynare

% The original replication codes used Dynare version 3.065, which allowed
% matrix notation in the model block. This, in addition to changes to
% stoch_simul.m and dynare_estimation.m going from Dynare 3 to Dynare 6,
% means that the original replication files are incompatible with modern
% Dynare.

% Notice: probably due to Kalman filter in model section, cannot use names 
% a and k 
% so pty = a and kk = k

% Declare variables
var xh      ${x}$             (long_name='Permanent TFP') 
    xhh     ${x^L}$           (long_name='Permanent TFP (lagged)')
    zh      ${z}$             (long_name='Transitory TFP')  
    pty     ${a}$             (long_name='Productivity') 
    s       ${s}$             (long_name='TFP Signal') 
    da      ${\Delta a}$      (long_name='Productivity Growth')
    d       ${\Delta^a}$      (long_name='IST')
    q       ${q}$             (long_name='Monetary Policy Shock')
    m_p     ${m^p}$           (long_name='Price Markups')
    epsma_p ${\varepsilon^p}$ (long_name='Price Markup Shock')
    m_w     ${m^w}$           (long_name='Wage Markups')
    epsma_w ${\varepsilon^w}$ (long_name='Wage Markup Shock')
    g       ${g}$             (long_name='Government Spending')
    lam     ${\lambda}$       (long_name='MU Consumption')
    phi     ${\phi}$          (long_name='LM (capital constraint)') % This was wrongly written as 'psi' in the original code
    c       ${c}$             (long_name='Consumption')
    i       ${inv}$           (long_name='Investment')            
    kbar    ${\bar{k}}$       (long_name='Capital Stock')
    kk      ${k}$             (long_name='Capital Services')
    u       ${u}$             (long_name='Capital Utilisation')
    y       ${y}$             (long_name='Output')
    n       ${n}$             (long_name='Labour')
    rk      ${r^k}$           (long_name='Capital Rental Rate')
    w       ${w}$             (long_name='Wages')
    pi      ${\pi}$           (long_name='Inflation')
    mc      ${mc}$            (long_name='Marginal Cost') 
    mc_w    ${mc^w}$          (long_name='Wage Marginal Cost') 
    r       ${r}$             (long_name='Real Interest Rate') 
    ri      ${i}$             (long_name='Nominal Interest Rate')   
    ca      ${c^a}$           (long_name='Consumption')
    ia      ${i^a}$           (long_name='Investment')
    ya      ${y^a}$           (long_name='Output')
    ka      ${k^a}$           (long_name='Capital') 
    wa      ${w^a}$           (long_name='Wages') 
    dca     ${\Delta c^a}$    (long_name='Consumption Growth') 
    diad    ${\Delta i^a}$    (long_name='Investment Growth') 
    dya     ${\Delta y^a}$    (long_name='Output Growth') 
    dn      ${\Delta n}$      (long_name='Labour Growth')
    dwa     ${\Delta w^a}$    (long_name='Wages Growth') 
    dr      ${\Delta r}$      (long_name='Real Interest Rate Change') 
    dpi     ${\Delta \pi}$    (long_name='Inflation Change') 
    dy      ${\Delta y}$      (long_name='Output Growth (detrended)')
;

% Declare shocks
varexo e_1    ${\epsilon}$      (long_name='Permanent Technology Shock')
       e_2    ${\eta}$          (long_name='Temporary Technology Shock')
       eps_d  ${\varepsilon^d}$ (long_name='IST Shock')
       eps_p  ${\varepsilon^p}$ (long_name='Price Markup Shock') 
       eps_w  ${\varepsilon^w}$ (long_name='Wage Markup Shock') 
       eps_q  ${\varepsilon^q}$ (long_name='Monetary Policy Shock') 
       eps_g  ${\varepsilon^g}$ (long_name='Government Spending Shock')
;   

% Declare parameters
parameters rho    ${\rho}$       (long_name='Persistence of Permanent TFP')
           rho_d  ${\rho_d}$     (long_name='Persistence of IST Shock')
           rho_q  ${\rho_q}$     (long_name='Persistence of Monetary Policy Shock') 
           rho_p  ${\rho_p}$     (long_name='Persistence of Price Markup Shock') 
           rho_w  ${\rho_w}$     (long_name='Persistence of Wage Markup Shock')
           rho_g  ${\rho_g}$     (long_name='Persistence of Government Spending Shock')
           sig_u  ${\sigma_u}$   (long_name='Standard Deviation of Permanent TFP Shock')
           sig_nu ${\sigma_\nu}$ (long_name='Standard Deviation of Noise Shock') % this was sig_s in the original code
           sig_d  ${\sigma_d}$   (long_name='Standard Deviation of IST Shock') 
           sig_q  ${\sigma_q}$   (long_name='Standard Deviation of Monetary Policy Shock')
           sig_p  ${\sigma_p}$   (long_name='Standard Deviation of Price Markup Shock')
           sig_w  ${\sigma_w}$   (long_name='Standard Deviation of Wage Markup Shock')
           sig_g  ${\sigma_g}$   (long_name='Standard Deviation of Government Spending Shock')
           zet    ${\zeta}$      (long_name='Inverse-Frisch Elasticity') % this and psi were wrongly swapped around in the original code
           h      ${h}$          (long_name='Habit Formation')
           del    ${\delta}$     (long_name='Depreciation Rate')
           bet    ${\beta}$      (long_name='Discount Factor')
           alp    ${\alpha}$     (long_name='Capital Share')
           chi    ${\chi}$       (long_name='Capital Adjustment Cost')
           xi     ${\xi}$        (long_name='Capital Utilisation Adjustment Cost') 
           iot    ${\iota}$      (long_name='Inflation Target')
           iot_w  ${\iota_w}$    (long_name='Wage Inflation Target')
           psi    ${\psi}$       (long_name='Steady state government spending ratio')      
           cal    ${\theta}$     (long_name='Calvo Parameter') 
           kap    ${\kappa}$     (long_name='NKPC Slope')
           cal_w  ${\theta_w}$   (long_name='Calvo Parameter for Wages')
           kap_w  ${\kappa_w}$   (long_name='NKPC Slope for Wages') 
           mu_p   ${\mu_p}$      (long_name='Demand Elasticity')
           psi_p  ${\psi_p}$     (long_name='Lagged Price Markup Shock Coefficient') % this was nu_p in original code  
           mu_w   ${\mu_w}$      (long_name='Wage Elasticity')
           psi_w  ${\psi_w}$     (long_name='Lagged Wage Markup Shock Coefficient') % this was nu_w in original code 
           rho_r  ${\rho_r}$     (long_name='Taylor Rule Persistence')
           gam_pi ${\gamma_\pi}$ (long_name='Taylor Rule Inflation Coefficient')
           gam_y  ${\gamma_y}$   (long_name='Taylor Rule Output Coefficient')
           BM_11  ${b_{11}}$     (long_name='Exogenous State Vector Coefficient B Matrix Element')
           BM_12  ${b_{12}}$     (long_name='Exogenous State Vector Coefficient B Matrix Element')
           BM_21  ${b_{21}}$     (long_name='Exogenous State Vector Coefficient B Matrix Element')
           BM_22  ${b_{22}}$     (long_name='Exogenous State Vector Coefficient B Matrix Element')
           BM_31  ${b_{31}}$     (long_name='Exogenous State Vector Coefficient B Matrix Element')
           BM_32  ${b_{32}}$     (long_name='Exogenous State Vector Coefficient B Matrix Element')
           FA_11  ${f_{11}}$     (long_name='Endogenous State Vector Coefficient F Matrix Element')
           FA_12  ${f_{12}}$     (long_name='Endogenous State Vector Coefficient F Matrix Element')
           FA_13  ${f_{13}}$     (long_name='Endogenous State Vector Coefficient F Matrix Element')
           FA_21  ${f_{21}}$     (long_name='Endogenous State Vector Coefficient F Matrix Element')
           FA_22  ${f_{22}}$     (long_name='Endogenous State Vector Coefficient F Matrix Element')
           FA_23  ${f_{23}}$     (long_name='Endogenous State Vector Coefficient F Matrix Element')
           CM_11  ${c_{11}}$     (long_name='Exogenous State Vector Coefficient C Matrix Element')
           CM_12  ${c_{12}}$     (long_name='Exogenous State Vector Coefficient C Matrix Element')
           CM_21  ${c_{21}}$     (long_name='Exogenous State Vector Coefficient C Matrix Element')     
           CM_22  ${c_{22}}$     (long_name='Exogenous State Vector Coefficient C Matrix Element')
;

% Parameters are initialised in the steady state file
% Need to convert matrix notation of Kalman filter in model block in 
% original code to scalars. This is also done in the steady state file.

model(linear);
[name='Permament TFP']
xh = (1+rho)*xh(-1) - rho*xhh(-1) + BM_11*e_1 + BM_12*e_2;

[name='Lagged permament TFP']
xhh = xh(-1) + BM_21*e_1 + BM_22*e_2;

[name='Temporary TFP']
zh = rho*zh(-1) + BM_31*e_1 + BM_32*e_2;

[name='Productivity']
pty = FA_11*xh(-1) + FA_12*xhh(-1) + FA_13*zh(-1) + CM_11*e_1 + CM_12*e_2;

[name='TFP signal']
s = FA_21*xh(-1) + FA_22*xhh(-1) + FA_23*zh(-1) + CM_21*e_1 + CM_22*e_2;

[name='Productivity growth']
da = pty - pty(-1);

[name='Productivity growth process']
d = rho_d*d(-1) + sig_d*eps_d;

[name='Monetary policy shock']
q = rho_q*q(-1) + sig_q*eps_q;

[name='Price markups']
m_p = rho_p*m_p(-1) + sig_p*(epsma_p - psi_p*epsma_p(-1));

[name='Price markups shock']
eps_p = epsma_p;

[name='Wage markups']
m_w = rho_w*m_w(-1) + sig_w*(epsma_w - psi_w*epsma_w(-1));

[name='Wage markups shock']
eps_w = epsma_w;

[name='Government LoM']
g = rho_g*g(-1) + sig_g*eps_g;

% lam = -c;
[name='Euler equation']
lam = (h*bet/((1-h*bet)*(1-h)))*c(+1) - ((1+h^2*bet)/((1-h*bet)*(1-h)))*c + (h/((1-h*bet)*(1-h)))*c(-1) + (h*bet/((1-h*bet)*(1-h)))*da(+1) - (h/((1-h*bet)*(1-h)))*da;

[name='Euler equation (capital returns)']
lam = lam(+1) - da(+1) + ri;

[name='Fisher equation']
ri = r - pi(+1);

[name='Capital optimality condition']
phi = (1-del)*bet*(phi(+1)-da(+1)) + (1-(1-del)*bet)*(lam(+1)-da(+1)+rk(+1));

[name='Investment optimality condition']
lam = phi + d - chi*(i-i(-1)+da) + bet*chi*(i(+1)-i+da(+1));

[name='Optimal capital utilisation']
rk = xi*u;

[name='Capital utilisation']
kk = u + kbar(-1) - da;

[name='Capital accumulation']
kbar = (1-del)*(kbar(-1)-da) + del*(d+i);

[name='Production']
y = alp*kk + (1-alp)*n;

[name='Inflation LoM']
pi = (iot/(1+iot*bet))*pi(-1) + (bet/(1+iot*bet))*pi(+1) + (((1-cal*bet)*(1-cal))/(cal*(1+iot*bet)))*mc + m_p;

[name='Marginal cost']
mc = alp*rk + (1-alp)*w;

[name='Capital-labour ratio condition']
kk - n = w - rk;

[name='Wage LoM']
w = (1/(1+bet))*w(-1) + (bet/(1+bet))*w(+1) - (1/(1+bet))*(pi + da) + bet/(1+bet)*(pi(+1) + da(+1)) - kap_w*mc_w + kap_w*m_w;

[name='Wage Marginal Cost']
mc_w = w - zet*n + lam;

[name='Taylor rule']
r = rho_r*r(-1) + (1-rho_r)*(gam_pi*pi + gam_y*y) + q;

[name='Output growth']
dy = y - y(-1);

% Other equations
#RkP = 1/bet - (1-del);
#WAP = ((alp^alp*(1-alp)^(1-alp)/(1+mu_p))*RkP^(-alp))^(1/(1-alp));
#KAN = (WAP/RkP)*(alp/(1-alp));
#IY = del*(KAN)^(1-alp);
#RkKPY = RkP*KAN^(1-alp);

[name='Resource constraint']
(1/psi - IY - RkKPY)*c + (IY)*i + RkKPY*u + (1/psi)*g = (1/psi)*y;

[name='Consumption + TFP']
ca = c + pty;

[name='Investment + TFP']
ia = i + pty;

[name='Output + TFP']
ya = y + pty;

[name='Capital + TFP']
ka = kk + pty;

[name='Wages + TFP']
wa = w + pty;

[name='Consumption growth']
dca = ca - ca(-1);

[name='Investment growth']
diad = ia - ia(-1);

[name='Output growth']
dya = ya - ya(-1);

[name='Labour growth']
dn = n - n(-1);

[name='Wages growth']
dwa = wa - wa(-1);

[name='Real interest rate change']
dr = r - r(-1);

[name='Inflation change']
dpi = pi - pi(-1);

end;


steady;

shocks;
var e_1; stderr 1;
var e_2; stderr 1;
var eps_d; stderr 1;
var eps_q; stderr 1;
var eps_p; stderr 1;
var eps_w; stderr 1;
var eps_g; stderr 1;
end;

stoch_simul(periods=2000,IRF=100,nograph);

% Need to create our own IRF_matrix since we are using Dynare newer than version 3.065
% Define the number of periods, shocks, and variables
T = 100; % Number of periods
shock_names = cellstr(M_.exo_names); % Get the shock names from M_.exo_names
var_names = cellstr(M_.endo_names); % Get the variable names from M_.endo_names
n_shock = length(shock_names); % Number of shocks
n_var = length(var_names); % Number of variables

% Initialise the IRF_matrix
IRF_matrix = zeros(T, n_shock, n_var);

% Fill the IRF_matrix with the IRFs from oo_.irfs
for j = 1:n_var
    for i = 1:n_shock
        var_name = var_names{j};
        shock_name = shock_names{i};
        field_name = [var_name '_' shock_name];
        if isfield(oo_.irfs, field_name)
            IRF_matrix(:, i, j) = oo_.irfs.(field_name);
        else
            warning(['Field ' field_name ' not found in oo_.irfs']);
        end
    end
end

% Save the IRF_matrix to a .mat file
save('IRF_matrix.mat', 'IRF_matrix');