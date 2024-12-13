function [ys,params,check] = simulate_steadystate(ys,exo,M_,options_) 
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impose restrictions on parameters)

% Read parameters to access by their name
for ii = 1:M_.param_nbr
    eval([ M_.param_names{ii} ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

%       [rho    rho_d   rho_q   rho_p   rho_w   rho_g   sig_u   sig_nu   sig_d   sig_q   sig_p   sig_w   sig_g   h   alp   zet     chi     xi      cal     cal_w   iot     iot_w   mu_p    psi_p   mu_w     psi_w   rho_r   gam_pi  gam_y ]
param = [0.6	0.60	0.4	    0.6	    0.6	    0.6	    0.5	    1	     0.15	 0.15	 0.15	 0.15	 0.5	 0.5 0.3   2       4	   5	   0.66	   0.66	   0       0	   0.3000  0.77	   0.0500   0.91	0.5	    1.5	    0.005];

rho = param(1); 
rho_d = param(2);
rho_q = param(3);
rho_p = param(4);
rho_w = param(5);
rho_g = param(6);

sig_u = param(7);
sig_nu = param(8);
sig_d = param(9);
sig_q = param(10);
sig_p = param(11);
sig_w = param(12);
sig_g = param(13);

psi = 0.2; % BLL cite JPT for this value and they set it to 1.28
h = param(14);
del = .025; 
bet = .99; 
alp = param(15);
zet = param(16);
chi = param(17);
xi = param(18);
cal = param(19);
cal_w = param(20);
iot = param(21);
iot_w = param(22);
mu_p = param(23);
psi_p = param(24);
mu_w = param(25);
psi_w = param(26);

rho_r = param(27);
gam_pi = param(28);
gam_y = param(29);

kap = (1-cal*bet)*(1-cal)/(cal*(1+iot*bet));
kap_w = ((1-cal_w*bet)*(1-cal_w))/(cal_w*(1+bet)*(1+zet*(1+1/mu_w)));

A = [ 1+rho -rho 0 ; 1 0 0 ; 0 0 rho ];
S_1 = [ (1-rho)*sig_u 0 0 ; 0 0 0 ; 0 0 (rho^.5)*sig_u ]^2;
F = [ 1 0 1 ; 1 0 0 ];
S_2 = [0 0; 0 sig_nu ]^2;
P = .01*eye(3,3);
K = P*F'/(F*P*F' + S_2);

for iter = 1:100000;
    P = A*P*A' - A*P*F'*K'*A' + S_1;
    Ki = P*F'/(F*P*F' + S_2);
    dif = max(max(abs(Ki-K))); 
    K = Ki;
    if dif < 1e-15
        break;
    end
    if iter == 99999
        display('conv not achieved !!!!');
    end;
end;

IKFA = (eye(3,3) - K*F)*A;
FA = F*A;
CM = real((F*P*F' + S_2)^.5);
BM = K*CM;

BM_11 = BM(1,1);
BM_12 = BM(1,2);
BM_21 = BM(2,1);
BM_22 = BM(2,2);
BM_31 = BM(3,1);
BM_32 = BM(3,2);

FA_11 = FA(1,1);
FA_12 = FA(1,2);
FA_13 = FA(1,3);
FA_21 = FA(2,1);
FA_22 = FA(2,2);
FA_23 = FA(2,3);

CM_11 = CM(1,1);
CM_12 = CM(1,2);
CM_21 = CM(2,1);
CM_22 = CM(2,2);

% Have to declare variables:
xh = 0;
xhh = 0;
zh = 0;
pty = 0;
s = 0;
da = 0;
d = 0;
q = 0;
m_p = 0;
epsma_p = 0;
m_w = 0;
epsma_w = 0;
g = 0;
lam = 0;
phi = 0;
c = 0;
i = 0;
kbar = 0;
kk = 0;
u = 0;
y = 0;
n = 0;
rk = 0;
w = 0;
pi = 0;
mc = 0;
mc_w = 0;
r = 0;
ri = 0;
ca = 0;
ia = 0;
ya = 0;
ka = 0;
wa = 0;
dca = 0;
diad = 0;
dya = 0;
dn = 0;
dwa = 0;
dr = 0;
dpi = 0;
dy = 0;

% End of the steady state model block

params=NaN(M_.param_nbr,1);
for iter = 1:M_.param_nbr %update parameters set in this file
    eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

for ii = 1:M_.orig_endo_nbr
    eval(['ys(' int2str(ii) ') = ' M_.endo_names{ii} ';']);
end