function [ys,params,check] = estimate_steadystate(ys,exo,M_,options_) 
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

% Initialise the matrices of the Kalman filter and then get scalars for modern Dynare
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