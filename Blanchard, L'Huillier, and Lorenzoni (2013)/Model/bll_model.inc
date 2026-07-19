% Shared scalar representation of the BLL model. The original pre-Dynare
% 4.5 files put matrices directly in the model block; BM, FA and CM are the
% scalar elements of those matrices and are updated by the steady-state file.
model(linear);
[name='Permanent TFP estimate']
xh = (1+rho)*xh(-1) - rho*xhh(-1) + BM_11*e_1 + BM_12*e_2;

[name='Lagged permanent TFP estimate']
xhh = xh(-1) + BM_21*e_1 + BM_22*e_2;

[name='Transitory TFP estimate']
zh = rho*zh(-1) + BM_31*e_1 + BM_32*e_2;

[name='Observed productivity']
pty = FA_11*xh(-1) + FA_12*xhh(-1) + FA_13*zh(-1) + CM_11*e_1 + CM_12*e_2;

[name='Observed TFP signal']
s = FA_21*xh(-1) + FA_22*xhh(-1) + FA_23*zh(-1) + CM_21*e_1 + CM_22*e_2;

[name='Productivity growth']
da = pty - pty(-1);

[name='Investment-specific technology process']
d = rho_d*d(-1) + sig_d*eps_d;

[name='Monetary policy shock process']
q = rho_q*q(-1) + sig_q*eps_q;

[name='Price markup process']
m_p = rho_p*m_p(-1) + sig_p*(epsma_p - psi_p*epsma_p(-1));

[name='Price markup innovation']
eps_p = epsma_p;

[name='Wage markup process']
m_w = rho_w*m_w(-1) + sig_w*(epsma_w - psi_w*epsma_w(-1));

[name='Wage markup innovation']
eps_w = epsma_w;

[name='Government spending process']
g = rho_g*g(-1) + sig_g*eps_g;

[name='Marginal utility of consumption']
lam = (h*bet/((1-h*bet)*(1-h)))*c(+1)
    - ((1+h^2*bet)/((1-h*bet)*(1-h)))*c
    + (h/((1-h*bet)*(1-h)))*c(-1)
    + (h*bet/((1-h*bet)*(1-h)))*da(+1)
    - (h/((1-h*bet)*(1-h)))*da;

[name='Euler equation']
lam = lam(+1) - da(+1) + ri;

[name='Fisher equation']
ri = r - pi(+1);

[name='Capital optimality condition']
phi = (1-del)*bet*(phi(+1)-da(+1))
    + (1-(1-del)*bet)*(lam(+1)-da(+1)+rk(+1));

[name='Investment optimality condition']
lam = phi + d - chi*(i-i(-1)+da) + bet*chi*(i(+1)-i+da(+1));

[name='Optimal capital utilisation']
rk = xi*u;

[name='Capital services']
kk = u + kbar(-1) - da;

[name='Capital accumulation']
kbar = (1-del)*(kbar(-1)-da) + del*(d+i);

[name='Production']
y = alp*kk + (1-alp)*n;

[name='Price Phillips curve']
pi = (iot/(1+iot*bet))*pi(-1)
    + (bet/(1+iot*bet))*pi(+1)
    + (((1-cal*bet)*(1-cal))/(cal*(1+iot*bet)))*mc
@#if PAPER_SPECIFICATION == 1
    + (((1-cal*bet)*(1-cal))/(cal*(1+iot*bet)))*m_p;
@#else
    + m_p;
@#endif

[name='Marginal cost']
mc = alp*rk + (1-alp)*w;

[name='Capital-labour ratio']
kk - n = w - rk;

% m_w is the reduced-form cost-push shock used in the authors' estimation
% code. Multiplying it by kappa_w would change the units of sig_w and would
% not reproduce Table 5. The slope itself is evaluated from the current
% values of cal_w and zet at every likelihood evaluation.
[name='Wage Phillips curve']
w = (1/(1+bet))*w(-1)
    + (bet/(1+bet))*w(+1)
    + (iot_w/(1+bet))*pi(-1)
    - ((1+bet*iot_w)/(1+bet))*pi
    + (bet/(1+bet))*pi(+1)
    + (iot_w/(1+bet))*da(-1)
    - ((1+bet*iot_w)/(1+bet))*da
    + (bet/(1+bet))*da(+1)
    - (((1-cal_w*bet)*(1-cal_w))/(cal_w*(1+bet)*(1+zet*(1+1/mu_w))))*mc_w
@#if PAPER_SPECIFICATION == 1
    + (((1-cal_w*bet)*(1-cal_w))/(cal_w*(1+bet)*(1+zet*(1+1/mu_w))))*m_w;
@#else
    + m_w;
@#endif

[name='Wage markup gap']
mc_w = w - zet*n + lam;

[name='Taylor rule']
r = rho_r*r(-1) + (1-rho_r)*(gam_pi*pi + gam_y*y) + q;

[name='Detrended output growth']
dy = y - y(-1);

#RkP = 1/bet - (1-del);
#WAP = ((alp^alp*(1-alp)^(1-alp)/(1+mu_p))*RkP^(-alp))^(1/(1-alp));
#KAN = (WAP/RkP)*(alp/(1-alp));
#IY = del*(KAN)^(1-alp);
#RkKPY = RkP*KAN^(1-alp);

@#if PAPER_SPECIFICATION == 1
#CY = 1-IY-psi;
[name='Resource constraint (Online Appendix equation 49)']
(1-psi)*y = CY*c + IY*i + RkKPY*u + g;
@#else
[name='Resource constraint (authors computational specification)']
(1/cy_ratio-IY-RkKPY)*c + IY*i + RkKPY*u + (1/cy_ratio)*g = (1/cy_ratio)*y;
@#endif

[name='Consumption with TFP trend']
ca = c + pty;

[name='Investment with TFP trend']
ia = i + pty;

[name='Output with TFP trend']
ya = y + pty;

[name='Capital with TFP trend']
ka = kk + pty;

[name='Real wage with TFP trend']
wa = w + pty;

[name='Consumption growth']
dca = ca - ca(-1);

[name='Investment growth']
diad = ia - ia(-1);

[name='Output growth']
dya = ya - ya(-1);

[name='Labour growth']
dn = n - n(-1);

[name='Real wage growth']
dwa = wa - wa(-1);

[name='Nominal interest rate change']
dr = r - r(-1);

[name='Inflation change']
dpi = pi - pi(-1);
end;
