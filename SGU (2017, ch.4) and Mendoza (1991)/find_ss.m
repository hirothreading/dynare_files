%Find steady state of SGU (2017, ch.4)/Mendoza (1991)
% Code written by David Murakami
 clear all;
 clc;
 
% Declare parameters
global betta sig del omega alphha rhoA sigA phii psii kap
global barh bark barRstar bary barnxy bard barc

betta = 0.96; 
sig = 2;
del = 0.1;
omega = 1.455;
alphha = 0.32;
rhoA = 0.42;
sigA = 1.29;
phii = 0.028;
psii = 0.000742;
kap = ((betta^-1 - 1 + del)/alphha)^(1/(alphha-1));
barh = ((1-alphha)*kap^alphha)^(1/(omega-1)); 
bark =kap*barh;
barRstar = 1/betta;
bary = ((1-alphha)*kap^(alphha*omega))^(1/(omega-1));
barnxy = 0.02;
bard = barnxy*bary/(barRstar-1);
barc = (1-barRstar)*bard + kap^alphha*barh - del*bark;

% Steady state model (annual)
lA = 1;
ls = 1;
h = barh;
lk = bark;
lR = barRstar;
ly = bary;
ld = bard;
lc = barc;
llam = (barc-barh^omega/omega)^(-sig);
lw = h^(omega-1);
li = lk - (1-del)*lk;
rk = alphha*lA*lk^(alphha-1)*h^(1-alphha);
lpi = rk*lk - li - phii/2*(lk-lk)^2;
lps = (betta*llam*lpi)/(llam - betta*llam);
lnx = ly - lc -li - phii/2*(lk-lk)^2;
ca = 0;

fprintf('\n//Steady state values:\n');
fprintf('\nld = %10.9f;\n',ld);
fprintf('\nlc = %10.9f;\n',lc);
fprintf('\nlps = %10.9f;\n',lps);
fprintf('\nls = %10.9f;\n',ls);
fprintf('\nlR = %10.9f;\n',lR);
fprintf('\nlpi = %10.9f;\n',lpi);
fprintf('\nh = %10.9f;\n',h);
fprintf('\nlw = %10.9f;\n',lw);
fprintf('\nllam = %10.9f;\n',llam);
fprintf('\nly = %10.9f;\n',ly);
fprintf('\nrk = %10.9f;\n',rk);
fprintf('\nlk = %10.9f;\n',lk);
fprintf('\nli = %10.9f;\n',li);
fprintf('\nlA = %10.9f;\n',lA);
fprintf('\nlnx = %10.9f;\n',lnx);
fprintf('\nca = %10.9f;\n',ca);
