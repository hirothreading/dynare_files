% Find the steady state for the simple SAM model
% David Murakami (University of Oxford)

clear all;
clc; 

% Declare parameters 
global betta s mu alppha U_ss q_ss y_ss rho_y siggma_y    
betta    = 0.99;
s        = 0.04;
mu       = 0.5;
alppha   = 2/3;
U_ss     = 0.06;
q_ss     = 0.7;
y_ss     = 1;
rho_y    = 0.95;
siggma_y = 0.01;

% Find ss values for: y w J V U M thetta f q
y = y_ss;
U = U_ss;
q = q_ss;
f = s/U - s;
thetta = f/q;
M = f*U;
V = M/q;
w = alppha*y;
J = (y - w)/(1-betta*(1-s));

% Print values
fprintf('\n//Steady state values:\n');
fprintf('\ny = %10.9f;\n',y);
fprintf('\nU = %10.9f;\n',U);
fprintf('\nq = %10.9f;\n',q);
fprintf('\nf = %10.9f;\n',f);
fprintf('\nthetta = %10.9f;\n',thetta);
fprintf('\nM = %10.9f;\n',M);
fprintf('\nV = %10.9f;\n',V);
fprintf('\nw = %10.9f;\n',w);
fprintf('\nJ = %10.9f;\n',J);