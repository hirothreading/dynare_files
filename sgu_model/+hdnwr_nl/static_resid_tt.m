function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 18);

T(1) = (1+params(9))^params(8);
T(2) = T(1)*(params(6)+params(7)*y(1));
T(3) = (1+params(9))^(params(8)*(1-params(5)))/(params(7)*(2-params(5)));
T(4) = (params(6)+params(7))/(params(6)+params(7)*y(1));
T(5) = (params(6)+params(7)*y(1))/(params(7)*(1+1/params(3)));
T(6) = T(4)^(1+1/params(3))-1;
T(7) = T(5)*T(6);
T(8) = T(4)^(1-params(5))-1;
T(9) = y(3)^params(4);
T(10) = y(2)^(-params(2));
T(11) = y(3)^(params(4)-1);
T(12) = (1+params(9))/params(1)*((1+y(7))/(1+params(9)))^params(11);
T(13) = (y(2)/(y(2)))^params(12);
T(14) = y(2)^params(2);
T(15) = (T(2)/(1+y(8)))^(-params(5));
T(16) = y(3)/(1-params(10))*T(15);
T(17) = T(16)^params(3);
T(18) = T(2)^(1-params(5));

end
