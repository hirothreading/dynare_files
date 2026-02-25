function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
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

assert(length(T) >= 22);

T = hdnwr_nl.static_resid_tt(T, y, x, params);

T(19) = getPowerDeriv(T(2)/(1+y(8)),(-params(5)),1);
T(20) = getPowerDeriv(T(16),params(3),1);
T(21) = T(6)*params(7)/(params(7)*(1+1/params(3)))+T(5)*(-(params(7)*(params(6)+params(7))))/((params(6)+params(7)*y(1))*(params(6)+params(7)*y(1)))*getPowerDeriv(T(4),1+1/params(3),1);
T(22) = getPowerDeriv(y(2),(-params(2)),1);

end
