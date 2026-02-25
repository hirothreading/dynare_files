function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = hdnwr_nl.static_g1_tt(T, y, x, params);
end
g1 = zeros(10, 10);
g1(1,2)=1;
g1(1,3)=(-(y(10)*getPowerDeriv(y(3),params(4),1)));
g1(1,10)=(-T(9));
g1(2,2)=T(22)-params(1)*(1+y(6))*T(22)/(1+y(7));
g1(2,6)=(-(T(10)*params(1)/(1+y(7))));
g1(2,7)=(-((-(T(10)*params(1)*(1+y(6))))/((1+y(7))*(1+y(7)))));
g1(3,3)=y(10)*params(4)*getPowerDeriv(y(3),params(4)-1,1);
g1(3,5)=(-1);
g1(3,10)=params(4)*T(11);
g1(4,2)=(-(y(9)*T(12)*((y(2))-y(2))/((y(2))*(y(2)))*getPowerDeriv(y(2)/(y(2)),params(12),1)));
g1(4,6)=1;
g1(4,7)=(-(y(9)*T(13)*(1+params(9))/params(1)*1/(1+params(9))*getPowerDeriv((1+y(7))/(1+params(9)),params(11),1)));
g1(4,9)=(-(T(12)*T(13)));
g1(5,7)=(-1);
g1(5,8)=1;
g1(6,1)=T(14)*y(3)/(1-params(10))*T(1)*params(7)/(1+y(8))*T(19)*T(20)-y(5)*T(1)*params(7)/(1+y(7));
g1(6,2)=T(17)*getPowerDeriv(y(2),params(2),1);
g1(6,3)=T(14)*T(20)*T(15)*1/(1-params(10));
g1(6,5)=(-(T(2)/(1+y(7))));
g1(6,7)=(-((-(T(2)*y(5)))/((1+y(7))*(1+y(7)))));
g1(6,8)=T(14)*T(20)*y(3)/(1-params(10))*T(19)*(-T(2))/((1+y(8))*(1+y(8)));
g1(7,1)=(-(T(3)*(-(params(7)*getPowerDeriv(params(6)+params(7)*y(1),2-params(5),1)))+T(18)+y(1)*T(1)*params(7)*getPowerDeriv(T(2),1-params(5),1)));
g1(7,8)=getPowerDeriv(1+y(8),1-params(5),1);
g1(8,1)=(-(((y(1)+T(7))*(1-params(10))*(T(21)-(T(8)*params(7)/(params(7)*(1-params(5)))+(params(6)+params(7)*y(1))/(params(7)*(1-params(5)))*(-(params(7)*(params(6)+params(7))))/((params(6)+params(7)*y(1))*(params(6)+params(7)*y(1)))*getPowerDeriv(T(4),1-params(5),1)))-(1-params(10))*(T(7)-(params(6)+params(7)*y(1))/(params(7)*(1-params(5)))*T(8))*(1+T(21)))/((y(1)+T(7))*(y(1)+T(7)))));
g1(8,4)=1;
g1(9,9)=1/y(9)-params(13)*1/y(9);
g1(10,10)=1/y(10)-params(14)*1/y(10);

end
