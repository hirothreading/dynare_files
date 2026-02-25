function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = hdnwr_nl.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(10, 17);
g1(1,5)=1;
g1(1,6)=(-(y(13)*T(38)));
g1(1,13)=(-T(5));
g1(2,5)=getPowerDeriv(y(5),(-params(2)),1);
g1(2,14)=(-(params(1)*(1+y(9))*T(37)/(1+y(15))));
g1(2,9)=(-(params(1)*T(6)/(1+y(15))));
g1(2,15)=(-((-(params(1)*(1+y(9))*T(6)))/((1+y(15))*(1+y(15)))));
g1(3,6)=y(13)*params(4)*T(39);
g1(3,8)=(-1);
g1(3,13)=params(4)*T(7);
g1(4,5)=(-(y(12)*T(8)*T(35)));
g1(4,9)=1;
g1(4,10)=(-(y(12)*T(9)*T(41)));
g1(4,12)=(-(T(8)*T(9)));
g1(5,1)=(-((1+y(10))*(-y(8))/(y(1)*y(1))));
g1(5,8)=(-((1+y(10))*1/y(1)));
g1(5,10)=(-(y(8)/y(1)));
g1(5,11)=1;
g1(6,4)=T(11)*T(24)*T(25)-y(1)*T(1)*params(7)/(1+y(10));
g1(6,5)=T(14)*T(36);
g1(6,6)=T(11)*T(25)*T(40);
g1(6,1)=(-(T(2)/(1+y(10))));
g1(6,10)=(-((-(T(2)*y(1)))/((1+y(10))*(1+y(10)))));
g1(6,11)=T(11)*T(25)*T(43);
g1(7,4)=(-(T(15)+y(4)*T(26)+T(3)*(-(params(7)*getPowerDeriv(params(6)+params(7)*y(4),2-params(5),1)))));
g1(7,11)=getPowerDeriv(1+y(11),1-params(5),1);
g1(8,4)=(-(((y(4)+T(19))*T(34)-(1-params(10))*(T(19)-T(20)*T(21))*(1+T(31)))/((y(4)+T(19))*(y(4)+T(19)))));
g1(8,7)=1;
g1(9,2)=(-(params(13)*1/y(2)));
g1(9,12)=1/y(12);
g1(9,16)=(-1);
g1(10,3)=(-(params(14)*1/y(3)));
g1(10,13)=1/y(13);
g1(10,17)=(-1);

end
