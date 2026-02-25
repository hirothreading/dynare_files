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
    T = hdnwr.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(10, 17);
g1(1,5)=1;
g1(1,6)=(-(y(13)*getPowerDeriv(y(6),params(4),1)));
g1(1,13)=(-T(4));
g1(2,5)=getPowerDeriv(y(5),(-params(2)),1);
g1(2,14)=(-(params(1)*(1+y(9))*getPowerDeriv(y(14),(-params(2)),1)/(1+y(15))));
g1(2,9)=(-(params(1)*T(5)/(1+y(15))));
g1(2,15)=(-((-(params(1)*(1+y(9))*T(5)))/((1+y(15))*(1+y(15)))));
g1(3,6)=y(13)*params(4)*getPowerDeriv(y(6),params(4)-1,1);
g1(3,8)=(-1);
g1(3,13)=params(4)*T(6);
g1(4,5)=(-(y(12)*T(7)*1/(steady_state(2))*getPowerDeriv(y(5)/(steady_state(2)),params(12),1)));
g1(4,9)=1;
g1(4,10)=(-(y(12)*T(8)*(1+params(9))/params(1)*1/(1+params(9))*getPowerDeriv((1+y(10))/(1+params(9)),params(11),1)));
g1(4,12)=(-(T(7)*T(8)));
g1(5,1)=(-((1+y(10))*(-y(8))/(y(1)*y(1))));
g1(5,8)=(-((1+y(10))*1/y(1)));
g1(5,10)=(-(y(8)/y(1)));
g1(5,11)=1;
g1(6,4)=T(9)*y(6)/(1-params(10))*T(1)*params(7)/(1+y(11))*T(19)*T(20)-y(1)*T(1)*params(7)/(1+y(10));
g1(6,5)=T(12)*getPowerDeriv(y(5),params(2),1);
g1(6,6)=T(9)*T(20)*T(10)*1/(1-params(10));
g1(6,1)=(-(T(2)/(1+y(10))));
g1(6,10)=(-((-(T(2)*y(1)))/((1+y(10))*(1+y(10)))));
g1(6,11)=T(9)*T(20)*y(6)/(1-params(10))*T(19)*(-T(2))/((1+y(11))*(1+y(11)));
g1(7,4)=(-(T(13)+y(4)*T(1)*params(7)*getPowerDeriv(T(2),1-params(5),1)+T(3)*(-(params(7)*getPowerDeriv(params(6)+params(7)*y(4),2-params(5),1)))));
g1(7,11)=getPowerDeriv(1+y(11),1-params(5),1);
g1(8,4)=(-(((y(4)+T(17))*(1-params(10))*(T(21)-(T(18)*params(7)/(params(7)*(1-params(5)))+(params(6)+params(7)*y(4))/(params(7)*(1-params(5)))*(-(params(7)*(params(6)+params(7))))/((params(6)+params(7)*y(4))*(params(6)+params(7)*y(4)))*getPowerDeriv(T(14),1-params(5),1)))-(1-params(10))*(T(17)-(params(6)+params(7)*y(4))/(params(7)*(1-params(5)))*T(18))*(1+T(21)))/((y(4)+T(17))*(y(4)+T(17)))));
g1(8,7)=1;
g1(9,2)=(-(params(13)*1/y(2)));
g1(9,12)=1/y(12);
g1(9,16)=(-1);
g1(10,3)=(-(params(14)*1/y(3)));
g1(10,13)=1/y(13);
g1(10,17)=(-1);

end
