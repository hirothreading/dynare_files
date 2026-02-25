function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 43);

T = hdnwr_nl.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(22) = T(1)*params(7)/(1+y(11));
T(23) = getPowerDeriv(T(2)/(1+y(11)),(-params(5)),1);
T(24) = T(10)*T(22)*T(23);
T(25) = getPowerDeriv(T(13),params(3),1);
T(26) = T(1)*params(7)*getPowerDeriv(T(2),1-params(5),1);
T(27) = params(7)/(params(7)*T(4));
T(28) = (-(params(7)*(params(6)+params(7))))/((params(6)+params(7)*y(4))*(params(6)+params(7)*y(4)));
T(29) = getPowerDeriv(T(16),T(4),1);
T(30) = T(28)*T(29);
T(31) = T(18)*T(27)+T(17)*T(30);
T(32) = getPowerDeriv(T(16),1-params(5),1);
T(33) = T(28)*T(32);
T(34) = (1-params(10))*(T(31)-(T(21)*params(7)/(params(7)*(1-params(5)))+T(20)*T(33)));
T(35) = 1/(steady_state(2))*getPowerDeriv(y(5)/(steady_state(2)),params(12),1);
T(36) = getPowerDeriv(y(5),params(2),1);
T(37) = getPowerDeriv(y(14),(-params(2)),1);
T(38) = getPowerDeriv(y(6),params(4),1);
T(39) = getPowerDeriv(y(6),params(4)-1,1);
T(40) = T(12)*1/(1-params(10));
T(41) = (1+params(9))/params(1)*1/(1+params(9))*getPowerDeriv((1+y(10))/(1+params(9)),params(11),1);
T(42) = (-T(2))/((1+y(11))*(1+y(11)));
T(43) = T(10)*T(23)*T(42);

end
