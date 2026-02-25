function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 21);

T(1) = (1+params(9))^params(8);
T(2) = T(1)*(params(6)+params(7)*y(4));
T(3) = (1+params(9))^(params(8)*(1-params(5)))/(params(7)*(2-params(5)));
T(4) = 1+1/params(3);
T(5) = y(6)^params(4);
T(6) = y(14)^(-params(2));
T(7) = y(6)^(params(4)-1);
T(8) = (1+params(9))/params(1)*((1+y(10))/(1+params(9)))^params(11);
T(9) = (y(5)/(steady_state(2)))^params(12);
T(10) = y(6)/(1-params(10));
T(11) = y(5)^params(2);
T(12) = (T(2)/(1+y(11)))^(-params(5));
T(13) = T(10)*T(12);
T(14) = T(13)^params(3);
T(15) = T(2)^(1-params(5));
T(16) = (params(6)+params(7))/(params(6)+params(7)*y(4));
T(17) = (params(6)+params(7)*y(4))/(params(7)*T(4));
T(18) = T(16)^T(4)-1;
T(19) = T(17)*T(18);
T(20) = (params(6)+params(7)*y(4))/(params(7)*(1-params(5)));
T(21) = T(16)^(1-params(5))-1;

end
