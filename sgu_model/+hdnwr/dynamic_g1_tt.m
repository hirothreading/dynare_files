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

assert(length(T) >= 21);

T = hdnwr.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(19) = getPowerDeriv(T(2)/(1+y(11)),(-params(5)),1);
T(20) = getPowerDeriv(T(11),params(3),1);
T(21) = T(16)*params(7)/(params(7)*(1+1/params(3)))+T(15)*(-(params(7)*(params(6)+params(7))))/((params(6)+params(7)*y(4))*(params(6)+params(7)*y(4)))*getPowerDeriv(T(14),1+1/params(3),1);

end
