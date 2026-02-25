function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
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
%   residual
%

if T_flag
    T = hdnwr_nl.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(10, 1);
    residual(1) = (y(5)) - (y(13)*T(5));
    residual(2) = (y(5)^(-params(2))) - (params(1)*(1+y(9))*T(6)/(1+y(15)));
    residual(3) = (y(13)*params(4)*T(7)) - (y(8));
    residual(4) = (1+y(9)) - (T(8)*T(9)*y(12));
    residual(5) = (1+y(11)) - ((1+y(10))*y(8)/y(1));
    residual(6) = (T(11)*T(14)) - (T(2)*y(1)/(1+y(10)));
    residual(7) = ((1+y(11))^(1-params(5))) - (y(4)*T(15)+T(3)*((params(6)+params(7))^(2-params(5))-(params(6)+params(7)*y(4))^(2-params(5))));
    residual(8) = (y(7)) - (params(10)+(1-params(10))*(T(19)-T(20)*T(21))/(y(4)+T(19)));
    residual(9) = (log(y(12))) - (params(13)*log(y(2))+x(it_, 1));
    residual(10) = (log(y(13))) - (params(14)*log(y(3))+x(it_, 2));

end
