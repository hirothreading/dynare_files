function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
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
%   residual
%

if T_flag
    T = hdnwr.static_resid_tt(T, y, x, params);
end
residual = zeros(10, 1);
    residual(1) = (y(2)) - (y(10)*T(9));
    residual(2) = (T(10)) - (T(10)*params(1)*(1+y(6))/(1+y(7)));
    residual(3) = (y(10)*params(4)*T(11)) - (y(5));
    residual(4) = (1+y(6)) - (T(12)*T(13)*y(9));
    residual(5) = (1+y(8)) - (1+y(7));
    residual(6) = (T(14)*T(17)) - (T(2)*y(5)/(1+y(7)));
    residual(7) = ((1+y(8))^(1-params(5))) - (T(3)*((params(6)+params(7))^(2-params(5))-(params(6)+params(7)*y(1))^(2-params(5)))+y(1)*T(18));
    residual(8) = (y(4)) - (params(10)+(1-params(10))*(T(7)-(params(6)+params(7)*y(1))/(params(7)*(1-params(5)))*T(8))/(y(1)+T(7)));
    residual(9) = (log(y(9))) - (log(y(9))*params(13)+x(1));
    residual(10) = (log(y(10))) - (log(y(10))*params(14)+x(2));

end
