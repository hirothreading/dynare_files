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
    T = be_simplified.static_resid_tt(T, y, x, params);
end
residual = zeros(7, 1);
    residual(1) = (y(1)-y(4)) - (y(1)-y(4)-1/params(1)*(y(3)-y(2)));
residual(2) = y(2)-(params(4)*(y(1)+params(10)*y(7))+params(5)*y(5)+y(2)*params(2));
    residual(3) = (y(3)) - (y(2)*params(3)+y(6));
    residual(4) = (y(4)) - (y(4)*params(11)+x(1));
    residual(5) = (y(5)) - (y(5)*params(12)+x(2));
    residual(6) = (y(6)) - (y(6)*params(13)+x(3));
    residual(7) = (y(7)) - (y(7)*params(14)+x(4));

end
