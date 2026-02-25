function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(18, 1);
end
[T_order, T] = hdnwr.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(10, 1);
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
