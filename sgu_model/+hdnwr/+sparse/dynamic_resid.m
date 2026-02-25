function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(18, 1);
end
[T_order, T] = hdnwr.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(10, 1);
    residual(1) = (y(12)) - (y(20)*T(4));
    residual(2) = (y(12)^(-params(2))) - (params(1)*(1+y(16))*T(5)/(1+y(27)));
    residual(3) = (y(20)*params(4)*T(6)) - (y(15));
    residual(4) = (1+y(16)) - (T(7)*T(8)*y(19));
    residual(5) = (1+y(18)) - ((1+y(17))*y(15)/y(5));
    residual(6) = (T(9)*T(12)) - (T(2)*y(5)/(1+y(17)));
    residual(7) = ((1+y(18))^(1-params(5))) - (y(11)*T(13)+T(3)*((params(6)+params(7))^(2-params(5))-(params(6)+params(7)*y(11))^(2-params(5))));
    residual(8) = (y(14)) - (params(10)+(1-params(10))*(T(17)-(params(6)+params(7)*y(11))/(params(7)*(1-params(5)))*T(18))/(y(11)+T(17)));
    residual(9) = (log(y(19))) - (params(13)*log(y(9))+x(1));
    residual(10) = (log(y(20))) - (params(14)*log(y(10))+x(2));
end
