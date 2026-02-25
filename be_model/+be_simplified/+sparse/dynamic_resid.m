function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = be_simplified.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(7, 1);
    residual(1) = (y(8)-y(11)) - (y(15)-y(18)-1/params(1)*(y(10)-y(16)));
residual(2) = (y(9)-(params(4)*(y(8)+params(10)*y(14))+params(5)*y(12)+y(16)*params(2)))*(1-params(15))+params(15)*(y(9)-(y(16)*params(2)+(y(8)+params(10)*y(14))*params(6)-params(8)+y(12)*params(7)));
    residual(3) = (y(10)) - (y(9)*params(3)+y(13));
    residual(4) = (y(11)) - (params(11)*y(4)+x(1));
    residual(5) = (y(12)) - (params(12)*y(5)+x(2));
    residual(6) = (y(13)) - (params(13)*y(6)+x(3));
    residual(7) = (y(14)) - (params(14)*y(7)+x(4));
end
