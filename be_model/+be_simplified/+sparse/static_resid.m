function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = be_simplified.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(7, 1);
    residual(1) = (y(1)-y(4)) - (y(1)-y(4)-1/params(1)*(y(3)-y(2)));
residual(2) = y(2)-(params(4)*(y(1)+params(10)*y(7))+params(5)*y(5)+y(2)*params(2));
    residual(3) = (y(3)) - (y(2)*params(3)+y(6));
    residual(4) = (y(4)) - (y(4)*params(11)+x(1));
    residual(5) = (y(5)) - (y(5)*params(12)+x(2));
    residual(6) = (y(6)) - (y(6)*params(13)+x(3));
    residual(7) = (y(7)) - (y(7)*params(14)+x(4));
end
