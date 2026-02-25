function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 18
    T = [T; NaN(18 - size(T, 1), 1)];
end
T(1) = (1+params(9))^params(8);
T(2) = T(1)*(params(6)+params(7)*y(1));
T(3) = (1+params(9))^(params(8)*(1-params(5)))/(params(7)*(2-params(5)));
T(4) = (params(6)+params(7))/(params(6)+params(7)*y(1));
T(5) = (params(6)+params(7)*y(1))/(params(7)*(1+1/params(3)));
T(6) = T(4)^(1+1/params(3))-1;
T(7) = T(5)*T(6);
T(8) = T(4)^(1-params(5))-1;
T(9) = y(3)^params(4);
T(10) = y(2)^(-params(2));
T(11) = y(3)^(params(4)-1);
T(12) = (1+params(9))/params(1)*((1+y(7))/(1+params(9)))^params(11);
T(13) = (y(2)/(y(2)))^params(12);
T(14) = y(2)^params(2);
T(15) = (T(2)/(1+y(8)))^(-params(5));
T(16) = y(3)/(1-params(10))*T(15);
T(17) = T(16)^params(3);
T(18) = T(2)^(1-params(5));
end
