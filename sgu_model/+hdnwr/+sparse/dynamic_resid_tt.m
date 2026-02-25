function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 18
    T = [T; NaN(18 - size(T, 1), 1)];
end
T(1) = (1+params(9))^params(8);
T(2) = T(1)*(params(6)+params(7)*y(11));
T(3) = (1+params(9))^(params(8)*(1-params(5)))/(params(7)*(2-params(5)));
T(4) = y(13)^params(4);
T(5) = y(22)^(-params(2));
T(6) = y(13)^(params(4)-1);
T(7) = (1+params(9))/params(1)*((1+y(17))/(1+params(9)))^params(11);
T(8) = (y(12)/(steady_state(2)))^params(12);
T(9) = y(12)^params(2);
T(10) = (T(2)/(1+y(18)))^(-params(5));
T(11) = y(13)/(1-params(10))*T(10);
T(12) = T(11)^params(3);
T(13) = T(2)^(1-params(5));
T(14) = (params(6)+params(7))/(params(6)+params(7)*y(11));
T(15) = (params(6)+params(7)*y(11))/(params(7)*(1+1/params(3)));
T(16) = T(14)^(1+1/params(3))-1;
T(17) = T(15)*T(16);
T(18) = T(14)^(1-params(5))-1;
end
