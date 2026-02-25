function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 21
    T = [T; NaN(21 - size(T, 1), 1)];
end
T(1) = (1+params(9))^params(8);
T(2) = T(1)*(params(6)+params(7)*y(11));
T(3) = (1+params(9))^(params(8)*(1-params(5)))/(params(7)*(2-params(5)));
T(4) = 1+1/params(3);
T(5) = y(13)^params(4);
T(6) = y(22)^(-params(2));
T(7) = y(13)^(params(4)-1);
T(8) = (1+params(9))/params(1)*((1+y(17))/(1+params(9)))^params(11);
T(9) = (y(12)/(steady_state(2)))^params(12);
T(10) = y(13)/(1-params(10));
T(11) = y(12)^params(2);
T(12) = (T(2)/(1+y(18)))^(-params(5));
T(13) = T(10)*T(12);
T(14) = T(13)^params(3);
T(15) = T(2)^(1-params(5));
T(16) = (params(6)+params(7))/(params(6)+params(7)*y(11));
T(17) = (params(6)+params(7)*y(11))/(params(7)*T(4));
T(18) = T(16)^T(4)-1;
T(19) = T(17)*T(18);
T(20) = (params(6)+params(7)*y(11))/(params(7)*(1-params(5)));
T(21) = T(16)^(1-params(5))-1;
end
