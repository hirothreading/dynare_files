function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = hdnwr_nl.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 43
    T = [T; NaN(43 - size(T, 1), 1)];
end
T(22) = T(1)*params(7)/(1+y(18));
T(23) = getPowerDeriv(T(2)/(1+y(18)),(-params(5)),1);
T(24) = T(10)*T(22)*T(23);
T(25) = getPowerDeriv(T(13),params(3),1);
T(26) = T(1)*params(7)*getPowerDeriv(T(2),1-params(5),1);
T(27) = params(7)/(params(7)*T(4));
T(28) = (-(params(7)*(params(6)+params(7))))/((params(6)+params(7)*y(11))*(params(6)+params(7)*y(11)));
T(29) = getPowerDeriv(T(16),T(4),1);
T(30) = T(28)*T(29);
T(31) = T(18)*T(27)+T(17)*T(30);
T(32) = getPowerDeriv(T(16),1-params(5),1);
T(33) = T(28)*T(32);
T(34) = (1-params(10))*(T(31)-(T(21)*params(7)/(params(7)*(1-params(5)))+T(20)*T(33)));
T(35) = 1/(steady_state(2))*getPowerDeriv(y(12)/(steady_state(2)),params(12),1);
T(36) = getPowerDeriv(y(12),params(2),1);
T(37) = getPowerDeriv(y(22),(-params(2)),1);
T(38) = getPowerDeriv(y(13),params(4),1);
T(39) = getPowerDeriv(y(13),params(4)-1,1);
T(40) = T(12)*1/(1-params(10));
T(41) = (1+params(9))/params(1)*1/(1+params(9))*getPowerDeriv((1+y(17))/(1+params(9)),params(11),1);
T(42) = (-T(2))/((1+y(18))*(1+y(18)));
T(43) = T(10)*T(23)*T(42);
end
