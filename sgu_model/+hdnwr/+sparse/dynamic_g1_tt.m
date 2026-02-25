function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = hdnwr.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 21
    T = [T; NaN(21 - size(T, 1), 1)];
end
T(19) = getPowerDeriv(T(2)/(1+y(18)),(-params(5)),1);
T(20) = getPowerDeriv(T(11),params(3),1);
T(21) = T(16)*params(7)/(params(7)*(1+1/params(3)))+T(15)*(-(params(7)*(params(6)+params(7))))/((params(6)+params(7)*y(11))*(params(6)+params(7)*y(11)))*getPowerDeriv(T(14),1+1/params(3),1);
end
