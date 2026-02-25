function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = hdnwr.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 22
    T = [T; NaN(22 - size(T, 1), 1)];
end
T(19) = getPowerDeriv(T(2)/(1+y(8)),(-params(5)),1);
T(20) = getPowerDeriv(T(16),params(3),1);
T(21) = T(6)*params(7)/(params(7)*(1+1/params(3)))+T(5)*(-(params(7)*(params(6)+params(7))))/((params(6)+params(7)*y(1))*(params(6)+params(7)*y(1)))*getPowerDeriv(T(4),1+1/params(3),1);
T(22) = getPowerDeriv(y(2),(-params(2)),1);
end
