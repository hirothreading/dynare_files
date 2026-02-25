function [T_order, T] = dynamic_g2_tt(y, x, params, steady_state, T_order, T)
if T_order >= 2
    return
end
[T_order, T] = hdnwr_nl.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
T_order = 2;
if size(T, 1) < 47
    T = [T; NaN(47 - size(T, 1), 1)];
end
T(44) = getPowerDeriv(T(2)/(1+y(18)),(-params(5)),2);
T(45) = getPowerDeriv(T(13),params(3),2);
T(46) = (-((-(params(7)*(params(6)+params(7))))*(params(7)*(params(6)+params(7)*y(11))+params(7)*(params(6)+params(7)*y(11)))))/((params(6)+params(7)*y(11))*(params(6)+params(7)*y(11))*(params(6)+params(7)*y(11))*(params(6)+params(7)*y(11)));
T(47) = T(27)*T(30)+T(27)*T(30)+T(17)*(T(29)*T(46)+T(28)*T(28)*getPowerDeriv(T(16),T(4),2));
end
