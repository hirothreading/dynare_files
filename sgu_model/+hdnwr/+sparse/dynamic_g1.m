function [g1, T_order, T] = dynamic_g1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 9
    T_order = -1;
    T = NaN(21, 1);
end
[T_order, T] = hdnwr.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
g1_v = NaN(34, 1);
g1_v(1)=(-((1+y(17))*(-y(15))/(y(5)*y(5))));
g1_v(2)=(-(T(2)/(1+y(17))));
g1_v(3)=(-(params(13)*1/y(9)));
g1_v(4)=(-(params(14)*1/y(10)));
g1_v(5)=T(9)*y(13)/(1-params(10))*T(1)*params(7)/(1+y(18))*T(19)*T(20)-y(5)*T(1)*params(7)/(1+y(17));
g1_v(6)=(-(T(13)+y(11)*T(1)*params(7)*getPowerDeriv(T(2),1-params(5),1)+T(3)*(-(params(7)*getPowerDeriv(params(6)+params(7)*y(11),2-params(5),1)))));
g1_v(7)=(-(((y(11)+T(17))*(1-params(10))*(T(21)-(T(18)*params(7)/(params(7)*(1-params(5)))+(params(6)+params(7)*y(11))/(params(7)*(1-params(5)))*(-(params(7)*(params(6)+params(7))))/((params(6)+params(7)*y(11))*(params(6)+params(7)*y(11)))*getPowerDeriv(T(14),1-params(5),1)))-(1-params(10))*(T(17)-(params(6)+params(7)*y(11))/(params(7)*(1-params(5)))*T(18))*(1+T(21)))/((y(11)+T(17))*(y(11)+T(17)))));
g1_v(8)=1;
g1_v(9)=getPowerDeriv(y(12),(-params(2)),1);
g1_v(10)=(-(y(19)*T(7)*1/(steady_state(2))*getPowerDeriv(y(12)/(steady_state(2)),params(12),1)));
g1_v(11)=T(12)*getPowerDeriv(y(12),params(2),1);
g1_v(12)=(-(y(20)*getPowerDeriv(y(13),params(4),1)));
g1_v(13)=y(20)*params(4)*getPowerDeriv(y(13),params(4)-1,1);
g1_v(14)=T(9)*T(20)*T(10)*1/(1-params(10));
g1_v(15)=1;
g1_v(16)=(-1);
g1_v(17)=(-((1+y(17))*1/y(5)));
g1_v(18)=(-(params(1)*T(5)/(1+y(27))));
g1_v(19)=1;
g1_v(20)=(-(y(19)*T(8)*(1+params(9))/params(1)*1/(1+params(9))*getPowerDeriv((1+y(17))/(1+params(9)),params(11),1)));
g1_v(21)=(-(y(15)/y(5)));
g1_v(22)=(-((-(T(2)*y(5)))/((1+y(17))*(1+y(17)))));
g1_v(23)=1;
g1_v(24)=T(9)*T(20)*y(13)/(1-params(10))*T(19)*(-T(2))/((1+y(18))*(1+y(18)));
g1_v(25)=getPowerDeriv(1+y(18),1-params(5),1);
g1_v(26)=(-(T(7)*T(8)));
g1_v(27)=1/y(19);
g1_v(28)=(-T(4));
g1_v(29)=params(4)*T(6);
g1_v(30)=1/y(20);
g1_v(31)=(-(params(1)*(1+y(16))*getPowerDeriv(y(22),(-params(2)),1)/(1+y(27))));
g1_v(32)=(-((-(params(1)*(1+y(16))*T(5)))/((1+y(27))*(1+y(27)))));
g1_v(33)=(-1);
g1_v(34)=(-1);
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 10, 32);
end
