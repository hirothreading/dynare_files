function [g1, T_order, T] = static_g1(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 8
    T_order = -1;
    T = NaN(22, 1);
end
[T_order, T] = hdnwr_nl.sparse.static_g1_tt(y, x, params, T_order, T);
g1_v = NaN(27, 1);
g1_v(1)=T(14)*y(3)/(1-params(10))*T(1)*params(7)/(1+y(8))*T(19)*T(20)-y(5)*T(1)*params(7)/(1+y(7));
g1_v(2)=(-(T(3)*(-(params(7)*getPowerDeriv(params(6)+params(7)*y(1),2-params(5),1)))+T(18)+y(1)*T(1)*params(7)*getPowerDeriv(T(2),1-params(5),1)));
g1_v(3)=(-(((y(1)+T(7))*(1-params(10))*(T(21)-(T(8)*params(7)/(params(7)*(1-params(5)))+(params(6)+params(7)*y(1))/(params(7)*(1-params(5)))*(-(params(7)*(params(6)+params(7))))/((params(6)+params(7)*y(1))*(params(6)+params(7)*y(1)))*getPowerDeriv(T(4),1-params(5),1)))-(1-params(10))*(T(7)-(params(6)+params(7)*y(1))/(params(7)*(1-params(5)))*T(8))*(1+T(21)))/((y(1)+T(7))*(y(1)+T(7)))));
g1_v(4)=1;
g1_v(5)=T(22)-params(1)*(1+y(6))*T(22)/(1+y(7));
g1_v(6)=(-(y(9)*T(12)*((y(2))-y(2))/((y(2))*(y(2)))*getPowerDeriv(y(2)/(y(2)),params(12),1)));
g1_v(7)=T(17)*getPowerDeriv(y(2),params(2),1);
g1_v(8)=(-(y(10)*getPowerDeriv(y(3),params(4),1)));
g1_v(9)=y(10)*params(4)*getPowerDeriv(y(3),params(4)-1,1);
g1_v(10)=T(14)*T(20)*T(15)*1/(1-params(10));
g1_v(11)=1;
g1_v(12)=(-1);
g1_v(13)=(-(T(2)/(1+y(7))));
g1_v(14)=(-(T(10)*params(1)/(1+y(7))));
g1_v(15)=1;
g1_v(16)=(-((-(T(10)*params(1)*(1+y(6))))/((1+y(7))*(1+y(7)))));
g1_v(17)=(-(y(9)*T(13)*(1+params(9))/params(1)*1/(1+params(9))*getPowerDeriv((1+y(7))/(1+params(9)),params(11),1)));
g1_v(18)=(-1);
g1_v(19)=(-((-(T(2)*y(5)))/((1+y(7))*(1+y(7)))));
g1_v(20)=1;
g1_v(21)=T(14)*T(20)*y(3)/(1-params(10))*T(19)*(-T(2))/((1+y(8))*(1+y(8)));
g1_v(22)=getPowerDeriv(1+y(8),1-params(5),1);
g1_v(23)=(-(T(12)*T(13)));
g1_v(24)=1/y(9)-params(13)*1/y(9);
g1_v(25)=(-T(9));
g1_v(26)=params(4)*T(11);
g1_v(27)=1/y(10)-params(14)*1/y(10);
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 10, 10);
end
