function [y, T, residual, g1] = static_3(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(7, 1);
  residual(1)=(y(10)*params(4)*y(3)^(params(4)-1))-(y(5));
  T(3)=(1+params(9))/params(1)*((1+y(7))/(1+params(9)))^params(11);
  T(4)=(y(2)/(y(2)))^params(12);
  residual(2)=(1+y(6))-(T(3)*T(4)*y(9));
  residual(3)=(1+y(8))-(1+y(7));
  T(5)=(1+params(9))^params(8);
  T(6)=T(5)*(params(6)+params(7)*y(1));
  T(7)=y(2)^params(2);
  T(8)=(T(6)/(1+y(8)))^(-params(5));
  T(9)=y(3)/(1-params(10))*T(8);
  T(10)=T(9)^params(3);
  residual(4)=(T(7)*T(10))-(T(6)*y(5)/(1+y(7)));
  T(11)=(1+params(9))^(params(8)*(1-params(5)))/(params(7)*(2-params(5)));
  T(12)=T(6)^(1-params(5));
  residual(5)=((1+y(8))^(1-params(5)))-(T(11)*((params(6)+params(7))^(2-params(5))-(params(6)+params(7)*y(1))^(2-params(5)))+y(1)*T(12));
  T(13)=y(2)^(-params(2));
  residual(6)=(T(13))-(T(13)*params(1)*(1+y(6))/(1+y(7)));
  residual(7)=(y(2))-(y(10)*y(3)^params(4));
  T(14)=getPowerDeriv(T(6)/(1+y(8)),(-params(5)),1);
  T(15)=getPowerDeriv(T(9),params(3),1);
  T(16)=getPowerDeriv(y(2),(-params(2)),1);
if nargout > 3
    g1_v = NaN(20, 1);
g1_v(1)=(-1);
g1_v(2)=(-(T(6)/(1+y(7))));
g1_v(3)=(-(y(9)*T(3)*((y(2))-y(2))/((y(2))*(y(2)))*getPowerDeriv(y(2)/(y(2)),params(12),1)));
g1_v(4)=T(10)*getPowerDeriv(y(2),params(2),1);
g1_v(5)=T(16)-params(1)*(1+y(6))*T(16)/(1+y(7));
g1_v(6)=1;
g1_v(7)=(-(y(9)*T(4)*(1+params(9))/params(1)*1/(1+params(9))*getPowerDeriv((1+y(7))/(1+params(9)),params(11),1)));
g1_v(8)=(-1);
g1_v(9)=(-((-(T(6)*y(5)))/((1+y(7))*(1+y(7)))));
g1_v(10)=(-((-(T(13)*params(1)*(1+y(6))))/((1+y(7))*(1+y(7)))));
g1_v(11)=T(7)*y(3)/(1-params(10))*T(5)*params(7)/(1+y(8))*T(14)*T(15)-y(5)*T(5)*params(7)/(1+y(7));
g1_v(12)=(-(T(11)*(-(params(7)*getPowerDeriv(params(6)+params(7)*y(1),2-params(5),1)))+T(12)+y(1)*T(5)*params(7)*getPowerDeriv(T(6),1-params(5),1)));
g1_v(13)=1;
g1_v(14)=T(7)*T(15)*y(3)/(1-params(10))*T(14)*(-T(6))/((1+y(8))*(1+y(8)));
g1_v(15)=getPowerDeriv(1+y(8),1-params(5),1);
g1_v(16)=1;
g1_v(17)=(-(T(13)*params(1)/(1+y(7))));
g1_v(18)=y(10)*params(4)*getPowerDeriv(y(3),params(4)-1,1);
g1_v(19)=T(7)*T(15)*T(8)*1/(1-params(10));
g1_v(20)=(-(y(10)*getPowerDeriv(y(3),params(4),1)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 7, 7);
end
end
