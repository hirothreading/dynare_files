function [y, T, residual, g1] = dynamic_3(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(7, 1);
  T(1)=(1+params(9))/params(1)*((1+y(17))/(1+params(9)))^params(11);
  T(2)=(y(12)/(steady_state(2)))^params(12);
  residual(1)=(1+y(16))-(T(1)*T(2)*y(19));
  residual(2)=(1+y(18))-((1+y(17))*y(15)/y(5));
  T(3)=(1+params(9))^params(8);
  T(4)=T(3)*(params(6)+params(7)*y(11));
  T(5)=(1+params(9))^(params(8)*(1-params(5)))/(params(7)*(2-params(5)));
  T(6)=T(4)^(1-params(5));
  residual(3)=((1+y(18))^(1-params(5)))-(y(11)*T(6)+T(5)*((params(6)+params(7))^(2-params(5))-(params(6)+params(7)*y(11))^(2-params(5))));
  residual(4)=(y(12))-(y(20)*y(13)^params(4));
  residual(5)=(y(20)*params(4)*y(13)^(params(4)-1))-(y(15));
  T(7)=y(12)^params(2);
  T(8)=(T(4)/(1+y(18)))^(-params(5));
  T(9)=y(13)/(1-params(10))*T(8);
  T(10)=T(9)^params(3);
  residual(6)=(T(7)*T(10))-(T(4)*y(5)/(1+y(17)));
  T(11)=y(22)^(-params(2));
  residual(7)=(y(12)^(-params(2)))-(params(1)*(1+y(16))*T(11)/(1+y(27)));
  T(12)=getPowerDeriv(T(4)/(1+y(18)),(-params(5)),1);
  T(13)=getPowerDeriv(T(9),params(3),1);
if nargout > 3
    g1_v = NaN(23, 1);
g1_v(1)=(-((1+y(17))*(-y(15))/(y(5)*y(5))));
g1_v(2)=(-(T(4)/(1+y(17))));
g1_v(3)=1;
g1_v(4)=(-(params(1)*T(11)/(1+y(27))));
g1_v(5)=1;
g1_v(6)=getPowerDeriv(1+y(18),1-params(5),1);
g1_v(7)=T(7)*T(13)*y(13)/(1-params(10))*T(12)*(-T(4))/((1+y(18))*(1+y(18)));
g1_v(8)=(-(T(6)+y(11)*T(3)*params(7)*getPowerDeriv(T(4),1-params(5),1)+T(5)*(-(params(7)*getPowerDeriv(params(6)+params(7)*y(11),2-params(5),1)))));
g1_v(9)=T(7)*y(13)/(1-params(10))*T(3)*params(7)/(1+y(18))*T(12)*T(13)-y(5)*T(3)*params(7)/(1+y(17));
g1_v(10)=(-(y(20)*getPowerDeriv(y(13),params(4),1)));
g1_v(11)=y(20)*params(4)*getPowerDeriv(y(13),params(4)-1,1);
g1_v(12)=T(7)*T(13)*T(8)*1/(1-params(10));
g1_v(13)=(-((1+y(17))*1/y(5)));
g1_v(14)=(-1);
g1_v(15)=(-(y(19)*T(2)*(1+params(9))/params(1)*1/(1+params(9))*getPowerDeriv((1+y(17))/(1+params(9)),params(11),1)));
g1_v(16)=(-(y(15)/y(5)));
g1_v(17)=(-((-(T(4)*y(5)))/((1+y(17))*(1+y(17)))));
g1_v(18)=(-(y(19)*T(1)*1/(steady_state(2))*getPowerDeriv(y(12)/(steady_state(2)),params(12),1)));
g1_v(19)=1;
g1_v(20)=T(10)*getPowerDeriv(y(12),params(2),1);
g1_v(21)=getPowerDeriv(y(12),(-params(2)),1);
g1_v(22)=(-((-(params(1)*(1+y(16))*T(11)))/((1+y(27))*(1+y(27)))));
g1_v(23)=(-(params(1)*(1+y(16))*getPowerDeriv(y(22),(-params(2)),1)/(1+y(27))));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 7, 21);
end
end
