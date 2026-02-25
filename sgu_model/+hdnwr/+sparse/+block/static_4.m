function [y, T] = static_4(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  T(17)=(params(6)+params(7)*y(1))/(params(7)*(1+1/params(3)))*(((params(6)+params(7))/(params(6)+params(7)*y(1)))^(1+1/params(3))-1);
  y(4)=params(10)+(1-params(10))*(T(17)-(params(6)+params(7)*y(1))/(params(7)*(1-params(5)))*(((params(6)+params(7))/(params(6)+params(7)*y(1)))^(1-params(5))-1))/(y(1)+T(17));
end
