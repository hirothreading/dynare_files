function [y, T] = dynamic_4(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  T(14)=(params(6)+params(7)*y(11))/(params(7)*(1+1/params(3)))*(((params(6)+params(7))/(params(6)+params(7)*y(11)))^(1+1/params(3))-1);
  y(14)=params(10)+(1-params(10))*(T(14)-(params(6)+params(7)*y(11))/(params(7)*(1-params(5)))*(((params(6)+params(7))/(params(6)+params(7)*y(11)))^(1-params(5))-1))/(y(11)+T(14));
end
