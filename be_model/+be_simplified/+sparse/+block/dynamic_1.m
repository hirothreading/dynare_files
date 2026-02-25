function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(11)=params(11)*y(4)+x(1);
  y(12)=params(12)*y(5)+x(2);
  y(13)=params(13)*y(6)+x(3);
  y(14)=params(14)*y(7)+x(4);
end
