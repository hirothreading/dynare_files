function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(2, 1);
  y(10)=y(9)*params(3)+y(13);
  residual(1)=((y(9)-(params(4)*(y(8)+params(10)*y(14))+params(5)*y(12)+y(16)*params(2)))*(1-params(15))+params(15)*(y(9)-(y(16)*params(2)+(y(8)+params(10)*y(14))*params(6)-params(8)+y(12)*params(7))))-(0);
  residual(2)=(y(8)-y(11))-(y(15)-y(18)-1/params(1)*(y(10)-y(16)));
if nargout > 3
    g1_v = NaN(4, 1);
g1_v(1)=1;
g1_v(2)=1/params(1)*params(3);
g1_v(3)=(1-params(15))*(-params(4))+params(15)*(-params(6));
g1_v(4)=1;
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 2, 2);
end
end
