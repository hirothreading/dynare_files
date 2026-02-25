function [y, T, residual, g1] = static_2(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(1, 1);
  T(2)=log(y(10));
  residual(1)=(T(2))-(T(2)*params(14)+x(2));
if nargout > 3
    g1_v = NaN(1, 1);
g1_v(1)=1/y(10)-params(14)*1/y(10);
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 1, 1);
end
end
