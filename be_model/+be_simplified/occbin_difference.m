function [binding, relax, err] = occbin_difference(zdatalinear, params, steady_state)
binding.constraint_1 = zdatalinear(:,1)>params(9);
relax.constraint_1 = zdatalinear(:,1)<=params(9);
err.binding_constraint_1 = abs((zdatalinear(:,1))-(params(9)));
err.relax_constraint_1 = abs((zdatalinear(:,1))-(params(9)));
end
