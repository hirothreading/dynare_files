function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 0
    T = [T; NaN(0 - size(T, 1), 1)];
end
end
