function [ys, params, check] = simulate_steadystate(ys, exo, M_, options_)
%SIMULATE_STEADYSTATE Use the same parameter-safe logic as estimation.
[ys, params, check] = estimate_steadystate(ys, exo, M_, options_);
end
