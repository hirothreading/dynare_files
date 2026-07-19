function [ys, params, check] = estimate_steadystate(ys, ~, M_, ~)
%ESTIMATE_STEADYSTATE Update only signal-extraction matrix coefficients.
%   Estimated structural parameters enter through M_.params and are copied
%   back unchanged. BM, FA and CM are derived from the current candidate
%   values of rho, sig_u and sig_nu at each likelihood evaluation.

params = M_.params;
ys(:) = 0;
check = 0;

try
    rho = getParameter(M_, 'rho');
    sigU = getParameter(M_, 'sig_u');
    sigNu = getParameter(M_, 'sig_nu');
    [BM, FA, CM] = bll_signal_extraction(rho, sigU, sigNu);

    derivedNames = { ...
        'BM_11','BM_12','BM_21','BM_22','BM_31','BM_32', ...
        'FA_11','FA_12','FA_13','FA_21','FA_22','FA_23', ...
        'CM_11','CM_12','CM_21','CM_22'};
    derivedValues = [ ...
        BM(1,1),BM(1,2),BM(2,1),BM(2,2),BM(3,1),BM(3,2), ...
        FA(1,1),FA(1,2),FA(1,3),FA(2,1),FA(2,2),FA(2,3), ...
        CM(1,1),CM(1,2),CM(2,1),CM(2,2)];

    for index = 1:numel(derivedNames)
        parameterIndex = find(strcmp(derivedNames{index}, M_.param_names), 1);
        if isempty(parameterIndex)
            error('BLL:MissingParameter', ...
                'Missing derived parameter %s.', derivedNames{index});
        end
        params(parameterIndex) = derivedValues(index);
    end
catch exception
    if startsWith(exception.identifier, 'BLL:') || ...
            startsWith(exception.identifier, 'MATLAB:validators:')
        check = 1;
        return
    end
    rethrow(exception)
end
end

function value = getParameter(M_, name)
index = find(strcmp(name, M_.param_names), 1);
if isempty(index)
    error('BLL:MissingParameter', 'Missing structural parameter %s.', name);
end
value = M_.params(index);
end
