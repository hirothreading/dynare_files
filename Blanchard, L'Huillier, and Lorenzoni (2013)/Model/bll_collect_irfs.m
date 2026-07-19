function reducedIRF = bll_collect_irfs(M_, oo_, periods)
%BLL_COLLECT_IRFS Convert oo_.irfs into a dense period-shock-variable array.

arguments
    M_ (1,1) struct
    oo_ (1,1) struct
    periods (1,1) double {mustBeInteger, mustBePositive}
end

variableNames = cellstr(M_.endo_names);
shockNames = cellstr(M_.exo_names);
reducedIRF = zeros(periods, numel(shockNames), numel(variableNames));

for variable = 1:numel(variableNames)
    for shock = 1:numel(shockNames)
        fieldName = [variableNames{variable}, '_', shockNames{shock}];
        if isfield(oo_.irfs, fieldName)
            values = oo_.irfs.(fieldName);
            count = min(periods, numel(values));
            reducedIRF(1:count,shock,variable) = values(1:count);
        end
    end
end
end
