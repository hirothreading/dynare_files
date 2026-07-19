%% Structural impulse responses at the published posterior means
clearvars;
close all;
clc;

dynare simulate.mod;

periods = 100;
reducedIRF = bll_collect_irfs(M_, oo_, periods);
variableNames = cellstr(M_.endo_names);
shockNames = {'Permanent technology','Transitory technology','Noise', ...
    'Investment-specific','Price markup','Wage markup','Monetary','Fiscal'};
structuralIRF = zeros(periods, numel(shockNames), numel(variableNames));

parameterNames = cellstr(M_.param_names);
rho = M_.params(strcmp('rho', parameterNames));
sigU = M_.params(strcmp('sig_u', parameterNames));
sigNu = M_.params(strcmp('sig_nu', parameterNames));

for variable = 1:numel(variableNames)
    structuralIRF(:,1:3,variable) = bll_transform_irfs( ...
        reducedIRF(:,1:2,variable), rho, sigU, sigNu);
    structuralIRF(:,4:8,variable) = reducedIRF(:,3:7,variable);
end

save('simulation_irfs.mat', 'structuralIRF', 'reducedIRF', ...
    'variableNames', 'shockNames', 'M_');

plotBLLFigure(structuralIRF, variableNames, {'ca','ia','ya','n'}, ...
    shockNames(1:3), 20, 'BLL_Figure5_posterior_mean.pdf');
plotBLLFigure(structuralIRF, variableNames, {'r','pi','wa'}, ...
    shockNames(1:3), 20, 'BLL_Figure6_posterior_mean.pdf');

function plotBLLFigure(irfs, allVariables, selectedVariables, rowNames, periods, outputFile)
figureHandle = figure('Color', 'w');
layout = tiledlayout(numel(rowNames), numel(selectedVariables), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for shock = 1:numel(rowNames)
    for variable = 1:numel(selectedVariables)
        index = find(strcmp(selectedVariables{variable}, allVariables), 1);
        nexttile;
        plot(1:periods, irfs(1:periods,shock,index), 'b-', 'LineWidth', 1.2);
        yline(0, 'k:');
        if shock == 1
            title(selectedVariables{variable}, 'Interpreter', 'none');
        end
        if variable == 1
            ylabel(rowNames{shock});
        end
    end
end
title(layout, 'BLL (2013): structural impulse responses');
exportgraphics(figureHandle, outputFile, 'ContentType', 'vector');
end
