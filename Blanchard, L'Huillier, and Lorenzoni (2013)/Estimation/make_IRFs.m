%% Draw-by-draw structural IRFs and variance decompositions
% Run estimate.mod with MH_REPLIC>0 first. Dynare stores paired parameter
% draws and reduced-form IRFs in Estimation/estimate/metropolis.
clearvars;
close all;
clc;

if exist('dynare_config', 'file') == 2
    dynare_config;
end
if exist('get_posterior_irf', 'file') ~= 2
    error('BLL:DynarePathMissing', ...
        'Add the Dynare matlab directory to the MATLAB path before running this script.');
end

if ~isfile('estimate_results.mat')
    error('BLL:MissingEstimationResults', ...
        'Run estimate.mod with posterior sampling before make_IRFs.m.');
end
% Dynare's posterior-IRF API reads these standard session globals.
global M_ oo_ options_ estim_params_ %#ok<GVMIS>
results = load('estimate_results.mat', 'M_', 'oo_', 'options_', 'estim_params_');
M_ = results.M_;
oo_ = results.oo_;
options_ = results.options_;
estim_params_ = results.estim_params_;
addpath('../Model');

periods = options_.irf;
coverage = options_.mh_conf_sig;
variableNames = {'ca','ia','ya','n','r','pi','wa'};
reducedShockNames = cellstr(M_.exo_names);
structuralShockNames = {'Permanent technology','Transitory technology','Noise', ...
    'Investment-specific','Price markup','Wage markup','Monetary','Fiscal'};

if ~isequal(reducedShockNames(:)', ...
        {'e_1','e_2','eps_d','eps_p','eps_w','eps_q','eps_g'})
    error('BLL:ShockOrderChanged', ...
        'The structural transformation requires the documented shock order.');
end

probe = get_posterior_irf(variableNames{1}, 'e_1');
drawCount = numel(probe);
if drawCount < 2
    error('BLL:InsufficientDraws', 'Posterior IRF files contain fewer than two draws.');
end

rhoDraws = zeros(drawCount,1);
sigUDraws = zeros(drawCount,1);
sigNuDraws = zeros(drawCount,1);
parameterNames = cellstr(M_.param_names);
for draw = 1:drawCount
    drawModel = set_parameters_locally(M_, probe(draw).draw);
    rhoDraws(draw) = drawModel.params(strcmp('rho', parameterNames));
    sigUDraws(draw) = drawModel.params(strcmp('sig_u', parameterNames));
    sigNuDraws(draw) = drawModel.params(strcmp('sig_nu', parameterNames));
end

structuralIRFDraws = zeros(periods, 8, numel(variableNames), drawCount);
for variable = 1:numel(variableNames)
    firstEFI = get_posterior_irf(variableNames{variable}, 'e_1');
    secondEFI = get_posterior_irf(variableNames{variable}, 'e_2');
    validateDrawCount(firstEFI, drawCount, variableNames{variable}, 'e_1');
    validateDrawCount(secondEFI, drawCount, variableNames{variable}, 'e_2');
    for draw = 1:drawCount
        reducedTechnologyIRF = [firstEFI(draw).irf, secondEFI(draw).irf];
        structuralIRFDraws(:,1:3,variable,draw) = bll_transform_irfs( ...
            reducedTechnologyIRF, rhoDraws(draw), sigUDraws(draw), sigNuDraws(draw));
    end

    for reducedShock = 3:7
        draws = get_posterior_irf(variableNames{variable}, ...
            reducedShockNames{reducedShock});
        validateDrawCount(draws, drawCount, variableNames{variable}, ...
            reducedShockNames{reducedShock});
        for draw = 1:drawCount
            structuralIRFDraws(:,reducedShock+1,variable,draw) = draws(draw).irf;
        end
    end
end

IRF_mean = mean(structuralIRFDraws, 4);
IRF_median = median(structuralIRFDraws, 4);
[IRF_low, IRF_high] = bll_hpd_bounds(structuralIRFDraws, coverage);
if any(IRF_low > IRF_high, 'all')
    error('BLL:InvalidPosteriorBands', 'Computed HPD bounds are not ordered.');
end

% Table 6: forecast-error variance shares are computed for every posterior
% draw after replacing the two EFI shocks with the three structural shocks.
horizons = [1, 4, 8, 12];
decompositionVariables = {'ca','ia','ya'};
varianceShareDraws = zeros(numel(horizons), 8, ...
    numel(decompositionVariables), drawCount);
for variable = 1:numel(decompositionVariables)
    variableIndex = find(strcmp(decompositionVariables{variable}, variableNames), 1);
    for draw = 1:drawCount
        for horizon = 1:numel(horizons)
            contributions = sum(structuralIRFDraws( ...
                1:horizons(horizon),:,variableIndex,draw).^2, 1);
            varianceShareDraws(horizon,:,variable,draw) = ...
                contributions/sum(contributions);
        end
    end
end
varianceShareMean = mean(varianceShareDraws, 4);
varianceShareMedian = median(varianceShareDraws, 4);
[varianceShareLow, varianceShareHigh] = ...
    bll_hpd_bounds(varianceShareDraws, coverage);

save('computed_IRFs.mat', 'structuralIRFDraws', 'IRF_mean', 'IRF_median', ...
    'IRF_low', 'IRF_high', 'varianceShareDraws', 'varianceShareMean', ...
    'varianceShareMedian', 'varianceShareLow', 'varianceShareHigh', ...
    'variableNames', 'structuralShockNames', 'decompositionVariables', ...
    'horizons', 'coverage', '-v7.3');

writeVarianceTable(varianceShareMean, horizons, decompositionVariables, ...
    structuralShockNames, 'BLL_Table6_posterior_mean.csv');
plotBLLPosteriorFigure(IRF_mean, IRF_low, IRF_high, variableNames, ...
    {'ca','ia','ya','n'}, structuralShockNames(1:3), 20, ...
    coverage, 'BLL_Figure5_posterior.pdf');
plotBLLPosteriorFigure(IRF_mean, IRF_low, IRF_high, variableNames, ...
    {'r','pi','wa'}, structuralShockNames(1:3), 20, ...
    coverage, 'BLL_Figure6_posterior.pdf');

function validateDrawCount(draws, expected, variable, shock)
if numel(draws) ~= expected
    error('BLL:PosteriorDrawMismatch', ...
        'Expected %d draws for %s/%s but found %d.', ...
        expected, variable, shock, numel(draws));
end
end

function writeVarianceTable(shares, horizons, variables, shocks, outputFile)
rowCount = numel(horizons)*numel(variables);
variableColumn = strings(rowCount,1);
horizonColumn = zeros(rowCount,1);
values = zeros(rowCount,numel(shocks));
row = 0;
for variable = 1:numel(variables)
    for horizon = 1:numel(horizons)
        row = row+1;
        variableColumn(row) = variables{variable};
        horizonColumn(row) = horizons(horizon);
        values(row,:) = shares(horizon,:,variable);
    end
end
result = table(variableColumn, horizonColumn, 'VariableNames', {'Variable','Quarter'});
for shock = 1:numel(shocks)
    columnName = matlab.lang.makeValidName(shocks{shock});
    result.(columnName) = values(:,shock);
end
writetable(result, outputFile);
end

function plotBLLPosteriorFigure(point, lower, upper, allVariables, ...
        selectedVariables, rowNames, periods, coverage, outputFile)
figureHandle = figure('Color', 'w');
layout = tiledlayout(numel(rowNames), numel(selectedVariables), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for shock = 1:numel(rowNames)
    for variable = 1:numel(selectedVariables)
        index = find(strcmp(selectedVariables{variable}, allVariables), 1);
        nexttile;
        plot(1:periods, lower(1:periods,shock,index), 'r--');
        hold on;
        plot(1:periods, point(1:periods,shock,index), 'b-', 'LineWidth', 1.2);
        plot(1:periods, upper(1:periods,shock,index), 'r--');
        yline(0, 'k:');
        if shock == 1
            title(selectedVariables{variable}, 'Interpreter', 'none');
        end
        if variable == 1
            ylabel(rowNames{shock});
        end
    end
end
title(layout, sprintf('BLL (2013): %.0f%% posterior HPD intervals', 100*coverage));
exportgraphics(figureHandle, outputFile, 'ContentType', 'vector');
end
