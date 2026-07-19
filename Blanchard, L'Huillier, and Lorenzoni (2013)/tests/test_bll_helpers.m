function tests = test_bll_helpers
tests = functiontests(localfunctions);
end

function testSignalExtractionMatchesLegacyAlgorithm(testCase)
parameterSets = [0.6, 0.5, 1; 0.9426, 1.1977, 1.4738; 0.85, 2, 0.2];
for row = 1:size(parameterSets,1)
    rho = parameterSets(row,1);
    sigU = parameterSets(row,2);
    sigNu = parameterSets(row,3);
    [actualBM, actualFA, actualCM] = bll_signal_extraction(rho, sigU, sigNu);
    [legacyBM, legacyFA, legacyCM] = legacySignalExtraction(rho, sigU, sigNu);
    verifyEqual(testCase, actualBM, legacyBM, 'AbsTol', 2e-11);
    verifyEqual(testCase, actualFA, legacyFA, 'AbsTol', 1e-14);
    verifyEqual(testCase, actualCM, legacyCM, 'AbsTol', 2e-11);
end
end

function testInnovationFactorIsValid(testCase)
[~, ~, CM, diagnostics] = bll_signal_extraction(0.9426, 1.1977, 1.4738);
verifyEqual(testCase, CM*CM', diagnostics.innovationCovariance, ...
    'AbsTol', 1e-11);
verifyLessThan(testCase, diagnostics.residual, 1e-13);
verifyGreaterThan(testCase, diagnostics.iterations, 1);
end

function testSteadyStatePreservesCandidateParameters(testCase)
structuralNames = {'rho','sig_u','sig_nu','cal_w','zet','sig_w'};
derivedNames = {'BM_11','BM_12','BM_21','BM_22','BM_31','BM_32', ...
    'FA_11','FA_12','FA_13','FA_21','FA_22','FA_23', ...
    'CM_11','CM_12','CM_21','CM_22'};
M_ = struct;
M_.param_names = [structuralNames, derivedNames]';
M_.params = [0.7; 0.8; 1.1; 0.82; 2.7; 0.3; zeros(numel(derivedNames),1)];
candidate = M_.params(1:numel(structuralNames));

[ys, params, check] = estimate_steadystate(ones(42,1), [], M_, struct);
verifyEqual(testCase, check, 0);
verifyEqual(testCase, ys, zeros(42,1));
verifyEqual(testCase, params(1:numel(structuralNames)), candidate, 'AbsTol', 0);
verifyNotEqual(testCase, params(numel(structuralNames)+1:end), ...
    zeros(numel(derivedNames),1));

M_.params = params;
M_.params(strcmp('rho', M_.param_names)) = 0.8;
[~, changedParams, changedCheck] = estimate_steadystate(zeros(42,1), [], M_, struct);
verifyEqual(testCase, changedCheck, 0);
verifyNotEqual(testCase, changedParams(strcmp('BM_11', M_.param_names)), ...
    params(strcmp('BM_11', M_.param_names)));
verifyEqual(testCase, changedParams(strcmp('cal_w', M_.param_names)), 0.82);
verifyEqual(testCase, changedParams(strcmp('zet', M_.param_names)), 2.7);
end

function testInvalidCandidateIsRejected(testCase)
names = {'rho','sig_u','sig_nu','BM_11','BM_12','BM_21','BM_22', ...
    'BM_31','BM_32','FA_11','FA_12','FA_13','FA_21','FA_22','FA_23', ...
    'CM_11','CM_12','CM_21','CM_22'};
M_ = struct('param_names', {names'}, 'params', [1.2; 1; 1; zeros(16,1)]);
[~, returned, check] = estimate_steadystate(zeros(2,1), [], M_, struct);
verifyEqual(testCase, check, 1);
verifyEqual(testCase, returned, M_.params);
end

function testStructuralIRFTransformAndBands(testCase)
reducedIRF = [1, -0.5; zeros(11,2)];
structuralIRF = bll_transform_irfs(reducedIRF, 0.8, 1.2, 0.7);
verifySize(testCase, structuralIRF, [12,3]);
verifyTrue(testCase, all(isfinite(structuralIRF), 'all'));

draws = randn(12,3,2,100);
[lower, upper] = bll_hpd_bounds(draws, 0.9);
verifySize(testCase, lower, [12,3,2]);
verifyTrue(testCase, all(lower <= upper, 'all'));
end

function testDynamicWageSlope(testCase)
bet = 0.99;
muW = 0.05;
slope = @(calW, zet) ((1-calW*bet)*(1-calW)) ...
    /(calW*(1+bet)*(1+zet*(1+1/muW)));
verifyNotEqual(testCase, slope(0.66,2), slope(0.86,2));
verifyNotEqual(testCase, slope(0.86,2), slope(0.86,3));
end

function [BM, FA, CM] = legacySignalExtraction(rho, sigU, sigNu)
A = [1+rho, -rho, 0; 1, 0, 0; 0, 0, rho];
S1 = diag([(1-rho)*sigU, 0, sqrt(rho)*sigU])^2;
F = [1, 0, 1; 1, 0, 0];
S2 = diag([0, sigNu])^2;
P = 0.01*eye(3);
K = P*F'/(F*P*F'+S2);
for iteration = 1:100000
    P = A*P*A' - A*P*F'*K'*A' + S1;
    nextK = P*F'/(F*P*F'+S2);
    difference = norm(nextK-K, Inf);
    K = nextK;
    if difference < 1e-15
        break
    end
end
FA = F*A;
CM = real((F*P*F'+S2)^0.5);
BM = K*CM;
end
